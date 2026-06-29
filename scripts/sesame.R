#!/usr/bin/env Rscript

# Description: Generate beta-value matrices and auditable probe-level QC reports
# from Illumina IDAT files using SeSAMe.

suppressPackageStartupMessages({
    library(logger)
    library(sesame)
    library(data.table)
    library(dplyr)
    library(BiocParallel)
    library(ExperimentHub)
    library(magrittr)
})

BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

source(file.path(this.path::this.dir(), "array_utils.R"))

sesame_checkVersion()

################################################################################
# Helpers
################################################################################

parse_bool <- function(value) {
    tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

parse_cli_args <- function(args) {
    opts <- list(
        idat_dir = NULL,
        out_path = NULL,
        manifest_file = NULL,
        pval_threshold = 0.05,
        threads = 1L,
        array_type = "EPICv2",
        sesame_data = NULL,
        skip_annotation = TRUE,
        save_raw_beta = FALSE,
        save_noncollapsed_beta = FALSE,
        save_detection_pvals = FALSE,
        save_intensity = FALSE
    )

    if (length(args) == 0L) {
        stop("No arguments supplied. Use --help for usage.", call. = FALSE)
    }

    i <- 1L
    while (i <= length(args)) {
        arg <- args[[i]]

        if (arg == "--help") {
            cat(
                paste0(
                    "Usage: sesame.R \\\n",
                    "  --idat_dir <directory> \\\n",
                    "  --out_path <directory> \\\n",
                    "  --array_type <EPICv2|MSA|RhelixaCustom|...> \\\n",
                    "  [--manifest_file <csv>] \\\n",
                    "  [--pval_threshold <number>] \\\n",
                    "  [--threads <integer>] \\\n",
                    "  [--sesame_data <directory>] \\\n",
                    "  [--skip_annotation <true|false>] \\n",
                    "  [--save_raw_beta <true|false>] \\n",
                    "  [--save_noncollapsed_beta <true|false>] \\n",
                    "  [--save_detection_pvals <true|false>] \\n",
                    "  [--save_intensity <true|false>]\n"
                )
            )
            quit(status = 0L)
        }

        if (i == length(args)) {
            stop("Missing value after argument: ", arg, call. = FALSE)
        }
        value <- args[[i + 1L]]

        if (arg %in% c("--idat_dir", "-i")) {
            opts$idat_dir <- value
        } else if (arg %in% c("--out_path", "-o")) {
            opts$out_path <- value
        } else if (arg == "--manifest_file") {
            opts$manifest_file <- value
        } else if (arg == "--pval_threshold") {
            opts$pval_threshold <- as.numeric(value)
        } else if (arg == "--threads") {
            # The Docker workflow intentionally forces serial execution.
            opts$threads <- 1L
        } else if (arg == "--array_type") {
            opts$array_type <- value
        } else if (arg == "--sesame_data") {
            opts$sesame_data <- value
        } else if (arg == "--skip_annotation") {
            opts$skip_annotation <- parse_bool(value)
        } else if (arg == "--save_raw_beta") {
            opts$save_raw_beta <- parse_bool(value)
        } else if (arg == "--save_noncollapsed_beta") {
            opts$save_noncollapsed_beta <- parse_bool(value)
        } else if (arg == "--save_detection_pvals") {
            opts$save_detection_pvals <- parse_bool(value)
        } else if (arg == "--save_intensity") {
            opts$save_intensity <- parse_bool(value)
        } else {
            stop("Unknown argument: ", arg, call. = FALSE)
        }

        i <- i + 2L
    }

    if (is.null(opts$idat_dir) || !nzchar(opts$idat_dir)) {
        stop("--idat_dir is required.", call. = FALSE)
    }
    if (is.null(opts$out_path) || !nzchar(opts$out_path)) {
        stop("--out_path is required.", call. = FALSE)
    }
    if (!dir.exists(opts$idat_dir)) {
        stop("IDAT directory does not exist: ", opts$idat_dir, call. = FALSE)
    }
    if (!is.finite(opts$pval_threshold) || opts$pval_threshold < 0 || opts$pval_threshold > 1) {
        stop("--pval_threshold must be between 0 and 1.", call. = FALSE)
    }

    opts
}

matrix_from_named_vectors <- function(vectors, sample_names, label) {
    if (length(vectors) == 0L) {
        stop("No vectors were produced for ", label, ".", call. = FALSE)
    }
    if (length(vectors) != length(sample_names)) {
        stop("Vector/sample count mismatch for ", label, ".", call. = FALSE)
    }

    reference_ids <- names(vectors[[1L]])
    if (is.null(reference_ids)) {
        stop("The first vector for ", label, " has no probe names.", call. = FALSE)
    }

    aligned <- lapply(seq_along(vectors), function(i) {
        x <- vectors[[i]]
        ids <- names(x)
        if (is.null(ids)) {
            stop("A vector for ", label, " has no probe names.", call. = FALSE)
        }
        if (!identical(ids, reference_ids)) {
            if (!setequal(ids, reference_ids)) {
                stop(
                    "Probe sets differ between samples while constructing ",
                    label,
                    ".",
                    call. = FALSE
                )
            }
            x <- x[reference_ids]
        }
        as.numeric(x)
    })

    result <- do.call(cbind, aligned)
    if (is.null(dim(result))) {
        result <- matrix(result, ncol = 1L)
    }
    rownames(result) <- reference_ids
    colnames(result) <- sample_names
    result
}

write_matrix_tsv <- function(mat, path) {
    dt <- data.table::as.data.table(mat, keep.rownames = "ProbeID")
    data.table::fwrite(
        dt,
        path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
    invisible(dt)
}

write_gzip_lines <- function(values, path) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(sort(unique(as.character(values))), con = con)
}

accumulate_pre_qc_collapsed_counts <- function(sdf_list, pval_threshold) {
    counts <- integer(0)

    for (i in seq_along(sdf_list)) {
        # Apply the same QCDPB signal-processing sequence as the final matrix,
        # but ignore its mask when extracting the pre-QC beta vector. This
        # isolates value loss caused by SeSAMe's accumulated Q/P masks rather
        # than conflating it with channel, dye-bias or noob transformations.
        processed <- sesame::prepSesame(
            sdf_list[[i]],
            prep = "QCDPB",
            prep_args = list(P = list(pval.threshold = pval_threshold))
        )
        unmasked_collapsed <- sesame::getBetas(
            processed,
            mask = FALSE,
            collapseToPfx = TRUE,
            collapseMethod = "mean"
        )

        ids <- names(unmasked_collapsed)
        keep <- !is.na(ids) & grepl("^cg[0-9]+$", ids)
        ids <- ids[keep]
        values_present <- !is.na(unmasked_collapsed[keep])

        new_ids <- setdiff(ids, names(counts))
        if (length(new_ids) > 0L) {
            counts <- c(counts, stats::setNames(integer(length(new_ids)), new_ids))
        }
        counts[ids] <- counts[ids] + as.integer(values_present)
    }

    counts
}

write_qc_probe_reports <- function(
    sdf_list,
    collapsed_betas_matrix,
    sample_names,
    out_path,
    array_type,
    pval_threshold
) {
    logger::log_info("Summarizing CpGs affected by QCDPB QC/preprocessing...")

    pre_counts <- accumulate_pre_qc_collapsed_counts(
        sdf_list,
        pval_threshold = pval_threshold
    )

    post_ids <- rownames(collapsed_betas_matrix)
    post_keep <- !is.na(post_ids) & grepl("^cg[0-9]+$", post_ids)
    post_counts <- stats::setNames(
        rowSums(!is.na(collapsed_betas_matrix[post_keep, , drop = FALSE])),
        post_ids[post_keep]
    )

    all_ids <- union(names(pre_counts), names(post_counts))
    pre_aligned <- stats::setNames(integer(length(all_ids)), all_ids)
    post_aligned <- stats::setNames(integer(length(all_ids)), all_ids)
    pre_aligned[names(pre_counts)] <- pre_counts
    post_aligned[names(post_counts)] <- post_counts

    unexpected_gain <- post_aligned > pre_aligned
    if (any(unexpected_gain)) {
        logger::log_warn(
            paste0(
                sum(unexpected_gain),
                " CpGs had more non-missing post-QC than pre-QC values; ",
                "their loss count was truncated at zero."
            )
        )
    }

    failed_counts <- pmax(pre_aligned - post_aligned, 0L)
    evaluable <- pre_aligned > 0L
    affected <- evaluable & failed_counts > 0L
    completely_lost <- affected & post_aligned == 0L

    n_affected <- sum(affected)
    batch_name <- basename(normalizePath(out_path, mustWork = FALSE))
    qc_probe_failures <- data.table::data.table(
        dataset = rep(array_type, n_affected),
        batch = rep(batch_name, n_affected),
        ProbeID = all_ids[affected],
        total_samples = rep(length(sample_names), n_affected),
        pre_qc_nonmissing_samples = as.integer(pre_aligned[affected]),
        post_qc_nonmissing_samples = as.integer(post_aligned[affected]),
        qc_failed_samples = as.integer(failed_counts[affected])
    )

    if (nrow(qc_probe_failures) > 0L) {
        qc_probe_failures[
            , qc_failed_fraction_of_evaluable :=
                qc_failed_samples / pre_qc_nonmissing_samples
        ]
        qc_probe_failures[
            , qc_failed_fraction_of_all_samples :=
                qc_failed_samples / total_samples
        ]
        qc_probe_failures[
            , qc_failed_in_any_sample := qc_failed_samples > 0L
        ]
        qc_probe_failures[
            , qc_failed_in_all_evaluable_samples :=
                qc_failed_samples == pre_qc_nonmissing_samples
        ]
        qc_probe_failures[
            , qc_completely_lost_in_batch :=
                post_qc_nonmissing_samples == 0L
        ]
        data.table::setorder(qc_probe_failures, -qc_failed_samples, ProbeID)
    } else {
        qc_probe_failures[
            , `:=`(
                qc_failed_fraction_of_evaluable = numeric(),
                qc_failed_fraction_of_all_samples = numeric(),
                qc_failed_in_any_sample = logical(),
                qc_failed_in_all_evaluable_samples = logical(),
                qc_completely_lost_in_batch = logical()
            )
        ]
    }

    data.table::fwrite(
        qc_probe_failures,
        file.path(out_path, "qc_probe_failures.tsv.gz"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA",
        compress = "gzip"
    )

    pre_available_ids <- names(pre_aligned)[pre_aligned > 0L]
    post_available_ids <- names(post_aligned)[post_aligned > 0L]

    write_gzip_lines(
        pre_available_ids,
        file.path(out_path, "pre_qc_collapsed_probe_ids.txt.gz")
    )
    write_gzip_lines(
        post_available_ids,
        file.path(out_path, "post_qc_collapsed_probe_ids.txt.gz")
    )

    qc_batch_summary <- data.table::data.table(
        dataset = array_type,
        batch = basename(normalizePath(out_path, mustWork = FALSE)),
        n_samples = length(sample_names),
        pval_threshold = pval_threshold,
        n_pre_qc_cpgs_with_any_value = sum(pre_aligned > 0L),
        n_post_qc_cpgs_with_any_value = sum(post_aligned > 0L),
        n_qc_affected_cpgs = sum(affected),
        n_qc_completely_lost_cpgs = sum(completely_lost),
        n_qc_partially_affected_cpgs = sum(affected & !completely_lost),
        n_qc_failed_sample_cpg_pairs = sum(failed_counts)
    )

    data.table::fwrite(
        qc_batch_summary,
        file.path(out_path, "qc_batch_summary.csv"),
        quote = TRUE,
        row.names = FALSE
    )

    logger::log_info(
        paste0(
            "QC report: ",
            sum(affected),
            " CpGs affected in at least one sample; ",
            sum(completely_lost),
            " CpGs completely lost in this batch."
        )
    )
}

################################################################################
# Main
################################################################################

main <- function() {
    opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))

    idat_dir <- opts$idat_dir
    out_path <- opts$out_path
    manifest_file <- opts$manifest_file
    pval_threshold <- opts$pval_threshold
    threads <- 1L
    array_type <- opts$array_type
    sesame_data <- opts$sesame_data
    skip_annotation <- opts$skip_annotation
    save_raw_beta <- opts$save_raw_beta
    save_noncollapsed_beta <- opts$save_noncollapsed_beta
    save_detection_pvals <- opts$save_detection_pvals
    save_intensity <- opts$save_intensity

    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    logger::log_info(paste0("IDAT directory: ", idat_dir))
    logger::log_info(paste0("Output directory: ", out_path))
    logger::log_info(paste0("Manifest file: ", manifest_file))
    logger::log_info(paste0("Array type: ", array_type))
    logger::log_info(paste0("p-value threshold: ", pval_threshold))
    logger::log_info(paste0("Number of threads: ", threads))
    logger::log_info(paste0("sesame_data: ", sesame_data))
    logger::log_info(paste0("skip_annotation: ", skip_annotation))
    logger::log_info(paste0("save_raw_beta: ", save_raw_beta))
    logger::log_info(paste0("save_noncollapsed_beta: ", save_noncollapsed_beta))
    logger::log_info(paste0("save_detection_pvals: ", save_detection_pvals))
    logger::log_info(paste0("save_intensity: ", save_intensity))

    bp <- BiocParallel::SerialParam()

    if (is.null(sesame_data)) {
        sesame_data <- array_utils.sesame_data_path()
    }
    if (!dir.exists(sesame_data)) {
        stop("SeSAMe cache directory does not exist: ", sesame_data, call. = FALSE)
    }
    ExperimentHub::setExperimentHubOption("CACHE", sesame_data)

    if (array_type == "RhelixaCustom") {
        if (is.null(manifest_file)) {
            manifest_file <- array_utils.manifest_path(array_type)
        }
        manifest <- array_utils.read_manifest(manifest_file, threads = threads)
    } else {
        manifest <- NULL
    }

    logger::log_info("Reading IDAT files...")
    idat_prefixes <- sesame::searchIDATprefixes(idat_dir)
    if (length(idat_prefixes) == 0L) {
        stop("No paired IDAT prefixes were found in: ", idat_dir, call. = FALSE)
    }
    sample_names <- basename(idat_prefixes)
    if (anyDuplicated(sample_names)) {
        stop("Duplicate sample prefixes were detected in the IDAT directory.", call. = FALSE)
    }

    sdf_list <- BiocParallel::bplapply(
        idat_prefixes,
        function(pfx) sesame::readIDATpair(pfx, manifest = manifest),
        BPPARAM = bp
    )
    names(sdf_list) <- sample_names

    ########################################
    # Per-sample QC summary
    ########################################
    logger::log_info("Calculating QC metrics...")
    qcs <- do.call(
        rbind,
        BiocParallel::bplapply(
            sdf_list,
            function(sdf) {
                as.data.frame(
                    sesame::sesameQC_calcStats(sdf, c("numProbes", "detection")),
                    stringsAsFactors = FALSE
                )
            },
            BPPARAM = bp
        )
    )
    rownames(qcs) <- sample_names
    qcs$frac_dt <- qcs$num_dt / qcs$num_probes
    qcs <- qcs[
        c(
            "num_probes", "num_dt", "frac_dt",
            "num_probes_cg", "num_dt_cg", "frac_dt_cg",
            "num_probes_ch", "num_dt_ch", "frac_dt_ch"
        )
    ]
    colnames(qcs) <- c(
        "N. Probes", "Detected Probes", "Detection Rate",
        "N. Probes(CG)", "Detected Probes(CG)", "Detection Rate(CG)",
        "N. Probes (CH)", "Detected Probes(CH)", "Detection Rate(CH)"
    )
    qcs <- data.table::as.data.table(qcs, keep.rownames = "SampleName")
    data.table::fwrite(
        qcs,
        file.path(out_path, "qc_summary.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
    rm(qcs)
    gc(reset = TRUE)

    ########################################
    # Optional raw beta matrix
    ########################################
    if (save_raw_beta) {
        logger::log_info("Creating the pre-filtering raw beta matrix...")
        raw_beta_vectors <- BiocParallel::bplapply(
            sdf_list,
            function(sdf) sesame::getBetas(sdf),
            BPPARAM = bp
        )
        raw_betas_matrix <- matrix_from_named_vectors(
            raw_beta_vectors,
            sample_names,
            "raw beta matrix"
        )
        write_matrix_tsv(raw_betas_matrix, file.path(out_path, "raw_beta_matrix.txt"))
        rm(raw_beta_vectors, raw_betas_matrix)
        gc(reset = TRUE)
    } else {
        logger::log_info("Skipping raw beta matrix creation because save_raw_beta=FALSE.")
    }

    ########################################
    # Optional non-collapsed post-QCDPB beta matrix
    ########################################
    betas <- NULL
    manifest_ann <- NULL

    if (save_noncollapsed_beta || !skip_annotation) {
        logger::log_info("Creating the post-filtering non-collapsed beta matrix...")
        beta_vectors <- BiocParallel::bplapply(
            sdf_list,
            function(sdf) {
                processed <- sesame::prepSesame(
                    sdf,
                    prep = "QCDPB",
                    prep_args = list(P = list(pval.threshold = pval_threshold))
                )
                sesame::getBetas(processed, mask = TRUE)
            },
            BPPARAM = bp
        )
        betas_matrix <- matrix_from_named_vectors(
            beta_vectors,
            sample_names,
            "post-filtering beta matrix"
        )

        if (save_noncollapsed_beta) {
            betas <- write_matrix_tsv(
                betas_matrix,
                file.path(out_path, "beta_matrix.txt")
            )
        } else {
            betas <- data.table::as.data.table(betas_matrix, keep.rownames = "ProbeID")
        }

        rm(beta_vectors, betas_matrix)
        gc(reset = TRUE)
    } else {
        logger::log_info(
            "Skipping non-collapsed post-filtering beta matrix creation because save_noncollapsed_beta=FALSE and skip_annotation=TRUE."
        )
    }

    ########################################
    # Optional annotation for non-collapsed beta matrix
    ########################################
    if (!skip_annotation) {
        logger::log_info("Creating the annotated beta matrix...")

        if (array_type == "RhelixaCustom") {
            manifest_ann <- manifest[c("Probe_ID", "CHR", "MAPINFO", "Strand_FR")] %>%
                tidyr::drop_na(CHR) %>%
                dplyr::mutate(
                    CHR = paste0("chr", CHR),
                    Strand_FR = ifelse(Strand_FR == "F", "+", "-"),
                    End = MAPINFO + 1
                ) %>%
                dplyr::rename(
                    ProbeID = Probe_ID,
                    Chromosome = CHR,
                    Start = MAPINFO,
                    Strand = Strand_FR
                ) %>%
                dplyr::filter(Chromosome != "chr0") %>%
                dplyr::mutate(
                    ProbeName = ifelse(
                        grepl("_", ProbeID),
                        gsub("_(?:[^_]*)$", "", ProbeID),
                        ProbeID
                    )
                ) %>%
                dplyr::select(ProbeID, ProbeName, Chromosome, Start, End, Strand)
        } else {
            manifest_ann <- sesameData::sesameData_getManifestGRanges(array_type) %>%
                as.data.frame(stringsAsFactors = FALSE) %>%
                dplyr::mutate(
                    ProbeID = rownames(.),
                    seqnames = as.character(seqnames),
                    strand = as.character(strand)
                ) %>%
                dplyr::mutate(
                    ProbeName = ifelse(
                        grepl("_", ProbeID),
                        gsub("_(?:[^_]*)$", "", ProbeID),
                        ProbeID
                    )
                ) %>%
                dplyr::rename(
                    Chromosome = seqnames,
                    Start = start,
                    End = end,
                    Strand = strand
                ) %>%
                dplyr::filter(Chromosome != "*") %>%
                dplyr::select(ProbeID, ProbeName, Chromosome, Start, End, Strand)
        }

        manifest_ann <- tryCatch(
            array_utils.add_ewas_atlas_annotation(
                manifest_ann,
                probeid_col = "ProbeName"
            ),
            error = function(e) {
                logger::log_warn(
                    paste0("EWAS Atlas annotation was skipped: ", conditionMessage(e))
                )
                manifest_ann
            }
        )
        manifest_ann <- tryCatch(
            array_utils.add_closest_gene_annotation(
                manifest_ann,
                chr = "Chromosome",
                start = "Start",
                end = "End",
                threads = threads
            ),
            error = function(e) {
                logger::log_warn(
                    paste0("Closest-gene annotation was skipped: ", conditionMessage(e))
                )
                manifest_ann
            }
        )
        manifest_ann <- tryCatch(
            array_utils.add_gene_detail_annotation(
                manifest_ann,
                geneid_col = "GeneID"
            ),
            error = function(e) {
                logger::log_warn(
                    paste0("Gene-detail annotation was skipped: ", conditionMessage(e))
                )
                manifest_ann
            }
        )
        manifest_ann <- manifest_ann %>%
            dplyr::select(-End) %>%
            dplyr::rename(Position = Start)

        betas$row_num <- seq_len(nrow(betas))
        betas_ann <- dplyr::right_join(manifest_ann, betas, by = "ProbeID") %>%
            dplyr::arrange(row_num) %>%
            dplyr::select(-row_num) %>%
            dplyr::mutate(
                dplyr::across(
                    dplyr::all_of(colnames(manifest_ann)),
                    ~ ifelse(is.na(.), "", .)
                )
            ) %>%
            dplyr::mutate(
                ProbeName = ifelse(
                    grepl("_", ProbeID),
                    gsub("_(?:[^_]*)$", "", ProbeID),
                    ProbeID
                )
            )

        data.table::fwrite(
            betas_ann,
            file.path(out_path, "beta_matrix_ann.txt"),
            sep = "	",
            quote = FALSE,
            row.names = FALSE,
            na = "NA"
        )
        rm(betas, betas_ann)
        gc(reset = TRUE)
    } else {
        logger::log_info("Skipping annotated beta matrix creation because skip_annotation=TRUE.")
        if (!is.null(betas)) {
            rm(betas)
            gc(reset = TRUE)
        }
    }

    ########################################
    # Suffix-collapsed post-QCDPB beta matrix and QC probe-loss report
    ########################################
    logger::log_info("Creating the suffix-collapsed beta matrix...")
    collapsed_beta_vectors <- BiocParallel::bplapply(
        sdf_list,
        function(sdf) {
            processed <- sesame::prepSesame(
                sdf,
                prep = "QCDPB",
                prep_args = list(P = list(pval.threshold = pval_threshold))
            )
            sesame::getBetas(
                processed,
                mask = TRUE,
                collapseToPfx = TRUE,
                collapseMethod = "mean"
            )
        },
        BPPARAM = bp
    )
    collapsed_betas_matrix <- matrix_from_named_vectors(
        collapsed_beta_vectors,
        sample_names,
        "suffix-collapsed beta matrix"
    )
    collapsed_betas_matrix <- collapsed_betas_matrix[
        !grepl("^ctl|^cgBK", rownames(collapsed_betas_matrix)),
        ,
        drop = FALSE
    ]

    write_qc_probe_reports(
        sdf_list = sdf_list,
        collapsed_betas_matrix = collapsed_betas_matrix,
        sample_names = sample_names,
        out_path = out_path,
        array_type = array_type,
        pval_threshold = pval_threshold
    )

    collapsed_betas <- write_matrix_tsv(
        collapsed_betas_matrix,
        file.path(out_path, "collapsed_beta_matrix.txt")
    )
    rm(collapsed_beta_vectors, collapsed_betas_matrix)
    gc(reset = TRUE)

    ########################################
    # Optional annotation for suffix-collapsed beta matrix
    ########################################
    if (!skip_annotation) {
        logger::log_info("Creating the annotated suffix-collapsed beta matrix...")
        collapsed_manifest_ann <- manifest_ann %>%
            dplyr::distinct(ProbeName, .keep_all = TRUE)
        collapsed_betas$row_num <- seq_len(nrow(collapsed_betas))
        collapsed_betas_ann <- dplyr::right_join(
            collapsed_manifest_ann,
            collapsed_betas,
            by = c("ProbeName" = "ProbeID")
        ) %>%
            dplyr::arrange(row_num) %>%
            dplyr::select(-row_num) %>%
            dplyr::mutate(
                dplyr::across(
                    dplyr::all_of(colnames(collapsed_manifest_ann)),
                    ~ ifelse(is.na(.), "", .)
                )
            ) %>%
            dplyr::select(-ProbeID) %>%
            dplyr::rename(ProbeID = ProbeName)

        data.table::fwrite(
            collapsed_betas_ann,
            file.path(out_path, "collapsed_beta_matrix_ann.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            na = "NA"
        )
        rm(collapsed_betas, collapsed_betas_ann)
        gc(reset = TRUE)
    } else {
        logger::log_info(
            "Skipping annotated suffix-collapsed beta matrix creation because skip_annotation=TRUE."
        )
        rm(collapsed_betas)
        gc(reset = TRUE)
    }

    ########################################
    # Optional detection p-value matrix corresponding to the Q-C-D-P sequence
    ########################################
    if (save_detection_pvals) {
        logger::log_info("Creating the p-value matrix...")
        pval_vectors <- BiocParallel::bplapply(
            sdf_list,
            function(sdf) {
                qcd <- sesame::prepSesame(sdf, prep = "QCD")
                sesame::pOOBAH(qcd, return.pval = TRUE)
            },
            BPPARAM = bp
        )
        pvals_matrix <- matrix_from_named_vectors(
            pval_vectors,
            sample_names,
            "detection p-value matrix"
        )
        pvals_matrix <- pvals_matrix[
            !grepl("^ctl|^cgBK", rownames(pvals_matrix)),
            ,
            drop = FALSE
        ]
        write_matrix_tsv(pvals_matrix, file.path(out_path, "detection_pvals.txt"))
        rm(pval_vectors, pvals_matrix)
        gc(reset = TRUE)
    } else {
        logger::log_info("Skipping detection p-value matrix creation because save_detection_pvals=FALSE.")
    }

    ########################################
    # Optional methylated/unmethylated-intensity matrices
    ########################################
    if (save_intensity) {
        logger::log_info("Creating the methylated-intensity matrix...")
        meth_vectors <- BiocParallel::bplapply(
            sdf_list,
            function(sdf) {
                sig_mu <- sesame::signalMU(sdf)
                stats::setNames(sig_mu$M, sig_mu$Probe_ID)
            },
            BPPARAM = bp
        )
        meth_matrix <- matrix_from_named_vectors(
            meth_vectors,
            sample_names,
            "methylated-intensity matrix"
        )
        meth_matrix <- meth_matrix[
            !grepl("^ctl|^cgBK", rownames(meth_matrix)),
            ,
            drop = FALSE
        ]
        write_matrix_tsv(meth_matrix, file.path(out_path, "methylated_intensity.txt"))
        rm(meth_vectors, meth_matrix)
        gc(reset = TRUE)

        logger::log_info("Creating the unmethylated-intensity matrix...")
        unmeth_vectors <- BiocParallel::bplapply(
            sdf_list,
            function(sdf) {
                sig_mu <- sesame::signalMU(sdf)
                stats::setNames(sig_mu$U, sig_mu$Probe_ID)
            },
            BPPARAM = bp
        )
        unmeth_matrix <- matrix_from_named_vectors(
            unmeth_vectors,
            sample_names,
            "unmethylated-intensity matrix"
        )
        unmeth_matrix <- unmeth_matrix[
            !grepl("^ctl|^cgBK", rownames(unmeth_matrix)),
            ,
            drop = FALSE
        ]
        write_matrix_tsv(
            unmeth_matrix,
            file.path(out_path, "unmethylated_intensity.txt")
        )
        rm(unmeth_vectors, unmeth_matrix)
        gc(reset = TRUE)
    } else {
        logger::log_info("Skipping methylated/unmethylated-intensity matrix creation because save_intensity=FALSE.")
    }

    rm(sdf_list)
    gc(reset = TRUE)

    logger::log_info("Finish!")
}

if (sys.nframe() == 0L) {
    main()
}
