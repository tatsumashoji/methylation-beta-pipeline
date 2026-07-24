#!/usr/bin/env Rscript

# Dockerized SeSAMe IDAT-to-beta pipeline with per-sample collapsed-CpG QC loss logging.
# Main outputs per batch:
#   collapsed_beta_matrix.txt
#   qc_probe_failures_by_sample.tsv.gz
#   qc_probe_failures.tsv.gz
#   pre_qc_collapsed_probe_ids.txt.gz
#   post_qc_collapsed_probe_ids.txt.gz
#   qc_batch_summary.csv
# Optional large outputs are controlled by SAVE_* environment variables or CLI flags.

suppressPackageStartupMessages(library(logger))
suppressPackageStartupMessages(library(sesame))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(this.path))

source(paste0(this.path::this.dir(), "/array_utils.R"))

suppressWarnings(try(sesame_checkVersion(), silent = TRUE))
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

str_to_bool <- function(x, default = FALSE) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return(default)
    tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

collapse_probe_id <- function(x) {
    ifelse(grepl("_", x), sub("_(?:[^_]*)$", "", x), x)
}

write_lines_gz <- function(x, path) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(sort(unique(as.character(x))), con)
}

safe_colnames <- function(mat, sample_ids) {
    if (is.null(dim(mat))) {
        mat <- matrix(mat, ncol = 1)
    }
    if (is.null(colnames(mat)) || length(colnames(mat)) != length(sample_ids)) {
        colnames(mat) <- sample_ids
    }
    mat
}

write_matrix_tsv <- function(mat, path) {
    dt <- data.table::as.data.table(mat, keep.rownames = "ProbeID")
    data.table::fwrite(dt, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    invisible(NULL)
}

parse_args <- function(args) {
    cfg <- list(
        threads = 1,
        pval_threshold = 0.05,
        array_type = "EPICv2",
        manifest_file = NULL,
        sesame_data = NULL,
        skip_annotation = TRUE,
        save_raw_beta = str_to_bool(Sys.getenv("SAVE_RAW_BETA", "0")),
        save_noncollapsed_beta = str_to_bool(Sys.getenv("SAVE_NONCOLLAPSED_BETA", "0")),
        save_detection_pvals = str_to_bool(Sys.getenv("SAVE_DETECTION_PVALS", "0")),
        save_intensity = str_to_bool(Sys.getenv("SAVE_INTENSITY", "0"))
    )

    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        if (arg == "--help") {
            cat("Usage: sesame.R --idat_dir <dir> --out_path <dir> --array_type <EPICv2|MSA|RhelixaCustom> --manifest_file <file> --pval_threshold <float> --threads <int> --sesame_data <dir> --skip_annotation true|false [--save_raw_beta true|false] [--save_noncollapsed_beta true|false] [--save_detection_pvals true|false] [--save_intensity true|false]\n")
            q(status = 0)
        } else if (arg %in% c("--idat_dir", "-i")) {
            cfg$idat_dir <- args[i + 1]; i <- i + 1
        } else if (arg %in% c("--out_path", "-o")) {
            cfg$out_path <- args[i + 1]; i <- i + 1
        } else if (arg == "--manifest_file") {
            cfg$manifest_file <- args[i + 1]; i <- i + 1
        } else if (arg == "--pval_threshold") {
            cfg$pval_threshold <- as.numeric(args[i + 1]); i <- i + 1
        } else if (arg == "--threads") {
            cfg$threads <- 1; i <- i + 1
        } else if (arg == "--array_type") {
            cfg$array_type <- args[i + 1]; i <- i + 1
        } else if (arg == "--sesame_data") {
            cfg$sesame_data <- args[i + 1]; i <- i + 1
        } else if (arg == "--skip_annotation") {
            cfg$skip_annotation <- str_to_bool(args[i + 1], default = TRUE); i <- i + 1
        } else if (arg == "--save_raw_beta") {
            cfg$save_raw_beta <- str_to_bool(args[i + 1]); i <- i + 1
        } else if (arg == "--save_noncollapsed_beta") {
            cfg$save_noncollapsed_beta <- str_to_bool(args[i + 1]); i <- i + 1
        } else if (arg == "--save_detection_pvals") {
            cfg$save_detection_pvals <- str_to_bool(args[i + 1]); i <- i + 1
        } else if (arg == "--save_intensity") {
            cfg$save_intensity <- str_to_bool(args[i + 1]); i <- i + 1
        } else {
            stop(paste0("Unknown argument: ", arg))
        }
        i <- i + 1
    }

    if (is.null(cfg$idat_dir)) stop("--idat_dir is required")
    if (is.null(cfg$out_path)) stop("--out_path is required")
    cfg
}

main <- function() {
    cfg <- parse_args(commandArgs(trailingOnly = TRUE))

    dir.create(cfg$out_path, showWarnings = FALSE, recursive = TRUE)

    logger::log_info(paste0("IDAT directory: ", cfg$idat_dir))
    logger::log_info(paste0("Output directory: ", cfg$out_path))
    logger::log_info(paste0("Manifest file: ", cfg$manifest_file))
    logger::log_info(paste0("Array type: ", cfg$array_type))
    logger::log_info(paste0("p-value threshold: ", cfg$pval_threshold))
    logger::log_info("Number of threads: 1")
    logger::log_info(paste0("sesame_data: ", cfg$sesame_data))
    logger::log_info(paste0("skip_annotation: ", cfg$skip_annotation))
    logger::log_info(paste0("save_raw_beta: ", cfg$save_raw_beta))
    logger::log_info(paste0("save_noncollapsed_beta: ", cfg$save_noncollapsed_beta))
    logger::log_info(paste0("save_detection_pvals: ", cfg$save_detection_pvals))
    logger::log_info(paste0("save_intensity: ", cfg$save_intensity))

    bp <- BiocParallel::SerialParam()

    if (is.null(cfg$sesame_data)) cfg$sesame_data <- array_utils.sesame_data_path()
    ExperimentHub::setExperimentHubOption("CACHE", cfg$sesame_data)

    if (cfg$array_type == "RhelixaCustom") {
        if (is.null(cfg$manifest_file)) cfg$manifest_file <- array_utils.manifest_path(cfg$array_type)
        manifest <- array_utils.read_manifest(cfg$manifest_file, threads = 1)
    } else {
        manifest <- NULL
    }

    logger::log_info("Reading IDAT files...")
    idat_prefixes <- sesame::searchIDATprefixes(cfg$idat_dir)
    if (length(idat_prefixes) == 0) stop("No IDAT prefixes found")

    sdf_list <- BiocParallel::bplapply(
        idat_prefixes,
        function(pfx) sesame::readIDATpair(pfx, manifest = manifest),
        BPPARAM = bp
    )
    sample_ids <- basename(idat_prefixes)
    names(sdf_list) <- sample_ids

    logger::log_info("Calculating QC metrics...")
    qcs <- do.call(
        rbind,
        BiocParallel::bplapply(
            sdf_list,
            function(sdf) {
                as.data.frame(sesame::sesameQC_calcStats(sdf, c("numProbes", "detection")), stringsAsFactors = FALSE)
            },
            BPPARAM = bp
        )
    )
    qcs$frac_dt <- qcs$num_dt / qcs$num_probes
    qcs <- qcs[c("num_probes", "num_dt", "frac_dt", "num_probes_cg", "num_dt_cg", "frac_dt_cg", "num_probes_ch", "num_dt_ch", "frac_dt_ch")]
    colnames(qcs) <- c("N. Probes", "Detected Probes", "Detection Rate", "N. Probes(CG)", "Detected Probes(CG)", "Detection Rate(CG)", "N. Probes (CH)", "Detected Probes(CH)", "Detection Rate(CH)")
    qcs <- data.table::as.data.table(qcs, keep.rownames = "SampleName")
    qcs[, SampleID := sample_ids]
    data.table::fwrite(qcs, file.path(cfg$out_path, "qc_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    rm(qcs); gc(reset = TRUE)

    if (cfg$save_raw_beta) {
        logger::log_info("Creating optional raw beta matrix...")
        raw_betas <- do.call(cbind, BiocParallel::bplapply(sdf_list, sesame::getBetas, BPPARAM = bp))
        raw_betas <- safe_colnames(raw_betas, sample_ids)
        write_matrix_tsv(raw_betas, file.path(cfg$out_path, "raw_beta_matrix.txt"))
        rm(raw_betas); gc(reset = TRUE)
    }

    if (cfg$save_noncollapsed_beta) {
        logger::log_info("Creating optional non-collapsed post-QC beta matrix...")
        betas <- do.call(
            cbind,
            BiocParallel::bplapply(
                sdf_list,
                function(sdf) {
                    sesame::getBetas(
                        sesame::prepSesame(sdf, prep = "QCDPB", prep_args = list(P = list(pval.threshold = cfg$pval_threshold)))
                    )
                },
                BPPARAM = bp
            )
        )
        betas <- safe_colnames(betas, sample_ids)
        betas <- betas[!grepl("^ctl|^cgBK", rownames(betas)), , drop = FALSE]
        write_matrix_tsv(betas, file.path(cfg$out_path, "beta_matrix.txt"))
        rm(betas); gc(reset = TRUE)
    }

    logger::log_info("Creating post-QC collapsed beta matrix and per-sample QC-loss records...")
    pre_qc_list <- list()
    post_qc_list <- list()
    qc_records <- list()

    for (sample_id in sample_ids) {
        sdf <- sdf_list[[sample_id]]
        prepped <- sesame::prepSesame(
            sdf,
            prep = "QCDPB",
            prep_args = list(P = list(pval.threshold = cfg$pval_threshold))
        )

        pre <- sesame::getBetas(prepped, mask = FALSE, collapseToPfx = TRUE)
        post <- sesame::getBetas(prepped, mask = TRUE, collapseToPfx = TRUE)

        pre <- pre[!grepl("^ctl|^cgBK", names(pre))]
        post <- post[!grepl("^ctl|^cgBK", names(post))]

        pre_qc_list[[sample_id]] <- pre
        post_qc_list[[sample_id]] <- post

        all_ids <- union(names(pre), names(post))
        pre_aligned <- pre[all_ids]
        post_aligned <- post[all_ids]
        failed <- !is.na(pre_aligned) & is.na(post_aligned)
        failed_ids <- all_ids[failed]

        if (length(failed_ids) > 0) {
            qc_records[[length(qc_records) + 1L]] <- data.table::data.table(
                SampleID = sample_id,
                ProbeID = failed_ids,
                pre_qc_beta = as.numeric(pre_aligned[failed_ids]),
                post_qc_beta = as.numeric(post_aligned[failed_ids])
            )
        }

        rm(prepped, pre, post, pre_aligned, post_aligned); gc(reset = TRUE)
    }

    all_post_ids <- sort(unique(unlist(lapply(post_qc_list, names), use.names = FALSE)))
    collapsed_betas <- matrix(NA_real_, nrow = length(all_post_ids), ncol = length(sample_ids), dimnames = list(all_post_ids, sample_ids))
    for (sample_id in sample_ids) {
        v <- post_qc_list[[sample_id]]
        collapsed_betas[names(v), sample_id] <- as.numeric(v)
    }
    write_matrix_tsv(collapsed_betas, file.path(cfg$out_path, "collapsed_beta_matrix.txt"))

    pre_ids <- sort(unique(unlist(lapply(pre_qc_list, names), use.names = FALSE)))
    post_ids <- sort(unique(unlist(lapply(post_qc_list, names), use.names = FALSE)))
    write_lines_gz(pre_ids, file.path(cfg$out_path, "pre_qc_collapsed_probe_ids.txt.gz"))
    write_lines_gz(post_ids, file.path(cfg$out_path, "post_qc_collapsed_probe_ids.txt.gz"))

    if (length(qc_records) > 0) {
        qc_by_sample <- data.table::rbindlist(qc_records, use.names = TRUE, fill = TRUE)
    } else {
        qc_by_sample <- data.table::data.table(SampleID = character(), ProbeID = character(), pre_qc_beta = numeric(), post_qc_beta = numeric())
    }
    data.table::fwrite(qc_by_sample, file.path(cfg$out_path, "qc_probe_failures_by_sample.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", compress = "gzip")

    if (nrow(qc_by_sample) > 0) {
        qc_summary <- qc_by_sample[, .(qc_failed_samples = data.table::uniqueN(SampleID)), by = ProbeID]
        qc_summary[, n_samples := length(sample_ids)]
        qc_summary[, qc_failed_fraction := qc_failed_samples / n_samples]
        qc_summary[, present_before_qc := ProbeID %in% pre_ids]
        qc_summary[, present_after_qc := ProbeID %in% post_ids]
        qc_summary[, qc_completely_lost := !(ProbeID %in% post_ids)]
        qc_summary[, qc_partially_affected := (ProbeID %in% post_ids)]
    } else {
        qc_summary <- data.table::data.table(
            ProbeID = character(), qc_failed_samples = integer(), n_samples = integer(),
            qc_failed_fraction = numeric(), present_before_qc = logical(),
            present_after_qc = logical(), qc_completely_lost = logical(), qc_partially_affected = logical()
        )
    }
    data.table::fwrite(qc_summary, file.path(cfg$out_path, "qc_probe_failures.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", compress = "gzip")

    batch_summary <- data.table::data.table(
        dataset = cfg$array_type,
        batch = basename(normalizePath(cfg$out_path, mustWork = FALSE)),
        n_samples = length(sample_ids),
        pval_threshold = cfg$pval_threshold,
        n_pre_qc_collapsed_cpgs = length(pre_ids),
        n_post_qc_collapsed_cpgs = length(post_ids),
        n_qc_affected_cpgs = data.table::uniqueN(qc_summary$ProbeID),
        n_qc_completely_lost_cpgs = sum(qc_summary$qc_completely_lost),
        n_qc_partially_affected_cpgs = sum(qc_summary$qc_partially_affected),
        n_qc_failed_sample_cpg_pairs = nrow(qc_by_sample)
    )
    data.table::fwrite(batch_summary, file.path(cfg$out_path, "qc_batch_summary.csv"), quote = FALSE, row.names = FALSE)

    if (cfg$save_detection_pvals) {
        logger::log_info("Creating optional detection p-value matrix...")
        pvals <- do.call(
            cbind,
            BiocParallel::bplapply(
                sdf_list,
                function(sdf) sesame::pOOBAH(sesame::dyeBiasNL(sesame::inferInfiniumIChannel(sdf)), return.pval = TRUE),
                BPPARAM = bp
            )
        )
        pvals <- safe_colnames(pvals, sample_ids)
        pvals <- pvals[!grepl("^ctl|^cgBK", rownames(pvals)), , drop = FALSE]
        write_matrix_tsv(pvals, file.path(cfg$out_path, "detection_pvals.txt"))
        rm(pvals); gc(reset = TRUE)
    }

    if (cfg$save_intensity) {
        logger::log_info("Creating optional methylated/unmethylated intensity matrices...")
        meth_intensity <- do.call(cbind, BiocParallel::bplapply(sdf_list, function(sdf) {
            sig <- sesame::signalMU(sdf); setNames(sig$M, sig$Probe_ID)
        }, BPPARAM = bp))
        meth_intensity <- safe_colnames(meth_intensity, sample_ids)
        meth_intensity <- meth_intensity[!grepl("^ctl|^cgBK", rownames(meth_intensity)), , drop = FALSE]
        write_matrix_tsv(meth_intensity, file.path(cfg$out_path, "methylated_intensity.txt"))
        rm(meth_intensity); gc(reset = TRUE)

        unmeth_intensity <- do.call(cbind, BiocParallel::bplapply(sdf_list, function(sdf) {
            sig <- sesame::signalMU(sdf); setNames(sig$U, sig$Probe_ID)
        }, BPPARAM = bp))
        unmeth_intensity <- safe_colnames(unmeth_intensity, sample_ids)
        unmeth_intensity <- unmeth_intensity[!grepl("^ctl|^cgBK", rownames(unmeth_intensity)), , drop = FALSE]
        write_matrix_tsv(unmeth_intensity, file.path(cfg$out_path, "unmethylated_intensity.txt"))
        rm(unmeth_intensity); gc(reset = TRUE)
    }

    logger::log_info("Finish!")
}

main()
