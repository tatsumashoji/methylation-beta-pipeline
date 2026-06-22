#!/usr/bin/env Rscript

# Description: Script for generating beta-value matrices from IDAT files using SeSAMe
# Singularity: sesame_25.03.16.sif
# Usage: sesame.R \
# --idat_dir <idat_dir> \cd /data/share/project/analysis/analysis/KK-PJ20240701/rawdata_all
# --out_path <out_path> \
# --manifest_file <manifest_file> \
# --array_type <EPICv2|RhelixaCustom|MSA|...> \
# --pval_threshold <pval_threshold> \
# --threads <threads>

########################################################################################

library(logger) 
suppressMessages(library(sesame))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(BiocParallel))
suppressMessages(library(ExperimentHub))
suppressMessages(library(magrittr))

source(paste0(this.path::this.dir(),"/array_utils.R"))

sesame_checkVersion()

########################################################################################

main <- function(){

    args <- commandArgs(trailingOnly=T)
    # default values
    threads <- 8
    pval_threshold <- 0.05
    array_type <- "EPICv2"
    manifest_file <- NULL
    sesame_data <- NULL
    skip_annotation <- TRUE
    for (i in seq_along(args)) {
        arg <- args[i]
        if (arg == "--help") {
            cat("Usage: sesame.R --idat_dir <idat directory> --beta_out_path <beta_out_path> --pval_out_path <pval_out_path> --array_type <EPICv2|MSA|RhelixaCustom|...> --manifest_file <manifest_file> --control_file <control_file> --pval_threshold <pval_threshold> --threads <threads>\n")
            q()
        } else if (arg %in% c("--idat_dir", "-i")) {
            idat_dir <- args[i+1]
        } else if (arg %in% c("--out_path", "-o")) {
            out_path <- args[i+1]
        } else if (arg == "--manifest_file") {
            manifest_file <- args[i+1]
        } else if (arg == "--pval_threshold") {
            pval_threshold <- as.numeric(args[i+1])
        } else if (arg == "--threads") {
            threads <- as.integer(args[i+1])
        } else if (arg == "--array_type") {
            array_type <- args[i+1]
        } else if (arg == "--sesame_data") {
            sesame_data <- args[i+1]
        } else if (arg == "--skip_annotation") {
            skip_annotation <- tolower(args[i+1]) %in% c("true", "t", "1", "yes", "y")
        } 
    }
    logger::log_info(paste0("IDAT directory: ", idat_dir))
    logger::log_info(paste0("Output directory: ", out_path))
    logger::log_info(paste0("Manifest file: ", manifest_file))
    logger::log_info(paste0("Array type: ", array_type))
    logger::log_info(paste0("p-value threshold: ", pval_threshold))
    logger::log_info(paste0("Number of threads: ", threads))
    logger::log_info(paste0("sesame_data: ", sesame_data))
    logger::log_info(paste0("skip_annotation: ", skip_annotation))

    if (threads <= 1) {
        bp <- BiocParallel::SerialParam()
    } else {
        bp <- BiocParallel::SnowParam(
            workers = threads,
            type = "SOCK",
            progressbar = FALSE,
            stop.on.error = TRUE
        )
    }

    # Specify the directory where sesameDataCache is stored
    if (is.null(sesame_data)) {
        sesame_data <- array_utils.sesame_data_path()
    }
    ExperimentHub::setExperimentHubOption("CACHE", sesame_data)

    if (array_type == "RhelixaCustom") {
        if (is.null(manifest_file)) {
            manifest_file <- array_utils.manifest_path(array_type)
        }
        manifest <- array_utils.read_manifest(manifest_file, threads=threads)
    } else {
        manifest <- NULL # NOTE: For non-custom arrays, manifest can be NULL
    }

    ########################################
    # Read IDAT files
    ########################################
    logger::log_info("Reading IDAT files...")
    idat_prefixes = searchIDATprefixes(idat_dir)
    sdf_list <- BiocParallel::bplapply(
        idat_prefixes, function(pfx) {
            sdf <- sesame::readIDATpair(pfx, manifest=manifest)
            return(sdf)
        },
        BPPARAM = bp
    )

    ########################################
    # Calculate QC metrics
    ########################################
    logger::log_info("Calculating QC metrics...")
    qcs <- do.call(
        rbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
                as.data.frame(sesame::sesameQC_calcStats(sdf, c("numProbes", "detection")), stringsAsFactors=FALSE)
            },
            BPPARAM = bp
        )
    )
    # NOTE: sesameQC_calcStats calculates mfrac_dt as num_dt/(num_probes - num_dtna); recalculate it as num_dt/num_probes.
    qcs$frac_dt <- qcs$num_dt / qcs$num_probes
    qcs <- qcs[c("num_probes", "num_dt", "frac_dt", "num_probes_cg", "num_dt_cg", "frac_dt_cg", "num_probes_ch", "num_dt_ch", "frac_dt_ch")]
    colnames(qcs) <- c("N. Probes", "Detected Probes", "Detection Rate", "N. Probes(CG)", "Detected Probes(CG)", "Detection Rate(CG)", "N. Probes (CH)", "Detected Probes(CH)", "Detection Rate(CH)")
    qcs <- data.table::as.data.table(qcs, keep.rownames = "SampleName")
    data.table::fwrite(qcs, paste0(out_path,"/qc_summary.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(qcs)
    gc(reset=TRUE)

    ########################################
    # Create the raw beta matrix, including ctl and cgBK probes
    ########################################  
    logger::log_info("Creating the pre-filtering beta matrix...")
    raw_betas <- do.call(
        cbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
                sesame::getBetas(sdf)
            },
            BPPARAM = bp
        )
    )
    # raw_betas <- raw_betas[!grepl("^ctl|^cgBK", rownames(raw_betas)),]
    raw_betas <- as.data.table(raw_betas, keep.rownames = "ProbeID")
    data.table::fwrite(raw_betas, paste0(out_path,"/raw_beta_matrix.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(raw_betas)
    gc(reset=TRUE)

    ########################################
    # Create the beta matrix after correction and p-value filtering, excluding ctl and cgBK probes
    ########################################  
    logger::log_info("Creating the post-filtering beta matrix...")
    betas <- do.call(
        cbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
                sesame::getBetas(prepSesame(sdf, prep = "QCDPB", prep_args=list(P=list(pval.threshold=pval_threshold))))
            }, 
            BPPARAM = bp
        )
    )
    # NOTE: If there is only one sample, column names may be dropped; restore them.
    sample_name <- colnames(betas)
    betas <- as.data.table(betas, keep.rownames = "ProbeID")
    colnames(betas) <- c("ProbeID", sample_name)
    data.table::fwrite(betas, paste0(out_path,"/beta_matrix.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")

    ########################################
    # Add annotation
    ########################################  
    if (!skip_annotation) {
    # NOTE: Only RhelixaCustom uses genomic positions from the manifest. Other arrays use SeSAMe functions to obtain genomic positions.
    logger::log_info("Creating the annotated beta matrix...")
    if(array_type=="RhelixaCustom"){
        manifest_ann <- manifest[c("Probe_ID", "CHR","MAPINFO","Strand_FR")] %>%
            tidyr::drop_na(CHR) %>%
            dplyr::mutate(
                CHR=paste0("chr", CHR), 
                Strand_FR=ifelse(Strand_FR=="F", "+", "-"), 
                End=MAPINFO+1
            ) %>%
            dplyr::rename(ProbeID=Probe_ID, Chromosome=CHR, Start=MAPINFO, Strand=Strand_FR) %>%
            dplyr::filter(Chromosome != "chr0") %>%
            dplyr::mutate(ProbeName = ifelse(grepl("_", ProbeID), gsub("_(?:[^_]*)$", "", ProbeID), ProbeID)) %>%
            dplyr::select(ProbeID, ProbeName, Chromosome, Start, End, Strand)
    } else {
        manifest_ann <- sesameData::sesameData_getManifestGRanges(array_type) %>% as.data.frame(stringsAsFactors=FALSE) %>%
            dplyr::mutate(
                ProbeID = rownames(.),
                seqnames = as.character(seqnames),
                strand = as.character(strand),
            ) %>%
            dplyr::mutate(ProbeName = ifelse(grepl("_", ProbeID), gsub("_(?:[^_]*)$", "", ProbeID), ProbeID)) %>%
            dplyr::rename(Chromosome=seqnames, Start=start, End=end, Strand=strand) %>%
            dplyr::filter(Chromosome != "*") %>%
            dplyr::select(ProbeID, ProbeName, Chromosome, Start, End, Strand)
    }
    manifest_ann <- tryCatch({
        array_utils.add_ewas_atlas_annotation(manifest_ann, probeid_col="ProbeName")
    }, error = function(e) {
        logger::log_warn(paste0("EWAS Atlas annotation was skipped: ", conditionMessage(e)))
        manifest_ann
    })
    manifest_ann <- tryCatch({
        array_utils.add_closest_gene_annotation(manifest_ann, chr="Chromosome", start="Start", end="End", threads=threads)
    }, error = function(e) {
        logger::log_warn(paste0("Closest-gene annotation was skipped: ", conditionMessage(e)))
        manifest_ann
    })
    manifest_ann <- tryCatch({
        array_utils.add_gene_detail_annotation(manifest_ann, geneid_col="GeneID")
    }, error = function(e) {
        logger::log_warn(paste0("Gene-detail annotation was skipped: ", conditionMessage(e)))
        manifest_ann
    })
    manifest_ann <- manifest_ann %>% dplyr::select(-End) %>% dplyr::rename(Position=Start)
    betas$row_num <- 1:nrow(betas)
    betas_ann <- dplyr::right_join(manifest_ann, betas, by="ProbeID") %>%
        dplyr::arrange(row_num) %>%
        dplyr::select(-row_num) %>%
        dplyr::mutate(across(all_of(colnames(manifest_ann)), ~ ifelse(is.na(.), "", .))) %>% # Fill missing values in columns derived only from the manifest with empty strings
        dplyr::mutate(ProbeName = ifelse(grepl("_", ProbeID), gsub("_(?:[^_]*)$", "", ProbeID), ProbeID))
    data.table::fwrite(betas_ann, paste0(out_path,"/beta_matrix_ann.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(betas,betas_ann)
    gc(reset=TRUE)

    } else {
        logger::log_info("Skipping annotated beta matrix creation because skip_annotation=TRUE.")
        rm(betas)
        gc(reset=TRUE)
    }

    ########################################
    # Create the suffix-collapsed beta matrix by removing suffixes such as _BC11 and averaging probes with the same base ID
    ########################################  
    logger::log_info("Creating the suffix-collapsed beta matrix...")
    collapsed_betas <- do.call(
        cbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
                sesame::getBetas(prepSesame(sdf, prep = "QCDPB", prep_args=list(P=list(pval.threshold=pval_threshold))), collapseToPfx = TRUE)
            }, 
            BPPARAM = bp
        )
    )
    sample_name <- colnames(collapsed_betas)
    collapsed_betas <- collapsed_betas[!grepl("^ctl|^cgBK", rownames(collapsed_betas)),]
    collapsed_betas <- as.data.table(collapsed_betas, keep.rownames = "ProbeID")
    colnames(collapsed_betas) <- c("ProbeID", sample_name)
    data.table::fwrite(collapsed_betas, paste0(out_path,"/collapsed_beta_matrix.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")

    ########################################
    # Add annotation to the suffix-collapsed beta matrix
    ########################################  
    if (!skip_annotation) {
    logger::log_info("Creating the annotated suffix-collapsed beta matrix...")
    # Remove duplicated rows by ProbeName and keep the first occurrence
    collapsed_manifest_ann <- manifest_ann %>% dplyr::distinct(ProbeName, .keep_all=TRUE)
    collapsed_betas$row_num <- 1:nrow(collapsed_betas)
    collapsed_betas_ann <- dplyr::right_join(collapsed_manifest_ann, collapsed_betas, by=c("ProbeName"="ProbeID")) %>%
        dplyr::arrange(row_num) %>%
        dplyr::select(-row_num) %>%
        dplyr::mutate(across(all_of(colnames(collapsed_manifest_ann)), ~ ifelse(is.na(.), "", .))) %>% # Fill missing values in columns derived only from the manifest with empty strings
        dplyr::select(-ProbeID) %>%
        dplyr::rename(ProbeID=ProbeName)
    data.table::fwrite(collapsed_betas_ann, paste0(out_path,"/collapsed_beta_matrix_ann.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(collapsed_betas, collapsed_betas_ann)
    gc(reset=TRUE)
    
    } else {
        logger::log_info("Skipping annotated suffix-collapsed beta matrix creation because skip_annotation=TRUE.")
        rm(collapsed_betas)
        gc(reset=TRUE)
    }

    ########################################
    # Create the detection p-value matrix
    ######################################## 
    logger::log_info("Creating the p-value matrix...")
    pvals <- do.call(
        cbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
            sesame::pOOBAH(dyeBiasNL(inferInfiniumIChannel(sdf)), return.pval=TRUE)
            }, 
            BPPARAM = bp
        )
    )
    sample_name <- colnames(pvals)
    pvals <- pvals[!grepl("^ctl|^cgBK", rownames(pvals)),]
    pvals <- data.table::as.data.table(pvals, keep.rownames = "ProbeID")
    colnames(pvals) <- c("ProbeID", sample_name)
    data.table::fwrite(pvals, paste0(out_path,"/detection_pvals.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(pvals)
    gc(reset=TRUE)

    ########################################
    # Create the methylated-intensity matrix
    ######################################## 
    logger::log_info("Creating the methylated-intensity matrix...")
    meth_intensity = do.call(
        cbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
                sig_MU <- sesame::signalMU(sdf)
                setNames(sig_MU$M, sig_MU$Probe_ID)
            }, 
            BPPARAM = bp
        )
    )
    sample_name <- colnames(meth_intensity)
    meth_intensity <- meth_intensity[!grepl("^ctl|^cgBK", rownames(meth_intensity)),]
    meth_intensity <- as.data.table(meth_intensity, keep.rownames = "ProbeID")
    colnames(meth_intensity) <- c("ProbeID", sample_name)
    data.table::fwrite(meth_intensity, paste0(out_path,"/methylated_intensity.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(meth_intensity)
    gc(reset=TRUE)

    ########################################
    # Create the unmethylated-intensity matrix
    ######################################## 
    logger::log_info("Creating the unmethylated-intensity matrix...")
    unmeth_intensity = do.call(
        cbind, BiocParallel::bplapply(
            sdf_list, function(sdf) {
                sig_MU <- sesame::signalMU(sdf)
                setNames(sig_MU$U, sig_MU$Probe_ID)
            }, 
            BPPARAM = bp
        )
    )
    sample_name <- colnames(unmeth_intensity)
    unmeth_intensity <- unmeth_intensity[!grepl("^ctl|^cgBK", rownames(unmeth_intensity)),]
    unmeth_intensity <- as.data.table(unmeth_intensity, keep.rownames = "ProbeID")
    colnames(unmeth_intensity) <- c("ProbeID", sample_name)
    data.table::fwrite(unmeth_intensity, paste0(out_path,"/unmethylated_intensity.txt"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    rm(unmeth_intensity)
    gc(reset=TRUE)

    logger::log_info("Finish!")
}

########################################################################################

if (length(sys.calls()) == 0) {
	main()
}
