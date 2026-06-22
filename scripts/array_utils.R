#!/usr/bin/env Rscript

# Description: Utility functions for methylation array processing
# Function names should be prefixed with array_utils.

####################################################################

# Get the database path
array_utils.database_path <- function(){
	suppressMessages(library(this.path))
	suppressMessages(library(RcppTOML))
	config_path <- paste0(this.path::this.dir(),"/../utils/config.toml")
	database_path <- RcppTOML::parseTOML(config_path)$database$path
	return(database_path)
}

# Get the manifest file path
array_utils.manifest_path <- function(array_type){
	database_path <- array_utils.database_path()
	manifest_dir_path <- paste0(database_path,"/MethylationArrayManifest/manifest_v25.01.13")
	if (array_type == "EPICv2") {
		manifest_path <- paste0(manifest_dir_path,"/EPIC-8v2-0_A1.csv")
	} else if (array_type == "MSA") {
		manifest_path <- paste0(manifest_dir_path,"/MSA-48v1-0_20102838_A1.csv")
	} else if (array_type == "RhelixaCustom") {
		manifest_path <- paste0(manifest_dir_path,"/RhelixaMethylEpiclockV1_GS_20042798X391986_A1_All_renamed.csv")
	} else {
		stop("array_type must be one of EPICv2, MSA, or RhelixaCustom")
	}
	return(manifest_path)
}

# Get the SeSAMeData path
array_utils.sesame_data_path <- function(){
	database_path <- array_utils.database_path()
	sesame_data_path <- paste0(database_path,"/SesameData/SesameData_25.02.06")
	return(sesame_data_path)
}

# Read a manifest file, convert it to a format compatible with SeSAMe, and return it
# This assumes that the manifest contains [Assay] and [Controls] sections
array_utils.read_manifest <- function(manifest_path, threads) {

	suppressMessages(library(data.table))

    data.table::setDTthreads(threads)

    lines <- readLines(manifest_path)
    section_starts <- grep("^\\[", lines)
    assay_start <- grep("^\\[Assay\\]", lines)
    control_start <- grep("^\\[Controls\\]", lines)
    assay_end <- control_start - 1
    control_end <- length(lines)

    manifest_df <- data.table::fread(manifest_path, skip=assay_start, nrows=(assay_end - assay_start - 1), na.strings = c("", "NA"))
    cols <- colnames(manifest_df)
    newcols <- replace(cols, cols == "IlmnID", "Probe_ID")
    newcols <- replace(newcols, newcols == "AddressA_ID", "U")
    newcols <- replace(newcols, newcols == "AddressB_ID", "M")
    colnames(manifest_df) <- newcols

    control_df <- data.table::fread(manifest_path, skip=control_start, nrows=(control_end - control_start), header=FALSE, na.strings = c("", "NA"))
    colnames(control_df) <- colnames(manifest_df)
    control_df$U <- control_df$Probe_ID
	control_df$AlleleA_ProbeSeq <- gsub("[^A-Za-z0-9_]", "_", control_df$AlleleA_ProbeSeq)
	control_df$AlleleA_ProbeSeq <- gsub("_$", "", control_df$AlleleA_ProbeSeq)
	control_df$AlleleA_ProbeSeq <- gsub("__+", "_", control_df$AlleleA_ProbeSeq)
	control_df$AlleleA_ProbeSeq <- paste0("ctl_", control_df$AlleleA_ProbeSeq)
	control_df$Probe_ID <- control_df$AlleleA_ProbeSeq
	control_df$Name <- NA
	control_df$AlleleA_ProbeSeq <- NA

    manifest_df <- as.data.frame(rbind(manifest_df, control_df))

    return(manifest_df)
}

# Add nearest-gene annotation
array_utils.add_closest_gene_annotation <- function(df, chr, start, end, threads=4){
	
	suppressMessages(library(ChIPpeakAnno))
	suppressMessages(library(org.Hs.eg.db))
	suppressMessages(library(dplyr))
	suppressMessages(library(BiocParallel))
	data("TSS.human.GRCh38")

	df_ann <- do.call(
		dplyr::bind_rows, BiocParallel::bplapply(
			split(df, cut(seq_len(nrow(df)), breaks=threads, labels=FALSE)), function(df_split) {
				if ("Tag" %in% colnames(df_split)) {
					stop("The column name 'Tag' cannot be used")
				}	
				df_split <- dplyr::mutate(df_split, Tag=rownames(df_split))
				df_bed <- df_split[c(chr, start, end, "Tag")] %>% dplyr::rename(seqnames = chr, start = start, end = end)
				df_bed <- df_bed[complete.cases(df_bed),] # Remove rows containing missing values
				df_GR <- suppressMessages(ChIPpeakAnno::toGRanges(df_bed, format="bed"))
				df_GR_ann <- ChIPpeakAnno::annotatePeakInBatch(df_GR, AnnotationData=TSS.human.GRCh38)
				# Reorder columns and standardize column names
				df_ann <- tibble::as_tibble(df_GR_ann) %>%
					dplyr::distinct(Tag, .keep_all=TRUE) %>%
					dplyr::select(-seqnames, -start, -end, -width, -strand, -peak, -shortestDistance, -fromOverlappingOrNearest) %>%
					dplyr::rename(GeneID = feature, Gene.Start = start_position, Gene.End = end_position, Gene.Strand = feature_strand, Gene.Inside = insideFeature, Gene.Distance = distancetoFeature)
				# Join the annotation back to the original table
				df_ann_split <- dplyr::left_join(df_split, df_ann, by="Tag") %>% dplyr::select(-Tag)
				return(df_ann_split)
			}, 
			BPPARAM = BiocParallel::MulticoreParam(threads)
		)
	)
	return(df_ann)
}

# Add detailed gene annotation information; this is copied from the RNA-seq annotation workflow
array_utils.add_gene_detail_annotation <- function(df, geneid_col = "GeneID") {
	
	suppressMessages(library(dplyr))
	suppressMessages(library(readr))
	
	ann_txt_path <- paste0(this.path::this.dir(),"/gene_annotation.txt")
	ann <- readr::read_tsv(ann_txt_path, show_col_types = FALSE)
	ann <- ann[c("GeneID", "Symbol", "Description", "Cellular_Component_ID", "Cellular_Component_Term", "Molecular_Function_ID", "Molecular_Function_Term", "Biological_Process_ID", "Biological_Process_Term")] %>%
		dplyr::rename(Gene.Symbol = Symbol, Gene.Description = Description)
	df_ann <- dplyr::left_join(df, ann, by = setNames(geneid_col, "GeneID"))

	return(df_ann)
}

# Add EWAS Atlas trait information
array_utils.add_ewas_atlas_annotation <- function(df, probeid_col = "ProbeID") {

	suppressMessages(library(dplyr))
	suppressMessages(library(readr))

	database_path <- array_utils.database_path()
	ewas_traits_tsv_path <- paste0(database_path, "/EWAS_Atlas/EWAS_Atlas_25.03.17/EWAS_Atlas_traits.tsv")
	ewas_traits <- readr::read_tsv(ewas_traits_tsv_path, show_col_types = FALSE) %>%
		dplyr::group_by(ProbeID) %>%
		dplyr::summarize(Traits = paste(Traits, collapse = " | "))
	df_ann <- df %>%
		dplyr::mutate(ProbeID_no_suffix = ifelse(grepl("_", .data[[probeid_col]]), gsub("_(?:[^_]*)$", "", .data[[probeid_col]]), .data[[probeid_col]])) %>%
		dplyr::left_join(ewas_traits, by = c("ProbeID_no_suffix" = "ProbeID")) %>%
		dplyr::select(-ProbeID_no_suffix) %>%
		dplyr::rename(EWAS_Atlas.Traits = Traits)
	return(df_ann)
}