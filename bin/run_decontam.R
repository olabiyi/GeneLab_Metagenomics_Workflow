#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to perform feature decontamination with decontam.
# E-mail: obadbotanist@yahoo.com
# Created: January 2026
# example: Rscript run_decontam.R \
#                  --feature-table 'kaiju_species_table_GLlbnMetag.tsv' \
#                  --feature-column 'Species' \
#                  --metadata-table 'mapping/metadata.csv' \
#                  --samples-column 'Sample_ID' \
#                  --prevalence-column 'NTC' \
#                  --frequency-column 'concentration' \
#                  --ntc-name  'true' \
#                  --threshold 0.1 \
#                  --classification-method 'kaiju' \
#                  --assay-suffix '_GLlbnMetag'
###############################################################################

library(optparse)

######## -------- Get input variables from the command line ----##############

version <- 1.0

# Input options
option_list <- list(

  make_option(c("-i", "--feature-table"), type = "character", default = NULL,
              help = "path to a tab or comma separated samples feature table 
                      i.e. species/functions table with species/functions as 
                      the first column and samples as other columns.",
              metavar = "path"),

  make_option(c("-f", "--feature-column"), type = "character", default = NULL,
              help = "Feature column name in feature table allowed values: 
                      ['Species', 'species', 'KO_ID', 'Pathway', 'KO', 'Uniref90']. 
                      Default: empty string",
              metavar = "Feature_Column"),

  make_option(c("-m", "--metadata-table"), type = "character", default = NULL,
              help = "path to a comma separated samples metadata file with the 
                      --prevalence-column and/or --frequency-column to be analyzed.",
              metavar = "path"),

  make_option(c("-s", "--samples-column"), type = "character", default = "Sample Name",
              help = "Column in metadata containing the sample names in the feature table. \
                      Default: 'Sample Name' ",
              metavar = "Sample_Name"),

  make_option(c("-p", "--prevalence-column"), type = "character", default = NULL,
              help = "Column in metadata to be use for prevalence based analysis. 
                      Default: 'Sample_or_Control' ",
              metavar = "prevalence_column"),

  make_option(c("-n", "--ntc-name"), type = "character", default = "Control_Sample",
              help = "Name of NTC in prevalence column. \
                      Default: 'Control_Sample' ",
              metavar = "NTC_Name"),

  make_option(c("-F", "--frequency-column"), type = "character", default = NULL,
              help = "Column in metadata to be use for frequency based analysis. \
                      Default: 'concentration' ",
              metavar = "frequency_column"),

  make_option(c("-t", "--threshold"), type = "numeric", default = 0.5,
              help = "Decontam's threshold (0-1) for both prevalence and frequency based analyses. \
                      Default: 0.1 ",
              metavar = "threshold"),

  make_option(c("-M", "--classification-method"), type = "character", default = "",
              help = "Taxonomy or functional method used to generate the input 
              feature table. The supplied string will be added to output file names
              ['kaiju', 'kraken2', 'metaphlan', 'contig-taxonomy',
               'gene-taxonomy', 'gene-function','Pathway-abundances', 
	       'Gene-families-KO', 'Gene-families-uniref' ]. Default: empty string.",
              metavar = ""),

  make_option(c("-o", "--output-prefix"), type = "character", default = "",
              help = "Unique name to tag onto output files. Default: empty string.",
              metavar = ""),

  make_option(c("-a", "--assay-suffix"), type = "character", default = "_GLMetagenomics",
              help = "Genelab assay suffix or empty string. One of ['_GLMetagenomics', 'GLlbsMetag', 'GLlblMetag', '']",
              metavar = "GLMetagenomics"),

  make_option(c("--version"), action = "store_true", type = "logical",
              default = FALSE,
              help = "Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list = option_list,
  usage = "Rscript %prog \\
                  --metadata-table 'mapping/metadata.csv' \\
                  --feature-table 'kaiju_species_table_GLlbnMetag' \\
                  --feature-column 'Species' \\ # ['Species', 'species', 'KO']
                  --prevalence-column 'NTC' \\
                  --frequency-column 'concentration' \\
                  --classification-method 'taxonomy' \\
                  --threshold 0.1 " ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to perform feature decontamination with decontam.",
                      "\nIt outputs decontam's primary output.",
                      "along with a decontaminated feature table",
                      sep = "")
)


opt <- parse_args(opt_parser)


if (opt$version) {
  cat("run_decontam.R version: ", version, "\n")
  options_tmp <- options(show.error.messages = FALSE)
  on.exit(options(options_tmp))
  stop()
}



if (is.null(opt[["metadata-table"]])) {
  stop("Path to a metadata file must be set.")
}

if (is.null(opt[["feature-table"]])) {
  stop("Path to a feature table e.g. species/functions table file must be set.")
}

if (opt[["samples-column"]] == "Sample Name") {
  message("I will assume that the sample names are in a column named 'Sample Name' \n")
}


library(decontam)
library(phyloseq)
library(glue)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(magrittr)

# Feature table decontamination with decontam
run_decontam <- function(feature_table, metadata, contam_threshold = 0.1,
                         prev_col = NULL, freq_col = NULL, ntc_name = "true") {

  # Retain metadata for only the samples present in the input feature table
  sub_metadata <- metadata[colnames(feature_table), ]
  # Modify NTC concentration
  # Often times the user may set the NTC concentration to zero because they think nothing
  # should be in the negative control but decontam fails if the value is set to zero.
  # To prevent decontam from failing, we replace zero with a very small concentration value
  # 0.0000001
  if (!is.null(freq_col)) {

    sub_metadata <- sub_metadata %>%
      mutate(!!freq_col := map_dbl(!!sym(freq_col), .f = function(conc) {
            if (conc == 0) return(0.0000001) else return(conc)
          }
        )
      )
    sub_metadata[, freq_col] <- as.numeric(sub_metadata[, freq_col])

  }
  
  # Create phyloseq object
  ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE), sample_data(sub_metadata))

  # In our phyloseq object, `prev_col` is the sample variable that holds the negative
  # control information. We'll summarize the data as a logical variable, with TRUE for control
  # samples, as that is the form required by isContaminant.
  sd <- as.data.frame(sample_data(ps)) # Extract sample metadata
  sd[, "is.neg"] <- 0 # Initialize
  sd[, "is.neg"] <- sample_data(ps)[[prev_col]] == ntc_name # Assign boolean value
  sample_data(ps) <- sd

  # Run Decontam
  if (!is.null(freq_col) && !is.null(prev_col)) {

    # Run decontam in both prevalence and frequency modes
    contamdf <- isContaminant(ps, neg = "is.neg", conc = freq_col, threshold = contam_threshold)

  } else if (!is.null(freq_col)) {

    # Run decontam in frequency mode
    contamdf <- isContaminant(ps, conc = freq_col, threshold = contam_threshold)

  } else if(!is.null(prev_col)){

    # Run decontam in prevalence mode
    contamdf <- isContaminant(ps, neg = "is.neg", threshold = contam_threshold)

  } else {

    cat("Both freq_col and prev_col cannot be set to NULL.\n")
    cat("Please supply either one or both column names in your metadata")
    cat("for frequency and prevalence based analysis, respectively\n")
    stop()

  }

  return(contamdf)
}


# Set input variables
feature_table_file <- opt[["feature-table"]] # 'kaiju_species_table_GLlbnMetag.csv'
metadata_file <- opt[["metadata-table"]] # "metadata.csv"
samples_column <-  opt[["samples-column"]] # 'Sample_ID'
freq_col <- opt[["frequency-column"]] # "input_conc_ng"
prev_col <- opt[["prevalence-column"]] # "NTC"
threshold <- opt[["threshold"]] # 0.5
ntc_name <- opt[["ntc-name"]] # "true"
# "kaiju", "kraken2", "metaphlan", "contig-taxonomy", "gene-taxonomy",
# "gene-function", "Pathway-abundances", "Gene-families-KO", "Gene-families-uniref"
method <- opt[["classification-method"]] # 'kaiju'
feature_column <- opt[["feature-column"]] # 'Species'
prefix <- opt[["output-prefix"]]
suffix <- opt[["assay-suffix"]] # GLlbnMetag



# Prepare feature table
feature_table <- read_delim(feature_table_file) %>%  as.data.frame()
row.names(feature_table) <- feature_table[[1]]
feature_table <- feature_table[, -1] %>% as.matrix()
colnames(feature_table) <-  colnames(feature_table) %>% str_remove_all("barcode") %>% str_remove_all(suffix)

# Prepare metadata
metadata <- read_delim(metadata_file)

# cast prevalence column (if present) to character and convert to lowercase
if (!is.null(prev_col) && prev_col %in% colnames(metadata)) {
  metadata <- metadata %>%
              mutate(!!prev_col := as.character(!!sym(prev_col)) %>%
              tolower()) %>%
              as.data.frame()
  # convert ntc_name to lowercase to match prevalence column setting above
  ntc_name <- tolower(ntc_name)
} else {
  metadata <- metadata %>% as.data.frame()
}
row.names(metadata) <- metadata[, samples_column]

# Subset metadata and feature table  to contain the samples
samples <- intersect(colnames(feature_table), rownames(metadata))
metadata <- metadata[samples, ]
feature_table <- feature_table[, samples]

if (method == "gene-function")  {

  name <- "Combined-gene-level-KO-function"

} else if (method == "gene-taxonomy") {

  name <- "Combined-gene-level-taxonomy"

} else if (method == "contig-taxonomy") {

  name <- "Combined-contig-level-taxonomy"

} else {

  name <- method

}

# Run decontam
# Assign prev and freq column names to NULL if the values in the supplied columns aren't unique
if (length(unique(metadata[,prev_col])) == 1) prev_col <- NULL
if (length(unique(metadata[,freq_col])) == 1) freq_col <- NULL

# Error if values in both prevalence and frequency columns are not different between samples within each column 
#i.e no difference between negative control(s) and other samples
if (is.null(freq_col) && is.null(prev_col)) {

  text2write <- "Values in both NTC and concentration columns are not unique between samples within each column.\nTherefore, feature decontamination with decontam cannot be performed."

  file_name <- glue("{prefix}{name}_decontam_failure.txt")

  cat(text2write, file = file_name)

  contamdf <- data.frame(x = rownames(feature_table), freq = NA,
        prev = NA, p.freq = NA, p.prev = NA, p = NA, contaminant = FALSE)
  colnames(contamdf)[1] <- feature_column

} else {

  contamdf <- run_decontam(feature_table, metadata, threshold, prev_col, freq_col, ntc_name) 
  contamdf <- as.data.frame(contamdf) %>% rownames_to_column(feature_column)

}

# Write decontaminated feature table and decontam's primary results
outfile <- glue("{prefix}{name}_decontam_results{suffix}.tsv")
write_tsv(x = contamdf, file = outfile)

taxonomy_methods <- c("kaiju", "kraken2", "metaphlan", "gene-taxonomy", "contig-taxonomy")

# Get the list of contaminants identified by decontam
contaminants <- contamdf %>%
                   filter(contaminant == TRUE) %>%
                   pull(!!sym(feature_column))

# Drop contaminants(s) if detected
if (length(contaminants) > 0) {

  # Drop contaminant features identified by decontam
  decontaminated_table <- feature_table %>%
    as.data.frame() %>%
    rownames_to_column(feature_column) %>%
    filter(str_detect(!!sym(feature_column),
                      pattern = str_c(contaminants,
                                      collapse = "|"),
                      negate = TRUE))

  rownames(decontaminated_table) <- decontaminated_table[[feature_column]]

  # Add _species string to output file name if it is a taxonomy method
  if (any(grepl(pattern = method, x = taxonomy_methods))) {
    outfile <- glue("{prefix}{name}_decontam_species_table{suffix}.tsv")
  } else {
    outfile <- glue("{prefix}{name}_decontam_table{suffix}.tsv")
  }

  write_tsv(x = decontaminated_table, file = outfile)

} else {
  message("No contaminant was detected by Decontam")
}
