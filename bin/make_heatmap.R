#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to create a heatmap.
# E-mail: obadbotanist@yahoo.com
# Created: November 2025
# example: Rscript make_heatmap.R \
#                  --metadata-table 'mapping/metadata.txt' \
#                  --feature-table 'Gene_function_table_filtered_GLMetagenomics.csv' \
#                  --group-column 'Description' \
#                  --samples-column 'sample_id' \
#                  --output-prefix  'Combined-gene-level-KO-function' \
#                  --assay-suffix   '_GLMetagenomics'

###############################################################################

library(optparse)

######## -------- Get input variables from the command line ----##############

version <- 1.0 

# Input options 
option_list <- list(
  
  make_option(c("-m", "--metadata-table"), type="character", default=NULL, 
              help="path to a tab or comma separated samples metadata file with the group to be plotted.",
              metavar="path"),
  
  make_option(c("-i", "--feature-table"), type="character", default=NULL, 
              help="path to a tab separated samples feature table 
              i.e. species/functions table with species/functions as the first column 
              and samples as other columns.",
              metavar="path"),

  make_option(c("-s", "--samples-column"), type="character", default="Sample Name", 
              help="Column in metadata containing the sample names in the feature table. \
                    Deafault: 'Sample Name' ",
              metavar="Sample_Name"),

  make_option(c("-g", "--group-column"), type="character", default=NULL, 
              help="Column in metadata to be use to facet/group plot. 
              Default: NULL ",
              metavar="group_column"),
  
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar="prefix"),
  
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLMetagenomics", 
              help="Genelab assay suffix.", metavar="GLMetagenomics"),
  
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)




opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --metadata-table 'metadata.csv' \\
                  --feature-table 'Gene_function_table_filtered_GLMetagenomics.csv' \\
                  --group-column 'Description' \\
                  --output-prefix  'Combined-gene-level-KO-function'" ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to create a heatmap.",
                      "\nIt outputs a static heatmap with features as rows and samples as columns.",
                      sep="")
)


opt <- parse_args(opt_parser)


if (opt$version) {
  cat("make_heatmap.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["metadata-table"]])) {
  stop("Path to a metadata file must be set.")
}

if(is.null(opt[["feature-table"]])) {
  stop("Path to a feature table e.g. species/functions table file must be set.")
}

if(is.null(opt[["group-column"]])) {
  stop("Please supply a string specifying the samples grouping column in your provided metadata.")
}

library(glue)
library(pheatmap)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(magrittr)


# Remove white colors based in there RGB values
cols <- colors()
rgb_vals <- col2rgb(cols) / 255  # normalize to 0–1
is_whitish <- colMeans(rgb_vals) > 0.8 # define "whitish"
# keep only NON-whitish colors
filtered_cols <- cols[!is_whitish]
colors2use  <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00",
                 "#CAB2D6","#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                 "#666666","#1B9E77", "#E6AB02","#00FFFFFF", "#B2182B","#FDDBC7","#D1E5F0","#CC0033",
		 "#FF00CC","#330033", "#999933","#FF9933","#FFFAFAFF", filtered_cols)

colours <- colorRampPalette(c('white','red'))(255)


feature_table_file <- opt[["feature-table"]] # 'kaiju_species_table_GLlbnMetag.csv'
metdata_file <-  opt[["metadata-table"]] # "metadata.csv"
samples_column <-  opt[["samples-column"]] # 'Sample_ID', 'Barcode'
group_column <-   opt[["group-column"]] # 'Description'
prefix <-  opt[["output-prefix"]] # "filtered-kaiju_species"
suffix <- opt[["assay-suffix"]]  # "GLlbnMetag"

#Prepare feature table
feature_table <- read_delim(feature_table_file) %>%  as.data.frame()
rownames(feature_table) <- feature_table[[1]]
feature_table <- feature_table[,-1] %>% as.matrix()
colnames(feature_table) <-  colnames(feature_table) %>% str_remove_all("barcode") %>% str_remove_all(suffix)

# Prepare metadata
metadata <- read_delim(metdata_file) %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column] %>% str_remove_all("barcode") %>% str_remove_all(suffix)

# GFet common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(feature_table), rownames(metadata))
feature_table <- feature_table[, common_samples]
metadata <- metadata[common_samples,]
metadata <- metadata %>% arrange(!!sym(group_column))

# Create column annotation
col_annotation <- as.data.frame(metadata)[,group_column, drop=FALSE]
rownames(col_annotation) <-  rownames(col_annotation)

# Calculate output plot width and height
number_of_samples <- ncol(feature_table)
width <- 1 * number_of_samples
if(width < 10) { width <- 10} # Set the mininum width to 10 inches
if(width > 100) { width <- 100} # Set the maximum width to 100 inches
number_of_features <- nrow(feature_table)
height <- 0.2 * number_of_features

if(height < 10) { height <- 10} # Set the mininum height to 10 inches
if(height > 100) { height <- 100} # Set the maximum height to 100 inches (highest that won't generate an error)

# Set colors by group
groups <- metadata[[group_column]] %>%  unique()
number_of_groups <-  length(groups)
my_colors <- colors2use[1:number_of_groups]
names(my_colors) <- groups
annotation_colors  <- list(my_colors)
names(annotation_colors) <- group_column

# Create heatmap
# All features may be difficult to visualize
png(filename = glue("{prefix}_heatmap{suffix}.png"), width = width,
    height = height, units = "in", res=300)
pheatmap(mat = feature_table[,rownames(col_annotation)],
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 90, 
         display_numbers = TRUE,
         fontsize = 12, 
         annotation_col = col_annotation,
         annotation_colors = annotation_colors ,
         number_format = "%.0f")
dev.off()


sorted_features <- rowSums(feature_table) %>% sort(decreasing = TRUE)

# Plot only top 50 features as it is often difficult to visualize all features at once
if(length(sorted_features) >= 50) { 
  
  top50 <- sorted_features[1:50]

png(filename = glue("{prefix}_top_50_heatmap{suffix}.png"), width = width,
    height = 12, units = "in", res=300)
pheatmap(mat = feature_table[names(top50), rownames(col_annotation)],
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 90, 
         display_numbers = TRUE,
         fontsize = 12, 
         annotation_col = col_annotation,
         annotation_colors = annotation_colors ,
         number_format = "%.0f")
dev.off()

}
