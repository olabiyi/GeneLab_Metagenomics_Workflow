#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to process contig or gene or function tables generated from assembly based processing.
# E-mail: obadbotanist@yahoo.com
# Created: January 2026
# example: Rscript process_assembly_table.R \
#                  --assembly-table 'Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv' \
#                  --assembly-summary 'assembly-summaries_GLmetagenomics.tsv' \
#                  --level  'Contig' \
#                  --type 'taxonomy'

###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 


option_list <- list(
  
  make_option(c("-t", "--assembly-table"), type="character",
              default=NULL, 
              help="path to a tab separated contig or gene level taxonomy or functions(KO) table. 
              i.e Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv or
              Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv or
              Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics",
              metavar="path"),
  make_option(c("-s", "--assembly-summary"), type="character",
              default=NULL, 
              help="path to assemblies summary file. assembly-summaries_GLmetagenomics.tsv or
              Assembly-based-processing-overview_GLmetagenomics.tsv",
              metavar="path"),
  make_option(c("-l", "--level"), type="character", default="Contig", 
              help="The assembly level at which taxonomy assignment was performed.
              i.e Contig or Gene. Default: Contig.",
              metavar="level"),
  make_option(c("-T", "--type"), type="character", default="taxonomy", 
              help="The type of assembly annotation performed in the input table.
              either taxonomy or KO. Default: taxonomy.",
              metavar="type"),
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLMetagenomics", 
              help="Genelab assay suffix.", metavar="GLMetagenomics"),
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)

library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(magrittr)
library(glue)

opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --taxonomy-table 'Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv'
                  --assembly-summary 'assembly-summaries_GLmetagenomics.tsv' 
                  --level  'Contig' 
                  --type 'taxonomy'" ,
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n  A script to reformat and dereplicate a contig or gene level taxonomy or functions table.",
                      "\nIt outputs a dereplicated contig or gene level species or functions table.",
                      sep="")
)


opt <- parse_args(opt_parser)

# print(opt)
# stop()


if (opt$version) {
  cat("process_assembly_table.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["assembly-table"]])) {
  stop("Path to a tab separated contig or gene level taxonomy or functions(KO) table file must be set.")
}

if(is.null(opt[["assembly-summary"]])) {
  stop("Path to a tab separated assemblies summary table file must be set.")
}


process_taxonomy <- function(taxonomy, prefix='\\w__') {
  # Function to process a taxonopmy assignment table
  #1. ~ taxonomy is a string specifying the taxonomic assignment file name
  #2 prefix ~ is a regular expression specifying the characters to remove
  # from the taxon names
  
  
  taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
  

  for (rank in colnames(taxonomy)) {
    #delete the taxonomy prefix
    taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],
                            replacement = '')
    indices <- which(is.na(taxonomy[,rank]))
    taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
    #replace empty cell
    indices <- which(taxonomy[,rank] == "")
    taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
  }
  taxonomy <- apply(X = taxonomy,MARGIN = 2,
                    FUN =  gsub,pattern = "_",replacement = " ") %>% 
    as.data.frame(stringAsfactor=F)
  return(taxonomy)
}

# Function to format a taxonomy assignment table by appending a suffix
# to a known name
format_taxonomy_table <- function(taxonomy=taxonomy.m,stringToReplace="Other",
                                  suffix=";Other") {
  
  for (taxa_index in seq_along(taxonomy)) {
    #indices <- which(taxonomy[,taxa_index] == stringToReplace)
    
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }
  return(taxonomy)
}

fix_names<- function(taxonomy,stringToReplace,suffix){
  #1~ taxonomy is a taxonomy dataframe with taxonomy ranks as column names
  #2~ stringToReplace is a vector of regex strings specifying what to replace
  #3~ suffix is a string specifying the replacement value
  
  
  for(index in seq_along(stringToReplace)){
    taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                      stringToReplace=stringToReplace[index], 
                                      suffix=suffix[index])
  }
  return(taxonomy)
}


# Function to read an input table into a dataframe
read_input_table <- function(file_name){
  
  # Get depth from file name
  df <- read_delim(file = file_name, comment = "#")
  return(df)
  
}


get_samples <- function(df, sample_names, end_col='species') {

  # Get common samples 
  cols <- colnames(df)
  start <- grep(end_col, cols)+1
  end <- length(cols)
  df_samples <- cols[start:end]
  sample_names <- intersect(df_samples, sample_names)

  return(sample_names)

}

# Read Assembly-based contig or gene taxonomy annotation table
read_taxonomy_table <- function(df, sample_names){
  
  taxonomy_table <- df %>%
    select(domain:species) %>%
    mutate(domain=replace_na(domain, "Unclassified"))
 
  sample_names <- get_samples(df, sample_names) 
  counts_table <- df %>% select(!!sample_names)
  
  taxonomy_table  <- process_taxonomy(taxonomy_table)
  taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")
  
  df <- bind_cols(taxonomy_table, counts_table)
  
  return(df)
}

assembly_summary <- opt[['assembly-summary']]
assembly_table <- opt[['assembly-table']]
level <- opt[['level']]
type <- opt[['type']]
prefix <- opt[["output-prefix"]]
suffix <- opt[["assay-suffix"]] 


# Read in assembly summary table
overview_table <-  read_input_table(assembly_summary) %>%
  select(
    where( ~all(!is.na(.)) )
  )

col_names <- names(overview_table) %>% str_remove_all("-assembly")

if(any(str_detect(col_names, "gene_calls_identified"))){
 # Input file is an Assembly passed processing overview file
  # Get the sample names for which genes were called
 sample_order <- overview_table %>%
                   filter(gene_calls_identified == "Yes") %>%
                   pull(Sample_ID) %>% sort()
    
}else {
  
  # Input is an assembly summary file
  # Get the sample names for which assemblies were generated
  sample_order <- col_names[-1] %>% sort()

}


if(type == "KO"){

  df <- read_input_table(assembly_table)
  # Get common sample ids
  sample_order <- get_samples(df, sample_order, "KO_function")
  
  table2write <- df  %>% select(KO_ID, !!sample_order)
  
}else{

# Deduplicate rows by summing together species values

  df <- read_input_table(assembly_table)
  # Get common sample ids
  sample_order <- get_samples(df, sample_order)
  table2write <- read_taxonomy_table(df, sample_order) %>%
               select(species, !!sample_order) %>%
               group_by(species) %>%
               summarise(across(everything(), sum)) %>% # Dereplicate
               filter(species != "Unclassified;_;_;_;_;_;_") %>% # Drop unclassified
               as.data.frame()

}

colnames(table2write) <-  colnames(table2write) %>% str_remove_all(suffix)
# type - taxonomy or KO
# level - Contig or Gene

if(type == "KO"){ type <- "KO-function" }
write_tsv(x = table2write, file = glue("{prefix}Combined-{tolower(level)}-level-{type}{suffix}.tsv"))
