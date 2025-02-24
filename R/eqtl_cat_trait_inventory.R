.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

# Load required libraries
library(data.table)
library(dplyr)
library(stringr)

setDTthreads(1)

setwd("/gpfs3/well/travis-prostate/projects/Coloc_pipeline/data/eQTL_catalogue")

# Get list of files
file_list <- list.files(pattern = "^QTD.*\\.cc\\.tsv\\.gz$")

gene_count = data.frame(study = character(), distinct_genes = double())

# Loop through each file and count distinct genes
for (file in sort(file_list)) {
  # Read the gzipped TSV file
  df <- fread(cmd = paste("zcat", file), select = "gene_id")
  
  # Check if 'gene_id' column exists
  if ("gene_id" %in% names(df)) {
    unique_genes <- uniqueN(df$gene_id)
    cat(file, ": ", unique_genes, " distinct genes\n")
    gene_count <- rbind(gene_count, data.frame(study = str_split_fixed(file,"\\.",2)[1], 
                                               distinct_genes = unique_genes))
  } else {
    cat(file, ": 'gene_id' column not found\n")
  }
  rm(df)       # Explicitly remove large object
  gc()         # Force garbage collection
}

write.table(gene_count %>% arrange(desc(distinct_genes)),
            file = "study_coverages.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")
