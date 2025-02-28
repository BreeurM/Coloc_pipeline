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


.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

library(dplyr)

setwd("/gpfs3/well/travis-prostate/projects/Coloc_pipeline")

batch_size <- 250000

cov <- read_table("data/eQTL_catalogue/study_coverages.txt", col_names = FALSE)
colnames(cov) <- c("study", "genes")

# Exclude the last study with no data at all

cov <- cov[1:nrow(cov)-1,]

cov <- cov %>% mutate(batch = ceiling(cumsum(genes)/batch_size))


for (b in 1:max(cov$batch)){
  filename <- paste0("data/eQTL_catalogue/studies_batch_", as.character(b),".txt")
  write.table(cov %>% filter(batch == b), filename, 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}

