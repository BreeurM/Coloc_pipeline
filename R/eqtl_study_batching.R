.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

args <- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(readr)

setwd("/gpfs3/well/travis-prostate/projects/Coloc_pipeline")

batch_size <- as.double(args[1])

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






