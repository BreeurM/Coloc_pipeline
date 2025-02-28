.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))

# Get study ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if study ID is provided
if (length(args) < 1) {
  stop("Please provide a study ID as an argument (e.g., QTD000261)", call. = FALSE)
}

study_id <- args[1]

# Construct file paths using study ID
input_file <- file.path("data/eQTL_catalogue", paste0(study_id, ".lbf_variable.txt.gz"))

# Read the dataset
data <- fread(input_file, select = c("molecular_trait_id", "chromosome"))

# Extract unique combinations of molecular_trait_id and chromosome
trait_ids <- unique(data, by = c("molecular_trait_id", "chromosome")) %>%
  mutate(chromosome = ifelse(chromosome == "X", 23, chromosome))

# Split into groups by chromosome
chromosome_groups <- split(trait_ids, trait_ids$chromosome)

# Write one file per chromosome
for (chr in names(chromosome_groups)) {
  output_file <- file.path("data/eQTL_catalogue/trait_ids/", paste0(study_id, "_chr", chr, "_trait_ids.txt"))
  write.table(chromosome_groups[[chr]],
              file = output_file,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t")
}
cat(study_id, "processed", "\n")