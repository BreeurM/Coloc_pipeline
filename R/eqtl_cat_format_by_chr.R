.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(tictoc)))

# Get study ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)

setwd("/gpfs3/well/travis-prostate/projects/Coloc_pipeline")

# Check if study ID is provided
# Check if both study ID and chromosome are provided
if (length(args) < 2) {
  stop("Please provide a study ID and chromosome as arguments (e.g., QTD000261 6)", call. = FALSE)
}


study_id <- args[1]
chr <- as.character(args[2])

tic("Extracting traits id for given study and chromosome")


# Construct file paths using study ID
input_file <- file.path("data/eQTL_catalogue", paste0(study_id, ".lbf_variable.txt.gz"))
output_file <- file.path("data/eQTL_catalogue/trait_ids/", paste0(study_id,"_chr", chr, "_trait_ids.txt"))

# Read the dataset
data <- fread(input_file)

data <- data %>% filter(chromosome == chr)

# Extract unique combinations of molecular_trait_id and chromosome
trait_ids <- data %>%
  select(molecular_trait_id, chromosome) %>%
  distinct()


# Write results to file
write.table(trait_ids,
            file = output_file,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")

cat("Trait IDs have been written to", output_file, "\n")

toc()