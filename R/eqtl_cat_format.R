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
output_file <- file.path("data/eQTL_catalogue", paste0(study_id, "_trait_ids.txt"))

# Read the dataset
data <- fread(input_file)

# Extract unique combinations of molecular_trait_id and chromosome
trait_ids <- data %>%
  select(molecular_trait_id, chromosome) %>%
  distinct()

# Convert chromosome X to 23
trait_ids <- trait_ids %>% mutate(
  chromosome = ifelse(chromosome == "X", 23, chromosome)
)

# Write results to file
write.table(trait_ids,
            file = output_file,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")

cat("Trait IDs have been written to", output_file, "\n")