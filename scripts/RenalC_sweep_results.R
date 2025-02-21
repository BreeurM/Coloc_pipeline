.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))

# Get study ID from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if study ID is provided
if (length(args) < 1) {
  stop("Please provide a study ID as an argument (e.g., QTD000261)", call. = FALSE)
}

study_id <- args[1]

setwd("/gpfs3/well/travis-prostate/projects/Coloc_pipeline")

# Construct results path using study ID
respath <- file.path("results", paste0(study_id, "_RenalC"))

# Get all chromosome directories
chr_dirs <- list.dirs(respath, full.names = TRUE, recursive = FALSE)
chr_dirs <- grep("chr_\\d+", chr_dirs, value = TRUE)  # Select only chromosome directories

# Initialize list to store results
all_coloc <- list()

for (chr_dir in chr_dirs) {
  
  print(paste("Processing", chr_dir))
  # Get chromosome number from directory name
  chr <- sub(".*chr_(\\d+)", "\\1", chr_dir)
  
  # Get all relevant .rds files in chromosome directory
  rds_files <- list.files(
    path = chr_dir,
    pattern = "^fullres_.*_RenalC\\.rds$",
    full.names = TRUE
  )
  
  if (length(rds_files) > 0) {
    for (file in rds_files) {
      # Extract trait/gene ID from filename
      trait_id <- sub("^fullres_(.*)_RenalC\\.rds$", "\\1", basename(file))
      
      # Read RDS file
      res <- readRDS(file)
      
      # Extract coloc results and add identifiers
      coloc_df <- res$coloc %>% dplyr::select(-c(!!sym(paste0("hit_", trait_id))))
      coloc_df$gene <- trait_id
      coloc_df$chromosome <- chr
      
      # Store in list
      all_coloc[[length(all_coloc) + 1]] <- coloc_df
    }
  }
}

# Combine all results into one data frame
combined_coloc <- dplyr::bind_rows(all_coloc)

# Save combined results
output_path <- file.path(respath, "combined_coloc_results.rds")
saveRDS(combined_coloc, output_path)

cat("Combined coloc results saved to:", output_path, "\n")