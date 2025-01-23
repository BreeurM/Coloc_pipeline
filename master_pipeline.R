rm(list = ls())

library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)
library(httr)
library(glue)
library(jsonlite)
library(ggrepel)
library(readr)
library(qqman)
library(tictoc)
library(GenomicRanges)

setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/scripts/pipeline_utils.R")


################################################################################
# List the data available from the eQTL catalogue
################################################################################


# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md
# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/coloc.susie/coloc_susie.md

# max_pulled_rows <- 1000
# URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")
# r <- GET(URL, accept_json())
# cont <- content(r, "text", encoding = "UTF-8")
# datasets <- fromJSON(cont)


################################################################################
# Extract LBF for kidney study
################################################################################


# kid_genexp_cs <- readr::read_tsv("eQTL_catalogue/QTD000261.credible_sets.tsv.gz", show_col_types = FALSE)
kid_genexp_lbf <- fread("data/eQTL_catalogue/QTD000261.lbf_variable.txt.gz")

length(table(kid_genexp_lbf$molecular_trait_id))
# 424 genes whose expression has been quantified here.

# regions <- names(table(kid_genexp_lbf$region))
# grouped_regions <- group_overlapping_regions(regions, tol = 100000)

## Pick a trait

trait <- "ENSG00000250508"
trait_lbf <- kid_genexp_lbf %>% filter(molecular_trait_id == trait)

path_to_out <- "N:/EPIC_genetics/Cancer_sumstats/RENAL_multi_ancestry/EURO_ONLY.tsv"


################################################################################
# Find out region of interest
################################################################################


## Subset the corresponding region in the cancer outcome data

lead_chr <- unique(trait_lbf$chr)
pos_interval <- str_split_fixed(unique(trait_lbf$region), ":", 2)[, 2]
pos_low <- as.numeric(str_split(pos_interval, "-")[[1]][1])
pos_high <- as.numeric(str_split(pos_interval, "-")[[1]][2])


## Query summary statistics in two batches because window is too large

dataset_id <- "QTD000261"

associations1 <- request_associations(dataset_id, pos_low, (pos_low + pos_high) / 2, as.numeric(lead_chr), trait)
associations2 <- request_associations(dataset_id, (pos_low + pos_high) / 2 + 1, pos_high, as.numeric(lead_chr), trait)

sumstats <- rbind(associations1, associations2)
rm(associations1, associations2)


## Merge rsid into trait_lbf using the variant column and add chr and pos information

trait_lbf <- trait_lbf %>%
  left_join(
    select(
      sumstats,
      variant,
      rsid,
      molecular_trait_id,
      pvalue
    ),
    by = c("variant", "molecular_trait_id")
  ) ## This will be the input into main_coloc_bf_susie


## Find out the 500kB region around min pvalue and whether it is covered

lead_pos <- trait_lbf %>%
  filter(pvalue == min(pvalue, na.rm = TRUE)) %>%
  pull(position) %>%
  mean()

required_start <- lead_pos - 500000
required_end <- lead_pos + 500000


################################################################################
# Run finemapping of renal cancer and get Bayes factors if needed
################################################################################


## Find out which regions are already covered

# Set the directory containing the files and list them
lbf_directory <- "data/renal_cancer_lbf"
lbf_file_list <- list.files(lbf_directory, full.names = FALSE)

# Extract chromosome and positions
coverage <- as.data.frame(
  do.call(rbind, lapply(lbf_file_list, function(file) {
    match <- regmatches(file, regexec("lbf_chr(\\d+)_([0-9]+)-([0-9]+)", file))
    if (length(match[[1]]) == 4) {
      return(data.frame(
        chromosome = as.integer(match[[1]][2]),
        start_position = as.integer(match[[1]][3]),
        end_position = as.integer(match[[1]][4])
      ))
    } else {
      return(NULL)
    }
  }))
)

# See if required region is covered
filtered_coverage <- subset(coverage, chromosome == lead_chr)
position_covered <- filtered_coverage[(filtered_coverage$start_position <= required_start) &
  (filtered_coverage$end_position >= required_end), ]



if (nrow(position_covered) > 0) {
   
  message(paste0("Region already fine-mapped. Laoding the BFs from ", lbf_directory))
  out_lbf <- readRDS(paste0(
    lbf_directory, "/lbf_chr",
    as.character(position_covered$chromosome),
    "_",
    as.character(position_covered$start_position),
    "-",
    as.character(position_covered$end_position),
    ".rds"
  ))
} else {
  # finemapping needs to be ran
  
  message(paste0("Region not covered, fine-mapping region", pos_interval, " on chr ", lead_chr))

  out_raw <- fread(path_to_out)
  # plot_can <- manhattan(out_raw %>% filter(p_value < 1e-5), chr = "chromosome", bp = "base_pair_location", snp = "rsid", p = "p_value")

  # can_lead_var <- out_raw[which.min(out_raw$p_value), variant_id]
  # # [1] "11_69422827_C_T"
  #

  out_trait_region <- out_raw %>%
    filter(chromosome == lead_chr) %>%
    filter(between(base_pair_location, pos_low, pos_high))


  out_trait_region <- TwoSampleMR::format_data(data.frame(out_trait_region),
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    pval_col = "p_value",
    log_pval = FALSE,
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )

  ## Find out if anything passes the conventional significance threshold

  if (any(out_trait_region$pval.exposure < 5e-2)) { 
    ## Run the finemapping

    tic("Running SuSiE")
    out_trait_susie <- finemap_susie(
      exp_data = out_trait_region,
      N_exp = 780000,
      exp_type = "cc",
      exp_sd = 0.034,
      LD_matrix = NULL,
      plink_loc = "plink",
      bfile_loc = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
      max_iter = 1000,
      repeat_until_convergence = TRUE,
      run_checks = FALSE
    )
    toc()

    out_lbf <- as.data.frame(t(out_trait_susie$lbf_variable)) %>%
      setNames(paste0("lbf_variable", 1:10)) %>%
      rownames_to_column("SNP") %>%
      left_join(out_trait_region %>%
        transmute(variant = paste0("chr", chr.exposure, "_", pos.exposure), SNP), by = "SNP")

    saveRDS(out_lbf, file = paste0(lbf_directory, "/lbf_chr", lead_chr, "_", pos_interval, ".rds"))
  }
}

################################################################################
# Run coloc based on BF
################################################################################


trait_mat <- as.matrix(dplyr::select(trait_lbf, lbf_variable1:lbf_variable10))
row.names(trait_mat) <- sub("(_[ACGT]+_[ACGT]+)$", "", trait_lbf$variant)
trait_mat <- t(trait_mat)


out_mat <- as.matrix(dplyr::select(out_lbf, lbf_variable1:lbf_variable10))
row.names(out_mat) <- out_lbf$variant
out_mat <- t(out_mat)


res_coloc <- coloc::coloc.bf_bf(trait_mat, out_mat)
temp <- res_coloc$summary %>% filter(PP.H4.abf > 0.5)
if(nrow(temp)>0){
  print(temp)
}else{
  print("No PP.H4 was greater than 0.5")
}


