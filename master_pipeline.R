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


plink_loc <- "plink"
bfile_loc <- "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered"

temp_dir_path <- "Temp"

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
sumstats <- fread("data/eQTL_catalogue/QTD000261.cc.tsv.gz")

length(table(kid_genexp_lbf$molecular_trait_id))
# 424 genes whose expression has been quantified here.

# regions <- names(table(kid_genexp_lbf$region))
# grouped_regions <- group_overlapping_regions(regions, tol = 100000)

## Pick a trait

trait <- "ENSG00000250508"
trait_data <- kid_genexp_lbf %>% filter(molecular_trait_id == trait)

region <- unique(trait_data$region)

path_to_out <- "N:/EPIC_genetics/Cancer_sumstats/RENAL_multi_ancestry/EURO_ONLY.tsv"


################################################################################
# Format and find out region of interest
################################################################################


## Merge rsid into trait using the variant column and add chr and pos information

trait_data <- trait_data %>%
  left_join(
    select(
      sumstats,
      variant,
      rsid,
      molecular_trait_id,
      pvalue,
      beta,
      se,
      ref,
      alt,
      maf
    ),
    by = c("variant", "molecular_trait_id")
  ) 

rm(sumstats, kid_genexp_lbf)

trait <- format_data(trait_data, snp_col = "rsid", eaf_col = "maf",
            effect_allele_col = "ref", other_allele_col = "alt",
            pval_col = "pvalue", chr_col = "chromosome",
            pos_col = "position") ## This will be the input into main_coloc_bf_susie


## Find out the 500kB region around min pvalue and whether it is covered

lead_chr <- unique(trait$chr)
pos_interval <- str_split_fixed(region, ":", 2)[, 2]
pos_low <- as.numeric(str_split(pos_interval, "-")[[1]][1])
pos_high <- as.numeric(str_split(pos_interval, "-")[[1]][2])

lead_pos <- trait %>%
  filter(pval == min(pvalue, na.rm = TRUE)) %>%
  pull(pos) %>%
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
  message("Loading outcome summary stats.")

  out_raw <- fread(path_to_out)

  out <- out_raw %>%
    filter(chromosome == lead_chr) %>%
    filter(between(base_pair_location, pos_low, pos_high))

  rm(out_raw)

  out <- format_data(data.frame(out),
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

  message("Loading outcome summary stats.")

  out_raw <- fread(path_to_out)
  # plot_can <- manhattan(out_raw %>% filter(p_value < 1e-5), chr = "chromosome", bp = "base_pair_location", snp = "rsid", p = "p_value")

  # can_lead_var <- out_raw[which.min(out_raw$p_value), variant_id]
  # # [1] "11_69422827_C_T"
  #

  out <- out_raw %>%
    filter(chromosome == lead_chr) %>%
    filter(between(base_pair_location, pos_low, pos_high))

  rm(out_raw)

  out <- format_data(data.frame(out),
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

  message(paste0("Region not covered, fine-mapping region ", pos_interval, " on chr ", lead_chr))

  if (any(out$pval < 5e-2)) {
    ## Run the finemapping

    tic("Running SuSiE")
    out_trait_susie <- finemap_susie(
      exp_data = out,
      N_exp = 780000,
      exp_type = "cc",
      exp_sd = 0.034,
      LD_matrix = NULL,
      plink_loc = plink_loc,
      bfile_loc = bfile_loc,
      max_iter = 1000,
      repeat_until_convergence = TRUE,
      run_checks = FALSE
    )
    toc()

    out_lbf <- as.data.frame(t(out_trait_susie$lbf_variable)) %>%
      setNames(paste0("lbf_variable", 1:10)) %>%
      rownames_to_column("SNP") %>%
      left_join(out %>%
        transmute(variant = paste0("chr", chr, "_", pos), SNP), by = "SNP")

    saveRDS(out_lbf, file = paste0(lbf_directory, "/lbf_chr", lead_chr, "_", pos_interval, ".rds"))
  } else {
    message("No variant crosses the significant threshold. Returning BF = 0")

    out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
      setNames(paste0("lbf_variable", 1:10)) %>%
      mutate(
        variant = paste0("chr", out$chr, "_", out$pos),
        SNP = out$SNP
      )
  }
}

################################################################################
# Run coloc based on BF
################################################################################


trait_mat <- as.matrix(dplyr::select(trait, lbf_variable1:lbf_variable10))
row.names(trait_mat) <- sub("(_[ACGT]+_[ACGT]+)$", "", trait$variant)
trait_mat <- t(trait_mat)


out_mat <- as.matrix(dplyr::select(out_lbf, lbf_variable1:lbf_variable10))
row.names(out_mat) <- out_lbf$variant
out_mat <- t(out_mat)


res_coloc <- coloc::coloc.bf_bf(trait_mat, out_mat)
temp <- res_coloc$summary %>% filter(PP.H4.abf > 0.5)
if (nrow(temp) > 0) {
  print(temp)
} else {
  print("No PP.H4 was greater than 0.5")
}

rm(trait_mat, out_mat, temp)


################################################################################
# Vanilla coloc from the summary stats
################################################################################


# Flip alleles for exp and out if not consistent.
# Hold on, needed? lbfs won't change...
# Will be needed for the MR and maybe the plots but we can opt out for now


# Run coloc based on lbfs computed in the format data function

trait_mat <- as.matrix(dplyr::select(trait, lbf)) # %>%
  # filter(!is.na(lbf)))
row.names(trait_mat) <- sub("(_[ACGT]+_[ACGT]+)$", "", trait$variant[!is.na(trait$lbf)])
trait_mat <- t(trait_mat)


out_mat <- as.matrix(dplyr::select(out, lbf))
row.names(out_mat) <- out$variant
out_mat <- t(out_mat)

res_abf_coloc <- coloc::coloc.bf_bf(trait_mat, out_mat)
res_abf_coloc$summary$PP.H4.abf > 0.5

rm(trait_mat, out_mat)

################################################################################
# ZZ plot
################################################################################


## Query LD matrix for desired region

SNP_list <- trait$SNP[between(trait$pos, required_start, required_end)]

LD_matrix <- get_ld_matrix_from_bim(SNP_list,
  plink_loc = plink_loc,
  bfile_loc = bfile_loc,
  with_alleles = T,
  temp_dir_path = temp_dir_path
)

## Extract lead snp from trait

lead_snp <- unique(trait$SNP[trait$pos == lead_pos])

## Flip z scores if needed, according to LD_mat

harm_trait <- align_to_LD(trait, LD_matrix)
harm_out   <- align_to_LD(out, LD_matrix)

harm_dat <- merge(harm_trait, harm_out)

## Into zz_plot fn

# Is there a coloc snp?

if(res_abf_coloc$summary$PP.H4.abf > 0.5){
  # Assign a value to coloc_snp, leave it NULL otherwise
}


zz_plot(as.data.frame(LD_matrix), lead_snp, harm_dat)


