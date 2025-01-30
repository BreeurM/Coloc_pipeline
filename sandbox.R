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

path_lbf      <- "data/eQTL_catalogue/QTD000261.lbf_variable.txt.gz"
path_sumstats <- "data/eQTL_catalogue/QTD000261.cc.tsv.gz"

kid_genexp_lbf <- fread(path_lbf)
traits <- names(table(kid_genexp_lbf$molecular_trait_id))
length(traits)
# 424 genes whose expression has been quantified here.

trait_id <- "ENSG00000184227"
trait <- format_eqtl_cat_trait(trait_id, path_lbf, path_sumstats)


################################################################################
# Run finemapping of renal cancer and get Bayes factors if needed
################################################################################


path_to_out <- "N:/EPIC_genetics/Cancer_sumstats/RENAL_multi_ancestry/EURO_ONLY.tsv"

## Find out which regions are already covered
region <- unique(trait$region)
lead_pos <- trait %>%
  filter(pval == min(pvalue, na.rm = TRUE)) %>%
  pull(pos) %>% mean
# Set the directory containing the files
lbf_directory <- "data/renal_cancer_lbf"

region_utils <- region_utils(region, lead_pos, lbf_directory)


if (region_utils$is_pos_covered) {
  message(paste0("Region already fine-mapped. Laoding the BFs from ", lbf_directory))
  out_lbf <- readRDS(paste0(
    lbf_directory, "/lbf_chr",
    as.character(region_utils$position_covered$chromosome),
    "_",
    as.character(region_utils$position_covered$start_position),
    "-",
    as.character(region_utils$position_covered$end_position),
    ".rds"
  ))
  
  
  message("Loading outcome summary stats.")

  out_raw <- fread(path_to_out)

  out <- out_raw %>%
    filter(SNP %in% out_lbf$SNP)

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

  # ggplot(out %>% filter(p_value < 1e-5), aes(x = pos, y= -log10(pval))) + geom_point()

  ## Find out if anything passes the conventional significance threshold

  message(paste0("Region not covered, fine-mapping region ", pos_interval, " on chr ", lead_chr))

  if (any(out$pval < 5e-2)) {
    ## Run the finemapping

    tic("Running SuSiE")
    out_trait_susie <- tryCatch(expr = {
      finemap_susie(
        exp_data = out,
        N_exp = 780000,
        exp_type = "cc",
        exp_sd = 0.034,
        LD_matrix = NULL,
        plink_loc = plink_loc,
        bfile_loc = bfile_loc,
        max_iter = 500,
        repeat_until_convergence = FALSE,
        run_checks = FALSE
      )
    }, error = function(e) {
      warning("SuSiE did not converge, returning BF = 0")
      NULL
    })
    toc()

    if (!is.null(out_trait_susie)) {
      out_lbf <- as.data.frame(t(out_trait_susie$lbf_variable)) %>%
        setNames(paste0("lbf_variable", 1:10)) %>%
        rownames_to_column("SNP") %>%
        left_join(out %>%
          transmute(variant = paste0("chr", chr, "_", pos), SNP), by = "SNP")
    } else {
      out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
        setNames(paste0("lbf_variable", 1:10)) %>%
        mutate(
          variant = paste0("chr", out$chr, "_", out$pos),
          SNP = out$SNP
        )
    }

    saveRDS(out_lbf, file = paste0(lbf_directory, "/lbf_chr", lead_chr, "_", pos_interval, ".rds"))
  } else {
    message("No variant crosses the significance threshold. Returning BF = 0")

    out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
      setNames(paste0("lbf_variable", 1:10)) %>%
      mutate(
        variant = paste0("chr", out$chr, "_", out$pos),
        SNP = out$SNP
      )
  }
}

out <- merge(out, out_lbf)

rm(position_covered)

################################################################################
# Run coloc based on BF
################################################################################


trait_mat <- as.matrix(dplyr::select(trait, lbf_variable1:lbf_variable10))
row.names(trait_mat) <- sub("(_[ACGT]+_[ACGT]+)$", "", trait$variant)
trait_mat <- t(trait_mat)
# 
# out_mat <- as.matrix(dplyr::select(out, lbf_variable1:lbf_variable10))
# row.names(out_mat) <- out$variant
# out_mat <- t(out_mat)

out_mat <- as.matrix(dplyr::select(out, lbf))
row.names(out_mat) <- out$variant
out_mat <- t(out_mat)

res_coloc <- coloc::coloc.bf_bf(trait_mat, out_mat)
coloc_snps <- res_coloc$summary %>% filter(PP.H4.abf > 0.5)
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
res_abf_coloc$summary

rm(trait_mat, out_mat)

################################################################################
# ZZ and locus plot
################################################################################


## Query LD matrix for desired region

SNP_list <- trait$SNP[between(trait$pos, required_start, required_end)]

LD_matrix <- get_ld_matrix_from_bim(SNP_list,
  plink_loc = plink_loc,
  bfile_loc = bfile_loc,
  with_alleles = T,
  temp_dir_path = temp_dir_path
)
rm(SNP_list)

## Extract lead snp from trait

lead_snp <- unique(trait$SNP[which.min(abs(trait$pos - lead_pos))])

## Flip z scores if needed, according to LD_mat

harm_dat <- merge(
  align_to_LD(trait, LD_matrix),
  align_to_LD(out, LD_matrix)
)


## Is there a coloc snp to highlight in the zz plot?

if (!is.na(res_abf_coloc$summary$PP.H4.abf) & res_abf_coloc$summary$PP.H4.abf > 0.5) {
  # Assign a value to coloc_snp, leave it NULL otherwise
}

## Make locus and zz plots, store them in list

plot_list <- locus_plot(as.data.frame(LD_matrix),
  as.data.frame(harm_dat),
  lead_SNP = lead_snp
)

plot_list$zzplot <- zz_plot(as.data.frame(LD_matrix), lead_snp, harm_dat)


################################################################################
# Investigate causality with MR
################################################################################


trait_mr <- TwoSampleMR::format_data(as.data.frame(trait)) %>%
  filter(pval.exposure < 1e-5)
out_mr <- TwoSampleMR::format_data(as.data.frame(out),
  snps = trait_mr$SNP,
  type = "outcome"
)
out_mr$samplesize.outcome <- 780000

harm_dat_mr <- TwoSampleMR::harmonise_data(trait_mr, out_mr)

res_mr_single <- mr_singlesnp(harm_dat_mr)
res_mr <- mr(harm_dat_mr, method_list = c("mr_wald_ratio", "mr_ivw", "mr_weighted_median"))
