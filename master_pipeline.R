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

# Trait paths

trait_id            <- "ENSG00000184227" 
trait_lbf_path      <- "data/eQTL_catalogue/QTD000261.lbf_variable.txt.gz"
trait_sumstats_path <- "data/eQTL_catalogue/QTD000261.cc.tsv.gz"

# Outcome paths and parameters

out_sumstats_path   <- "N:/EPIC_genetics/Cancer_sumstats/RENAL_multi_ancestry/EURO_ONLY.tsv"
out_lbf_dir_path    <- "data/renal_cancer_lbf"
N_out               <- 780000
out_type            <- "cc"
out_sd              <- 0.034

# Utility paths

# All the .wrapper functions have these paths as default:
# plink_path          <- "plink"
# bfile_path          <- "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered"
# temp_dir_path       <- "Temp"
# They can be changed when calling the .wrappers functions

################################################################################
# Extract LBF for kidney study
################################################################################


trait <- format_eqtl_cat_trait(trait_id, trait_lbf_path, trait_sumstats_path)


################################################################################
# Loading outcome summary stats
# Run finemapping and get Bayes factors if needed
################################################################################


## Loading and formatting

out_raw <- fread(out_sumstats_path)

out <- out_raw %>% filter((chromosome == unique(trait$chr)) &
  (between(base_pair_location, min(trait$pos), max(trait$pos))))

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


out <- finemap.wrapper(out, out_lbf_dir_path, N_out, out_type, out_sd,trait)


################################################################################
# Run coloc based on  (SuSiE and vanilla)
################################################################################


res_coloc <- coloc.wrapper(out, trait)


################################################################################
# ZZ and locus plot
################################################################################


plot_list <- plot.wrapper(trait, out, res_coloc)


################################################################################
# Investigate causality with MR
################################################################################


res_mr <- mr.wrapper(trait, out, N_out, res_coloc)






