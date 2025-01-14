rm(list = ls())

library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)

setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/coloc_susie_utils.R")

########## Digging deeper into the runsusie outputs


sum_stat_path <- "N:/EPIC_genetics/Annotated_windows/Annotated_windows/"
file_list <- list.files(sum_stat_path, full.names = T, recursive = T)

#### Protein first

prot <- "MSMB"
restricted_file_list <- file_list[grepl(prot, file_list)]
file_name <- "N:/EPIC_genetics/Annotated_windows/Annotated_windows/MSMB_P08118_OID20275/MSMB_P08118_OID20275_rs10993994.rsids.csv"
  exp_raw <- read_csv(file_name)
# Formatting to ensure that column names are consistent
exp_raw <- format_data(exp_raw,
  snp_col = "RSID", beta_col = "BETA",
  se_col = "SE", log_pval = T, pval_col = "LOG10P",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pos_col = "GENPOS", chr_col = "CHROM", eaf_col = "A1FREQ", min_pval = NA
)
# Find lead variant and extract region around it

# Find lead variants
temp <- exp_raw$pos.exposure[exp_raw$pval.exposure == min(exp_raw$pval.exposure)]
cat("Potential lead variants span a range of", max(temp) - min(temp), "kb \n")

ggplot(exp_raw, aes(x = pos.exposure, y = -log10(pval.exposure))) +
  geom_point()

# Set broad manual window
# exp_data <- exp_raw %>% filter(pos.exposure > 4.25e7 & pos.exposure < 4.3e7)

# Another option if the lead variant is known
lead_var <- "rs10993994"
lead_pos <- exp_raw$pos.exposure[exp_raw$SNP == lead_var]
width <- 250000
exp_data <- exp_raw %>% filter(between(pos.exposure, lead_pos - width, lead_pos + width))


exp_susie <- finemap_susie(
  exp_data = exp_data,
  N_exp = 34000,
  exp_type = "quant",
  LD_matrix = NULL, 
  plink_loc = "plink", 
  bfile_loc = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
  max_iter = 1000
)

###### Cancer second

rm(list = setdiff(ls(), "exp_susie"))
# to make sure the two parts are independent

library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)

setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/coloc_susie_utils.R")


out_raw <- fread("N:/EPIC_genetics/Cancer_sumstats/ELLIPSE_V2_META_EUROPEAN_Prostate_012121.txt")
out_raw <- out_raw %>% filter(!is.na(SNP_Id))

lead_var <- "rs10993994"
lead_pos <- out_raw$Position[out_raw$SNP_Id == lead_var]
width <- 25000
out_raw <- out_raw %>% filter(between(Position, lead_pos - width, lead_pos + width))

sum(out_raw$EA != out_raw$Allele_1)
# [1] 0

out_data <- TwoSampleMR::format_data(data.frame(out_raw),
                                     chr_col = "Chromosome",
                                     pos_col = "Position",
                                     snp_col = "SNP_Id",
                                     beta_col = "Estimate_Effect",
                                     se_col = "SE",
                                     pval_col = "P_value",
                                     log_pval = FALSE,
                                     eaf_col = "EAF_Control",
                                     effect_allele_col = "Allele_1",
                                     other_allele_col = "Allele_2")


# Remove snps for which eaf is missing

out_data <- out_data %>% filter(!is.na(eaf.exposure))


out_susie <- finemap_susie(
  exp_data = out_data,
  N_exp = 180000,
  exp_type = "cc",
  exp_sd = 0.48,
  LD_matrix = NULL, 
  plink_loc = "plink", 
  bfile_loc = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
  max_iter = 1000
)
