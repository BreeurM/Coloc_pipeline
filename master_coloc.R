rm(list = ls())

library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)
library(readr)
library(ggplot2)
library(ggrepel)

setwd("~/Code/Coloc_pipeline")
source("coloc_susie_utils.R")


################################################################################ 
# Load and format protein data
################################################################################ 


# exp_raw <- read_csv("PanC_coloc/RET_P07949_OID21346_rs2744077.rsids.csv")
exp_raw <- read_csv("Polycythaemia_vera/PDGFRA_P16234_OID20133_rs35597368.rsids.csv")
# Formatting to ensure that column names are consistent
exp_raw <- format_data(exp_raw,
                           snp_col = "RSID", beta_col = "BETA",
                           se_col = "SE", log_pval = T, pval_col = "LOG10P",
                           effect_allele_col = "ALLELE1",
                           other_allele_col = "ALLELE0",
                           pos_col = "GENPOS", chr_col = "CHROM", 
                           eaf_col = "A1FREQ", min_pval = NA)

# Find lead variant and extract region around it

# Find lead variants
temp = exp_raw$pos.exposure[exp_raw$pval.exposure == min(exp_raw$pval.exposure)]
cat("Potential lead variants span a range of", max(temp) - min(temp), "kb \n")

ggplot(exp_raw, aes(x = pos.exposure, y= -log10(pval.exposure))) + geom_point()

# Set broad manual window 
# exp_data <- exp_raw %>% filter(pos.exposure > 4.25e7 & pos.exposure < 4.3e7)

# Another option if the lead variant is known
lead_var <- "rs35597368"
lead_pos <- exp_raw$pos.exposure[exp_raw$SNP == lead_var]
width <- 250000
exp_data <- exp_raw %>% filter(between(pos.exposure,lead_pos - width, lead_pos + width))

################################################################################ 
# Load and format PanC data
################################################################################ 

out_raw <- fread("Polycythaemia_vera/Polycythaemia_vera.tsv")
# out_data <- TwoSampleMR::format_data(data.frame(out_raw), 
#                                        snps = exp_data$SNP,
#                                        type = "outcome",
#                                        chr_col = "CHR",
#                                        pos_col = "PO",
#                                        snp_col = "rsid",
#                                        beta_col = "Effect",
#                                        se_col = "StdErr",
#                                        pval_col = "P-value",
#                                        log_pval = FALSE,
#                                        eaf_col = "Freq1",
#                                        effect_allele_col = "ALLELE1", 
#                                        other_allele_col = "ALLELE2")
# Careful, this is a plink outcome, idk for sure whether one should use A1/A2 or
# REF/ALT. To be checked.

out_data <- TwoSampleMR::format_data(data.frame(out_raw), 
                                     snps = exp_data$SNP,
                                     type = "outcome",
                                     chr_col = "`#CHR`",
                                     pos_col = "POS",
                                     snp_col = "rsid",
                                     beta_col = "all_inv_var_meta_beta",
                                     se_col = "all_inv_var_meta_sebeta",
                                     pval_col = "all_inv_var_meta_p",
                                     log_pval = FALSE,
                                     eaf_col = "FINNGEN_af_alt",
                                     effect_allele_col = "ALT",
                                     other_allele_col = "REF")

rm(out_raw) # To spare memory

# Remove snps for which eaf is missing

out_data <- out_data %>% filter(!is.na(eaf.outcome))

################################################################################ 
# Run colocalisation
################################################################################ 

# Run coloc with default parameters
res <- main_coloc(exp_data = exp_data, N_exp = 36000, exp_type = "quant", exp_sd = 1,
               out_data = out_data, N_out = 725000, out_type = "cc", out_sd = .002,
               LD_matrix = NULL, plink_loc = "plink", 
               bfile_loc = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
               zz_plot = TRUE, coloc_snp = lead_var)


# To provide custom LD matrix, make sure that the column names are in the format
# "rsID_Allele1_Allele2"
# If not possible, make sure that the alleles and of your exposure data, outcome
# data, and LD matrix are already all aligned.


################################################################################ 
# Examine results
################################################################################ 


res$coloc.res$susie.res
# nsnps
# <lgcl>
#   1:     NA

# Looks like we did not get any result here.

res$finemapping.res$out_susie$sets
# $cs
# NULL
# 
# $coverage
# NULL
# 
# $requested_coverage
# [1] 0.95

# No credible set was found for outcome, which is in line with the high p-vals in out_data
min(out_data$pval.outcome)
# [1] 0.000283

# Tried with smaller coverage (setting out_coverage = .1 when calling main_coloc)
# Still no credible set. Very likely that there is no signal at all there for PanC


res$coloc.res$abf.res$summary
# nsnps     PP.H0.abf     PP.H1.abf     PP.H2.abf     PP.H3.abf     PP.H4.abf 
# 1.241000e+03 5.289982e-279  2.002386e-01 1.543288e-279  5.767510e-02  7.420863e-01

# Yet, vanilla ABF coloc supports colocalisation
# H0: no asso
# H1: asso with exp only <-- Highest posterior proba, this hypothesis is supported
# H2: asso with out only
# H3: asso with both, but not colocalized
# H4: asso with both, colocalized

# SuSiE fails to derive a credible set because the signal is very weak. 


################################################################################ 
# Format results into dataframe
################################################################################ 


df_res <- format_main_coloc_results(res)
#Add two colums to keep track of traits
df_res <- df_res %>% mutate(trait_1 = "P16234", trait_2 = "Polycythaemia_vera")
