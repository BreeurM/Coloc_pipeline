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

setwd("~/Code/Coloc_pipeline")
source("coloc_susie_utils.R")


################################################################################ 
# Load and format protein data
################################################################################ 


exp_raw <- read_csv("PanC_coloc/RET_P07949_OID21346_rs2744077.rsids.csv")
# Formatting to ensure that column names are consistent
exp_raw <- format_data(exp_raw,
                           snp_col = "RSID", beta_col = "BETA",
                           se_col = "SE", log_pval = T, pval_col = "LOG10P",
                           effect_allele_col = "ALLELE1",
                           other_allele_col = "ALLELE0",
                           pos_col = "GENPOS", chr_col = "CHROM", 
                           eaf_col = "A1FREQ", min_pval = NA)

# Find lead variant and extract region around it

# Pb: too many lead variants here:
temp = exp_raw$pos.exposure[exp_raw$pval.exposure == min(exp_raw$pval.exposure)]
cat("Potential lead variants span a range of", max(temp) - min(temp), "kb \n")

ggplot(exp_raw, aes(x = pos.exposure, y= -log10(pval.exposure))) + geom_point()

# Set broad manual window 
exp_data <- exp_raw %>% filter(pos.exposure > 4.25e7 & pos.exposure < 4.3e7)



################################################################################ 
# Load and format PanC data
################################################################################ 

out_raw <- fread("PanC_coloc/Panc_Rsid_Annot.txt") 
out_data <- TwoSampleMR::format_data(data.frame(out_raw), 
                                       snps = exp_data$SNP,
                                       type = "outcome",
                                       chr_col = "CHR",
                                       pos_col = "PO",
                                       snp_col = "rsid",
                                       beta_col = "Effect",
                                       se_col = "StdErr",
                                       pval_col = "P-value",
                                       log_pval = FALSE,
                                       eaf_col = "Freq1",
                                       effect_allele_col = "ALLELE1", 
                                       other_allele_col = "ALLELE2")
# Careful, this is a plink outcome, idk for sure whether one should use A1/A2 or
# REF/ALT. To be checked.

rm(out_raw) # To spare memory


################################################################################ 
# Run colocalisation
################################################################################ 

# Run coloc with default parameters
res <- main_coloc(exp_data = exp_data, N_exp = 36000, exp_type = "quant", exp_sd = 1,
               out_data = out_data, N_out = 15000, out_type = "cc", out_sd = .51,
               LD_matrix = NULL, 
               bfile_loc = "N:/EPIC_genetics/1000G_EUR/1000G_EUR/QC_1000G_P3")



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

# No credible set was found for PanC, which is in line with the low p-vals in out_data
min(out_data$pval.outcome)
# [1] 0.01166011

# Tried with smaller coverage (setting out_coverage = .1 when calling main_coloc)
# Still no credible set. Very likely that there is no signal at all there for PanC


res$coloc.res$abf.res$summary
# nsnps    PP.H0.abf    PP.H1.abf    PP.H2.abf    PP.H3.abf    PP.H4.abf 
# 1.815000e+03 0.000000e+00 8.304898e-01 0.000000e+00 6.108796e-02 1.084222e-01 

# Vanilla ABF coloc supports this as well, and finds that the snps are likely only 
# associated with REF 
# H0: no asso
# H1: asso with exp only <-- Highest posterior proba, this hypothesis is supported
# H2: asso with out only
# H3: asso with both, but not colocalized
# H4: asso with both, colocalized



res$flags$exp_alignment_check$alignment_check
# [1] 0.6085741
res$flags$out_alignment_check$alignment_check
# [1] 0.6657118

# Both alignment flags were close calls
# Look at the kriging plots. Suspected misalignments should be flagged in red:
res$flags$exp_alignment_check$kriging$plot
res$flags$out_alignment_check$kriging$plot
# No misalignment.


