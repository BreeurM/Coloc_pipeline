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

setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/pipeline_utils.R")


################################################################################
# List the data available from the eQTL catalogue
################################################################################


# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md
# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/coloc.susie/coloc_susie.md

max_pulled_rows <- 1000 
URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")
r <- GET(URL, accept_json())
cont <- content(r, "text", encoding = "UTF-8")
datasets <- fromJSON(cont)


################################################################################
# Extract LBF for kidney study
################################################################################


kid_genexp_cs <- readr::read_tsv("eQTL_catalogue/QTD000261.credible_sets.tsv.gz", show_col_types = FALSE)
kid_genexp_lbf <- fread("eQTL_catalogue/QTD000261.lbf_variable.txt.gz")

# Merge rsid into kid_genexp_lbf using the variant column and add chr and pos information
kid_genexp_lbf <- kid_genexp_lbf %>%
  left_join(select(kid_genexp_cs, variant, rsid, molecular_trait_id), 
            by = c("variant", "molecular_trait_id"))%>%
  separate(variant, into = c("chr", "bp", "extra1", "extra2"), sep = "_", remove = FALSE) %>%
  mutate(
    chr = str_remove(chr, "chr"),  # Remove "chr" prefix
    chr = ifelse(chr == "X", "23", chr),  # Replace "X" with "23"
    chr = as.numeric(chr),
    bp = as.numeric(bp)
  )

length(table(kid_genexp_lbf$molecular_trait_id))
# 424 genes whose expression has been quantified here.


################################################################################
# Load renal cancer data
################################################################################


out_raw <- fread("N:/EPIC_genetics/Cancer_sumstats/RENAL_multi_ancestry/EURO_ONLY.tsv")
plot_can <- manhattan(out_raw %>% filter(p_value <1e-5), chr="chromosome", bp = "base_pair_location", snp = "rsid", p = "p_value")

can_lead_var <- out_raw[which.min(out_raw$p_value),variant_id]
# [1] "11_69422827_C_T"

## Keep region around variant 

lead_pos <- out_raw[which.min(out_raw$p_value),base_pair_location]
width <- 100000
out_trait_region <- out_raw %>% filter(between(base_pair_location,lead_pos - width, lead_pos + width))
rm(out_raw)


################################################################################
# Run coloc for a chosen trait
################################################################################


## Select gene expr. where lead variant appears in credible set

candidate_traits <- kid_genexp_lbf$molecular_trait_id[grepl(can_lead_var, kid_genexp_lbf$variant)]
# [1] "ENSG00000132740" "ENSG00000250508" "ENSG00000197345"

trait <- "ENSG00000132740"
trait_lfb <- kid_genexp_lbf %>% filter(molecular_trait_id == trait)

## Subset the corresponding region in the cancer outcome data
# specific to eQTL, will have to be remade for the general case
# 
# pos_interval <- str_split_fixed(unique(trait_lfb$region), ":", 2)[,2] 
# # [1] "67903863-69903863"
# pos_low  <- as.numeric(str_split(pos_interval, "-")[[1]][1])
# pos_high <- as.numeric(str_split(pos_interval, "-")[[1]][2])


## Format the cancer data

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
                                     other_allele_col = "other_allele")

## Run the finemapping

out_trait_susie <- finemap_susie(
  exp_data = out_trait_region,
  N_exp = 780000,
  exp_type = "cc",
  exp_sd = 0.034,
  LD_matrix = NULL, 
  plink_loc = "plink", 
  bfile_loc = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
  max_iter = 1000
)












