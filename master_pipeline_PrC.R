rm(list = ls())


library(renv)
library(httr)
library(glue)
library(jsonlite)
library(ggrepel)
library(ggpubr)
library(readr)
library(tictoc)
# library(GenomicRanges)
# library(rtracklayer)
# library(IRanges)
library(liftOver)
library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)


# Install packages to infer genome build if necessary
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")



setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/scripts/pipeline_utils.R")

# Trait paths

trait_id            <- "MSMB"
trait_lbf_path      <- "data/PrC_coloc/MSMB"
trait_sumstats_path <- "data/PrC_coloc/MSMB/MSMB_P08118_OID20275_rs10993994.rsids.csv"
N_exp               <- 35000
exp_type            <- "quant"
exp_sd              <- 1

# Outcome paths and parameters

out_id              <- "Prostate_cancer"
out_sumstats_path   <- "data/PrC_coloc/ELLIPSE_V2_META_EUROPEAN_Results_012121.txt"
out_lbf_dir_path    <- "data/PrC_coloc/prostate_cancer_lbf"
N_out               <- 178000
out_type            <- "cc"
out_sd              <- 0.48

# Utility paths

# All the .wrapper functions have these paths as default:
# plink_path          <- "plink"
# bfile_path          <- "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered"
# temp_dir_path       <- "Temp"
# They can be changed when calling the .wrappers functions
chain_path  <- "data/hg19ToHg38.over.chain"
rsid38_path <- "data/hg38_rsid_map_file.txt"

# Where to store the results
respath <- "results"

################################################################################
# Extract LBF for kidney study
################################################################################


trait <- fread(trait_sumstats_path)

trait <- format_data(data.frame(trait),
                   chr_col = "CHROM",
                   pos_col = "GENPOS",
                   snp_col = "RSID",
                   beta_col = "BETA",
                   se_col = "SE",
                   pval_col = "LOG10P",
                   log_pval = TRUE,
                   eaf_col = "A1FREQ",
                   effect_allele_col = "ALLELE1",
                   other_allele_col = "ALLELE0"
)

trait <- finemap.wrapper(out = trait, out_lbf_dir_path = trait_lbf_path, 
                N_out = N_exp, out_type = exp_type, out_sd = exp_sd)

################################################################################
# Loading outcome summary stats
# Run finemapping and get Bayes factors if needed
################################################################################


## Loading and formatting

out_raw <- fread(out_sumstats_path)

out_raw <- out_raw %>% 
  filter((Chromosome == unique(trait$chr)))%>% 
  mutate(chr = Chromosome, pos = Position, SNP = SNP_Id)

build <- get_genome_build_local(out_raw, sampled_snps = 100, path_to_38 = rsid38_path)
if (build == "GRCH37") {
  out_raw <- lift_coordinates(out_raw, chain_path)
}

out <- out_raw %>% filter(between(pos, min(trait$pos), max(trait$pos)))

rm(out_raw)

out <- format_data(data.frame(out),
  chr_col = "chr",
  pos_col = "pos",
  snp_col = "SNP_Id",
  beta_col = "Estimate_Effect",
  se_col = "SE",
  pval_col = "P_value",
  log_pval = FALSE,
  eaf_col = "EAF_Control",
  effect_allele_col = "Allele_1",
  other_allele_col = "Allele_2"
)


out <- finemap.wrapper(out, out_lbf_dir_path, N_out, out_type, out_sd, trait)


################################################################################
# Run coloc based on  (SuSiE and vanilla)
################################################################################


res_coloc <- coloc.wrapper(trait, out, trait_id, out_id)


################################################################################
# Investigate causality with MR
################################################################################


res_mr <- mr.wrapper(trait, out, N_out, res_coloc, trait_id, out_id)


################################################################################
# ZZ and locus plot
################################################################################


plots <- plot.wrapper(trait, out, res_coloc, res_mr, trait_id, out_id, out_type)


################################################################################
# Store the results
################################################################################


res <- list(coloc = res_coloc,
            mr = res_mr,
            plots = plots)

resfile <- paste0(respath, "/fullres_", 
                  str_replace(trait_id, " ", "_"), "_", 
                  str_replace(out_id, " ", "_"), ".rds")
saveRDS(res, resfile)


