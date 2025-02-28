rm(list = ls())


library(renv)
library(httr)
library(glue)
library(jsonlite)
library(ggrepel)
library(ggpubr)
library(readr)
library(tictoc)
library(liftOver)
library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)


################################################################################
# Set parameters and paths
################################################################################


# Main arguments ---------------------------------------------------------------

workdir <- "~/Code/Coloc_pipeline"
path_to_utils <- "scripts/pipeline_utils.R"
plink_path <- "plink"
bfile_path <- "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered"
chain_path <- "data/gen_utils/hg19ToHg38.over.chain"
rsid38_path <- "data/gen_utils/hg38_rsid_map_file.txt"
temp_dir_path <- "Temp"
respath <- "results"

# Trait parameters -------------------------------------------------------------

trait_id <- "MSMB"
trait_lbf_path <- "data/PrC_coloc/MSMB"
trait_sumstats_path <- "data/PrC_coloc/MSMB/MSMB_P08118_OID20275_rs10993994.rsids.csv"
trait_build <- 38 # or 37, or NA
N_trait <- 35000
trait_type <- "quant"
trait_sd <- 1
# Trait column mappings
trait_chr_col <- "CHROM" # Trait chromosome column name
trait_pos_col <- "GENPOS" # Trait position column name
trait_snp_col <- "RSID" # Trait SNP ID column name
trait_beta_col <- "BETA" # Trait beta column name
trait_se_col <- "SE" # Trait standard error column name
trait_pval_col <- "LOG10P" # Trait p-value column name
trait_log_pval <- TRUE # Trait p-values are log-transformed (TRUE/FALSE)
trait_eaf_col <- "A1FREQ" # Trait effect allele frequency column name
trait_effect_allele_col <- "ALLELE1" # Trait effect allele column name
trait_other_allele_col <- "ALLELE0" # Trait other allele column name

# Outcome parameters -----------------------------------------------------------

out_id <- "Prostate_cancer"
out_sumstats_path <- "data/PrC_coloc/ELLIPSE_V2_META_EUROPEAN_Results_012121.txt"
out_lbf_path <- "data/PrC_coloc/prostate_cancer_lbf"
out_build <- 37 # or 38, or NA
N_out <- 178000
out_type <- "cc"
out_sd <- 0.48
# Outcome column mappings
out_chr_col <- "Chromosome" # Out chromosome column name
out_pos_col <- "Position" # Out position column name
out_snp_col <- "SNP_Id" # Out SNP ID column name
out_beta_col <- "Estimate_Effect" # Out beta column name
out_se_col <- "SE" # Out standard error column name
out_pval_col <- "P_value" # Out p-value column name
out_log_pval <- FALSE # Out p-values are log-transformed (TRUE/FALSE)
out_eaf_col <- "EAF_Control" # Out effect allele frequency column name
out_effect_allele_col <- "Allele_1" # Out effect allele column name
out_other_allele_col <- "Allele_2" # Out other allele column name


################################################################################
# Check parameter integrity and create folders if needed
################################################################################

setwd(workdir)
source(path_to_utils)

if (!file.exists(rsid38_path) & (is.na(trait_build) | is.na(out_build))) {
  warning("No hg38 reference, cannot check genome builds.
          Assuming trait and outcome are in hg38.")
  trait_build <- out_build <- 38
}

if (!file.exists(chain_path)) {
  warning("No chain file found, assuming trait and outcome are in hg38.")
  trait_build <- out_build <- 38
}


if (!dir.exists(trait_lbf_path)) {
  dir.create(trait_lbf_path, recursive = TRUE)
}

if (!dir.exists(out_lbf_path)) {
  dir.create(out_lbf_path, recursive = TRUE)
}

if (!dir.exists(respath)) {
  dir.create(respath, recursive = TRUE)
}


################################################################################
# Extract LBF for trait
################################################################################

cat("Processing trait", trait_id, "\n")
cat("\n")


trait <- fread(trait_sumstats_path)

if (is.na(trait_build)) {
  build <- get_genome_build_local(trait,
    sampled_snps = 100, path_to_38 = rsid38_path,
    snp_col = trait_snp_col,
    chr_col = trait_chr_col,
    pos_col = trait_pos_col
  )
  trait_build <- ifelse(build == "GRCH37", 37, 38)
}

if (trait_build == 37) {
  trait <- lift_coordinates(trait, chain_path,
    snp_col = trait_snp_col,
    chr_col = trait_chr_col,
    pos_col = trait_pos_col
  )
  cat("Trait lifted to hg38.\n")
}


trait <- format_data(data.frame(trait),
  chr_col = trait_chr_col,
  pos_col = trait_pos_col,
  snp_col = trait_snp_col,
  beta_col = trait_beta_col,
  se_col = trait_se_col,
  pval_col = trait_pval_col,
  log_pval = as.logical(trait_log_pval),
  eaf_col = trait_eaf_col,
  effect_allele_col = trait_effect_allele_col,
  other_allele_col = trait_other_allele_col
)


trait <- finemap.wrapper(
  out = trait, out_lbf_dir_path = trait_lbf_path,
  N_out = N_trait, out_type = trait_type, out_sd = trait_sd,
  plink_path = plink_path,
  bfile_path = bfile_path,
  temp_dir_path = temp_dir_path
)

cat("Trait processed.\n")
cat("\n")


################################################################################
# Process outcome data
################################################################################


cat("Processing outcome", out_id, "\n")
cat("\n")


out_raw <- fread(out_sumstats_path)


if (is.na(out_build)) {
  out <- out_raw %>%
    filter(!!sym(out_chr_col) == unique(trait$chr))
  build <- get_genome_build_local(out,
    sampled_snps = 100, path_to_38 = rsid38_path,
    snp_col = out_snp_col,
    chr_col = out_chr_col,
    pos_col = out_pos_col
  )
  out_build <- ifelse(build == "GRCH37", 37, 38)
}

if (out_build == 37) {
  out <- out_raw %>%
    filter(!!sym(out_chr_col) == unique(trait$chr))
  
  out <- lift_coordinates(out, chain_path,
    snp_col = out_snp_col,
    chr_col = out_chr_col,
    pos_col = out_pos_col
  )
  cat("Outcome lifted to hg38.\n")
  
  out <- out %>% filter(between(!!sym(out_pos_col), min(trait$pos), max(trait$pos)))
  
}else{
  out <- out_raw %>% filter((!!sym(out_chr_col) == unique(trait$chr)) & 
                                  (between(!!sym(out_pos_col), min(trait$pos), max(trait$pos))))
  }


rm(out_raw)

out <- format_data(data.frame(out),
  chr_col = out_chr_col,
  pos_col = out_pos_col,
  snp_col = out_snp_col,
  beta_col = out_beta_col,
  se_col = out_se_col,
  pval_col = out_pval_col,
  log_pval = as.logical(out_log_pval),
  eaf_col = out_eaf_col,
  effect_allele_col = out_effect_allele_col,
  other_allele_col = out_other_allele_col
)


out <- finemap.wrapper(out, out_lbf_path, N_out, out_type, out_sd, trait,
  plink_path = plink_path,
  bfile_path = bfile_path,
  temp_dir_path = temp_dir_path
)

cat("Outcome processed.\n")
cat("\n")


################################################################################
# Run coloc based on  (SuSiE and vanilla)
################################################################################


res_coloc <- coloc.wrapper(trait, out, trait_id, out_id)

cat("Colocalisation performed.\n")
cat("\n")


################################################################################
# Investigate causality with MR
################################################################################


res_mr <- mr.wrapper(trait, out, N_out, res_coloc, trait_id, out_id)


################################################################################
# Store the results
################################################################################


res <- list(trait = trait,
            out = out, 
            coloc = res_coloc,
            mr = res_mr)

resfile <- paste0(respath, "/fullres_", 
                  str_replace(trait_id, " ", "_"), "_", 
                  str_replace(out_id, " ", "_"), ".rds")
saveRDS(res, resfile)

cat("Results written to", resfile, "\n")


################################################################################
# ZZ and locus plot
################################################################################


plots <- plot.wrapper(trait, out, res_coloc, res_mr, 
                      trait_id, out_id, out_type,
                      plink_path = plink_path,
                      bfile_path = bfile_path,
                      temp_dir_path = temp_dir_path)
plotfile <- paste0(respath, "/plots_", 
                   str_replace(trait_id, " ", "_"), "_", 
                   str_replace(out_id, " ", "_"), ".pdf")
ggsave(plotfile, plots, dpi = 300, width = 15, height = 10)


cat("Plots written to", plotfile, "\n")
