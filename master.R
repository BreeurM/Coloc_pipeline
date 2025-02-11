library(argparse)

# Create parser
parser <- ArgumentParser(description = 'Coloc analysis pipeline with column mappings')

# Main arguments ---------------------------------------------------------------
parser$add_argument('--workdir', type='character', default='/gpfs3/well/travis-prostate/projects/Coloc_pipeline',
                    help='Working directory')
parser$add_argument('--plink_path', type='character', default='/well/travis-prostate/projects/plink',
                    help='Path to PLINK executable')
parser$add_argument('--bfile_path', type='character', default='data/LD_mat/LD_REF_DAT_MAF_MAC_Filtered',
                    help='Path to LD reference files')
parser$add_argument('--respath', type='character', default='results',
                    help='Path to store results')
parser$add_argument('--chain_path', type='character', default='data/hg19ToHg38.over.chain',
                    help='Path to chain file for liftover')
parser$add_argument('--rsid38_path', type='character', default='data/hg38_rsid_map_file.txt',
                    help='Path to GRCh38 rsID map file')

# Trait parameters -------------------------------------------------------------
parser$add_argument('--trait_id', type='character',
                    help='Trait identifier')
parser$add_argument('--is_trait_eqtl_cat', type='logical', default=FALSE,
                    help='Is trait data downloaded from the eQTL catalogue')
parser$add_argument('--trait_sumstats_path', type='character', 
                    help='Path to trait summary statistics')
parser$add_argument('--trait_lbf_path', type='character', 
                    help='Path to trait LBF directory')
parser$add_argument('--N_trait', type='integer',
                    help='Sample size for exposure')
parser$add_argument('--trait_type', type='character', default='quant',
                    help='Exposure type (quant/cc)')
parser$add_argument('--trait_sd', type='double', default=1.0,
                    help='Exposure standard deviation')

# Trait column mappings --------------------------------------------------------
parser$add_argument('--trait_chr_col', type='character', default='chr',
                    help='Trait chromosome column name')
parser$add_argument('--trait_pos_col', type='character', default='pos',
                    help='Trait position column name')
parser$add_argument('--trait_snp_col', type='character', default='SNP',
                    help='Trait SNP ID column name')
parser$add_argument('--trait_beta_col', type='character', default='beta',
                    help='Trait beta column name')
parser$add_argument('--trait_se_col', type='character', default='se',
                    help='Trait standard error column name')
parser$add_argument('--trait_pval_col', type='character', default='pval',
                    help='Trait p-value column name')
parser$add_argument('--trait_log_pval', type='character', default='FALSE',
                    help='Trait p-values are log-transformed (TRUE/FALSE)')
parser$add_argument('--trait_eaf_col', type='character', default='eaf',
                    help='Trait effect allele frequency column name')
parser$add_argument('--trait_effect_allele_col', type='character', default='effect_allele',
                    help='Trait effect allele column name')
parser$add_argument('--trait_other_allele_col', type='character', default='other_allele',
                    help='Trait other allele column name')

# Outcome parameters -----------------------------------------------------------
parser$add_argument('--out_id', type='character',
                    help='Outcome identifier')
parser$add_argument('--is_out_eqtl_cat', type='logical', default=FALSE,
                    help='Is outcome data downloaded from the eQTL catalogue')
parser$add_argument('--out_sumstats_path', type='character', 
                    help='Path to outcome summary statistics')
parser$add_argument('--out_lbf_dir_path', type='character', 
                    help='Path to outcome LBF directory')
parser$add_argument('--N_out', type='integer',
                    help='Sample size for outcome')
parser$add_argument('--out_type', type='character', default='quant',
                    help='Outcome type (quant/cc)')
parser$add_argument('--out_sd', type='double', default=1.0,
                    help='Outcome standard deviation')

# Outcome column mappings ------------------------------------------------------
parser$add_argument('--out_chr_col', type='character', default='chr',
                    help='Outcome formatted chromosome column name')
parser$add_argument('--out_pos_col', type='character', default='pos',
                    help='Outcome formatted position column name')
parser$add_argument('--out_snp_col', type='character', default='SNP',
                    help='Outcome formatted SNP ID column name')
parser$add_argument('--out_beta_col', type='character', default='beta',
                    help='Outcome beta column name')
parser$add_argument('--out_se_col', type='character', default='se',
                    help='Outcome standard error column name')
parser$add_argument('--out_pval_col', type='character', default='pval',
                    help='Outcome p-value column name')
parser$add_argument('--out_log_pval', type='character', default='FALSE',
                    help='Outcome p-values are log-transformed (TRUE/FALSE)')
parser$add_argument('--out_eaf_col', type='character', default='eaf',
                    help='Outcome effect allele frequency column name')
parser$add_argument('--out_effect_allele_col', type='character', default='effect_allele',
                    help='Outcome effect allele column name')
parser$add_argument('--out_other_allele_col', type='character', default='other_allele',
                    help='Outcome other allele column name')


# Parse arguments
args <- parser$parse_args()

# Set working directory first
setwd(args$workdir)

# Rest of script remains mostly the same with variable substitutions
.libPaths("/well/travis-prostate/users/cbe235/miniconda3/envs/renvironment/lib/R/library")

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

source("scripts/pipeline_utils.R")


################################################################################
# Extract LBF for trait
################################################################################

if (args$is_trait_eqtl_cat){
  trait <- format_eqtl_cat_trait(args$trait_id, args$trait_lbf_path, args$trait_sumstats_path)
}else{
  trait <- fread(args$trait_sumstats_path)
  
  build <- get_genome_build_local(trait, sampled_snps = 100, path_to_38 = args$rsid38_path,
                                  snp_col = args$trait_snp_col,
                                  chr_col = args$trait_chr_col,
                                  pos_col = args$trait_pos_col)
  if (build == "GRCH37") {
    trait <- lift_coordinates(trait, args$chain_path,
                              snp_col = args$trait_snp_col,
                              chr_col = args$trait_chr_col,
                              pos_col = args$trait_pos_col)
  }
  
  
  trait <- format_data(data.frame(trait),
                       chr_col = args$trait_chr_col,
                       pos_col = args$trait_pos_col,
                       snp_col = args$trait_snp_col,
                       beta_col = args$trait_beta_col,
                       se_col = args$trait_se_col,
                       pval_col = args$trait_pval_col,
                       log_pval = as.logical(args$trait_log_pval),
                       eaf_col = args$trait_eaf_col,
                       effect_allele_col = args$trait_effect_allele_col,
                       other_allele_col = args$trait_other_allele_col
  )
}

trait <- finemap.wrapper(out = trait, out_lbf_dir_path = args$trait_lbf_path, 
                         N_out = args$N_trait, out_type = args$trait_type, out_sd = args$trait_sd,
                         plink_path = args$plink_path,
                         bfile_path = args$bfile_path)

print(colnames(trait))

################################################################################
# Process outcome data
################################################################################


if(args$is_out_eqtl_cat){
  out <- format_eqtl_cat_trait(args$out_id, args$out_lbf_path, args$out_sumstats_path)
}else{
  out_raw <- fread(args$out_sumstats_path)
  
  print(colnames(out_raw))
  # Filter and rename using raw column names
  out_raw <- out_raw %>% 
    filter(!!sym(args$out_chr_col) == unique(trait$chr)) 
  
  build <- get_genome_build_local(out_raw, sampled_snps = 100, path_to_38 = args$rsid38_path,
                                  snp_col = args$out_snp_col,
                                  chr_col = args$out_chr_col,
                                  pos_col = args$out_pos_col)
  if (build == "GRCH37") {
    out_raw <- lift_coordinates(out_raw, args$chain_path,
                                snp_col = args$out_snp_col,
                                chr_col = args$out_chr_col,
                                pos_col = args$out_pos_col)
  }
  
  out <- out_raw %>% filter(between(!!sym(args$out_pos_col), min(trait$pos), max(trait$pos)))
  rm(out_raw)
  
  out <- format_data(data.frame(out),
                     chr_col = args$out_chr_col,
                     pos_col = args$out_pos_col,
                     snp_col = args$out_snp_col,
                     beta_col = args$out_beta_col,
                     se_col = args$out_se_col,
                     pval_col = args$out_pval_col,
                     log_pval = as.logical(args$out_log_pval),
                     eaf_col = args$out_eaf_col,
                     effect_allele_col = args$out_effect_allele_col,
                     other_allele_col = args$out_other_allele_col
  )
}


out <- finemap.wrapper(out, args$out_lbf_dir_path, args$N_out, args$out_type, args$out_sd, trait,
                       plink_path = args$plink_path,
                       bfile_path = args$bfile_path)



################################################################################
# Run coloc based on  (SuSiE and vanilla)
################################################################################


res_coloc <- coloc.wrapper(trait, out, args$trait_id, args$out_id)


################################################################################
# ZZ and locus plot
################################################################################


plots <- plot.wrapper(trait, out, res_coloc, args$trait_id, args$out_id,
                      plink_path = args$plink_path,
                      bfile_path = args$bfile_path)
plotfile <- paste0(args$respath, "/plots_", 
                   str_replace(args$trait_id, " ", "_"), "_", 
                   str_replace(args$out_id, " ", "_"), ".pdf")
ggsave(plotfile, plots, dpi = 300, width = 15, height = 10)


################################################################################
# Investigate causality with MR
################################################################################


res_mr <- mr.wrapper(trait, out, args$N_out, res_coloc, args$trait_id, args$out_id)


################################################################################
# Store the results
################################################################################


res <- list(coloc = res_coloc,
            mr = res_mr,
            plots = plots)

resfile <- paste0(args$respath, "/fullres_", 
                  str_replace(args$trait_id, " ", "_"), "_", 
                  str_replace(args$out_id, " ", "_"), ".rds")
saveRDS(res, resfile)



