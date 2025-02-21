#!/bin/bash 

#SBATCH -J msmb_prc
#SBATCH -c 4

#SBATCH -o logs/msmb_prc-%j.out 
#SBATCH -e logs/msmb_prc-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 

Rscript master.R \
  --trait_id "MSBM" \
  --trait_lbf_path "data/PrC/MSMB/"\
  --trait_sumstats_path "data/PrC/MSMB/MSMB_P08118_OID20275_rs10993994.rsids.csv" \
  --N_trait 35000 \
  --trait_chr_col "CHROM" \
  --trait_pos_col "GENPOS" \
  --trait_snp_col "RSID" \
  --trait_beta_col "BETA" \
  --trait_se_col "SE" \
  --trait_pval_col "LOG10P" \
  --trait_log_pval TRUE \
  --trait_eaf_col "A1FREQ" \
  --trait_effect_allele_col "ALLELE1" \
  --trait_other_allele_col "ALLELE0" \
  --out_id "PrC" \
  --out_lbf_dir_path "data/PrC/prostate_cancer_lbf" \
  --out_sumstats_path "data/PrC/ELLIPSE_V2_META_EUROPEAN_Results_012121.txt" \
  --N_out 178000 \
  --out_type "cc" \
  --out_sd 0.48 \
  --out_chr_col "Chromosome" \
  --out_pos_col "Position" \
  --out_snp_col "SNP_Id" \
  --out_beta_col "Estimate_Effect" \
  --out_se_col "SE" \
  --out_pval_col "P_value" \
  --out_eaf_col "EAF_Control"\
  --out_effect_allele_col "Allele_1" \
  --out_other_allele_col "Allele_2" \

  
