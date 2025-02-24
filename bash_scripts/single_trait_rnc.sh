#!/bin/bash 

#SBATCH -J eqtl_rnc
#SBATCH -c 4

#SBATCH -o logs/eqtl_rnc-%j.out 
#SBATCH -e logs/eqtl_rnc-%j.err 
#SBATCH -p short 
 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 

Rscript scripts/master.R \
  --temp_dir_path "temps/ENSG00000143537.grp_2.contained.ENST00000360674_temp" \
  --respath "results/QTD000264_RenalC/chr_2" \
  --trait_id "ENSG00000143537.grp_2.contained.ENST00000360674" \
  --is_trait_eqtl_cat TRUE \
  --trait_lbf_path "data/eQTL_catalogue/QTD000264.lbf_variable.txt.gz"\
  --trait_sumstats_path "data/eQTL_catalogue/QTD000264.cc.tsv.gz" \
  --N_trait 35000 \
  --out_id "RenalC" \
  --out_lbf_dir_path "data/RenalC/renal_cancer_lbf" \
  --out_sumstats_path "data/RenalC/EURO_ONLY.tsv" \
  --N_out 780000 \
  --out_type "cc" \
  --out_sd 0.034 \
  --out_chr_col "chromosome" \
  --out_pos_col "base_pair_location" \
  --out_snp_col "rsid" \
  --out_beta_col "beta" \
  --out_se_col "standard_error" \
  --out_pval_col "p_value" \
  --out_eaf_col "effect_allele_frequency"\
  --out_effect_allele_col "effect_allele" \
  --out_other_allele_col "other_allele" \
