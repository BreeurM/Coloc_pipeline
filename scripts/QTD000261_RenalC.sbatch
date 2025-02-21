#!/bin/bash
#SBATCH -J qtd000261_rnc
#SBATCH -c 4
#SBATCH -o logs/qtd000261_rnc-%A_%a.out
#SBATCH -e logs/qtd000261_rnc-%A_%a.err
#SBATCH -p short
#SBATCH --array=1-22


echo "------------------------------------------------" 
echo "Processing Chromosome $SLURM_ARRAY_TASK_ID" 
echo "Run on host: $(hostname)" 
echo "Started at: $(date)" 
echo "------------------------------------------------" 


# File that lists trait IDs and their chromosome numbers.
TRAITS_FILE="data/eQTL_catalogue/QTD000261_trait_ids.txt"

# Get current chromosome from array index
CHROMOSOME=$SLURM_ARRAY_TASK_ID

# Extract all trait_ids for this chromosome
TRAIT_IDS=$(awk -v chr="$CHROMOSOME" '$2 == chr {print $1}' $TRAITS_FILE)

# Loop over each trait_id for this chromosome
for trait_id in $TRAIT_IDS; do
  echo " "
  echo "Processing trait: $trait_id (Chromosome $CHROMOSOME)"
  
  Rscript scripts/master.R \
    --temp_dir_path "temps/${trait_id}_temp" \
    --respath "results/QTD000261_RenalC/chr_${CHROMOSOME}" \
    --trait_id "${trait_id}" \
    --is_trait_eqtl_cat TRUE \
    --trait_lbf_path "data/eQTL_catalogue/QTD000261.lbf_variable.txt.gz" \
    --trait_sumstats_path "data/eQTL_catalogue/QTD000261.cc.tsv.gz" \
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
    --out_eaf_col "effect_allele_frequency" \
    --out_effect_allele_col "effect_allele" \
    --out_other_allele_col "other_allele"
  
  echo "Finished processing trait: $trait_id"
  echo " "
done

echo "------------------------------------------------" 
echo "All traits processed"
echo "Finished at: $(date)" 
echo "------------------------------------------------"

