#!/bin/bash

#SBATCH -J ""$1"_rnc" 
#SBATCH -c 4
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH -p long
#SBATCH --array=1-100%50


# Check if study ID is provided
if [ -z "$1" ]; then
  echo "ERROR: Study ID not provided. Usage: $0 <study_id>"
  exit 1
fi
STUDY_ID="$1"

TRAITS_FILE="data/eQTL_catalogue/${STUDY_ID}_trait_ids.txt"

# Get the current line from the array task ID
LINE=$SLURM_ARRAY_TASK_ID

# Extract trait_id and chromosome from the corresponding line
trait_id=$(awk -v line="$LINE" 'NR == line {print $1}' "$TRAITS_FILE")
CHROMOSOME=$(awk -v line="$LINE" 'NR == line {print $2}' "$TRAITS_FILE")

echo "------------------------------------------------" 
echo "Processing Study: $STUDY_ID, Trait: $trait_id (Chromosome $CHROMOSOME)" 
echo "Run on host: $(hostname)" 
echo "Started at: $(date)" 
echo "------------------------------------------------" 

Rscript scripts/master.R \
    --temp_dir_path "temps/${STUDY_ID}_${trait_id}_temp" \
    --respath "results/${STUDY_ID}_RenalC/chr_${CHROMOSOME}" \
    --trait_id "${trait_id}" \
    --is_trait_eqtl_cat TRUE \
    --trait_lbf_path "data/eQTL_catalogue/${STUDY_ID}.lbf_variable.txt.gz" \
    --trait_sumstats_path "data/eQTL_catalogue/${STUDY_ID}.cc.tsv.gz" \
    --N_trait 35000 \
    --out_id "RenalC" \
    --out_lbf_dir_path "data/RenalC/rnc_lbf" \
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

echo "------------------------------------------------" 
echo "Finished processing trait: $trait_id"
echo "Finished at: $(date)" 
echo "------------------------------------------------"



