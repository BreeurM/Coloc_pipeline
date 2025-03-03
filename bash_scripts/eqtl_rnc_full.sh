#!/bin/bash

#SBATCH -J eqtl_rnc_full
#SBATCH -c 8
#SBATCH -o logs/%x-%A.out
#SBATCH -e logs/%x-%A.err
#SBATCH -p long

echo "------------------------------------------------" 
echo "Run on host: $(hostname)" 
echo "Started at: $(date)" 
echo "------------------------------------------------" 


# conda activate renvironment


# File that lists studies in batch
STUDY_FILE="data/eQTL_catalogue/study_coverages.txt"

# Extract all study ids
LIST_STUDIES=$(awk '{print $1}' "$STUDY_FILE")


echo -e "Extracting traits and chromosomes for unprocessed studies...\n" 


for STUDY_ID in $LIST_STUDIES; do
    # Check if no files start with STUDY_ID in the directory
    if ! compgen -G "data/eQTL_catalogue/trait_ids/${STUDY_ID}*" > /dev/null; then
        echo "Processing $STUDY_ID..."
        Rscript scripts/eqtl_study_format.R "${STUDY_ID}"
    fi
done

echo " "
echo "Finished at: $(date)"
echo " "

for BATCH in {1..6}; do

    echo "------------------------------------------------" 
    echo "Launching pipeline for batch $BATCH " 
    echo "------------------------------------------------" 
    
    bash_scripts/eqtl_rnc_batch.sh ${BATCH}
    
    echo -e "\nWaiting for the jobs to run...\n"
    
    # Wait for all the result files to exist, indicating that the batch has run through
    
    # Lists studies in batch
    
    STUDY_FILE_BATCH="data/eQTL_catalogue/studies_batch_${BATCH}.txt"
    LIST_STUDIES_BATCH=$(awk '{print $1}' "$STUDY_FILE_BATCH")
    
    # timeout=$((432000 / BATCH))
    
    timeout=864000
    
    while [[ $timeout -gt 0 ]]; do
        # Assume all files exist until proven otherwise
        all_files_exist=1
    
        # Check each STUDY_ID in the list
        for STUDY_ID in $LIST_STUDIES_BATCH; do
            file="results/${STUDY_ID}_RenalC/combined_coloc_results.rds"
            if [[ ! -f "$file" ]]; then
                # Flag that at least one file is missing
                all_files_exist=0 
                break
            fi
        done
    
        # Exit loop if all files exist
        if [[ $all_files_exist -eq 1 ]]; then
            break
        fi
    
        # Otherwise, wait and decrement timeout
        sleep 3600
        ((timeout -= 3600))
    done

    # Timeout handling
    if [[ $timeout -le 0 ]]; then
        echo "Timeout reached while waiting for batch $BATCH to finish"
        exit 1
    fi

done


