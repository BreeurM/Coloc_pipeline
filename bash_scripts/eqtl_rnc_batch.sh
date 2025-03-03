if [ -z "$1" ]; then
  echo "ERROR: Batch ID is not provided. Usage: $0 <batch_id> "
  exit 1
fi

BATCH="$1"


echo " " 
echo "Batch specific run: batch $BATCH"
echo "Run on host: $(hostname)" 
echo "Started at: $(date)" 
echo " " 


# File that lists studies in batch
STUDY_FILE="data/eQTL_catalogue/studies_batch_${BATCH}.txt"

# Extract all study ids
LIST_STUDIES=$(awk '{print $1}' "$STUDY_FILE")

# Associative array to track array job IDs per study
declare -A study_jobs

# Initialize study_jobs entries for each study
for STUDY_ID in $LIST_STUDIES; do
    study_jobs["$STUDY_ID"]=""
done


# Launch all the array jobs once the formatting has finished

echo " " 
echo "Submitting pipeline jobs" 
echo " " 

for CHR_ID in {1..22}; do
    echo "Chromosome $CHR_ID"
  
    # Reset dependency tracker for each chromosome
    previous_array_job_id=""

    for STUDY_ID in $LIST_STUDIES; do
    
        TRAITS_FILE="data/eQTL_catalogue/trait_ids/${STUDY_ID}_chr${CHR_ID}_trait_ids.txt"
        
        # Check if the trait file exists
        if [ ! -f "$TRAITS_FILE" ]; then
            echo "Trait file $TRAITS_FILE not found. Skipping study $STUDY_ID for chromosome $CHR_ID."
            # Skip rest of loop for this STUDY_ID
            continue
        fi
    
        if [ -n "$previous_array_job_id" ]; then
            array_dependency="afterany:$previous_array_job_id"
        else
            array_dependency=""
        fi
    
        # Extract the number of rows to create array
        total_traits=$(wc -l < "$TRAITS_FILE")
        
        # Submit the array job with dependency on previous study
        
        # Set the default array task limit
        array_limit=5

        # Check if the current chromosome is 6 and adjust the limit
        if [[ "${CHR_ID}" == "6" ]]; then
            array_limit=20
        fi
        
        mkdir -p "logs/${STUDY_ID}/${CHR_ID}"
        
        array_job_id=$(sbatch --parsable --dependency="$array_dependency" \
            -J "${STUDY_ID: -3}_${CHR_ID}_rnc" \
            -o "logs/${STUDY_ID}/${CHR_ID}/%x-%A_%a.out" \
            -e "logs/${STUDY_ID}/${CHR_ID}/%x-%A_%a.err" \
            --array=1-"${total_traits}"%${array_limit} \
            bash_scripts/eqtl_rnc_study_chr.sh "${STUDY_ID}" "${CHR_ID}")
        
        # Update dependency counters
        
        # Update dependency tracker for next chromosome
        previous_array_job_id="$array_job_id"
        
        # Store array job ID for the study
        study_jobs["$STUDY_ID"]+=" $array_job_id"
    
    done
done

# Submitting sweep jobs after all array jobs complete for each study

echo " " 
echo "Submitting sweep jobs" 
echo " " 

for STUDY_ID in $LIST_STUDIES; do

    # Process job IDs: trim leading space and replace spaces with colons
    job_ids="${study_jobs[$STUDY_ID]# }"
    
    if [ -n "$job_ids" ]; then
        dependency_list="${job_ids// /:}"
        capture=$(sbatch --parsable -J "${STUDY_ID: -3}_sweep" \
            --dependency=afterany:"${dependency_list}" -o logs/${STUDY_ID}/%x_%A.out\
            --wrap="Rscript scripts/eqtl_rnc_sweep_results.R ${STUDY_ID}")
    else
        echo "No jobs found for ${STUDY_ID}, skipping sweep job"
    fi

done


echo " " 
echo "Finished at: $(date)" 
echo " " 
