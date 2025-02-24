echo "------------------------------------------------" 
echo "Run on host: $(hostname)" 
echo "Started at: $(date)" 
echo "------------------------------------------------" 

# File that lists studies available

STUDY_FILE="data/eQTL_catalogue/test.txt"


# Extract all study ids

LIST_STUDIES=$(awk '{print $2}' "$STUDY_FILE")

for STUDY_ID in $LIST_STUDIES; do

  echo " "
  echo "Processing study: $STUDY_ID"
  
  which sbatch 
  
  # Step 1: Format eQTL data
  # Dependency on previous study's coloc job
  
  if [ -n "$previous_array_job_id" ]; then
    array_dependency="afterok:$previous_array_job_id"
  else
    array_dependency=""
  fi


  format_job_id=$(sbatch --parsable --dependency="$array_dependency" -p short -J "${STUDY_ID}_format" -o logs/%x_%A.out --wrap="Rscript scripts/eqtl_cat_format.R ${STUDY_ID}")
  
  TRAITS_FILE="data/eQTL_catalogue/${STUDY_ID}_trait_ids.txt"
  total_traits=$(wc -l < "$TRAITS_FILE")
    
  # Step 2: Submit the array job with dynamic study ID and array range, 
  # Capture the job ID for future iteration
  
  array_job_id=$(sbatch --parsable --dependency=afterok:$format_job_id -J "${STUDY_ID}_rnc" --array=1-"${total_traits}"%100  bash_scripts/eqtl_single_study_rnc.sh "${STUDY_ID}")
  
    
  # Step 3: Submit the final sweep job with dependency
  
  sbatch -J "${STUDY_ID}_sweep" --dependency=afterok:"${array_job_id}" -o logs/%x_%A.out --wrap="Rscript scripts/eqtl_rnc_sweep_results.R ${STUDY_ID}"


  # Update dependency tracker for next iteration
  previous_array_job_id="$array_job_id"

  echo "Done."
  echo " "
done
  
  