format_data:
	- outputs unnecessary columns, add a select at the end
	- original 2SMR read_exp_data function seemed to deal with cases where -log10 is reported automatically
	
get_genome_build/lift_coordinates:
  - merging gets rid of the NA rsID
  - chromosome X not in the ref file
	
finemap_wrapper
  - case where less than 10 cs is not accounted for?
	- Don't output BF when no cs
	- Handle the flags from finemap_susie
  - Only works if "trait" is in eQTL catalogue format, i.e. needs "region" to be defined
  - Does not echo the correct susie error
  
finemap-susie
  - grid search over L
  
master(_pipeline):
  - Figure out how to get the sample size (for now only one outcome of interest, can be set externally)
  - Detect whether lbf are in dataset + whether data comes from eqtl_cat
  - Prefix of lbf files is hard coded
  
Overall:
  - effect of the priors: estimate pi from empirical Bayes?
  - rn absence of cs is dealt with in a dirty way.
      - Trim eqtl data to get rid of 0s
      - Don't output BF when no cs in region_susie