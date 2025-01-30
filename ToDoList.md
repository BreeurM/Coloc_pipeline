format_data:
	- outputs unnecessary columns, add a select at the end
	
region_susie
  - case where less than 10 cs is not accounted for?
	- Don't output BF when no cs
	- Handle the flags from finemap_susie
	
master_pipeline:
  - Figure out how to get the sample size (for now only one outcome of interest, can be set externally)
  - Detect whether lbf are in dataset + whether data comes from eqtl_cat
  - Prefix of lbf files is hard coded
  
Overall:
  - effect of the priors: estimate pi from empirical Bayes?
  - rn absence of cs is dealt with in a dirty way.
      - Trim eqtl data to get rid of 0s
      - Don't output BF when no cs in region_susie