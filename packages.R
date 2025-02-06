# Install CRAN packages
install.packages("remotes")
remotes::install_version("Matrix", version = "1.6-1")
remotes::install_version("MASS", version = "7.3-60") 

install.packages(c(
  "renv", "httr", "glue", "jsonlite", "ggrepel", "ggpubr",
  "readr", "tictoc", "data.table", "tidyverse", "coloc", 
  "susieR", "TwoSampleMR"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("liftOver")


