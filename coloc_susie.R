rm(list = ls())

library(data.table)
library(tidyverse)
library(coloc)
library(susieR)
library(TwoSampleMR)


##' Input: 
##' @param path_to_exposure:  path to exposure data, formatted with TwoSampleMr, 
##'                           file type tbd
##' @param path_to_outcome:   path to outcome data, formatted with TwoSampleMr, 
##'                           file type tbd
##' @param path_to_LD_matrix: file type tbd, NULL by default. If NULL, LD matrix will be
##'                           extracted from 1000G panel using TwoSampleMR
##' @param exposure_BF_col:   name of the column containing the Bayes factors in exposure data
##'                           NULL by default. If NULL, BFs will be computed.
##' @param outcome_BF_col:    name of the column containing the Bayes factors in exposure data
##'                           NULL by default. If NULL, BFs will be computed.
##' @param window_size:       window_size around SNP of interest, default = 1000kb

main <- function(path_to_exposure,
                 path_to_outcome,
                 path_to_LD_matrix = NULL,
                 exposure_BF_column = NULL,
                 outcome_BF_column = NULL,
                 window_size = 1000){
  
  # Get snp rsid in the region of interest
  
  # either with plink and 1000G or just by position in exp/out data ?
  
  
  # For now assume exp/out data stored in csv
  
  exp_data <- read.csv(path_to_exposure)
  out_data <- read.csv(path_to_outcome)
  
  
}