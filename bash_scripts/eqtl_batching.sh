#!/bin/bash 

#SBATCH -J eqtl_inventory
#SBATCH -c 2

#SBATCH -o logs/eqtl_inventory-%j.out 
#SBATCH -e logs/eqtl_inventory-%j.err 
#SBATCH -p long

conda activate renvironment

Rscript scripts/eqtl_cat_trait_inventory.R
 
Rscript scripts/eqtl_study_batching.R 250000