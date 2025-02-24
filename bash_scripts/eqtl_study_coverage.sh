#!/bin/bash 

#SBATCH -J eqtl_inventory
#SBATCH -c 2

#SBATCH -o logs/eqtl_inventory-%j.out 
#SBATCH -e logs/eqtl_inventory-%j.err 
#SBATCH -p long

Rscript scripts/eqtl_cat_trait_inventory.R
 