import os 
import pandas as pd 
import numpy as np
import sys
import argparse

# Function to parse command-line arguments
def check_arg(args=None):
    """
    Parse and validate command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Extract GTEx eQTL information',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-c', '--chrom', help='Current chromosome', required=True)
    parser.add_argument('-e', '--ensembl', help='Ensembl Gene ID', required=True)
    parser.add_argument('-r', '--rsid', help='SNP rsid', required=True)
    parser.add_argument('-p', '--plink', help='Full path for PLINK file')
    parser.add_argument('-t', '--tissue', help='Full path to the GTEx file (allpairs)', required=True)
    parser.add_argument('-n', '--tissue_name', help='Tissue name for graph and file output', required=True)
    parser.add_argument('-ref', '--reference', help='Path to the directory containing the GTEx Map file', required=True)
    parser.add_argument('-d', '--dataset', help='AGs master dataset')
    parser.add_argument('-g', '--gene', help='Gene name for graph', required=True)
    
    results = parser.parse_args(args)
    return (
        results.chrom, results.ensembl, results.rsid, results.plink, 
        results.tissue, results.tissue_name, results.reference, 
        results.dataset, results.gene
    )


# Function to automate GTEx eQTL extraction and processing
def automate(ENSEMBLE, CHR, RSID, GENE_NAME, PLINK, TISSUE, TISSUE_NAME, PATH_REF, DATASET):
    """
    Extract and process eQTL data for the specified gene, chromosome, and SNP.

    Parameters:
    ENSEMBLE    : Ensembl Gene ID
    CHR         : Chromosome number
    RSID        : SNP rsid
    GENE_NAME   : Gene name
    PLINK       : Path to PLINK file (Optional)
    TISSUE      : Path to GTEx file
    TISSUE_NAME : Name of the tissue for output labeling
    PATH_REF    : Directory path to GTEx reference Map file
    DATASET     : AGs master dataset (Optional)

    Process:
    1. Extract relevant entries from the GTEx file for the given Ensembl ID.
    2. Merge extracted data with reference Map data.
    3. Save the processed data to an output file.

    Outputs:
    A text file containing processed eQTL information.
    """
    # Generate output filename for extracted data
    file_name_lung = RSID + ".find." + GENE_NAME + "_" + TISSUE_NAME + ".txt"

    # Extract lines matching the specified Ensembl Gene ID from the GTEx file
    command1 = "zgrep " + str(ENSEMBLE) + " " + TISSUE + " > " + file_name_lung
    os.system(command1)  # Execute the shell command

    # Define column headers for the extracted data
    header = ['gene_id', 'variant_id', 'tss_distance', 'ma_samples', 
              'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se']

    # Load the extracted data into a DataFrame
    lung_find = pd.read_csv(file_name_lung, sep="\t", header=None)
    lung_find.columns = header

    # Load the reference data for the specified chromosome
    file_name = PATH_REF + "ref_GTEX_hg38_chr" + str(CHR) + ".txt"
    chrom_ref = pd.read_csv(file_name, sep="\t")

    # Merge the extracted data with the reference data
    lung_find = lung_find.merge(chrom_ref, on="variant_id")

    # Save the processed data to the output file
    lung_find.to_csv(file_name_lung, sep="\t", index=None)



def main(chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset , gene ):
    automate(ensembl, chrom, rsid, gene, plink, tissue, tissue_name, reference, dataset)



if __name__ == '__main__':

    chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene  = check_arg(sys.argv[1:])
    main(chrom, ensembl , rsid , plink , tissue , tissue_name ,reference, dataset, gene )

