import h5py
import warnings
from tables import NaturalNameWarning
warnings.filterwarnings('ignore', category=NaturalNameWarning)  # Suppress warnings for natural naming in HDF5
import numpy as np
import pandas as pd
import pyranges as pr  # For working with genomic ranges
import sys
import argparse

# Function to parse and validate command-line arguments
def check_arg(args=None):
    """
    Parse command-line arguments for the colocalization test script.
    """
    parser = argparse.ArgumentParser(
        description='Colocalisation test between SNPs from 2 GWAS traits',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-i', '--input', help='Summary statistics from GWAS trait', required=True)
    parser.add_argument('-o', '--output', help='File name of the results', default='results.txt')
    parser.add_argument('-b', '--beta', help='Name of the beta column in GWAS1 if it is not named as beta', default='beta')
    parser.add_argument('-se', '--se', help='Name of the se column in GWAS1 if it is not named as se', default='se')
    parser.add_argument('-var', '--variant_col', help='Name of the rsid column in GWAS1 (default: SNP)', default='SNP')
    parser.add_argument('-c', '--cytoband', help='File with cytoband info (default: cytoband_hg19.txt)', default='cytoband_hg19.txt')
    parser.add_argument('-chr', '--chrom_col', help='Name of chromosome column', default='chr')
    parser.add_argument('-p', '--position_col', help='Name of position column', default='pos')
    parser.add_argument('-r', '--ref_col', help='Name of the reference allele column', default='ref')
    parser.add_argument('-a', '--alt_col', help='Name of the alternate allele column', default='alt')
    parser.add_argument('-t', '--test_col', help='Name of the test allele column', default='minor_allele')
    
    results = parser.parse_args(args)
    return (
        results.input, results.output, results.beta, results.se, results.variant_col, 
        results.cytoband, results.chrom_col, results.position_col, results.ref_col, 
        results.alt_col, results.test_col
    )

# Function to store results in HDF5 format
def store_hd5f(df, name):
    """
    Save DataFrame to an HDF5 file grouped by cytoband.

    Parameters:
    df  : DataFrame with processed GWAS data, including cytoband info.
    name: Name of the output HDF5 file.
    """
    for i in df.cytoband.unique():
        temp = df[df["cytoband"] == i]  # Subset by cytoband
        if not temp.empty:
            key_name = str(i)  # Use cytoband name as the key
            temp.to_hdf(name, key_name, format='t', complevel=4)  # Save to HDF5 with compression
        del temp  # Free memory for large datasets

# Main function to perform colocalization analysis
def main(input_file, output, beta, se, variant_col, cyto, chrom, pos, ref, alt, test):
    """
    Perform colocalization analysis on GWAS summary statistics.

    Steps:
    1. Load cytoband reference file.
    2. Process input GWAS data in chunks.
    3. Calculate Z-scores, R-scores, and log Bayes Factors (LBF).
    4. Join GWAS data with cytoband regions.
    5. Store the processed results in HDF5 format.

    Parameters:
    input_file : Path to GWAS summary statistics file.
    output     : Path to save processed results.
    beta       : Column name for effect size.
    se         : Column name for standard error.
    variant_col: Column name for SNP ID.
    cyto       : Path to cytoband file.
    chrom      : Column name for chromosome.
    pos        : Column name for position.
    ref        : Column name for reference allele.
    alt        : Column name for alternate allele.
    test       : Column name for test allele.
    """
    # Load cytoband reference data
    cyto_df = pd.read_csv(cyto, sep="\t")
    gr = pr.PyRanges(cyto_df)  # Convert to PyRanges object

    # Process input GWAS file in chunks for memory efficiency
    for gm_chunk in pd.read_csv(input_file, sep="\t", chunksize=2000000, low_memory=False):
        # Select relevant columns and rename for consistency
        temp = gm_chunk[[variant_col, beta, se, chrom, pos, ref, alt, test]].copy()
        temp.columns = ["SNP", "beta", "se", "Chromosome", "Start", "ref", "alt", "test_allele"]

        # Add End column for genomic range calculations
        temp["End"] = temp.Start + 1

        # Calculate additional metrics for colocalization analysis
        W = 0.2  # Weight parameter for LBF calculation
        temp["Z"] = temp["beta"] / temp["se"]  # Z-score
        temp["V"] = temp["se"] ** 2  # Variance
        temp["R"] = W ** 2 / (W ** 2 + temp["V"])  # R-score
        temp["lbf"] = 0.5 * (np.log(1 - temp["R"]) + (temp["R"] * temp["Z"] ** 2))  # Log Bayes Factor

        # Filter columns for further analysis
        temp_again = temp[["Chromosome", "Start", "End", "SNP", "beta", "se", "Z", "lbf", "ref", "alt", "test_allele"]].copy()
        gr2 = pr.PyRanges(temp_again)  # Convert to PyRanges for genomic operations

        # Join GWAS data with cytoband regions
        results = gr2.join(gr)
        results = results.df  # Convert joined PyRanges object to DataFrame

        # Select final columns for output
        results = results[["Chromosome", "Start", "SNP", "beta", "se", "cytoband", "Z", "lbf", "ref", "alt", "test_allele"]].copy()

        # Store results in HDF5 format
        store_hd5f(results, output)

# Entry point for script execution
if __name__ == '__main__':
    # Parse command-line arguments
    input_file, output, beta, se, variant_col, cyto, chrom_col, position_col, ref_col, alt_col, test_col = check_arg(sys.argv[1:])

    # Run the main analysis function
    main(input_file, output, beta, se, variant_col, cyto, chrom_col, position_col, ref_col, alt_col, test_col)
