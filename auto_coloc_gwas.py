"""
Used to perform coloc between 2 GWAS traits
"""

## Import necessary libraries
import pandas as pd  # Data manipulation
import math  # Mathematical operations
import numpy as np  # Numerical operations
import sys  # System-specific functions
import argparse  # Command-line argument parsing
import warnings  # Suppress warnings
warnings.filterwarnings('ignore')  # Ignore warnings for clean output
from itertools import repeat  # Repeating elements in loops
import sumstats  # For handling GWAS summary statistics
import h5py  # Reading HDF5 files
import random  # Random number generation (not used here)
from hashlib import md5  # Hashing (not used here)
import pyranges as pr  # Genomic range operations
import os.path  # File path handling
from tqdm import tqdm  # Progress bar for loops

# Parse command-line arguments
def check_arg(args=None):
    """
    Parse input arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description='Colocalisation test between SNPs from 2 GWAS traits',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('-g1', '--gwas1', help='HDF5 Summary statistics from GWAS trait 1', required='True')
    parser.add_argument('-g2', '--gwas2', help='HDF5 Summary statistics from GWAS trait 2', required='True')
    parser.add_argument('-l', '--list', help='Text file of SNPs to test - one SNP per line', default="list.txt")  
    parser.add_argument('-o', '--output', help='Output file for results', default='results.txt')  
    parser.add_argument('-kb', '--kb_window', help='Window size (in kb) around the SNP (default: 100000)', default=100)
    parser.add_argument('-c', '--cytoband', help='Cytoband file for genomic ranges', default='cytoband_hg19.txt')
    results = parser.parse_args(args)
    
    return (results.gwas1, results.gwas2, results.list, results.output, results.kb_window, results.cytoband)

# Perform colocalization analysis
def coloc1(
    trait1_lnbfs,  # Log-bayes factors (LBFs) for GWAS 1
    trait2_lnbfs,  # Log-bayes factors (LBFs) for GWAS 2
    prior1: float = 1e-4,  # Prior probability for GWAS 1
    prior2: float = 1e-4,  # Prior probability for GWAS 2
    prior12: float = 1e-5  # Prior probability for shared association
):
    """
    Compute posterior probabilities for colocalization.
    
    Parameters:
    trait1_lnbfs: Log-bayes factors (LBFs) for trait 1
    trait2_lnbfs: Log-bayes factors (LBFs) for trait 2
    prior1      : default 1e-4, Prior probability for trait 1
    prior2      : default 1e-4, Prior probability for trait 2
    prior12     : default 1e-5, Prior probability for shared association  
    """
    # Logarithmic probabilities for different scenarios
    log_numerators = (
        0,  # PP0: Neither trait associated
        math.log(prior1) + sumstats.log_sum(trait1_lnbfs),  # PP1: Only GWAS 1 associated
        math.log(prior2) + sumstats.log_sum(trait2_lnbfs),  # PP2: Only GWAS 2 associated
        math.log(prior1) + math.log(prior2) + sumstats.log_sum(  # PP3: Both traits associated but independent
            trait1_lnbf + trait2_lnbf
            for i, trait1_lnbf in enumerate(trait1_lnbfs)
            for j, trait2_lnbf in enumerate(trait2_lnbfs)
            if i != j
        ),
        math.log(prior12) + sumstats.log_sum(  # PP4: Shared association
            trait1_lnbf + trait2_lnbf
            for trait1_lnbf, trait2_lnbf in zip(trait1_lnbfs, trait2_lnbfs)
        )
    )
    # Normalize to compute posterior probabilities
    yield from (
        math.exp(log_numerator - sumstats.log_sum(log_numerators))
        for log_numerator in log_numerators
    )

# Retrieve PP4 value
def get_PP4(title: str, pp0: float, pp1: float, pp2: float, pp3: float, pp4: float):
    """
    Extract PP4 value for shared association.
    """
    return pp4


def get_PP(df, snp):
    """
    Compute the colocalization posterior probability (PP4) for a SNP.
    
    Args:
    df (pd.DataFrame): DataFrame containing the GWAS summary statistics, 
                        with columns for the log-bayes factors (LBF) for both traits.
    snp (str)        : The SNP identifier for which to compute the PP4.

    Returns:
    pd.DataFrame: A DataFrame containing the SNP and its computed PP4 value.
    """
    # Compute the colocalization posterior probabilities using the coloc1 function
    PP4 = get_PP4("snp", *coloc1(df.gwas_1_lbf.to_numpy(), df.gwas_2_lbf.to_numpy()))  # Compute PP4 using LBFs from both GWAS traits
    
    PP4_rounded = round(PP4, 2)  # Round PP4 to two decimal places for better readability
    
    # Prepare the output data with SNP identifier and corresponding PP4 value
    raw = [[snp, PP4_rounded]]  # Create a list containing the SNP and its rounded PP4 value
    
    # Create a DataFrame to store the result, with columns "SNP" and "PP4"
    df2 = pd.DataFrame(data=raw, columns=["SNP", "PP4"])
    
    # Return the DataFrame containing the SNP and its PP4 value
    return df2


def get_snp_window(kb_window, SNP, SNP_file, gwas1, gwas2, cyto):
    """
    Extract data within a specified window around a SNP.
    
    Args:
    kb_window (int)         : The size of the window (in kilobases) around the SNP.
    SNP (str)               : The SNP identifier for which to extract data.
    SNP_file (pd.DataFrame) : DataFrame containing SNP positions and identifiers.
    gwas1 (str)             : File path to trait 1 sumstat file (in HDF5 format).
    gwas2 (str)             : File path to trait 2 sumstat file (in HDF5 format).
    cyto (pyranges.PyRanges): PyRanges object containing cytoband information.
    
    Returns:
    pd.DataFrame: A DataFrame containing merged GWAS data for the SNP within the specified window.
    """
    # Find the row corresponding to the SNP and reset the index for further operations
    small = SNP_file[SNP_file["SNP"] == SNP].reset_index()  # Locate the SNP's position in the file
    
    position = small.at[0, 'Start']  # Extract the SNP's starting position
    upper = int(position) + kb_window  # Define the upper boundary of the window
    lower = int(position) - kb_window  # Define the lower boundary of the window
    chrom = small.Chromosome.astype(int)[0]  # Extract chromosome information for the SNP
    
    # Create a PyRanges object representing the genomic region of interest
    data = {'Chromosome': [chrom], 'Start': [lower], 'End': [upper]}
    gr3 = pr.from_dict(data)  # Convert the data dictionary into a genomic range
    
    # Join the genomic range with cytoband information
    cyto_temp = cyto.join(gr3).df  # Match cytoband data based on the genomic region
    cytoband_list = cyto_temp.cytoband.to_list()  # Extract the list of cytobands in the region

    # Read the relevant sections of the GWAS data files based on cytoband information
    gwas1_df = pd.concat([pd.read_hdf(gwas1, i) for i in cytoband_list], ignore_index=True)  # Load GWAS 1 data
    gwas2_df = pd.concat([pd.read_hdf(gwas2, i) for i in cytoband_list], ignore_index=True)  # Load GWAS 2 data

    # Filter GWAS data within the defined window (both lower and upper bounds)
    gwas1_clean = gwas1_df[(gwas1_df["Start"] >= lower) & (gwas1_df["Start"] <= upper)]
    gwas2_clean = gwas2_df[(gwas2_df["Start"] >= lower) & (gwas2_df["Start"] <= upper)]

    # Merge the filtered GWAS data from both traits by the SNP identifier
    merged = gwas1_clean.merge(gwas2_clean, on="SNP")  # Merge based on the SNP column

    # Return the relevant columns from the merged data for downstream analysis
    return merged[["SNP", "gwas_1_lbf", "gwas_2_lbf"]]  # Return a DataFrame with SNP and LBF values


def error_snps(SNP):
    """
    Return a placeholder result for problematic SNPs.
    """
    return pd.DataFrame([[SNP, np.nan]], columns=["SNP", "PP4"])


def process_snps(kb_window, SNP, SNP_file, gwas1, gwas2, cyto):
    """
    Process SNP data and compute posterior probabilities.
    
    Args:
    kb_window (int)         : The size of the window around the SNP in kilobases.
    SNP (str)               : The SNP identifier to process.
    SNP_file (pd.DataFrame) : DataFrame containing SNP positions and identifiers.
    gwas1 (str)             : File path to trait 1 sumstat file (in HDF5 format).
    gwas2 (str)             : File path to trait 2 sumstat file (in HDF5 format).
    cyto (pyranges.PyRanges): PyRanges object containing cytoband information.
    
    Returns:
    pd.DataFrame: A DataFrame containing the computed posterior probabilities for the SNP.
    """
    try:
        # Get SNP data for the specified window (using the provided kb_window)
        merged = get_snp_window(kb_window, SNP, SNP_file, gwas1, gwas2, cyto)
        
        # Calculate the posterior probabilities (PP4) for the SNP
        return get_PP(merged, SNP)
    except:
        # In case of any error, return a default response with NaN value for the PP
        return error_snps(SNP)


def main(gwas1, gwas2, list_file, output, kb_window, cyto):
    """
    Execute colocalization for all SNPs in the list.
    
    Args:
    gwas1 (str)    : File path to trait 1 sumstat file (in HDF5 format).
    gwas2 (str)    : File path to trait 2 sumstat file (in HDF5 format).
    list_file (str): File path to a text file containing SNPs to process (one SNP per line).
    output (str)   : File path to save the results.
    kb_window (int): The size of the window around each SNP in kilobases.
    cyto (str)     : File path to the cytoband file containing chromosome, start, end, and cytoband info.
    
    Returns:
    None: The function saves the results as a text file.
    """
    # Convert the specified window size from kilobases to base pairs (divide by 2 for the window)
    kb_window = round((int(kb_window) * 1000) / 2, 1)
    
    # Read cytoband file into a DataFrame
    cyto_df = pd.read_csv(cyto, sep="\t")
    
    # Convert cytoband data to PyRanges format for genomic range operations
    gr = pr.PyRanges(cyto_df)
    
    # Read SNP list (from file) and assign proper column names
    SNP_list_raw = pd.read_csv(list_file, sep="\t", header=None)
    SNP_list_raw.columns = ["Chromosome", "Start", "SNP"]  # Set column names for SNP list
    SNP_list_raw["End"] = SNP_list_raw.Start.astype(int) + 1  # Define 'End' column (start position + 1)
    
    # Join SNP list with cytoband data (for genomic position matching)
    SNP_list = pr.PyRanges(SNP_list_raw).join(gr).df
    
    # Extract SNPs as a flat list for processing
    flat_list = SNP_list.SNP.tolist()

    # Initialize an empty DataFrame to store the results
    results = pd.DataFrame()
    
    # Iterate through the SNP list with a progress bar
    for snp in tqdm(flat_list):
        # Process each SNP and append the result to the results DataFrame
        results = results.append(process_snps(kb_window, snp, SNP_list, gwas1, gwas2, gr))
    
    # Save the results DataFrame to a file (tab-separated format)
    results.to_csv(output, sep="\t", index=None)


if __name__ == '__main__':
    gwas1, gwas2, list_file, output, kb_window, cyto = check_arg(sys.argv[1:])
    main(gwas1, gwas2, list_file, output, kb_window, cyto)
