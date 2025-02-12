################################################################################
# Utility functions - memory and temporary files
################################################################################



#' Get Available Memory on the System
#'
#' This helper function determines the total available memory on the system,
#' which can be used to dynamically allocate resources for memory-intensive tasks.
#'
#' @return Numeric value representing the total memory in megabytes (MB)
get_available_memory <- function() {
  if (Sys.info()["sysname"] == "Windows") {
    # On Windows, estimate memory (simplified, may vary by system)
    # Define the input
    system_output <- system("wmic OS get TotalVisibleMemorySize /Value", intern = TRUE)
    as.numeric(sub(".*=(\\d+).*", "\\1", system_output[3])) / 1024
  } else {
    # On Unix-based systems, use 'free' command
    as.numeric(system("free -m | awk '/Mem:/ {print $2}'", intern = TRUE))
  }
}

# Functions for managing temporary files in a specified directory

#' Create an empty temporary directory
#'
#' @param path The path to the directory to create
#'
set.temp.dir <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE) # Remove existing directory and contents
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE) # Create a new empty directory
  return(path)
}

#' Create a temporary file in the specified directory
#'
#' @param name The name of the file to create
#' @param dir The directory where the file should be created
#'
#' @return The full path to the created file
#'
#' @throws Error if the temporary directory is not specified
#'
temp.file <- function(name, dir = NULL) {
  if (is.null(dir)) {
    stop("Temporary directory not specified. Use set.temp.dir first.")
  }
  file_path <- file.path(dir, name)
  file.create(file_path) # Create the file
  return(file_path) # Return the file path
}

#' Delete the temporary directory and its contents
#'
#' @param path The path to the directory to delete
#'
cleanup.temp.dir <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE) # Remove directory and all contents
  }
}



################################################################################
# Data formatting
################################################################################



#' Read and format exposure or outcome data
#' Checks and organises columns for use with MR or enrichment tests.
#' Infers p-values when possible from beta and se.
#'
#' @param dat Data frame. Must have header with at least SNP column present.
#' @param suffixe Is this the exposure or the outcome data that is being read in? The default is `"exposure"`.
#' @param snps SNPs to extract. If NULL then doesn't extract any and keeps all. The default is `NULL`.
#' @param snp_col Required name of column with SNP rs IDs. The default is `"SNP"`.
#' @param beta_col Required for MR. Name of column with effect sizes. The default is `"beta"`.
#' @param se_col Required for MR. Name of column with standard errors. The default is `"se"`.
#' @param eaf_col Required for MR. Name of column with effect allele frequency. The default is `"eaf"`.
#' @param effect_allele_col Required for MR. Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"effect_allele"`.
#' @param other_allele_col Required for MR. Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"other_allele"`.
#' @param pval_col Required for enrichment tests. Name of column with p-value. The default is `"pval"`.
#' @param units_col Optional column name for units. The default is `"units"`.
#' @param ncase_col Optional column name for number of cases. The default is `"ncase"`.
#' @param ncontrol_col Optional column name for number of controls. The default is `"ncontrol"`.
#' @param samplesize_col Optional column name for sample size. The default is `"samplesize"`.
#' @param gene_col Optional column name for gene name. The default is `"gene"`.
#' @param id_col The default is `"id"`.
#' @param min_pval Minimum allowed p-value. The default is `1e-200`.
#' @param z_col The default is `"z"`.
#' @param info_col The default is `"info_col"`.
#' @param chr_col The default is `"chr_col"`.
#' @param pos_col The default is `"pos"`.
#' @param log_pval The pval is -log10(P). The default is `FALSE`.
#'
#' @author Optimised from TwoSampleMR https://github.com/MRCIEU/TwoSampleMR/blob/master/R/read_data.R
format_data <- function(dat, header = TRUE, snp_col = "SNP",
                        beta_col = "beta", se_col = "se", eaf_col = "eaf",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele", pval_col = "pval",
                        z_col = "z", chr_col = "chr",
                        pos_col = "pos", log_pval = FALSE,
                        other_to_keep = NULL, w = .2) {
  # if (inherits(dat, "data.table")) {
  #   datname <- deparse(substitute(dat))
  #   stop(paste0(
  #     "Your ", datname, " data.frame is also of class 'data.table', ",
  #     "please reformat as simply a data.frame with ", datname, " <- data.frame(",
  #     datname, ") and then rerun your format_data() call."
  #   ))
  # }
  
  all_cols <- c(
    snp_col, beta_col, se_col, eaf_col,
    effect_allele_col, other_allele_col, pval_col,
    z_col, chr_col, pos_col, other_to_keep
  )
  
  dat <- dat %>%
    dplyr::select(any_of(all_cols))
  
  if (!(snp_col %in% names(dat))) {
    stop("SNP column not found")
  }
  
  dat <- dat %>%
    mutate(SNP = tolower(!!sym(snp_col)) %>% str_replace_all("[[:space:]]", "")) %>%
    filter(!is.na(SNP)) %>% # Might get rid of it when we switch to chr and pos
    dplyr::select(-!!sym(snp_col))
  
  dat <- dat %>%
    mutate(dup = duplicated(SNP)) %>%
    filter(!dup) %>%
    dplyr::select(-dup)
  
  dat <- dat %>%
    mutate(
      beta = as.numeric(!!sym(beta_col)),
      se = as.numeric(!!sym(se_col)),
      eaf = as.numeric(!!sym(eaf_col)),
      effect_allele = toupper(as.character(!!sym(effect_allele_col))),
      other_allele = toupper(as.character(!!sym(other_allele_col))),
      pval = as.numeric(!!sym(pval_col))
    )
  
  if (log_pval) {
    dat$pval <- 10^-dat$pval
  }
  
  dat <- dat %>%
    mutate(
      beta = ifelse(!is.finite(beta), NA, beta),
      se = ifelse(!is.finite(se) | se <= 0, NA, se),
      eaf = ifelse(!is.finite(eaf) | eaf <= 0 | eaf >= 1, NA, eaf),
      effect_allele = ifelse(!grepl("^[ACTG]+$", effect_allele) & !effect_allele %in% c("D", "I"), NA, effect_allele),
      other_allele = ifelse(!grepl("^[ACTG]+$", other_allele) & !other_allele %in% c("D", "I"), NA, other_allele),
      # pval = ifelse(!is.finite(pval) | pval <= 0 | pval > 1, NA, pval),
      # pval = ifelse(pval < min_pval, min_pval, pval)
    )
  
  if ("beta" %in% names(dat) && "se" %in% names(dat) && !"pval" %in% names(dat)) {
    dat <- dat %>%
      mutate(
        pval = pnorm(abs(beta) / se, lower.tail = FALSE) * 2,
        pval_origin = "inferred"
      )
  }
  
  if (z_col %in% names(dat)) {
    dat <- dat %>%
      mutate(z = !!sym(z_col))
  } else {
    dat <- dat %>% mutate(z = beta / se)
  }
  
  if (chr_col %in% names(dat)) {
    dat <- dat %>%
      mutate(chr = !!sym(chr_col))
  }
  
  if (pos_col %in% names(dat)) {
    dat <- dat %>%
      mutate(pos = !!sym(pos_col))
  }
  
  for (col in c("SNP", "beta", "se", "effect_allele", "other_allele", "eaf")) {
    if (!col %in% names(dat)) {
      dat[[col]] <- NA
    }
  }
  
  rownames(dat) <- NULL
  
  
  dat <- dat %>% mutate(
    variant = paste0("chr", chr, "_", pos),
    R = w^2 / (w^2 + se^2), # Calculate R-score MIGHT NOT BE CORRECT
    lbf_variable0 = 0.5 * (log(1 - R) + (R * z^2)) # Calculate log Bayes factor
  )
  
  return(dat)
}


#' Infers the genome build of the summary statistics file (GRCh37 or GRCh38)
#' from the data. Uses SNP (RSID) & CHR & BP to get genome build.
#'
#' @param data data table/data frame obj of the summary statistics file for
#' the GWAS ,or file path to summary statistics file.
#' @param nThread Number of threads to use for parallel processes.
#' @param sampled_snps Downsample the number of SNPs used when inferring genome
#' build to save time.
#' @param standardise_headers Run
#' @param standardise_headers Run
#' \code{standardise_data_column_headers_crossplatform}.
#' @param mapping_file \pkg{Mungedata} has a pre-defined
#' column-name mapping file
#' which should cover the most common column headers and their interpretations.
#' However, if a column header that is in your file is missing of the mapping we
#' give is incorrect you can supply your own mapping file. Must be a 2 column
#' dataframe with column names "Uncorrected" and "Corrected". See
#' \code{data(dataColHeaders)} for default mapping and necessary format.
#' @param dbSNP version of dbSNP to be used (144 or 155). Default is 155.
#' @param header_only Instead of reading in the entire \code{data} file,
#' only read in the first N rows where N=\code{sampled_snps}.
#' This should help speed up cases where you have to read in \code{data}
#' from disk each time.
#' @param allele_match_ref Instead of returning the genome_build this will
#' return the propotion of matches to each genome build for each allele (A1,A2).
#' @inheritParams format_data
#' @inheritParams get_genome_builds
#'
#' @return ref_genome the genome build of the data
#' @author M.Breeur, adapted from https://rdrr.io/github/neurogenomics/MungeSumstats/src/R/get_genome_build.R
get_genome_build <- function(data, sampled_snps = 500, dbSNP = 155) {
  seqnames <- chr <- SNP <- pos <- NULL
  
  # Do some filtering first to avoid errors
  data <- data[complete.cases(SNP, chr, pos)]
  # also remove common incorrect formatting of SNP
  data <- data[grepl("^rs", SNP), ]
  data <- data[SNP != ".", ]
  # also deal with common misformatting of CHR
  # if chromosome col has chr prefix remove it
  data[, chr := gsub("chr", "", chr)]
  
  # if removing erroneous cases leads to <min(10k,50% org dataset) will fail -
  # NOT ENOUGH DATA TO INFER
  nrow_clean <- nrow(data)
  
  # Downsample SNPs to save time
  if ((nrow(data) > sampled_snps) && !(is.null(sampled_snps))) {
    snps <- sample(data$SNP, sampled_snps)
  } else { # nrow(data)<10k
    snps <- data$SNP
  }
  
  data <- data[SNP %in% snps, ]
  
  snp_loc_data_37 <- MungeSumstats:::load_ref_genome_data(
    snps = snps,
    ref_genome = "GRCH37",
    dbSNP = dbSNP
  )
  snp_loc_data_38 <- MungeSumstats:::load_ref_genome_data(
    snps = snps,
    ref_genome = "GRCH38",
    dbSNP = dbSNP
  )
  # convert CHR filed in ref genomes to character not factor
  snp_loc_data_37[, seqnames := as.character(seqnames)]
  snp_loc_data_38[, seqnames := as.character(seqnames)]
  # convert CHR filed in data to character if not already
  data[, chr := as.character(chr)]
  # Now check which genome build has more matches to data
  num_37 <-
    nrow(snp_loc_data_37[data, ,
                         on = c("SNP" = "SNP", "pos" = "pos", "seqnames" = "chr"),
                         nomatch = FALSE
    ])
  num_38 <-
    nrow(snp_loc_data_38[data, ,
                         on = c("SNP" = "SNP", "pos" = "pos", "seqnames" = "chr"),
                         nomatch = FALSE
    ])
  # if no matches throw error
  if (num_37 == 0 && num_38 == 0) {
    msg_err <-
      paste0(
        "No matches found in either reference genome for your ",
        "SNPs.\nPlease check their formatting (SNP, CHR and BP",
        " columns) or supply the genome build."
      )
    stop(msg_err)
  }
  if (num_37 > num_38) {
    ref_gen_num <- num_37
    ref_genome <- "GRCH37"
  } else {
    ref_gen_num <- num_38
    ref_genome <- "GRCH38"
  }
  
  message("Inferred genome build: ", ref_genome)
  # add a warning if low proportion of matches found
  msg <- paste0(
    "WARNING: Less than 10% of your sampled SNPs matched that of ",
    "either reference genome, this may question the quality of ",
    "your summary statistics file."
  )
  if (ref_gen_num / length(snps) < 0.1) {
    message(msg)
  }
  return(ref_genome)
}



#' Infer Genome Build from Input Data by Comparing to GRCh38 Reference
#' 
#' This function determines whether the input data is aligned to the GRCh37 or GRCh38 genome build
#' by comparing SNP chromosomal positions against a GRCh38 reference. It assumes GRCh37 if there's 
#' insufficient matching with GRCh38.
#'
#' @param data A data.table containing GWAS summary statistics. Must include columns for SNP (rsID), chromosome, and position, which can be specified with snp_col, chr_col, pos_col.
#' @param sampled_snps Number of SNPs to sample for comparison (default: 500). Uses all SNPs if dataset is smaller.
#' @param path_to_38 Path to reference file for GRCh38 build (expected columns: SNP, chr, pos).
#' @param snp_col Name of the column in `data` containing SNP rsIDs. Default is "SNP".
#' @param chr_col Name of the column in `data` containing chromosome numbers. Default is "chr".
#' @param pos_col Name of the column in `data` containing base-pair positions. Default is "pos".
#' 
#' @return Character string indicating inferred genome build ("GRCH37" or "GRCH38")
#' 
#' @author M. Breeur
get_genome_build_local <- function(data, sampled_snps = 500, path_to_38,
                                   snp_col = "SNP", chr_col = "chr", pos_col = "pos") {
  
  data <- data %>% mutate(
    SNP = !!sym(snp_col),
    chr = !!sym(chr_col),
    pos = !!sym(pos_col)
  )
  # Remove rows with missing values in key columns
  data <- data[complete.cases(SNP, chr, pos)]
  
  # Filter out non-rsID SNPs and placeholder values
  data <- data[grepl("^rs", SNP), ]  # Keep only rsIDs
  data <- data[SNP != ".", ]         # Remove missing SNP codes
  
  # Standardize chromosome formatting by removing 'chr' prefix
  data[, chr := gsub("chr", "", chr)]
  
  # Downsample if dataset exceeds sampling threshold
  if ((nrow(data) > sampled_snps) && !(is.null(sampled_snps))) {
    snps <- sample(data$SNP, sampled_snps)
  } else {
    snps <- data$SNP
  }
  data <- data[SNP %in% snps, ]
  
  # Load GRCh38 reference data
  ref38 <- fread(path_to_38)
  
  # Ensure chromosome column type matches in both datasets
  ref38[, chr := gsub("chr", "", chr)]
  ref38[, chr := as.character(chr)]
  data[, chr := as.character(chr)]
  
  # Calculate match rate with GRCh38 positions
  num_38 <- nrow(
    ref38[data, , on = c("SNP" = "SNP", "pos" = "pos", "chr" = "chr"), nomatch = FALSE]
  )/sampled_snps
  
  # Determine genome build based on match proportion
  if (num_38 > 0.5) {
    ref_genome <- "GRCH38"
  } else {
    ref_genome <- "GRCH37"
  }
  
  message("Inferred genome build: ", ref_genome)
  
  # Optional low-match warning (currently commented out)
  # if (ref_gen_num/length(snps) < 0.1) {
  #   message("WARNING: Less than 10% of sampled SNPs matched reference genomes")
  # }
  
  return(ref_genome)
}




#' Lift genomic coordinates from hg19 to hg38 using chain file
#' 
#' This function converts chromosomal positions between genome builds (hg19 to hg38)
#' using UCSC chain files. It handles chromosome formatting and merges lifted positions
#' back with original data.
#'
#' @param sumstats A data.frame/data.table containing SNP coordinates. Must contain:
#'   - SNP: SNP identifier (rsID)
#'   - chr: Chromosome (format "chr1" or "1")
#'   - pos: Genomic position (hg19 coordinates)
#' @param chain_file Path to UCSC chain file for hg19-to-hg38 conversion
#' 
#' @return A data.frame with:
#'   - All original columns except chr/pos
#'   - New chr/pos columns in hg38 coordinates
#'   - Only SNPs that successfully lifted over
#' @note
#' Requires GenomicRanges and rtracklayer packages. The chain file can be downloaded
#' from UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
#' 
#' @details
#' Key steps:
#' 1. Standardizes chromosome formatting
#' 2. Creates genomic ranges object
#' 3. Applies liftOver using provided chain file
#' 4. Merges successful lifts with original data
#' 5. Returns dataframe with updated coordinates
#' 
#' @author M. Breeur
lift_coordinates <- function(sumstats, chain_file, snp_col = "SNP", chr_col = "chr", pos_col = "pos") {
  # Check required columns
  sumstats <- sumstats %>% mutate(
    SNP = !!sym(snp_col),
    chr = !!sym(chr_col),
    pos = !!sym(pos_col)
  )
  
  # Add 'chr' prefix if missing
  if (!any(grepl("^chr", sumstats$chr[1]))) {
    sumstats$chr <- paste0("chr", sumstats$chr)
  }
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = sumstats$chr,
    ranges = IRanges(start = sumstats$pos, end = sumstats$pos),
    strand = "*",
    SNP = sumstats$SNP
  )
  
  # Import chain file
  chain <- import.chain(chain_file)
  
  # Perform liftOver
  gr38 <- liftOver(gr, chain)
  gr38 <- unlist(gr38)
  
  # Convert back to dataframe
  sumstats38 <- as.data.frame(gr38)
  colnames(sumstats38) <- c("chr38", "pos38", "pos38end", "width", "strand", "SNP")
  
  # Clean chromosome names
  sumstats38$chr38 <- sub("^chr", "", sumstats38$chr38)
  sumstats38 <- sumstats38[complete.cases(sumstats38$SNP), ]
  
  # Merge with original data
  merged <- merge(sumstats, sumstats38[, c("SNP", "chr38", "pos38")], by = "SNP")
  
  # Remove original coordinates and rename new ones
  merged <- as.data.frame(merged) %>% 
    dplyr::select(-c(chr, pos)) 
  merged <- merged %>% mutate(
    !!sym(chr_col) := chr38,
    !!sym(pos_col) := pos38
  )%>% dplyr::select(-c(chr38, pos38)) 
  
  names(merged)[names(merged) %in% c("chr38", "pos38")] <- c("chr", "pos")
  
  return(merged)
}



#' Format eQTL Catalogue Trait Data
#'
#' This function reads and formats data downloaded from the eQTL Catalogue.
#'
#' @param trait A character string specifying the molecular trait ID of interest.
#' @param path_lbf The file path to the LBF data file (loaded using `fread`).
#' @param path_sumstats The file path to the summary statistics file (loaded using `fread`).
#'
#' @return A formatted data frame containing eQTL summary statistics with harmonized column names and additional variables retained.
#'
#' @import data.table
#' @import dplyr
#' @export
format_eqtl_cat_trait <- function(trait, path_lbf, path_sumstats) {
  # Load LBF and summary statistics data
  data_lbf <- fread(path_lbf)
  sumstats <- fread(path_sumstats)
  
  # Filter LBF data for the specified molecular trait
  trait_data <- data_lbf %>% filter(molecular_trait_id == trait)
  
  # Merge LBF data with summary statistics
  trait_data <- trait_data %>%
    left_join(
      select(
        sumstats,
        variant,   # Variant identifier
        rsid,      # Reference SNP ID
        molecular_trait_id,
        pvalue,    # P-value of association
        beta,      # Effect size estimate
        se,        # Standard error
        ref,       # Reference allele
        alt,       # Alternative allele
        maf        # Minor allele frequency
      ),
      by = c("variant", "molecular_trait_id")
    )
  
  # Format the data for downstream analysis
  trait <- format_data(
    trait_data,
    snp_col = "rsid",
    eaf_col = "maf",
    effect_allele_col = "ref",
    other_allele_col = "alt",
    pval_col = "pvalue",
    chr_col = "chromosome",
    pos_col = "position",
    other_to_keep = c("variant", "region", paste0("lbf_variable", 1:10))
  )
  
  # Remove allele information from the variant column
  trait$variant <- sub("(_[ACGT]+_[ACGT]+)$", "", trait$variant)
  
  return(trait)
}


################################################################################
# LD matrices utilites
################################################################################



#' Generate a Linkage Disequilibrium (LD) Matrix from 1000G
#'
#' This function calculates a linkage disequilibrium (LD) matrix for a given list of RSIDs
#' using https://mrcieu.github.io/TwoSampleMR/reference/ld_matrix.html.
#' It also formats the LD matrix and returns it in two forms: a cleaned matrix for analysis and a data frame for plotting.
#'
#' @param rsid_list    A character vector of RSIDs for which to compute the LD matrix.
#' @param with_alleles A logical value indicating whether allele information should be included
#' in the LD matrix. Default is `TRUE`.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{LD_Anal}: A cleaned LD matrix for analysis, with incomplete rows/columns removed.
#'   \item \code{LD_Plot}: A data frame version of the LD matrix, including RSIDs for plotting.
#' }
#'
#' @details
#' This function leverages the `ieugwasr::ld_matrix_local` function to compute the LD matrix
#' locally using PLINK. Any temporary files created during the computation are cleaned up
#' after the LD matrix is generated. Columns or rows with excessive missing values are
#' removed from the analytical matrix, ensuring the returned matrix is complete.
#'
#' @author K. Smith-Byrne, M. Breeur
get_ld_matrix_1000g <- function(rsid_list, with_alleles = T) {
  # Compute the LD matrix using the ieugwasr package
  LD_Full <- ieugwasr::ld_matrix(rsid_list, with_alleles = with_alleles)
  
  # Format the LD matrix for analysis
  LD <- LD_Full[, !colnames(LD_Full) %in% names(which(colSums(is.na(as.matrix(LD_Full))) > dim(LD_Full) - 1))]
  LD <- LD[complete.cases(LD), ]
  rownames(LD) <- colnames(LD)
  
  # Prepare a version for plotting
  LD_Full <- as.data.frame(LD_Full)
  LD_Full$RS_number <- rownames(LD)
  
  # Return the results
  return(list(LD_Anal = LD, LD_Plot = LD_Full))
}


#' Compute Linkage Disequilibrium (LD) Matrix from UK Biobank Data
#'
#' This function computes the LD matrix for a given list of SNPs (rsIDs) using PLINK and UK Biobank data.
#' It supports allele annotation for the LD matrix and automatically cleans up temporary files after computation.
#'
#' @param rsid_list    A character vector of SNP rsIDs for which the LD matrix is to be calculated.
#' @param plink_loc    A string specifying the path to the PLINK executable.
#' @param bfile_loc    A string specifying the path to the PLINK binary file prefix ("bfile").
#' @param plink_memory A integer with the memory in MB allocated to plink
#'                     Default is `NULL`, allocates 75% of the memory to plink.
#' @param with_alleles Logical (default = TRUE). If TRUE, row and column names in the resulting matrix will include allele information (rsID_allele1_allele2). If FALSE, only rsIDs are used.
#'
#' @return A square LD matrix (as a numeric matrix) where rows and columns correspond to the SNPs specified in `rsid_list`. The matrix contains pairwise LD values (R²).
#'
#' @details
#' The function follows these steps:
#'
#' 1. Identify Shell Type: Determines whether the system shell is "cmd" (Windows) or "sh" (Unix-based).
#' 2. Create Temporary File: Writes the `rsid_list` to a temporary file for input into PLINK.
#' 3. Generate BIM File:
#'    - Constructs a command to run PLINK with `--make-just-bim`, filtering the binary files for SNPs in `rsid_list`.
#'    - Executes the command and reads the resulting BIM file, which contains SNP metadata.
#' 4. Compute LD Matrix:
#'    - Builds another PLINK command to calculate the LD matrix using `--r square`.
#'    - Runs the command and reads the resulting LD matrix into R as a numeric matrix.
#' 5. Annotate Matrix:
#'    - If `with_alleles = TRUE`, sets matrix row and column names to include rsID and allele information (e.g., rsID_A1_A2).
#'    - Otherwise, uses only the rsID for row and column names.
#' 6. Clean Up: Deletes all temporary files created during the process.
#'
#' @examples
#' # Example usage:
#' rsid_list <- c("rs123", "rs456", "rs789")
#' plink_loc <- "/path/to/plink"
#' bfile_loc <- "/path/to/bfile"
#' ld_matrix <- get_ld_matrix_UKBB(rsid_list, plink_loc, bfile_loc, with_alleles = TRUE)
#'
#' @importFrom utils read.table write.table
#' @importFrom stats as.matrix
#' @author K. Smith-Byrne, M. Breeur
#' @export
get_ld_matrix_from_bim <- function(rsid_list, plink_loc, bfile_loc, plink_memory = NULL, with_alleles = TRUE, temp_dir_path = NULL) {
  # Determine shell type based on operating system
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
  
  # Create a temporary file to store rsID list
  set.temp.dir(temp_dir_path)
  fn <- temp.file("extracted_snp", temp_dir_path)
  write.table(data.frame(rsid_list),
              file = fn, row.names = FALSE,
              col.names = FALSE, quote = FALSE
  )
  
  # Calculate memory allocation for PLINK
  if (is.null(plink_memory)) {
    # Use 75% of available memory
    total_memory <- get_available_memory()
    plink_memory <- as.character(floor(total_memory * 0.75))
  } else {
    plink_memory <- as.character(plink_memory)
  }
  
  # Generate BIM file with PLINK
  fun1 <- paste0(
    shQuote(plink_loc, type = shell), " --bfile ",
    shQuote(bfile_loc, type = shell), " --extract ", shQuote(fn, type = shell),
    " --make-just-bim --keep-allele-order --out ", shQuote(fn, type = shell),
    " --memory ", plink_memory
  )
  system(fun1, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # Read the BIM file containing SNP metadata
  bim <- read.table(paste0(fn, ".bim"), stringsAsFactors = FALSE)
  
  # Compute LD matrix with PLINK
  fun2 <- paste0(
    shQuote(plink_loc, type = shell), " --bfile ",
    shQuote(bfile_loc, type = shell), " --extract ", shQuote(fn, type = shell),
    " --r square --keep-allele-order --out ", shQuote(fn, type = shell),
    " --memory ", plink_memory
  )
  system(fun2, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # Read the LD matrix into R as a numeric matrix
  res <- read.table(paste0(fn, ".ld"), header = FALSE) %>% as.matrix()
  
  # Annotate matrix with allele information if specified
  if (with_alleles) {
    rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep = "_")
  } else {
    rownames(res) <- colnames(res) <- bim$V2
  }
  
  # Clean up temporary files
  cleanup.temp.dir(temp_dir_path)
  
  return(res)
}




#' Align SNPs to LD Matrix
#'
#' This function aligns SNPs in a dataset to match the allele orientation of a given LD matrix.
#'
#' @param data A data frame containing SNPs with columns including `SNP`,`effect_allele`, `other_allele`, `beta`, `z`, and `eaf`.
#' @param LD_matrix A matrix with column names formatted as `SNP_A1_A2`, where `SNP` is the SNP ID, `A1` is the reference allele, and `A2` is the alternate allele.
#'
#' @return A data frame with SNPs aligned to the LD matrix, with `beta`, `z`, and `eaf` values adjusted accordingly. SNPs that cannot be aligned are excluded with a warning.
#'
#' @import dplyr
#' @import stringr
#' 
#' @author M.Breeur
align_to_LD <- function(data, LD_matrix) {
  # Extract SNP, reference allele (A1), and alternate allele (A2) from LD matrix column names
  LD_alignment <- data.frame(
    SNP = str_split_fixed(colnames(LD_matrix), "_", 3)[, 1],
    LD_A1 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 2],
    LD_A2 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 3]
  )
  
  # Merge input data with LD alignment information
  temp <- merge(data, LD_alignment) %>%
    mutate(
      temp_effect_allele = effect_allele,  # Store original effect allele
      temp_other_allele = other_allele     # Store original other allele
    )
  
  # Identify and flip alleles if necessary
  temp <- temp %>%
    mutate(
      flipped = effect_allele != LD_A1,  # Determine if flipping is required
      effect_allele = if_else(flipped, temp_other_allele, temp_effect_allele),
      other_allele = if_else(flipped, temp_effect_allele, temp_other_allele),
      beta = if_else(flipped, -beta, beta),  # Flip beta coefficient if alleles are flipped
      z = if_else(flipped, -z, z),          # Flip z-score accordingly
      eaf = if_else(flipped, 1 - eaf, eaf)  # Adjust effect allele frequency
    )
  
  # Check for inconsistencies in allele alignment and exclude problematic SNPs
  if (any(temp$LD_A1 != temp$effect_allele)) {
    n_excl <- sum(temp$LD_A1 != temp$effect_allele)  # Count excluded SNPs
    temp <- temp %>% filter(LD_A1 == effect_allele)  # Remove misaligned SNPs
    warning(paste0("Inconsistent alleles when matching to LD matrix. Excluding ", as.character(n_excl), " variants from the analysis."))
  }
  
  # Return cleaned dataset with unnecessary columns removed
  return(temp %>% dplyr::select(-c(temp_other_allele, temp_effect_allele)))
}




################################################################################
# Fine-mapping and wrapper
################################################################################


#' Perform fine mapping with SuSiE on one dataset
#'
#' This function implements the SuSiE method for a given dataset amd outputs credible sets and Bayes Factors.
#'
#' @param exp_data  Either a data frame or a file path to the exposure dataset.
#'                  The dataset must be formatted as required by TwoSampleMR
#'                  and include at least the columns: `SNP`, `beta`, `se`, and `eaf`.
#' @param N_exp     An integer specifying the sample size for the exposure dataset.
#' @param exp_type  A string specifying the type of exposure study.
#'                  Either `"quant"` (quantitative trait) or `"cc"` (case-control).
#' @param exp_sd    A numeric value specifying the std of the trait for quantitative studies or the proportion of cases for case-control studies. Default is `1`.
#'
#' @param LD_matrix Either `NULL`, a precomputed LD matrix, or a file path to a CSV file containing the LD matrix.
#'                  If `NULL`, the LD matrix is computed from 1000G using the `get_ld_matrix` function.
#'
#' @param max_iter  Max number of iterations to run susie_rss
#' @param ...       Parameters to be passed on to susie_rss
#'
#' @return A SuSiE object (`susie.res`) containing information on shared genetic signals between the exposure and outcome datasets.
#'         If SuSiE fails to converge, the function returns `NULL`.
#'
#' @author M. Breeur
#' @author Adapted from C. Wallace: https://github.com/chr1swallace/coloc/blob/main/R/susie.R
finemap_susie <- function(exp_data, N_exp, exp_type, exp_sd = 1,
                          LD_matrix = NULL, plink_loc = "plink", bfile_loc = NULL, temp_dir_path = "Temp",
                          repeat_until_convergence = FALSE, run_checks = TRUE, max_iter = 10000,
                          ...) {
  # Check input consistency
  # exp_type
  
  # Import exposure data if provided as a file path
  if (is.character(exp_data)) {
    exp_data <- read.csv(exp_data)
  }
  
  harm_data <- exp_data
  # Import LD_matrix if provided as file paths, or import from 1000G if NULL
  
  tic("Loading LD matrix")
  if (is.null(LD_matrix)) {
    if (is.null(bfile_loc) | is.null(plink_loc)) {
      warning("Either bfile_loc or plink_loc is missing. Getting LD matrix from 1000G in TwoSampleMR.")
      LD_matrix <- tryCatch(expr = {
        get_ld_matrix_1000g(harm_data$SNP, with_alleles = T)$LD_Anal
      }, error = function(e) {
        stop("Window too large, SNP list must be smaller than 500. Try reducing window or provide a local ld reference.")
        NULL
      })
    } else {
      LD_matrix <- get_ld_matrix_from_bim(harm_data$SNP,
                                          plink_loc = plink_loc,
                                          bfile_loc = bfile_loc,
                                          with_alleles = T,
                                          temp_dir_path = temp_dir_path
      )
    }
  } else {
    if (is.character(LD_matrix)) {
      LD_matrix <- read.csv(LD_matrix)
    }
  }
  toc()
  
  # Check if allele information is present in the LD matrix and align if necessary
  tic("Flipping alleles to match LD_mat")
  if (!(all(grepl("^[ACTG]+$", str_split_fixed(colnames(LD_matrix), "_", 3)[, 2])) &
        all(grepl("^[ACTG]+$", str_split_fixed(colnames(LD_matrix), "_", 3)[, 3])))) {
    # Handle case where allele information is missing or incorrect
    warning("LD matrix has no or incorrect allele information, allele alignment will be inferred from expected Z-scores vs. observed Z-scores")
    harm_data <- harm_data %>%
      dplyr::select(
        SNP, pos,
        beta,
        se,
        eaf
      )
  } else {
    # Flip alleles in harmonized data to align with the LD matrix
    LD_alignment <- data.frame(
      SNP = str_split_fixed(colnames(LD_matrix), "_", 3)[, 1],
      LD_A1 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 2],
      LD_A2 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 3]
    )
    temp <- merge(harm_data, LD_alignment) %>%
      mutate(
        temp_effect_allele = effect_allele,
        temp_other_allele = other_allele
      )
    temp <- temp %>%
      mutate(
        flipped = effect_allele != LD_A1,
        effect_allele = if_else(flipped, temp_other_allele, temp_effect_allele),
        other_allele = if_else(flipped, temp_effect_allele, temp_other_allele),
        beta = if_else(flipped, -beta, beta),
        eaf = if_else(flipped, 1 - eaf, eaf)
      )
    
    if (any(temp$LD_A1 != temp$effect_allele)) {
      # Stop execution if alleles cannot be aligned
      n_excl <- sum(temp$LD_A1 != temp$effect_allele)
      temp <- temp %>% filter(LD_A1 == effect_allele)
      warning(paste0("Inconsistent alleles when matching to LD matrix. Excluding ", as.character(n_excl), " variants from the analysis."))
    }
    
    harm_data <- temp %>%
      dplyr::select(
        SNP, pos,
        beta,
        se,
        eaf
      )
  }
  toc()
  
  # Remove alleles from LD matrix column names for cross ref with ext/out data
  colnames(LD_matrix) <- str_split_fixed(colnames(LD_matrix), "_", 2)[, 1]
  rownames(LD_matrix) <- colnames(LD_matrix)
  
  # Prepare data for SuSiE colocalization analysis
  tic("Formatting data for SuSiE")
  exp_for_coloc <- harm_data %>%
    dplyr::select(
      beta,
      se,
      SNP,
      eaf
    )
  exp_for_coloc$se <- exp_for_coloc$se^2
  colnames(exp_for_coloc) <- c("beta", "varbeta", "snp", "MAF")
  exp_for_coloc <- as.list(exp_for_coloc)
  exp_for_coloc$type <- exp_type
  if (exp_type == "quant") {
    exp_for_coloc$sdY <- exp_sd
  } else {
    exp_for_coloc$s <- exp_sd
  }
  exp_for_coloc$N <- N_exp
  exp_for_coloc$LD <- LD_matrix[exp_for_coloc$snp, exp_for_coloc$snp]
  toc()
  
  # Check alignment between exposure data and LD matrix
  # Flags for faulty allele alignment
  if (run_checks) {
    tic("Running checks for SuSiE input")
    exp_qc <- list(
      kriging = kriging_rss(exp_for_coloc$beta / sqrt(exp_for_coloc$varbeta),
                            exp_for_coloc$LD,
                            n = N_exp
      ),
      alignment_check = tryCatch(expr = {
        check_alignment(exp_for_coloc)
      }, error = function(e) {
        warning("Alignment could not be checked.")
        NULL
      })
    )
    print(exp_qc$alignment_check)
    if (!is.null(exp_qc$alignment_check)) {
      if (exp_qc$alignment_check < 0.7) {
        warning("Suspected alignment error.")
      }
    }
    toc()
  } else {
    exp_qc <- NULL
  }
  
  # Calculate Z-scores for SuSiE
  z <- harm_data$beta / harm_data$se
  
  snp <- harm_data$SNP
  names(z) <- snp
  LD <- LD_matrix[snp, snp]
  
  # Initialize convergence flag
  converged <- FALSE
  
  # Set defaults for SuSiE arguments
  susie_args <- list(...)
  
  # Include sample size parameter if not already specified
  susie_args <- c(list(n = N_exp), susie_args)
  
  # Set prior variance with regard to trait type and sdY
  if (exp_type == "quant") {
    susie_args <- c(list(prior_variance = (.15 / exp_sd)^2, estimate_prior_variance = FALSE), susie_args)
  } else {
    susie_args <- c(list(prior_variance = .2^2, estimate_prior_variance = FALSE), susie_args)
  }
  
  # Iteratively run SuSiE until convergence or maximum iterations are reached
  while (!converged) {
    message("Running max iterations: ", max_iter)
    susie.res <- do.call(
      susie_rss,
      c(list(z = z, R = LD, max_iter = max_iter), susie_args)
    )
    converged <- susie.res$converged
    message("\tConverged: ", converged)
    if (!converged && repeat_until_convergence == FALSE) {
      stop("susie_rss() did not converge in ", max_iter, " iterations. Try running with run_until_convergence=TRUE")
    }
    if (!converged) {
      max_iter <- max_iter * 10
    } # Increase iterations if not converged
  }
  
  return(susie.res)
}


region_utils <- function(region, lead_pos, lbf_directory) {
  # Find out region of interest
  
  lead_chr <- str_remove(str_split_fixed(region, ":", 2)[1], "chr")
  pos_interval <- str_split_fixed(region, ":", 2)[, 2]
  if (length(str_split(pos_interval, "-")[[1]]) == 2) {
    pos_low <- as.numeric(str_split(pos_interval, "-")[[1]][1])
    pos_high <- as.numeric(str_split(pos_interval, "-")[[1]][2])
  } else {
    # First entry negative
    pos_low <- -as.numeric(str_split(pos_interval, "-")[[1]][2])
    pos_high <- as.numeric(str_split(pos_interval, "-")[[1]][3])
  }
  
  required_start <- lead_pos - 500000
  required_end <- lead_pos + 500000
  
  # List files in lbf_directory
  
  lbf_file_list <- list.files(lbf_directory, full.names = FALSE)
  
  # Extract chromosome and positions
  coverage <- as.data.frame(
    do.call(rbind, lapply(lbf_file_list, function(file) {
      match <- regmatches(file, regexec("lbf_chr(\\d+)_([0-9]+)-([0-9]+)", file))
      if (length(match[[1]]) == 4) {
        return(data.frame(
          chromosome = as.integer(match[[1]][2]),
          start_position = as.integer(match[[1]][3]),
          end_position = as.integer(match[[1]][4])
        ))
      } else {
        return(NULL)
      }
    }))
  )
  
  
  # See if required region is covered
  if (nrow(coverage) > 0) {
    filtered_coverage <- subset(coverage, chromosome == lead_chr)
    position_covered <- filtered_coverage[(filtered_coverage$start_position <= required_start) &
                                            (filtered_coverage$end_position >= required_end), ]
  } else {
    position_covered <- coverage
  }
  
  return(list(
    region = list(
      chr = lead_chr,
      pos_high = pos_high,
      pos_low = pos_low
    ),
    is_pos_covered = nrow(position_covered) > 0,
    existing_coverage = position_covered
  ))
}




#' Conditional Fine-Mapping Wrapper Using SuSiE
#' 
#' Checks for existing fine-mapping results in a region and runs SuSiE if needed.
#' Handles LD matrix computation, result caching, and error recovery.
#'
#' @param out Data.frame for target trait containing:
#'   - SNP: Variant IDs
#'   - chr: Chromosome
#'   - pos: Genomic positions
#'   - pval: P-values
#' @param out_lbf_dir_path Directory path for storing/loading Bayes factor (BF) results
#' @param N_out Sample size for outcome trait
#' @param out_type Trait type ("quant" for quantitative, "cc" for case-control)
#' @param out_sd Trait standard deviation (required for quantitative traits)
#' @param trait Optional data.frame for secondary trait to define region (default: NULL)
#' @param plink_path Path to PLINK executable (default: "plink")
#' @param bfile_path Path to LD reference panel (default UKBB path shown)
#' @param temp_dir_path Temporary directory for intermediate files (default: "Temp")
#' 
#' @return Original dataframe merged with:
#'   - lbf_variable1-10: Bayes factors across SuSiE iterations
#'   - variant: Formatted chr_pos identifier
#' 
#' @note
#' Requires:
#' - PLINK installed and accessible via `plink_path`
#' - LD reference panel in PLINK binary format
#' - SuSiE and associated fine-mapping utilities
#' - Internal functions: lbf_in_dataframe, region_utils, finemap_susie
#' 
#' @details
#' Workflow:
#' 1. Checks for existing BFs in dataset
#' 2. Defines genomic region using either:
#'    - Secondary trait coordinates (if provided)
#'    - Outcome trait coordinates (default)
#' 3. Loads pre-computed BFs if region exists in storage
#' 4. Runs SuSiE fine-mapping if:
#'    - No existing results found
#'    - At least one variant with p < 5e-2
#' 5. Handles failed convergence by returning null BFs (-5)
#' 6. Caches new results to avoid recomputation
#' 
#' @author M.Breeur
finemap.wrapper <- function(out, out_lbf_dir_path, N_out, out_type, out_sd,
                            trait = NULL,
                            plink_path = "plink",
                            bfile_path = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
                            temp_dir_path = "Temp") {
  lbf_in_dataframe <- function(df) {
    required_cols <- paste0("lbf_variable", 1:10)
    return(any(required_cols %in% colnames(df)))
  }
  if (lbf_in_dataframe(out)) {
    message("BF are already in the dataset.")
    return(out)
  } else {
    if (!is.null(trait)) {
      lead_pos <- trait %>%
        filter(abs(z) == max(abs(z), na.rm = TRUE)) %>%
        pull(pos) %>%
        mean()
      region <- paste0(
        "chr", as.character(unique(trait$chr)), ":",
        as.character(min(trait$pos)), "-",
        as.character(max(trait$pos))
      )
      region_utils <- region_utils(region, lead_pos, out_lbf_dir_path)
    } else {
      lead_pos <- out %>%
        filter(abs(z) == max(abs(z), na.rm = TRUE)) %>%
        pull(pos) %>%
        mean()
      region <- paste0(
        "chr", as.character(unique(out$chr)), ":",
        as.character(min(out$pos)), "-",
        as.character(max(out$pos))
      )
      region_utils <- region_utils(region, lead_pos, out_lbf_dir_path)
    }
    
    
    if (region_utils$is_pos_covered) {
      message(paste0("Region already fine-mapped. Loading the BFs from ", out_lbf_dir_path))
      out_trait_susie <- NULL
      out_lbf <- readRDS(paste0(
        out_lbf_dir_path, "/lbf_chr",
        as.character(region_utils$existing_coverage$chromosome),
        "_",
        as.character(region_utils$existing_coverage$start_position),
        "-",
        as.character(region_utils$existing_coverage$end_position),
        ".rds"
      ))
    } else {
      message("Region not covered. Computing the BFs...")
      
      ## Find out if anything passes the conventional significance threshold
      
      if (any(out$pval < 5e-2)) {
        ## Run the finemapping
        
        tic("Running SuSiE")
        out_trait_susie <- tryCatch(expr = {
          finemap_susie(
            exp_data = out,
            N_exp = N_out,
            exp_type = out_type,
            exp_sd = out_sd,
            LD_matrix = NULL,
            plink_loc = plink_path,
            bfile_loc = bfile_path,
            temp_dir_path = temp_dir_path,
            max_iter = 10000,
            repeat_until_convergence = FALSE,
            run_checks = FALSE
          )
        }, error = function(e) {
          warning("SuSiE did not converge, returning BF = 0")
          NULL
        })
        toc()
        
        if (!is.null(out_trait_susie)) {
          out_lbf <- as.data.frame(t(out_trait_susie$lbf_variable)) %>%
            setNames(paste0("lbf_variable", 1:10)) %>%
            rownames_to_column("SNP") %>%
            left_join(out %>%
                        transmute(variant = paste0("chr", chr, "_", pos), SNP), by = "SNP")
        } else {
          out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
            setNames(paste0("lbf_variable", 1:10)) %>%
            mutate(
              variant = paste0("chr", out$chr, "_", out$pos),
              SNP = out$SNP
            )
        }
        
        saveRDS(out_lbf, file = paste0(
          out_lbf_dir_path, "/lbf_chr",
          region_utils$region$chr, "_",
          as.character(region_utils$region$pos_low), "-",
          as.character(region_utils$region$pos_high), ".rds"
        ))
      } else {
        message("No variant crosses the significance threshold. Returning BF = 0")
        out_trait_susie <- NULL
        
        out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
          setNames(paste0("lbf_variable", 1:10)) %>%
          mutate(
            variant = paste0("chr", out$chr, "_", out$pos),
            SNP = out$SNP
          )
      }
    }
    return(merge(out, out_lbf))
  }
}



################################################################################
# Coloc wrapper
################################################################################


#' Perform Colocalisation Analysis Using Bayes Factors from SuSiE
#'
#' This function performs colocalisation analysis on a pair of datasets containing Bayes factors 
#' for credible sets computed via SuSiE. It utilizes the `coloc.bf_bf` function from the `coloc` package.
#'
#' @param trait1 A data frame containing Bayes factors for credible sets of the first trait. Must include columns 
#'   prefixed with "lbf_variable" and a "variant" column.
#' @param trait2 A data frame containing Bayes factors for credible sets of the second trait. Must include columns 
#'   prefixed with "lbf_variable" and a "variant" column.
#' @param trait1_name A character string specifying the name of the first trait. Default is "exposure".
#' @param trait2_name A character string specifying the name of the second trait. Default is "outcome".
#'
#' @return A data frame summarizing the colocalisation results, including indices of credible sets and colocalisation method classification.
#'
#' @importFrom dplyr select mutate rename_with
#' @importFrom stringr str_replace
#' @importFrom coloc coloc.bf_bf
#' 
#' @author M. Breeur
coloc.wrapper <- function(trait1, trait2, trait1_name = "exposure", trait2_name = "outcome") {
  
  trait1_mat <- as.matrix(dplyr::select(trait1, starts_with("lbf_variable")))
  row.names(trait1_mat) <- trait1$variant
  # Reorder columns if necessary
  col_names <- colnames(trait1_mat)
  numbers <- as.numeric(gsub("lbf_variable", "", col_names))
  column_order <- order(numbers)
  trait1_mat_ordered <- t(trait1_mat[, column_order])
  
  trait2_mat <- as.matrix(dplyr::select(trait2, starts_with("lbf_variable")))
  row.names(trait2_mat) <- trait2$variant
  # Reorder columns if necessary
  col_names <- colnames(trait2_mat)
  numbers <- as.numeric(gsub("lbf_variable", "", col_names))
  column_order <- order(numbers)
  trait2_mat_ordered <- t(trait2_mat[, column_order])
  
  res <- coloc::coloc.bf_bf(trait2_mat_ordered, trait1_mat_ordered)$summary
  # inverting 1 and 2 because coloc_bf_bf takes entry 1 as hit2 and vice versa
  
  res$idx1 <- res$idx1 - 1
  res$idx2 <- res$idx2 - 1
  
  res <- res %>% mutate(method = ifelse(idx1 == 0 & idx2 == 0, "Vanilla", # coloc between the two non-finemapped signals
                                        ifelse(idx1 == 0 | idx2 == 0, "Hybrid", # coloc between one susie credible set and an original signal
                                               "SuSiE")))%>% # coloc between susie credible sets
    rename_with(~ c(paste0("hit_", str_replace(trait1_name," ", "_") ), 
                    paste0("hit_", str_replace(trait2_name," ", "_"))), 
                c(hit1, hit2))
  return(res)
}



################################################################################
# MR wrapper
################################################################################


#' Perform Mendelian Randomization (MR) If Traits Are Colocalised
#'
#' This function takes colocalisation results and performs Mendelian Randomization (MR) if the two traits are colocalised.
#' It uses the `TwoSampleMR` package to format and harmonise the data before conducting MR analysis.
#'
#' @param trait A data frame containing summary statistics for the exposure trait.
#' @param out A data frame containing summary statistics for the outcome trait.
#' @param N_out An integer specifying the sample size for the outcome trait.
#' @param res_coloc A data frame containing colocalisation results, including a column `PP.H4.abf` for colocalisation probability.
#' @param trait_name A character string specifying the name of the exposure trait. Default is "exposure".
#' @param out_name A character string specifying the name of the outcome trait. Default is "outcome".
#'
#' @return A data frame containing MR results from `mr_singlesnp`, or `NULL` if no colocalisation is detected.
#'
#' @importFrom dplyr filter
#' @importFrom TwoSampleMR format_data harmonise_data mr_singlesnp
#' 
#' @author M. Breeur
mr.wrapper <- function(trait, out, N_out, res_coloc, trait_name = "exposure", out_name = "outcome") {
  colocalised <- (nrow(res_coloc %>% filter(PP.H4.abf > 0.5)) > 0)
  if (colocalised) {
    trait$pheno <- trait_name
    trait_mr <- TwoSampleMR::format_data(as.data.frame(trait), phenotype_col = "pheno") %>%
      filter(pval.exposure < 1e-5)
    out$pheno <- out_name
    out_mr <- TwoSampleMR::format_data(as.data.frame(out),
                                       snps = trait_mr$SNP,
                                       type = "outcome",
                                       phenotype_col = "pheno"
    )
    out_mr$samplesize.outcome <- N_out
    
    harm_dat_mr <- TwoSampleMR::harmonise_data(trait_mr, out_mr)
    
    res_mr_single <- mr_singlesnp(harm_dat_mr, all_method = c())
    
    res_mr_single <- merge(res_mr_single,
                           trait %>% dplyr::select(c(SNP, variant)),
                           on = "SNP")
  } else {
    res_mr_single <- NULL
  }
  
  return(res_mr_single)
}




################################################################################
# Plot functions
################################################################################


#' zz_plot Function
#'
#' This function generates a Z-Z scatter plot for two traits (exposure and outcome), where SNPs (Single Nucleotide Polymorphisms)
#' are color-coded based on their linkage disequilibrium (LD) with a lead SNP. Additionally, the coloc SNP (a key SNP for colocalization analysis) is highlighted.
#'
#' @param LD_Mat Dataframe containing the linkage disequilibrium (LD) matrix with SNPs as row names.
#' @param lead_SNP Character string indicating the lead SNP to be labeled.
#' @param Harm_dat Dataframe containing harmonized data for SNPs, including beta and standard errors for exposure and outcome.
#' @param coloc_SNP Character string specifying the coloc SNP. Default is `NULL`
#' @param exposure_name Character string for the name of the exposure trait (default: unique value from Harm_dat's "exposure" column).
#' @param outcome_name Character string for the name of the outcome trait (default: unique value from Harm_dat's "outcome" column).
#'
#' @return A ggplot object representing the Z-Z plot where SNPs are color-coded based on their LD with the coloc SNP,
#' and the coloc and lead SNPs are highlighted.
#'
#' @examples
#' zz_plot(LD_Mat, lead_SNP = "rs123", Harm_dat, coloc_SNP = "rs456")
#'
#' @author K. Smith-Byrne
#'
zz_plot <- function(LD_Mat, lead_SNP, Harm_dat, coloc_SNP = NULL,
                    exp_name = "exposure",
                    out_name = "outcome") {
  ## Step 1: Prepare LD Data
  # Add a column 'RS_number' to the LD matrix that stores SNP identifiers
  colnames(LD_Mat) <- str_split_fixed(colnames(LD_Mat), "_", 3)[, 1]
  LD_Mat$RS_number <- str_split_fixed(rownames(LD_Mat), "_", 3)[, 1]
  
  # Subset the LD matrix to retain only the lead_SNP column for plotting
  LD_TEMP <- LD_Mat[, c("RS_number", lead_SNP)]
  
  ## Step 2: Merge LD information with harmonized data
  # Merge the harmonized data (Harm_dat) with the LD data, matching SNP identifiers
  temp_dat_format <- as.data.frame(merge(Harm_dat, LD_TEMP, by.x = "SNP", by.y = "RS_number", all.x = TRUE))
  
  # Square the LD values to calculate r^2 for the coloc SNP
  temp_dat_format[, lead_SNP] <- temp_dat_format[, lead_SNP]^2
  
  ## Step 3: Handle missing LD values and define LD color categories
  # Replace NA LD values with 0 for consistency
  temp_dat_format[, lead_SNP] <- ifelse(is.na(temp_dat_format[, lead_SNP]), 0, temp_dat_format[, lead_SNP])
  
  # Create a new 'LD' column categorising LD strength into bins
  temp_dat_format <- temp_dat_format %>%
    mutate(
      !!sym(lead_SNP) := replace_na(!!sym(lead_SNP), 0),
      LD = case_when(
        !!sym(lead_SNP) >= 0 & !!sym(lead_SNP) <= 0.2 ~ "LD < 0.2",
        !!sym(lead_SNP) > 0.2 & !!sym(lead_SNP) <= 0.4 ~ "0.2 < LD < 0.4",
        !!sym(lead_SNP) > 0.4 & !!sym(lead_SNP) <= 0.6 ~ "0.4 < LD < 0.6",
        !!sym(lead_SNP) > 0.6 & !!sym(lead_SNP) <= 0.8 ~ "0.6 < LD < 0.8",
        !!sym(lead_SNP) > 0.8 ~ "LD > 0.8",
        TRUE ~ "No LD"
      )
    )
  
  # Assign the "Lead SNP" label to the coloc SNP
  temp_dat_format$LD <- ifelse(temp_dat_format$SNP == lead_SNP, "Lead SNP", temp_dat_format$LD)
  
  # Create a flag to highlight the coloc SNP
  temp_dat_format <- temp_dat_format %>% mutate(coloc_flag = case_when(SNP == lead_SNP ~ "Yes", TRUE ~ "No"))
  
  ## Step 4: Calculate Z-scores for exposure and outcome
  # Compute Z-scores: beta divided by the standard error
  
  ## Step 5: Generate the Z-Z plot
  # Use jitter to slightly shift points for better visualization
  pos <- position_jitter(width = 0.5, seed = 1)
  
  # Plot the Z-scores with SNPs color-coded by LD category
  plot <- temp_dat_format %>%
    mutate(LD = fct_reorder(LD, get(lead_SNP))) %>% # Reorder LD levels for better display
    ggplot(aes(z.x, z.y, color = LD)) +
    geom_point(size = 2) + # Add scatter points
    theme_minimal() + # Use a clean, black-and-white theme
    xlab(paste(exp_name, " Z-score", sep = "")) +
    ylab(paste(out_name, " Z-score", sep = "")) +
    ggtitle(paste("Z-Z Locus Plot ", exp_name, " VS ", out_name, sep = "")) +
    theme(
      axis.text = element_text(hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position="top"
    ) +
    # Add labels to the lead SNP and coloc SNP
    geom_label_repel(
      size = 5,
      data = temp_dat_format %>% filter(SNP %in% c(lead_SNP, coloc_SNP)),
      aes(label = SNP), show.legend = FALSE
    ) +
    # Customize legend title and colors for LD categories
    labs(color = "LD with Lead SNP") +
    scale_color_manual(values = c(
      "No LD" = "#D3D3D3",
      "LD < 0.2" = "#225EA8",
      "0.2 < LD < 0.4" = "#41B6C4",
      "0.4 < LD < 0.6" = "#7FCDBB",
      "0.6 < LD < 0.8" = "#FE9929",
      "LD > 0.8" = "#8856A7",
      "Lead SNP" = "#F768A1"
    ))
  
  ## Step 6: Return the final plot
  return(plot)
}


#' Generate Locus Plots for GWAS Summary Statistics
#'
#' @param LD_Mat A linkage disequilibrium (LD) matrix containing SNP information.
#' @param harm_dat A data frame containing harmonized GWAS summary statistics.
#' @param coloc_SNP The lead SNP used for colocalization analysis.
#'
#' @return A list of three ggplot objects:
#'         - `pvalues_at_locus`: Scatter plot comparing exposure and outcome -log10 p-values.
#'         - `outcome_at_locus`: Manhattan plot for the outcome GWAS at the locus.
#'         - `exposure_at_locus`: Manhattan plot for the exposure GWAS at the locus.
#'
#' @author K. Smith-Byrne, M. Breeur
#'
#' @export
locus_plot <- function(LD_Mat, harm_dat, lead_SNP, coloc_SNP = NULL, exp_name = "exposure", out_name = "outcome") {
  # Compute -log10(p-value) for exposure and outcome
  harm_dat$log10p.exp <- -pchisq((harm_dat$beta.x / harm_dat$se.x)^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  harm_dat$log10p.out <- -pchisq((harm_dat$beta.y / harm_dat$se.y)^2, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  
  # Extract the lead SNP column for plotting
  colnames(LD_Mat) <- str_split_fixed(colnames(LD_Mat), "_", 3)[, 1]
  LD_Mat$RS_number <- str_split_fixed(rownames(LD_Mat), "_", 3)[, 1]
  LD_TEMP <- LD_Mat[, c("RS_number", lead_SNP)]
  harm_dat <- merge(harm_dat, LD_TEMP, by.x = "SNP", by.y = "RS_number", all.x = TRUE)
  harm_dat[, lead_SNP] <- harm_dat[, lead_SNP]^2 # Square the LD values
  
  # Assign LD category labels
  harm_dat <- harm_dat %>%
    mutate(
      !!lead_SNP := replace_na(!!sym(lead_SNP), 0),
      LD = case_when(
        !!sym(lead_SNP) > 0 & !!sym(lead_SNP) <= 0.2 ~ "LD < 0.2",
        !!sym(lead_SNP) > 0.2 & !!sym(lead_SNP) <= 0.4 ~ "0.2 < LD < 0.4",
        !!sym(lead_SNP) > 0.4 & !!sym(lead_SNP) <= 0.6 ~ "0.4 < LD < 0.6",
        !!sym(lead_SNP) > 0.6 & !!sym(lead_SNP) <= 0.8 ~ "0.6 < LD < 0.8",
        !!sym(lead_SNP) > 0.8 ~ "LD > 0.8",
        TRUE ~ "No LD"
      )
    )
  
  
  # Mark lead SNP in LD category
  harm_dat$LD <- ifelse(harm_dat$SNP == lead_SNP, "Lead SNP", harm_dat$LD)
  harm_dat <- harm_dat %>% mutate(lead_flag = case_when(SNP == lead_SNP ~ "Yes", TRUE ~ "No"))
  
  # Define factor levels for LD categories
  harm_dat$LD <- factor(harm_dat$LD, levels = c("No LD", "LD < 0.2", "0.2 < LD < 0.4", "0.4 < LD < 0.6", "0.6 < LD < 0.8", "LD > 0.8", "Lead SNP"))
  
  # Define LD color mapping
  ld_colors <- c(
    "No LD" = "#D3D3D3", "LD < 0.2" = "#225EA8", "0.2 < LD < 0.4" = "#41B6C4",
    "0.4 < LD < 0.6" = "#7FCDBB", "0.6 < LD < 0.8" = "#FE9929", "LD > 0.8" = "#8856A7", "Lead SNP" = "#F768A1"
  )
  
  # Generate scatter plot for p-values
  p_p_plot <- ggplot(harm_dat, aes(x = log10p.exp, y = log10p.out, color = LD)) +
    geom_point(size = 2) +
    scale_x_continuous(limits = c(0, max(harm_dat$log10p.exp) * 1.05)) +
    labs(
      x = paste0(exp_name, " -log10(p-value)"),
      y = paste0(out_name, " -log10(p-value)"),
      color = "LD with Lead SNP"
    ) +
    theme_minimal() +
    ggtitle(paste0("P-values ", exp_name, " vs ", out_name)) +
    scale_color_manual(values = ld_colors, guide = "none") +
    theme(
      axis.text = element_text(hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position="none"
    )+
    geom_label_repel(data = harm_dat %>% filter(SNP %in% c(lead_SNP, coloc_SNP)), aes(label = SNP), size = 5, show.legend = FALSE)
  
  # Generate Manhattan plots for exposure and outcome
  locus_plot_generator <- function(y_var, name, title_prefix) {
    ggplot(harm_dat, aes(x = pos.x, y = !!sym(y_var), color = LD)) +
      geom_point(size = 2) +
      labs(x = "BP Position", y = "-log10(p-value)", color = "LD with Lead SNP") +
      theme_minimal() +
      ggtitle(paste0(title_prefix, name)) +
      theme(plot.title = element_text(size = 10)) +
      scale_color_manual(values = ld_colors, guide = "none") +
      geom_label_repel(data = harm_dat %>% filter(SNP %in% c(lead_SNP, coloc_SNP)), aes(label = SNP), size = 3, show.legend = FALSE)
  }
  
  locus_out <- locus_plot_generator("log10p.exp", exp_name, "Locus Manhattan plot for ")
  locus_exp <- locus_plot_generator("log10p.out", out_name, "Locus Manhattan plot for ")
  
  # Return a list of plots
  return(list(pvalues_at_locus = p_p_plot, outcome_at_locus = locus_out, exposure_at_locus = locus_exp))
}


#' Generate a Colocalization and Mendelian Randomization Summary Table
#'
#' This function creates a summary table of colocalization and Mendelian Randomization (MR) results,
#' filtering colocalization results where the posterior probability (PP.H4.abf) exceeds 0.5.
#'
#' @param trait A data frame containing SNP-variant mappings for the exposure trait.
#' @param out A data frame containing SNP-variant mappings for the outcome trait.
#' @param res_coloc A data frame containing colocalization results, including `PP.H4.abf` and `method`.
#' @param res_mr A data frame containing Mendelian Randomization results, including `SNP` and `p` values.
#' @param trait_name A character string specifying the name of the exposure trait (default: "exposure").
#' @param out_name A character string specifying the name of the outcome trait (default: "outcome").
#'
#' @return A `ggtexttable` summary table displaying SNP associations with colocalization and MR results.
#' If no colocalization result passes the threshold, a message table is returned.
#' 
#' @import dplyr ggpubr stringr
#' @author M.Breeur
coloc_mr_table <- function(trait, out, res_coloc, res_mr, trait_name = "exposure", out_name = "outcome") {
  
  # Filter colocalization results where PP.H4.abf > 0.5 and select relevant columns
  temp <- res_coloc %>%
    filter(PP.H4.abf > 0.5) %>%
    select(
      !!sym(paste0("hit_", str_replace_all(trait_name, " ", "_"))), 
      !!sym(paste0("hit_", str_replace_all(out_name, " ", "_"))), 
      PP.H4.abf, method
    ) %>%
    distinct() %>%
    rename(
      Coloc_method = method,
      PPH4 = PP.H4.abf
    )
  
  # Join with trait SNP data
  temp <- temp %>%
    left_join(
      trait %>% select(SNP, variant), 
      by = setNames("variant", paste0("hit_", trait_name))
    ) %>%
    rename(!!sym(paste0("SNP_", trait_name)) := SNP)
  
  # Join with outcome SNP data
  temp <- temp %>%
    left_join(
      out %>% select(SNP, variant), 
      by = setNames("variant", paste0("hit_", out_name))
    ) %>%
    rename(!!sym(paste0("SNP_", out_name)) := SNP)
  
  # Remove original hit columns
  temp <- temp %>%
    select(-all_of(c(
      paste0("hit_", str_replace_all(trait_name, " ", "_")), 
      paste0("hit_", str_replace_all(out_name, " ", "_"))
    )))
  
  # Add MR results if colocalization data exists
  if (nrow(temp) > 0){
    temp <- temp %>%
      left_join(
        res_mr %>% select(SNP, p, b, se),
        by = setNames("SNP", paste0("SNP_", trait_name))
      ) %>% mutate(MR_beta = sprintf("%.2f (%.2f, %.2f)", 
                                     b, 
                                     b - 1.96 * se, 
                                     b + 1.96 * se)) %>%
      rename(MR_pvalue = p) %>% dplyr::select(-c(b, se))
    
    # Define column order
    cols <- c(
      paste0("SNP_", trait_name),
      paste0("SNP_", out_name),
      "Coloc_method",
      "PPH4",
      "MR_beta",
      "MR_pvalue"
    )
    
    temp <- temp %>% select(all_of(cols))
  }
  
  
  # Generate table output
  if (nrow(temp) > 0) {
    table <- ggtexttable(temp, rows = NULL, theme = ttheme("blank")) %>%
      tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
      tab_add_hline(at.row = nrow(temp) + 1, row.side = "bottom", linewidth = 2)
  } else {
    table <- ggtexttable(
      c("No colocalisation method crossed the 0.5 threshold"),
      theme = ttheme("blank")
    ) %>%
      tab_add_hline(at.row = 1, row.side = "top", linewidth = 2) %>%
      tab_add_hline(at.row = 1, row.side = "bottom", linewidth = 2)
  }
  
  return(table)
}


#' Generate Combined Z-Z Plot and Locus Plots for Trait Comparison
#' 
#' Creates a multi-panel visualization containing:
#' - Z-Z plot showing effect direction correlation
#' - Locus plots for both exposure and outcome traits
#' - Colocalization highlights when applicable
#' 
#' @param trait Data.frame for primary trait (exposure) containing:
#'   - SNP: Variant IDs
#'   - pos: Genomic positions
#'   - z: Z-scores
#' @param out Data.frame for secondary trait (outcome) with same columns as `trait`
#' @param res_coloc Data.frame from coloc analysis containing PP.H4.abf column
#' @param trait_name Name for exposure trait (default: "exposure")
#' @param out_name Name for outcome trait (default: "outcome")
#' @param plink_path Path to PLINK executable (default: "plink")
#' @param bfile_path Path to LD reference panel files (default UKBB path shown)
#' @param temp_dir_path Temporary directory for intermediate files (default: "Temp")
#' 
#' @return A ggarrange object containing:
#' - Top panel: Z-Z plot
#' - Bottom panel: 
#'   - Left: Combined p-value tracks
#'   - Right: Exposure/outcome locus plots
#' 
#' @note
#' Requires:
#' - PLINK installed and accessible via `plink_path`
#' - LD reference panel in PLINK binary format
#' - ggpubr, grid, and coloc packages
#' - Internal functions: get_ld_matrix_from_bim, align_to_LD, locus_plot, zz_plot
#' 
#' @details
#' Workflow:
#' 1. Identifies lead SNP (±500kb window) from exposure trait
#' 2. Extracts LD matrix using reference panel
#' 3. Aligns effect directions between traits using LD
#' 4. Highlights colocalized SNPs when PP.H4.abf > 0.5
#' 5. Arranges plots in publication-ready layout
#' 
#' @author M. Breeur
plot.wrapper <- function(trait, out, res_coloc, res_mr, trait_name = "exposure", out_name = "outcome",
                         plink_path = "plink",
                         bfile_path = "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered",
                         temp_dir_path = "Temp") {
  # Define plot window
  lead_pos <- trait %>%
    filter(abs(z) == max(abs(z), na.rm = TRUE)) %>%
    pull(pos) %>%
    mean()
  
  required_start <- lead_pos - 500000
  required_end <- lead_pos + 500000
  
  # Extract LD matrix
  SNP_list <- trait$SNP[between(trait$pos, required_start, required_end)]
  
  LD_matrix <- get_ld_matrix_from_bim(SNP_list,
                                      plink_loc = plink_path,
                                      bfile_loc = bfile_path,
                                      with_alleles = T,
                                      temp_dir_path = temp_dir_path
  )
  
  ## Extract lead snp from trait
  
  lead_snp <- unique(trait$SNP[which.min(abs(trait$pos - lead_pos))])
  
  ## Flip z scores if needed, according to LD_mat
  
  harm_dat <- merge(
    align_to_LD(trait, LD_matrix),
    align_to_LD(out, LD_matrix),
    by = c("variant", "SNP")
  )
  
  ## Make locus and zz plots, store them in list
  
  # Highlight coloc_snp if colocalised
  colocalised <- (nrow(res_coloc %>% filter(PP.H4.abf > 0.5)) > 0)
  if (colocalised) {
    # coloc_snp = hit in hit_out_name
    vars <- as.data.frame(res_coloc) %>% dplyr::select(!!sym(paste0("hit_", str_replace(out_name," ","_"))))
    coloc_snp <- out$SNP[out$variant == vars[which.max(res_coloc$PP.H4.abf)]]
  } else {
    coloc_snp <- NULL
  }
  
  plot_list <- locus_plot(as.data.frame(LD_matrix),
                          as.data.frame(harm_dat),
                          lead_SNP = lead_snp,
                          coloc_SNP = coloc_snp,
                          exp_name = trait_name,
                          out_name = out_name
  )
  
  
  plot_list$zzplot <- zz_plot(
    LD_Mat = as.data.frame(LD_matrix),
    lead_SNP = lead_snp,
    coloc_SNP = coloc_snp,
    Harm_dat = harm_dat,
    exp_name = trait_name,
    out_name = out_name
  )
  
  plot_list$table <- coloc_mr_table(
    trait, out, res_coloc, res_mr, 
    trait_name, out_name)
  
  return(ggarrange(ggarrange(plot_list$zzplot, plot_list$table, ncol = 2, widths = c(1, .8)),
                   ggarrange(plot_list$pvalues_at_locus,
                             ggarrange(plot_list$exposure_at_locus,
                                       plot_list$outcome_at_locus,
                                       nrow = 2),
                             ncol = 2, widths = c(1, .8)),
                   nrow = 2))
  
  return(plot_list)
}
