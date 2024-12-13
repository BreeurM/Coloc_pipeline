#' SNP Proxy and Colocalisation Analysis Toolkit
#'
#' This script provides a suite of functions for genomic analysis, including identifying proxy SNPs,
#' computing linkage disequilibrium (LD) matrices, extracting genomic regions for colocalisation,
#' and performing colocalisation analysis using the coloc + SuSiE algorithm. The toolkit leverages PLINK
#' TwoSampleMR, and coloc for processing of SNP data, LD computation, and colocalisation analysis
#'
#' ## Key Functions:
#' 1. `find_proxy_snps`: Identifies proxy SNPs within a specified genomic window based on LD thresholds.
#' 2. `get_ld_matrix`: Computes and formats an LD matrix for a set of SNPs using PLINK.
#' 3. `extract_regions_for_coloc`: Extracts genomic regions around SNPs for colocalisation analysis.
#' 4. `main_coloc`: Performs colocalisation analysis using SuSiE and LD fine-mapping.
#'
#' ## Requirements:
#' - PLINK installed and accessible via command line.
#' - `ieugwasr`, `TwoSampleMR`, and `coloc` R packages.
#'
#' ## Notes:
#' - Input data must be formatted correctly for PLINK and TwoSampleMR.
#' - Ensure compatibility of allele information across datasets to avoid errors during harmonization.
#' - Temporary files generated during execution are automatically cleaned up.
#'
#' ## Outputs:
#' - Proxy SNPs for a given target SNP.
#' - LD matrices for analysis and plotting.
#' - Colocalisation results with posterior probabilities and fine-mapping outputs.


#' Read and format exposure or outcome data
#' Checks and organises columns for use with MR or enrichment tests.
#' Infers p-values when possible from beta and se.
#'
#' @param dat Data frame. Must have header with at least SNP column present.
#' @param type Is this the exposure or the outcome data that is being read in? The default is `"exposure"`.
#' @param snps SNPs to extract. If NULL then doesn't extract any and keeps all. The default is `NULL`.
#' @param header The default is `TRUE`.
#' @param phenotype_col Optional column name for the column with phenotype name corresponding the the SNP. If not present then will be created with the value `"Outcome"`. The default is `"Phenotype"`.
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
format_dat <- function(dat, type = "exposure", snps = NA,
                        phenotype_col = NA, snp_col = NA,
                        beta_col = NA, se_col = NA, eaf_col = NA,
                        effect_allele_col = NA, other_allele_col = NA,
                        pval_col = NA, units_col = NA, ncase_col = NA,
                        ncontrol_col = NA, samplesize_col = NA,
                        gene_col = NA, id_col = NA, min_pval = 1e-200,
                        z_col = NA, info_col = NA, chr_col = NA,
                        pos_col = NA, log_pval = FALSE) {
  # Check for required SNP column
  if (!snp_col %in% names(dat)) {
    stop("SNP column not found")
  }

  # Select and standardize relevant columns
  all_cols <- c(
    phenotype_col, snp_col, beta_col, se_col, eaf_col,
    effect_allele_col, other_allele_col, pval_col, units_col,
    ncase_col, ncontrol_col, samplesize_col, gene_col, id_col,
    z_col, info_col, chr_col, pos_col
  )

  dat <- dat %>%
    select(any_of(all_cols)) %>%
    rename_with(~"SNP", all_of(snp_col)) %>%
    mutate(SNP = tolower(SNP) %>% str_replace_all("[[:space:]]", "")) %>%
    filter(!is.na(SNP))

  # Filter SNPs if provided
  if (!is.null(snps)) {
    dat <- dat %>% filter(SNP %in% snps)
  }

  # Add or rename phenotype column
  dat <- dat %>% mutate(!!type := if_else(phenotype_col %in% names(dat), phenotype_col, type))
  if (phenotype_col %in% names(dat) && phenotype_col != type) {
    dat <- dat %>% select(-all_of(phenotype_col))
  }

  # Convert log p-values if needed
  if (log_pval && pval_col %in% names(dat)) {
    dat <- dat %>% mutate(!!pval_col := 10^-(.data[[pval_col]]))
  }

  # Remove duplicated SNPs
  dat <- dat %>% distinct(SNP, .keep_all = TRUE)

  # Check and clean columns for MR
  required_cols <- c(beta_col, se_col, effect_allele_col)
  optional_cols <- c(other_allele_col, eaf_col)
  
  if (!all(required_cols %in% names(dat))) {
    warning("Missing required columns for MR analysis: ", paste(setdiff(required_cols, names(dat)), collapse = ", "))
    dat <- dat %>% mutate(mr_keep = FALSE)
  } else {
    dat <- dat %>% mutate(mr_keep = rowSums(across(all_of(required_cols), ~ !is.na(.))) == length(required_cols))
  }

  # Infer p-values if missing
  if (!pval_col %in% names(dat) && all(c(beta_col, se_col) %in% names(dat))) {
    dat <- dat %>% mutate(!!pval_col := 2 * pnorm(abs(.data[[beta_col]]) / .data[[se_col]], lower.tail = FALSE))
  }

  # Generate sample size if missing
  if (!samplesize_col %in% names(dat) && all(c(ncase_col, ncontrol_col) %in% names(dat))) {
    dat <- dat %>% mutate(!!samplesize_col := .data[[ncase_col]] + .data[[ncontrol_col]])
  }

  # Finalize column names for output
  dat <- dat %>% rename_with(~ paste0(.x, ".", type), setdiff(names(dat), c("SNP", type)))

  return(dat)
}




#' Identify Proxy SNPs within a Defined Genomic Window
#'
#' This function identifies proxy SNPs for a given target SNP within a specified
#' genomic window (in kilobases) using PLINK. It calculates linkage disequilibrium (LD)
#' for the SNPs in the window and selects potential proxies based on a specified LD threshold.
#'
#' @param plink_path     A string specifying the path to the PLINK executable. Default is `"plink"`.
#' @param bfile_prefix   A string specifying the prefix of the PLINK binary files (.bed, .bim, .fam).
#' @param snp_dat        A data frame containing information about the target SNP. Must include at least the column `SNP`.
#' @param window_size_kb An integer specifying the size of the genomic window around the target SNP in kilobases.
#' @param outcome_dat    A data frame containing outcome data for SNPs. Should include at least a `SNP` column.
#' @param file_list      A character vector of file paths containing additional data relevant to the analysis.
#'
#' @return A data frame containing information about the best proxy SNP for the target SNP,
#' including annotation and metadata. If no proxies are identified, returns an empty data frame.
#'
#' @details
#' The function uses PLINK to extract SNPs within the specified genomic window around the target SNP.
#' It computes the LD matrix using the `get_ld_matrix` function and filters SNPs with an LD value
#' greater than 0.8. Additional filters are applied to match SNPs against the provided outcome
#' and exposure datasets. The best proxy SNP is selected based on the highest LD value.
#'
#' Temporary files created during the analysis are automatically cleaned up. If the specified SNP
#' is not found in the data, the function returns an empty data frame with the same structure as `snp_dat`.
find_proxy_snps <- function(plink_path = "plink",
                            bfile_prefix = "N:/EPIC_genetics/1000G_EUR/1000G_EUR/QC_1000G_P3",
                            snp_dat,
                            window_size_kb,
                            outcome_dat,
                            file_list) {
  # Check if PLINK is installed
  if (system(paste(plink_path, "--version"), ignore.stdout = TRUE) != 0) {
    stop("PLINK not found at the specified path. Please provide the correct path to the PLINK executable.")
  }

  cat("Extracting SNPs in a", window_size_kb, "kb window around", snp_dat$SNP, "\n")

  # Execute PLINK command
  extracted_snps_file <- tempfile("extracted_snps", fileext = ".txt")
  cmd <- paste(
    plink_path, "--bfile", bfile_prefix,
    "--snp", snp_dat$SNP,
    "--window", window_size_kb,
    "--write-snplist",
    "--out", extracted_snps_file,
    "--allow-no-sex"
  )

  # Run the command and capture both stdout and stderr
  system(cmd, intern = TRUE, ignore.stderr = FALSE)

  # Read extracted SNPs list with tryCatch
  extracted_snps <- tryCatch(
    {
      read.table(paste(extracted_snps_file, ".snplist", sep = ""), header = FALSE, col.names = "SNP")
    },
    error = function(e) {
      message("An error occurred while reading the extracted SNPs file. Returning an empty data.frame.")
      data.frame()
    }
  )

  if (dim(extracted_snps)[1] > 0) {
    cat("Extracted", dim(extracted_snps)[1], "SNPs around", snp_dat$SNP, "\n")

    # Compute LD matrix
    LD_Return <- get_ld_matrix(extracted_snps$SNP, plink_loc = plink_path, bfile_loc = bfile_prefix)
    LD_Mat <- as.data.frame(abs(LD_Return$LD_Anal))

    # Identify proxy SNPs
    LD_Return_index <- tryCatch(
      expr = {
        LD_Mat %>%
          dplyr::select(contains(snp_dat$SNP)) %>%
          filter_at(1, all_vars(. > 0.8)) %>%
          mutate(SNP = row.names(.))
      },
      error = function(e) {
        message("Query SNP doesn't pass QC")
        data.frame()
      }
    )

    cat("There are", dim(LD_Return_index)[1], "potential proxy SNPs for", snp_dat$SNP, "\n")

    # Process and filter additional data
    if (length(paste0(file_list[grepl(snp_dat$SNP, file_list) &
      grepl(str_replace_all(paste0(snp_dat$gene.exposure, "_"), "-", "-"), file_list)])) > 0) {
      temp_prot <- TwoSampleMR::read_exposure_data(
        paste0(file_list[grepl(snp_dat$SNP, file_list) &
          grepl(str_replace_all(paste0(snp_dat$gene.exposure, "_"), "-", "-"), file_list)]),
        sep = ",",
        snp_col = "RSID", beta_col = "BETA", se_col = "SE", log_pval = TRUE, pval_col = "LOG10P",
        effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
        pos_col = "GENPOS", chr_col = "CHROM", eaf_col = "A1FREQ", min_pval = NA
      ) %>% filter(SNP %in% LD_Return_index$SNP)

      temp_prot$exposure <- snp_dat$Assay
      temp_out <- outcome_dat %>% filter(SNP %in% temp_prot$SNP)

      cat("There are", dim(temp_out)[1], "potential proxy SNPs in the outcome data for", snp_dat$SNP, "\n")
      if (dim(LD_Return_index)[1] > 0) {
        top_snp <- LD_Return_index %>%
          filter(SNP %in% temp_out$SNP) %>%
          arrange(desc(.[[1]])) %>%
          slice(1) %>%
          select(SNP) %>%
          as.character()

        cat(top_snp, "is the proxy SNP for", snp_dat$SNP, "\n")

        new_exp_dat <- temp_prot %>%
          filter(SNP == top_snp) %>%
          mutate(
            gene.exposure = snp_dat$gene.exposure,
            info.exposure = snp_dat$info.exposure,
            id.exposure = snp_dat$id.exposure,
            exposure = snp_dat$exposure,
            query_snp = snp_dat$SNP
          ) %>%
          select(all_of(c(colnames(snp_dat), "query_snp")))
        return(new_exp_dat)
      } else {
        cat("No proxies in our data: ", snp_dat$SNP, "\n")
        return(data.frame(matrix(ncol = length(colnames(snp_dat)), nrow = 0, dimnames = list(NULL, colnames(snp_dat)))))
      }
    }
  } else {
    cat("SNP not in our data: ", snp_dat$SNP, "\n")
    return(data.frame(matrix(ncol = length(colnames(snp_dat)), nrow = 0, dimnames = list(NULL, colnames(snp_dat)))))
  }
}



#' Generate a Linkage Disequilibrium (LD) Matrix
#'
#' This function calculates a linkage disequilibrium (LD) matrix for a given list of RSIDs
#' using PLINK and a provided bfile. It also formats the LD matrix and returns it in two forms:
#' a cleaned matrix for analysis and a data frame for plotting.
#'
#' @param rsid_list    A character vector of RSIDs for which to compute the LD matrix.
#' @param plink_loc    A string specifying the path to the PLINK binary.
#' @param bfile_loc    A string specifying the path to the PLINK bfile (prefix of binary files).
#' @param with_alleles A logical value indicating whether allele information should be included
#' in the LD matrix. Default is `FALSE`.
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
get_ld_matrix <- function(rsid_list, plink_loc, bfile_loc, with_alleles = F) {
  # Set the temporary directory
  temp_dir <- "Temp_dir/"
  tempdir(temp_dir)

  # Compute the LD matrix using the ieugwasr package
  LD_Full <- ieugwasr::ld_matrix_local(rsid_list,
    bfile = bfile_loc,
    plink_bin = plink_loc,
    with_alleles = with_alleles
  )

  # Clean up temporary files
  temp_dir <- tempdir()
  files <- list.files(temp_dir)
  if (length(files) > 0) {
    unlink(file.path(temp_dir, files), recursive = TRUE)
  }

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



##' Function extracting the desired windows to colocalise
##' @param exp_data: either region to colocalise formatted with TwoSampleMr,
##'                  or str - path to exposure data in csv
##'
##' @param out_data: either outcome data formatted with TwoSampleMr
##'                  suppose that the extracted snps correspond to exposure snps
##'                  or path to outcome data in csv
##'
##' @param window_size: window_size around SNP of interest, default = 1000kb

extract_regions_for_coloc <- function(exp_data, out_data, window_size = 1000) {
  # For now assume exp/out data stored in csv, TO BE CHANGED

  if (is.character(exp_data)) {
    exp_data <- read.csv(exp_data)
  }
  if (is.character(out_data)) {
    out_data <- read.csv(out_data)
  }

  # Keep region of interest

  ## Find set of leading SNPs in exp data
  ## Take all SNPs within window
  ## Exceptions: window too large (exceeds range in the file)
  ##             lead SNPs too distant (find out how to quantify)

  extracted_snps <- # list of rsID to keep


    exp_data <- exp_data %>% filter(SNP %in% extracted_snps)
  out_data <- out_data %>% filter(SNP %in% extracted_snps)

  return(list(exp_data = exp_data, out_data = out_data))
}



#' Perform Colocalization Analysis with SuSiE
#'
#' This function performs colocalization analysis using the SuSiE method for a given pair of exposure and outcome datasets.
#' It computes Bayes factors for each dataset, aligns the SNPs and alleles with an optional LD matrix, and integrates the results using SuSiE.
#'
#' @param exp_data      Either a data frame or a file path to the exposure dataset.
#'                      The dataset must be formatted as required by TwoSampleMR
#'                      and include at least the columns: `SNP`, `beta`, `se`, and `eaf`.
#' @param N_exp         An integer specifying the sample size for the exposure dataset.
#' @param exp_type      A string specifying the type of exposure study.
#'                      Either `"quant"` (quantitative trait) or `"cc"` (case-control).
#' @param exp_sd        A numeric value specifying the std of the trait for quantitative studies or the proportion of cases for case-control studies. Default is `1`.
#'
#' @param out_data      Either a data frame or a file path to the outcome dataset.
#'                      The dataset must be formatted as required by TwoSampleMR
#'                      and include at least the columns: `SNP`, `beta`, `se`, and `eaf`.
#' @param N_out         An integer specifying the sample size for the outcome dataset.
#' @param out_type      A string specifying the type of outcome study.
#'                      Either `"quant"` (quantitative trait) or `"cc"` (case-control).
#' @param out_sd        An numeric value specifying the std of the trait for quantitative studies or the proportion of cases for case-control studies. Default is `1`.
#'
#' @param LD_matrix     Either `NULL`, a precomputed LD matrix, or a file path to a CSV file containing the LD matrix.
#'                      If `NULL`, the LD matrix is computed from 1000G using the `get_ld_matrix` function.
#'
#' @param exp_BF_column An optional string specifying the column name containing Bayes factors in the exposure dataset.
#'                      If `NULL`, Bayes factors are computed within the function.
#' @param out_BF_column An optional string specifying the column name containing Bayes factors in the outcome dataset.
#'                      If `NULL`, Bayes factors are computed within the function.
#'
#' @return A SuSiE colocalization result object (`susie.res`) containing information on shared genetic signals between the exposure and outcome datasets.
#'         If SuSiE fails to converge for either dataset, the function returns `NULL`.
main_coloc <- function(exp_data, N_exp, exp_type, exp_sd = 1,
                       out_data, N_out, out_type, out_sd = 1,
                       LD_matrix = NULL,
                       exp_BF_column = NULL,
                       out_BF_column = NULL) {
  # Check input consistency

  # Import exposure and outcome data if provided as file paths
  if (is.character(exp_data)) {
    exp_data <- read.csv(exp_data)
  }
  if (is.character(out_data)) {
    out_data <- read.csv(out_data)
  }

  # Harmonise the exposure and outcome data to ensure consistent SNP and allele representation
  harm_data <- harmonise_data(exp_data, out_data)
  if (any(harm_data$effect_allele.exposure != harm_data$effect_allele.outcome) |
    any(harm_data$other_allele.exposure != harm_data$other_allele.outcome)) {
    stop("Inconsistent effect alleles for exposure and outcome after harmonisation.")
  }
  harm_data <- harm_data %>%
    mutate(
      effect_allele = effect_allele.exposure,
      other_allele = other_allele.exposure
    ) %>%
    select(-c(
      effect_allele.exposure,
      other_allele.exposure,
      effect_allele.outcome,
      other_allele.outcome
    ))

  # Import LD_matrix if provided as file paths, or import from 1000G if NULL
  if (is.null(LD_matrix)) {
    LD_matrix <- get_ld_matrix(harm_data$SNP, plink_loc = plink_loc, bfile_loc = bfile_loc, with_alleles = T)$LD_Anal
  } else {
    if (is.character(LD_matrix)) {
      LD_matrix <- read.csv(LD_matrix)
    }
  }

  # Check for allele information in the LD matrix, and align if possible
  if (all(str_split_fixed(colnames(LD_matrix), "_", 3)[, 2] %in% c("A", "C", "T", "G")) &
    all(str_split_fixed(colnames(LD_matrix), "_", 3)[, 3] %in% c("A", "C", "T", "G"))) {
    # SNP names in the LD_matrix do not contain allele information, or their format is wrong
    warning("LD matrix has no or incorrect allele information, allele alignment will be inferred from expected Zscores VS observed Zscores")
    harm_data <- harm_data %>%
      mutate(pos = pos.exposure) %>%
      select(
        SNP, pos,
        beta.exposure, beta.outcome,
        se.exposure, se.outcome,
        eaf.exposure, eaf.outcome
      )
  } else {
    # Flip alleles in harmonised data to align with LD matrix
    LD_alignment <- data.frame(
      SNP = str_split_fixed(colnames(LD_matrix), "_", 3)[, 1],
      LD_A1 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 2],
      LD_A2 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 3]
    )
    temp <- merge(harm_data, LD_alignment)
    temp <- temp %>%
      mutate(
        flipped = effect_allele != LD_A1,
        effect_allele = if_else(flipped, other_allele, effect_allele),
        effect_allele = if_else(flipped, effect_allele, other_allele),
        beta.exposure = if_else(flipped, -beta.exposure, beta.exposure),
        beta.outcome = if_else(flipped, -beta.outcome, beta.outcome),
        eaf.exposure = if_else(flipped, 1 - eaf.exposure, eaf.exposure),
        eaf.outcome = if_else(flipped, 1 - eaf.outcome, eaf.outcome)
      )

    if (any(temp$LD_A1 != temp$effect_allele)) {
      # Mismatched alleles after the flip, LD alleles and data alleles do not correspond
      stop("Could not flip the alleles to match LD matrix. Check that the allele info provided is correct.")
    }

    harm_data <- temp %>%
      mutate(pos = pos.exposure) %>%
      select(
        SNP, pos,
        beta.exposure, beta.outcome,
        se.exposure, se.outcome,
        eaf.exposure, eaf.outcome
      )
  }

  # Prepare data for SuSiE colocalization analysis

  # Remove alleles from LD matrix column names for cross ref with ext/out data
  colnames(LD_matrix) <- str_split_fixed(colnames(LD_matrix), "_", 2)[, 1]
  rownames(LD_matrix) <- colnames(LD_matrix)

  exp_for_coloc <- harm_data %>%
    select(
      beta.exposure,
      se.exposure,
      SNP,
      eaf.exposure
    )
  exp_for_coloc$se.exposure <- exp_for_coloc$se.exposure^2
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


  out_for_coloc <- harm_region %>%
    select(
      beta.outcome,
      se.outcome,
      SNP,
      eaf.outcome
    )
  out_for_coloc$se.outcome <- out_for_coloc$se.outcome^2
  colnames(out_for_coloc) <- c("beta", "varbeta", "snp", "MAF")
  out_for_coloc <- as.list(out_for_coloc)
  out_for_coloc$type <- out_type
  if (out_type == "quant") {
    out_for_coloc$sdY <- out_sd
  } else {
    out_for_coloc$s <- out_sd
  }
  out_for_coloc$N <- N_out
  out_for_coloc$LD <- LD_matrix[out_for_coloc$snp, out_for_coloc$snp]

  # Run vanilla colocalisation as a warm up
  ABF <- coloc.abf(exp_for_coloc, out_for_coloc)

  # Flags for faulty allele alignment, TBC
  exp_alignment_qc <- list(
    kriging = kriging_rss(exp_for_coloc$beta / sqrt(exp_for_coloc$varbeta),
      exp_for_coloc$LD,
      n = N_exp
    ),
    alignment_check = check_alignment(exp_for_coloc)
  )

  out_alignment_qc <- list(
    kriging = kriging_rss(out_for_coloc$beta / sqrt(out_for_coloc$varbeta),
      out_for_coloc$LD,
      n = N_out
    ),
    alignment_check = check.alignment(out_for_coloc)
  )

  # Run SuSiE for fine-mapping and colocalization
  exp_susie <- tryCatch(
    expr = {
      temp <- coloc::runsusie(exp_for_coloc,
        repeat_until_convergence = F,
        maxit = 1000
      )
    },
    error = function(e) {
      warning("SuSiE did not converge for exposure")
      NULL
    }
  )

  out_susie <- tryCatch(
    expr = {
      temp <- coloc::runsusie(out_for_coloc,
        repeat_until_convergence = F,
        maxit = 1000
      )
    },
    error = function(e) {
      warning("SuSiE did not converge for outcome")
      NULL
    }
  )
  if (!is.null(exp_susie) & !is.null(out_susie)) {
    susie.res <- coloc::coloc.susie(exp_susie, out_susie)
  } else {
    warning("Returning NULL SuSiE object")
    susie.res <- NULL
  }

  return(list(
    susie.res = susie.res,
    abf.res = ABF,
    exp_susie = exp_susie,
    out_susie = out_susie
  ))
}





#' Perform fine mapping with SuSiE on one dataset
#'
#' This function implements the SuSiE method for a given dataset amd outputs credible sets and Bayes Factors.
#'
#' @param exp_data      Either a data frame or a file path to the exposure dataset.
#'                      The dataset must be formatted as required by TwoSampleMR
#'                      and include at least the columns: `SNP`, `beta`, `se`, and `eaf`.
#' @param N_exp         An integer specifying the sample size for the exposure dataset.
#' @param exp_type      A string specifying the type of exposure study.
#'                      Either `"quant"` (quantitative trait) or `"cc"` (case-control).
#' @param exp_sd        A numeric value specifying the std of the trait for quantitative studies or the proportion of cases for case-control studies. Default is `1`.
#'
#' @param LD_matrix     Either `NULL`, a precomputed LD matrix, or a file path to a CSV file containing the LD matrix.
#'                      If `NULL`, the LD matrix is computed from 1000G using the `get_ld_matrix` function.
#' @return A SuSiE object (`susie.res`) containing information on shared genetic signals between the exposure and outcome datasets.
#'         If SuSiE fails to converge, the function returns `NULL`.
finemap_susie <- function(exp_data, N_exp, exp_type, exp_sd = 1,
                          LD_matrix = NULL) {
  # Check input consistency

  # Import exposure and outcome data if provided as file paths
  if (is.character(exp_data)) {
    exp_data <- read.csv(exp_data)
  }

  # Harmonise the exposure and outcome data to ensure consistent SNP and allele representation
  harm_data <- harm_data %>%
    mutate(
      effect_allele = effect_allele.exposure,
      other_allele = other_allele.exposure
    ) %>%
    select(-c(
      effect_allele.exposure,
      other_allele.exposure,
      effect_allele.outcome,
      other_allele.outcome
    ))

  # Import LD_matrix if provided as file paths, or import from 1000G if NULL
  if (is.null(LD_matrix)) {
    LD_matrix <- get_ld_matrix(harm_data$SNP, plink_loc = plink_loc, bfile_loc = bfile_loc, with_alleles = T)$LD_Anal
  } else {
    if (is.character(LD_matrix)) {
      LD_matrix <- read.csv(LD_matrix)
    }
  }

  # Check for allele information in the LD matrix, and align if possible
  if (all(str_split_fixed(colnames(LD_matrix), "_", 3)[, 2] %in% c("A", "C", "T", "G")) &
    all(str_split_fixed(colnames(LD_matrix), "_", 3)[, 3] %in% c("A", "C", "T", "G"))) {
    # SNP names in the LD_matrix do not contain allele information, or their format is wrong
    warning("LD matrix has no or incorrect allele information, allele alignment will be inferred from expected Zscores VS observed Zscores")
    harm_data <- harm_data %>%
      mutate(pos = pos.exposure) %>%
      select(
        SNP, pos,
        beta.exposure, beta.outcome,
        se.exposure, se.outcome,
        eaf.exposure, eaf.outcome
      )
  } else {
    # Flip alleles in harmonised data to align with LD matrix
    LD_alignment <- data.frame(
      SNP = str_split_fixed(colnames(LD_matrix), "_", 3)[, 1],
      LD_A1 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 2],
      LD_A2 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 3]
    )
    temp <- merge(harm_data, LD_alignment)
    temp <- temp %>%
      mutate(
        flipped = effect_allele != LD_A1,
        effect_allele = if_else(flipped, other_allele, effect_allele),
        effect_allele = if_else(flipped, effect_allele, other_allele),
        beta.exposure = if_else(flipped, -beta.exposure, beta.exposure),
        beta.outcome = if_else(flipped, -beta.outcome, beta.outcome),
        eaf.exposure = if_else(flipped, 1 - eaf.exposure, eaf.exposure),
        eaf.outcome = if_else(flipped, 1 - eaf.outcome, eaf.outcome)
      )

    if (any(temp$LD_A1 != temp$effect_allele)) {
      # Mismatched alleles after the flip, LD alleles and data alleles do not correspond
      stop("Could not flip the alleles to match LD matrix. Check that the allele info provided is correct.")
    }

    harm_data <- temp %>%
      mutate(pos = pos.exposure) %>%
      select(
        SNP, pos,
        beta.exposure, beta.outcome,
        se.exposure, se.outcome,
        eaf.exposure, eaf.outcome
      )
  }

  # Prepare data for SuSiE colocalization analysis

  # Remove alleles from LD matrix column names for cross ref with ext/out data
  colnames(LD_matrix) <- str_split_fixed(colnames(LD_matrix), "_", 2)[, 1]
  rownames(LD_matrix) <- colnames(LD_matrix)

  exp_for_coloc <- harm_data %>%
    select(
      beta.exposure,
      se.exposure,
      SNP,
      eaf.exposure
    )
  exp_for_coloc$se.exposure <- exp_for_coloc$se.exposure^2
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


  out_for_coloc <- harm_region %>%
    select(
      beta.outcome,
      se.outcome,
      SNP,
      eaf.outcome
    )
  out_for_coloc$se.outcome <- out_for_coloc$se.outcome^2
  colnames(out_for_coloc) <- c("beta", "varbeta", "snp", "MAF")
  out_for_coloc <- as.list(out_for_coloc)
  out_for_coloc$type <- out_type
  if (out_type == "quant") {
    out_for_coloc$sdY <- out_sd
  } else {
    out_for_coloc$s <- out_sd
  }
  out_for_coloc$N <- N_out
  out_for_coloc$LD <- LD_matrix[out_for_coloc$snp, out_for_coloc$snp]

  # Run vanilla colocalisation as a warm up
  ABF <- coloc.abf(exp_for_coloc, out_for_coloc)

  # Flags for faulty allele alignment, TBC
  exp_alignment_qc <- list(
    kriging = kriging_rss(exp_for_coloc$beta / sqrt(exp_for_coloc$varbeta),
      exp_for_coloc$LD,
      n = N_exp
    ),
    alignment_check = check_alignment(exp_for_coloc)
  )

  out_alignment_qc <- list(
    kriging = kriging_rss(out_for_coloc$beta / sqrt(out_for_coloc$varbeta),
      out_for_coloc$LD,
      n = N_out
    ),
    alignment_check = check.alignment(out_for_coloc)
  )

  # Run SuSiE for fine-mapping and colocalization
  exp_susie <- tryCatch(
    expr = {
      temp <- coloc::runsusie(exp_for_coloc,
        repeat_until_convergence = F,
        maxit = 1000
      )
    },
    error = function(e) {
      warning("SuSiE did not converge for exposure")
      NULL
    }
  )

  out_susie <- tryCatch(
    expr = {
      temp <- coloc::runsusie(out_for_coloc,
        repeat_until_convergence = F,
        maxit = 1000
      )
    },
    error = function(e) {
      warning("SuSiE did not converge for outcome")
      NULL
    }
  )
  if (!is.null(exp_susie) & !is.null(out_susie)) {
    susie.res <- coloc::coloc.susie(exp_susie, out_susie)
  } else {
    warning("Returning NULL SuSiE object")
    susie.res <- NULL
  }

  return(list(
    susie.res = susie.res,
    abf.res = ABF,
    exp_susie = exp_susie,
    out_susie = out_susie
  ))
}
