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
format_dat <- function(dat, suffixe = "", snps = NA,
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
  
  all_cols <- allcols[!is.na(all_cols)]
  
  dat <- dat %>%
    dplyr::select(any_of(all_cols)) %>%
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
    dat <- dat %>% dplyr::select(-all_of(phenotype_col))
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

#' Infers the genome build of the summary statistics file (GRCh37 or GRCh38)
#' from the data. Uses SNP (RSID) & CHR & BP to get genome build.
#'
#' @param sumstats data table/data frame obj of the summary statistics file for
#' the GWAS ,or file path to summary statistics file.
#' @param nThread Number of threads to use for parallel processes.
#' @param sampled_snps Downsample the number of SNPs used when inferring genome
#' build to save time.
#' @param standardise_headers Run
#' @param standardise_headers Run
#' \code{standardise_sumstats_column_headers_crossplatform}.
#' @param mapping_file \pkg{MungeSumstats} has a pre-defined
#' column-name mapping file
#' which should cover the most common column headers and their interpretations.
#' However, if a column header that is in your file is missing of the mapping we
#' give is incorrect you can supply your own mapping file. Must be a 2 column
#' dataframe with column names "Uncorrected" and "Corrected". See
#' \code{data(sumstatsColHeaders)} for default mapping and necessary format.
#' @param dbSNP version of dbSNP to be used (144 or 155). Default is 155.
#' @param header_only Instead of reading in the entire \code{sumstats} file,
#' only read in the first N rows where N=\code{sampled_snps}.
#' This should help speed up cases where you have to read in \code{sumstats}
#' from disk each time.
#' @param allele_match_ref Instead of returning the genome_build this will
#' return the propotion of matches to each genome build for each allele (A1,A2).
#' @inheritParams format_sumstats
#' @inheritParams get_genome_builds
#'
#' @return ref_genome the genome build of the data
#' @importFrom data.table setDT :=
#' @keywords internal
get_genome_build <- function(dataset, dbSNP = 155) {
  ### Add this to avoid confusing BiocCheck
  seqnames <- CHR <- SNP <- BP <- alt_alleles <- NULL
  sumstats <- data.table::setDT(sumstats)
  # need SNP ID column (RS ID) CHR and BP (POS) to infer build
  # - check these are present, considering all known names
  if (standardise_headers) {
    sumstats_return <-
      standardise_sumstats_column_headers_crossplatform(
        sumstats_dt = sumstats,
        mapping_file = mapping_file
      )
    sumstats <- sumstats_return$sumstats_dt
  }
  
  err_msg <-
    paste0(
      "SNP ID column (RS ID), CHR and BP (POSITION) columns are needed ",
      "to infer the genome build. These could not be\nfound in your ",
      "dataset. Please specify the genome build manually to run ",
      "format_sumstats()."
    )
  # Infer genome build using SNP & CHR & BP
  if (!all(c("SNP", "CHR", "BP") %in% colnames(sumstats))) {
    # want it returned rather than throwing an error
    if (isTRUE(allele_match_ref)) {
      return(err_msg)
    } else {
      stop(err_msg)
    }
  }
  
  #### Do some filtering first to avoid errors ####
  nrow_org <- nrow(sumstats)
  sumstats <- sumstats[complete.cases(SNP, BP, CHR)]
  err_msg2 <-
    paste0(
      "SNP ID column (RS ID), CHR and BP (POSITION)",
      "columns are needed to",
      " infer the genome build.",
      "These contain too many\nmissing values in",
      " your dataset to be used.",
      "Please specify the genome build manually",
      " to run format_sumstats()"
    )
  # also remove common incorrect formatting of SNP
  sumstats <- sumstats[grepl("^rs", SNP), ]
  sumstats <- sumstats[SNP != ".", ]
  # also deal with common misformatting of CHR
  # if chromosome col has chr prefix remove it
  sumstats[, CHR := gsub("chr", "", CHR)]
  
  # for internal testing - filter to specified chromosomes
  if (!is.null(chr_filt)) {
    sumstats <- sumstats[CHR %in% chr_filt]
  }
  
  # if removing erroneous cases leads to <min(10k,50% org dataset) will fail -
  # NOT ENOUGH DATA TO INFER
  nrow_clean <- nrow(sumstats)
  size_okay <- FALSE
  if (nrow_clean > sampled_snps || (nrow_clean != 0 &&
                                    (nrow_clean / nrow_org) > .5)) {
    size_okay <- TRUE
  }
  if (!size_okay) {
    # want it returned rather than throwing an error
    if (isTRUE(allele_match_ref)) {
      return(err_msg2)
    } else {
      stop(err_msg2)
    }
  }
  #### Downsample SNPs to save time ####
  if ((nrow(sumstats) > sampled_snps) && !(is.null(sampled_snps))) {
    snps <- sample(sumstats$SNP, sampled_snps)
  } else { # nrow(sumstats)<10k
    snps <- sumstats$SNP
  }
  
  sumstats <- sumstats[SNP %in% snps, ]
  
  # now split into functions two roles
  if (isTRUE(allele_match_ref)) {
    # 1. checking for matches for A1/A2 to ref genomes
    if (is.null(ref_genome)) {
      # have to check multiple
      snp_loc_data_37 <- load_ref_genome_data(
        snps = snps,
        ref_genome = "GRCH37",
        dbSNP = dbSNP
      )
      snp_loc_data_38 <- load_ref_genome_data(
        snps = snps,
        ref_genome = "GRCH38",
        dbSNP = dbSNP
      )
      # convert CHR filed in ref genomes to character not factor
      snp_loc_data_37[, seqnames := as.character(seqnames)]
      snp_loc_data_38[, seqnames := as.character(seqnames)]
      # convert CHR filed in data to character if not already
      sumstats[, CHR := as.character(CHR)]
      # Now check which genome build has more matches to data
      num_37 <-
        nrow(snp_loc_data_37[sumstats, ,
                             on = c("SNP" = "SNP", "pos" = "BP", "seqnames" = "CHR"),
                             nomatch = FALSE
        ])
      num_38 <-
        nrow(snp_loc_data_38[sumstats, ,
                             on = c("SNP" = "SNP", "pos" = "BP", "seqnames" = "CHR"),
                             nomatch = FALSE
        ])
      
      if (num_37 > num_38) {
        ref_gen_num <- num_37
        ref_genome <- "GRCH37"
        snp_loc_data <- snp_loc_data_37
      } else {
        ref_gen_num <- num_38
        ref_genome <- "GRCH38"
        snp_loc_data <- snp_loc_data_38
      }
    } else {
      # only check one chosen
      snp_loc_data <- load_ref_genome_data(
        snps = snps,
        ref_genome = ref_genome,
        dbSNP = dbSNP
      )
      # convert CHR filed in ref genomes to character not factor
      snp_loc_data[, seqnames := as.character(seqnames)]
      # convert CHR filed in data to character if not already
      sumstats[, CHR := as.character(CHR)]
    }
    # Now check which allele has more matches to data
    # want to match on ref and alt alleles too
    # need to take first alt from list to do this
    snp_loc_data[, alt_alleles := unlist(lapply(
      alt_alleles,
      function(x) x[[1]]
    ))]
    num_a1 <-
      nrow(snp_loc_data[sumstats, ,
                        on = c(
                          "SNP" = "SNP", "pos" = "BP", "seqnames" = "CHR",
                          "ref_allele" = "A1", "alt_alleles" = "A2"
                        ),
                        nomatch = FALSE
      ])
    num_a2 <-
      nrow(snp_loc_data[sumstats, ,
                        on = c(
                          "SNP" = "SNP", "pos" = "BP", "seqnames" = "CHR",
                          "ref_allele" = "A2", "alt_alleles" = "A1"
                        ),
                        nomatch = FALSE
      ])
    
    if (num_a1 >= num_a2) {
      message("Effect/frq column(s) relate to A2 in the inputted sumstats")
      # this is what MSS expects so no action required
      switch_req <- FALSE
    } else { # num_a1<num_a2
      message("Effect/frq column(s) relate to A1 in the inputted sumstats")
      switch_req <- TRUE
    }
    
    return(switch_req)
  } else {
    # 2. checking for ref genome
    # otherwise SNP, CHR, BP were all found and can infer
    snp_loc_data_37 <- load_ref_genome_data(
      snps = snps,
      ref_genome = "GRCH37",
      dbSNP = dbSNP,
      chr_filt = chr_filt
    )
    snp_loc_data_38 <- load_ref_genome_data(
      snps = snps,
      ref_genome = "GRCH38",
      dbSNP = dbSNP,
      chr_filt = chr_filt
    )
    # convert CHR filed in ref genomes to character not factor
    snp_loc_data_37[, seqnames := as.character(seqnames)]
    snp_loc_data_38[, seqnames := as.character(seqnames)]
    # convert CHR filed in data to character if not already
    sumstats[, CHR := as.character(CHR)]
    # Now check which genome build has more matches to data
    num_37 <-
      nrow(snp_loc_data_37[sumstats, ,
                           on = c("SNP" = "SNP", "pos" = "BP", "seqnames" = "CHR"),
                           nomatch = FALSE
      ])
    num_38 <-
      nrow(snp_loc_data_38[sumstats, ,
                           on = c("SNP" = "SNP", "pos" = "BP", "seqnames" = "CHR"),
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
}


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
#' 1. **Identify Shell Type**: Determines whether the system shell is "cmd" (Windows) or "sh" (Unix-based).
#' 2. **Create Temporary File**: Writes the `rsid_list` to a temporary file for input into PLINK.
#' 3. **Generate BIM File**:
#'    - Constructs a command to run PLINK with `--make-just-bim`, filtering the binary files for SNPs in `rsid_list`.
#'    - Executes the command and reads the resulting BIM file, which contains SNP metadata.
#' 4. **Compute LD Matrix**:
#'    - Builds another PLINK command to calculate the LD matrix using `--r square`.
#'    - Runs the command and reads the resulting LD matrix into R as a numeric matrix.
#' 5. **Annotate Matrix**:
#'    - If `with_alleles = TRUE`, sets matrix row and column names to include rsID and allele information (e.g., rsID_A1_A2).
#'    - Otherwise, uses only the rsID for row and column names.
#' 6. **Clean Up**: Deletes all temporary files created during the process.
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
get_ld_matrix_from_bim <- function(rsid_list, plink_loc, bfile_loc, plink_memory = NULL, with_alleles = TRUE) {
  # Determine shell type based on operating system
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")

  # Create a temporary file to store rsID list
  fn <- tempfile()
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
  tryCatch(expr = {
    file.remove(list.files(tempdir(), pattern = ".ld", recursive = TRUE, full.names = TRUE))
    file.remove(list.files(tempdir(), pattern = ".bim", recursive = TRUE, full.names = TRUE))
    file.remove(list.files(tempdir(), pattern = ".log", recursive = TRUE, full.names = TRUE))
    file.remove(list.files(tempdir(), pattern = ".nosex", recursive = TRUE, full.names = TRUE))
  }, error = function(e) {
    warning("Temporary files used for LD computation could not be removed.")
  })

  return(res)
}


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
                          LD_matrix = NULL, plink_loc = "plink", bfile_loc = NULL,
                          ...) {
  # Check input consistency
  # exp_type

  # Import exposure data if provided as a file path
  if (is.character(exp_data)) {
    exp_data <- read.csv(exp_data)
  }

  harm_data <- exp_data
  # Import LD_matrix if provided as file paths, or import from 1000G if NULL
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
      LD_matrix <- get_ld_matrix_from_bim(harm_data$SNP, plink_loc = plink_loc, bfile_loc = bfile_loc, with_alleles = T)
    }
  } else {
    if (is.character(LD_matrix)) {
      LD_matrix <- read.csv(LD_matrix)
    }
  }

  # Check if allele information is present in the LD matrix and align if necessary
  if (!(all(grepl("^[ACTG]+$", str_split_fixed(colnames(LD_matrix), "_", 3)[, 2])) &
    all(grepl("^[ACTG]+$", str_split_fixed(colnames(LD_matrix), "_", 3)[, 3])))) {
    # Handle case where allele information is missing or incorrect
    warning("LD matrix has no or incorrect allele information, allele alignment will be inferred from expected Z-scores vs. observed Z-scores")
    harm_data <- harm_data %>%
      mutate(pos = pos.exposure) %>%
      dplyr::select(
        SNP, pos,
        beta.exposure,
        se.exposure,
        eaf.exposure
      )
  } else {
    # Flip alleles in harmonized data to align with the LD matrix
    LD_alignment <- data.frame(
      SNP = str_split_fixed(colnames(LD_matrix), "_", 3)[, 1],
      LD_A1 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 2],
      LD_A2 = str_split_fixed(colnames(LD_matrix), "_", 3)[, 3]
    )
    temp <- merge(harm_data, LD_alignment)
    temp <- temp %>%
      mutate(
        flipped = effect_allele.exposure != LD_A1,
        effect_allele = if_else(flipped, other_allele.exposure, effect_allele.exposure),
        other_allele = if_else(flipped, effect_allele.exposure, other_allele.exposure),
        beta.exposure = if_else(flipped, -beta.exposure, beta.exposure),
        eaf.exposure = if_else(flipped, 1 - eaf.exposure, eaf.exposure)
      )

    if (any(temp$LD_A1 != temp$effect_allele)) {
      # Stop execution if alleles cannot be aligned
      n_excl <- sum(temp$LD_A1 != temp$effect_allele)
      temp <- temp %>% filter(LD_A1 == effect_allele)
      warning(paste0("Inconsistent alleles when matching to LD matrix. Excluding ", as.character(n_excl), " variants from the analysis."))
    }

    harm_data <- temp %>%
      mutate(pos = pos.exposure) %>%
      dplyr::select(
        SNP, pos,
        beta.exposure,
        se.exposure,
        eaf.exposure
      )
  }

  # Remove alleles from LD matrix column names for cross ref with ext/out data
  colnames(LD_matrix) <- str_split_fixed(colnames(LD_matrix), "_", 2)[, 1]
  rownames(LD_matrix) <- colnames(LD_matrix)

  # Prepare data for SuSiE colocalization analysis
  exp_for_coloc <- harm_data %>%
    dplyr::select(
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

  # Check alignment between exposure data and LD matrix
  # Flags for faulty allele alignment
  exp_qc <- list(
    kriging = kriging_rss(exp_for_coloc$beta / sqrt(exp_for_coloc$varbeta),
      exp_for_coloc$LD,
      n = N_exp
    ),
    alignment_check = check_alignment(exp_for_coloc)
  )
  if (exp_qc$alignment_check < 0.7) {
    warning("Suspected alignment error.")
  }

  # Calculate Z-scores for SuSiE
  z <- harm_data$beta.exposure / harm_data$se.exposure

  snp <- harm_data$SNP
  names(z) <- snp
  LD <- LD_matrix[snp, snp]

  # Initialize convergence flag
  converged <- FALSE

  # Set defaults for SuSiE arguments
  susie_args <- list(...)
  if ("max_iter" %in% names(susie_args)) {
    maxit <- susie_args$max_iter
    susie_args <- susie_args[setdiff(names(susie_args), "max_iter")]
  }

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
    message("Running max iterations: ", maxit)
    res <- do.call(
      susie_rss,
      c(list(z = z, R = LD, max_iter = maxit), susie_args)
    )
    converged <- res$converged
    message("\tConverged: ", converged)
    if (!converged && repeat_until_convergence == FALSE) {
      stop("susie_rss() did not converge in ", maxit, " iterations. Try running with run_until_convergence=TRUE")
    }
    if (!converged) {
      maxit <- maxit * 100
    } # Increase iterations if not converged
  }

  # Annotate SuSiE results
  susie.res <- annotate_susie(res, snp, LD)

  # Format the data

  LBF <- as.matrix(susie.res$lbf_variable)

  return(list(susie.res = susie.res))
}


#' Extract Credible Sets and Posterior Statistics from a SuSiE Object
#'
#' This function processes a SuSiE (Sum of Single Effects) object to extract credible sets,
#' their purity, and posterior statistics such as posterior means, standard deviations,
#' and log Bayes factors. It also identifies overlapping credible sets and computes
#' metrics related to their quality.
#'
#' @param susie_object A SuSiE object, typically the result of running `susieR::susie()`,
#'                     containing information about credible sets, posterior inclusion
#'                     probabilities (PIPs), and other statistical metrics.
#'
#' @return A list with three data frames:
#'   \item{cs_df}{A data frame containing credible set information, including size,
#'                log10 Bayes factors, and purity metrics.}
#'   \item{variant_df}{A data frame with variant-level statistics, including posterior
#'                     means, standard deviations, and inclusion in credible sets.}
#'   \item{lbf_df}{A data frame with log Bayes factors for each variable across credible sets.}
#'
#' @details
#' The function identifies credible sets, evaluates their purity based on correlation metrics,
#' and flags low-purity sets. Posterior statistics, including means, standard deviations,
#' and Bayes factors, are extracted for each variant. Overlapping credible sets are also
#' tracked and excluded from certain outputs.
#'
#' @importFrom dplyr as_tibble mutate filter select bind_cols left_join
#' @importFrom purrr map_df
#' @importFrom susieR susie_get_posterior_mean susie_get_posterior_sd
#'
#' @author Kaur Alasoo: https://github.com/eQTL-Catalogue/qtlmap/blob/master/bin/run_susie.R
#' @export
extractResults <- function(susie_object) {
  # Extract credible sets from the susie_object
  credible_sets <- susie_object$sets$cs
  cs_list <- list()

  # Initialize purity data and add additional columns for analysis
  susie_object$sets$purity <- dplyr::as_tibble(susie_object$sets$purity) %>%
    dplyr::mutate(
      cs_id = rownames(susie_object$sets$purity), # Add credible set ID
      cs_size = NA, # Initialize size of credible set
      cs_log10bf = NA, # Initialize log10 Bayes factor
      overlapped = NA # Initialize overlap flag
    )

  # Track variants already included in credible sets
  added_variants <- c()

  # Iterate over each credible set
  for (index in seq_along(credible_sets)) {
    cs_variants <- credible_sets[[index]] # Variants in the current credible set
    cs_id <- susie_object$sets$cs_index[[index]] # Index of the credible set

    # Check if any variants overlap with previously added sets
    is_overlapped <- any(cs_variants %in% added_variants)
    susie_object$sets$purity$overlapped[index] <- is_overlapped
    susie_object$sets$purity$cs_size[index] <- length(cs_variants) # Size of credible set
    susie_object$sets$purity$cs_log10bf[index] <- log10(exp(susie_object$lbf[cs_id])) # Log10 Bayes factor

    # If the set is not overlapping, add it to the results
    if (!is_overlapped) {
      cs_list[[index]] <- dplyr::tibble(
        cs_id = paste0("L", cs_id), # Format the credible set ID
        variant_id = susie_object$variant_id[cs_variants] # Variants in the credible set
      )
      added_variants <- append(added_variants, cs_variants) # Update added variants
    }
  }

  # Combine all credible set data into a single data frame
  df <- purrr::map_df(cs_list, identity)

  # Extract purity values for all sets
  purity_res <- susie_object$sets$purity

  # Check for empty purity results; if not empty, process further
  if (nrow(purity_res) > 0) {
    purity_df <- dplyr::as_tibble(purity_res) %>%
      dplyr::filter(!overlapped) %>% # Filter out overlapping sets
      dplyr::mutate(
        cs_avg_r2 = mean.abs.corr^2, # Average R^2
        cs_min_r2 = min.abs.corr^2, # Minimum R^2
        low_purity = min.abs.corr < 0.5 # Flag for low purity
      ) %>%
      dplyr::select(cs_id, cs_log10bf, cs_avg_r2, cs_min_r2, cs_size, low_purity) # Select relevant columns
  } else {
    purity_df <- dplyr::tibble() # Empty purity data frame
  }

  # Extract posterior means, standard deviations, and other related values
  mean_vec <- susieR::susie_get_posterior_mean(susie_object)
  sd_vec <- susieR::susie_get_posterior_sd(susie_object)

  # Extract alpha, mean, and standard deviation matrices
  alpha_mat <- t(susie_object$alpha)
  colnames(alpha_mat) <- paste0("alpha", seq(ncol(alpha_mat)))

  mean_mat <- t(susie_object$alpha * susie_object$mu) / susie_object$X_column_scale_factors
  colnames(mean_mat) <- paste0("mean", seq(ncol(mean_mat)))

  sd_mat <- sqrt(t(susie_object$alpha * susie_object$mu2 - (susie_object$alpha * susie_object$mu)^2)) /
    susie_object$X_column_scale_factors
  colnames(sd_mat) <- paste0("sd", seq(ncol(sd_mat)))

  # Extract log Bayes factors for variables
  lbf_variable_mat <- t(susie_object$lbf_variable)
  colnames(lbf_variable_mat) <- paste0("lbf_variable", seq(ncol(lbf_variable_mat)))

  # Create a data frame for posterior statistics
  posterior_df <- dplyr::tibble(
    variant_id = rownames(alpha_mat),
    pip = susie_object$pip, # Posterior inclusion probabilities
    z = susie_object$z, # Z-scores
    posterior_mean = mean_vec, # Posterior mean
    posterior_sd = sd_vec # Posterior standard deviation
  ) %>% dplyr::bind_cols(purrr::map(list(alpha_mat, mean_mat, sd_mat), dplyr::as_tibble))

  # Create a data frame for log Bayes factors
  lbf_df <- dplyr::tibble(variant_id = rownames(lbf_variable_mat)) %>%
    dplyr::bind_cols(dplyr::as_tibble(lbf_variable_mat))

  # Combine results into final data frames if data is non-empty
  if (nrow(df) > 0 & nrow(purity_df) > 0 & ncol(lbf_df) > 10) { # Check data validity
    cs_df <- purity_df # Credible set data frame
    variant_df <- dplyr::left_join(posterior_df, df, by = "variant_id") %>%
      dplyr::left_join(cs_df, by = "cs_id") # Combine variant data
  } else {
    cs_df <- NULL
    variant_df <- NULL
    lbf_df <- NULL
  }

  # Return a list of results
  return(list(cs_df = cs_df, variant_df = variant_df, lbf_df = lbf_df))
}
