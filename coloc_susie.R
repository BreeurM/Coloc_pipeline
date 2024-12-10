rm(list = ls())

library(data.table)
library(tidyverse)
library(coloc)
library(susieR)
library(TwoSampleMR)

setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/coloc_susie_utils.R")

########## Digging deeper into the runsusie outputs

setwd("~/Code/MR_EPIC")
can <- "Prostate" ####### TO CHANGE
file_path <- paste0("harmonised_data/", can, "/")
res_file <- paste0("res/", can, ".csv")
sum_stat_path <- "N:/EPIC_genetics/Annotated_windows/Annotated_windows/"
file_list <- list.files(sum_stat_path, full.names = T, recursive = T)

harm_dat <- fread(paste0(file_path, "All_SNPs.csv"))
res_mr <- read_csv(res_file)


prot <- "MSMB"

# Get the IV and info on the gene.exposure
res_prot <- res_mr %>% filter(Assay == prot & analysis == "Cis")
snp_dat <- harm_dat %>% filter(SNP == res_prot$method[which.min(res_prot$p)] &
                                 exposure == prot)
sum_stat_file <- file_list[grepl(snp_dat$SNP, file_list) &
                             grepl(str_replace_all(paste0(snp_dat$gene.exposure, "_"), "-", "-"), file_list)]

############ Define the region to keep
bfile_loc <- "N:/EPIC_genetics/1000G_EUR/1000G_EUR/QC_1000G_P3"
plink_loc <- "plink"

window_size_kb <- 1000

cat("Extracting SNPs in a", window_size_kb, "kb window around", snp_dat$SNP, "\n")
# Execute PLINK command
extracted_snps_file <- tempfile("extracted_snps", fileext = ".txt")
cmd <- paste(
  plink_loc, "--bfile", bfile_loc,
  "--snp", snp_dat$SNP,
  "--window", window_size_kb,
  "--write-snplist",
  "--out", extracted_snps_file,
  "--allow-no-sex"
)

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


###################### Load datasets for coloc

# Load summary statistics for the IV region
prot_region <- TwoSampleMR::read_exposure_data(sum_stat_file,
                                               sep = ",",
                                               snp_col = "RSID", beta_col = "BETA",
                                               se_col = "SE", log_pval = T, pval_col = "LOG10P",
                                               effect_allele_col = "ALLELE1",
                                               other_allele_col = "ALLELE0",
                                               pos_col = "GENPOS", chr_col = "CHROM", eaf_col = "A1FREQ", min_pval = NA
) %>% filter(SNP %in% extracted_snps$SNP)


harm_data <- prot_region

harm_data <- harm_data %>%
  mutate(
    effect_allele = effect_allele.exposure,
    other_allele = other_allele.exposure
  )

harm_data <- harm_data %>%
  select(-c(
    effect_allele.exposure,
    other_allele.exposure
  ))

LD_matrix <- get_ld_matrix(harm_data$SNP, plink_loc = plink_loc, bfile_loc = bfile_loc, with_alleles = T)$LD_Anal

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
    eaf.exposure = if_else(flipped, 1 - eaf.exposure, eaf.exposure)
  )

harm_data <- temp %>%
  mutate(pos = pos.exposure) %>%
  select(
    SNP, pos,
    beta.exposure,
    se.exposure,
    eaf.exposure
  )


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
exp_for_coloc$type <- "quant"
exp_for_coloc$sdY <- 1
exp_for_coloc$N <- 36000

colnames(LD_matrix) <- str_split_fixed(colnames(LD_matrix), "_", 2)[, 1]
rownames(LD_matrix) <- colnames(LD_matrix)
exp_for_coloc$LD <- LD_matrix[exp_for_coloc$snp, exp_for_coloc$snp]

check_alignment(exp_for_coloc)

s1 <- coloc::runsusie(exp_for_coloc,
                            repeat_until_convergence = F,
                            maxit = 1000)


cs1=s1$sets
idx1=cs1$cs_index
bf1=s1$lbf_variable[idx1,,drop=FALSE]



