rm(list = ls())

library(data.table)
library(tidyverse)
library(dplyr)
library(tidyr)
library(coloc)
library(susieR)
library(TwoSampleMR)
library(httr)
library(glue)
library(jsonlite)
library(ggrepel)
library(ggpubr)
library(readr)
library(qqman)
library(tictoc)
library(GenomicRanges)

setwd("~/Code/Coloc_pipeline")
source("~/Code/Coloc_pipeline/scripts/pipeline_utils.R")


plink_path <- "plink"
bfile_path <- "N:/EPIC_genetics/UKBB/LD_REF_FILES/LD_REF_DAT_MAF_MAC_Filtered"

temp_dir_path <- "Temp"

################################################################################
# List the data available from the eQTL catalogue
################################################################################


# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md
# https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/coloc.susie/coloc_susie.md

# max_pulled_rows <- 1000
# URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")
# r <- GET(URL, accept_json())
# cont <- content(r, "text", encoding = "UTF-8")
# datasets <- fromJSON(cont)


################################################################################
# Extract LBF for kidney study
################################################################################

# path_lbf      <- "data/eQTL_catalogue/QTD000264.lbf_variable.txt.gz"
# path_sumstats <- "data/eQTL_catalogue/QTD000264.cc.tsv.gz"

# kid_genexp_lbf <- fread(path_lbf)
# traits <- names(table(kid_genexp_lbf$molecular_trait_id))
# length(traits)


trait_id <- "ENSG00000143537.grp_2.contained.ENST00000360674"
fullres  <- readRDS(paste0("~/Code/Coloc_pipeline/results/fullres_", trait_id, "_RenalC.rds"))
trait    <- fullres$trait


################################################################################
# Run finemapping of renal cancer and get Bayes factors if needed
################################################################################


out <- fullres$out


# path_to_out <- "N:/EPIC_genetics/Cancer_sumstats/RENAL_multi_ancestry/EURO_ONLY.tsv"
# 
# ## Find out which regions are already covered
# region <- unique(trait$region)
# lead_pos <- trait %>%
#   filter(pval == min(pvalue, na.rm = TRUE)) %>%
#   pull(pos) %>% mean
# # Set the directory containing the files
# lbf_directory <- "data/renal_cancer_lbf"
# 
# region_utils <- region_utils(region, lead_pos, lbf_directory)
# 
# 
# if (region_utils$is_pos_covered) {
#   message(paste0("Region already fine-mapped. Laoding the BFs from ", lbf_directory))
#   out_lbf <- readRDS(paste0(
#     lbf_directory, "/lbf_chr",
#     as.character(region_utils$position_covered$chromosome),
#     "_",
#     as.character(region_utils$position_covered$start_position),
#     "-",
#     as.character(region_utils$position_covered$end_position),
#     ".rds"
#   ))
#   
#   
#   message("Loading outcome summary stats.")
# 
#   out_raw <- fread(path_to_out)
# 
#   out <- out_raw %>%
#     filter(SNP %in% out_lbf$SNP)
# 
#   rm(out_raw)
# 
#   out <- format_data(data.frame(out),
#     chr_col = "chromosome",
#     pos_col = "base_pair_location",
#     snp_col = "rsid",
#     beta_col = "beta",
#     se_col = "standard_error",
#     pval_col = "p_value",
#     log_pval = FALSE,
#     eaf_col = "effect_allele_frequency",
#     effect_allele_col = "effect_allele",
#     other_allele_col = "other_allele"
#   )
# 
# } else {
#   # finemapping needs to be ran
# 
#   message("Loading outcome summary stats.")
# 
#   out_raw <- fread(path_to_out)
#   # out_man <- out_raw %>% filter(p_value < 1e-5)
#   # snps <- out_man %>% filter(((chromosome == "1")& between(base_pair_location, 156213007, 156240042))|
#   #                              ((chromosome == "1")& between(base_pair_location, 155050566,155062775))|
#   #                              ((chromosome == "2")& between(base_pair_location, 8559833,8583792))|
#   #                              ((chromosome == "2")& between(base_pair_location, 127638426,127681786))|
#   #                              ((chromosome == "2")& between(base_pair_location, 201116164,201176687))|
#   #                              ((chromosome == "7")& between(base_pair_location, 128937032,128950038))|
#   #                              ((chromosome == "16")& between(base_pair_location, 28822999,28837232))|
#   #                              ((chromosome == "16")& between(base_pair_location, 79721312,79770532))|
#   #                              ((chromosome == "22")& between(base_pair_location, 37807934,37817183))) %>% 
#   #   select(rsid)
#   # plot_can <- manhattan(out_man, chr = "chromosome", bp = "base_pair_location", snp = "rsid", p = "p_value", highlight = c(snps$rsid))
# 
#   # can_lead_var <- out_raw[which.min(out_raw$p_value), variant_id]
#   # # [1] "11_69422827_C_T"
#   #
# 
#   out <- out_raw %>%
#     filter(chromosome == lead_chr) %>%
#     filter(between(base_pair_location, pos_low, pos_high))
# 
#   rm(out_raw)
# 
#   out <- format_data(data.frame(out),
#     chr_col = "chromosome",
#     pos_col = "base_pair_location",
#     snp_col = "rsid",
#     beta_col = "beta",
#     se_col = "standard_error",
#     pval_col = "p_value",
#     log_pval = FALSE,
#     eaf_col = "effect_allele_frequency",
#     effect_allele_col = "effect_allele",
#     other_allele_col = "other_allele"
#   )
# 
#   # ggplot(out %>% filter(p_value < 1e-5), aes(x = pos, y= -log10(pval))) + geom_point()
# 
#   ## Find out if anything passes the conventional significance threshold
# 
#   message(paste0("Region not covered, fine-mapping region ", pos_interval, " on chr ", lead_chr))
# 
#   if (any(out$pval < 5e-2)) {
#     ## Run the finemapping
# 
#     tic("Running SuSiE")
#     out_trait_susie <- tryCatch(expr = {
#       finemap_susie(
#         exp_data = out,
#         N_exp = 780000,
#         exp_type = "cc",
#         exp_sd = 0.034,
#         LD_matrix = NULL,
#         plink_loc = plink_loc,
#         bfile_loc = bfile_loc,
#         max_iter = 500,
#         repeat_until_convergence = FALSE,
#         run_checks = FALSE
#       )
#     }, error = function(e) {
#       warning("SuSiE did not converge, returning BF = 0")
#       NULL
#     })
#     toc()
# 
#     if (!is.null(out_trait_susie)) {
#       out_lbf <- as.data.frame(t(out_trait_susie$lbf_variable)) %>%
#         setNames(paste0("lbf_variable", 1:10)) %>%
#         rownames_to_column("SNP") %>%
#         left_join(out %>%
#           transmute(variant = paste0("chr", chr, "_", pos), SNP), by = "SNP")
#     } else {
#       out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
#         setNames(paste0("lbf_variable", 1:10)) %>%
#         mutate(
#           variant = paste0("chr", out$chr, "_", out$pos),
#           SNP = out$SNP
#         )
#     }
# 
#     saveRDS(out_lbf, file = paste0(lbf_directory, "/lbf_chr", lead_chr, "_", pos_interval, ".rds"))
#   } else {
#     message("No variant crosses the significance threshold. Returning BF = 0")
# 
#     out_lbf <- as.data.frame(matrix(-5, nrow = nrow(out), ncol = 10)) %>%
#       setNames(paste0("lbf_variable", 1:10)) %>%
#       mutate(
#         variant = paste0("chr", out$chr, "_", out$pos),
#         SNP = out$SNP
#       )
#   }
# }
# 
# out <- merge(out, out_lbf)
# 
# rm(position_covered)

################################################################################
# Run coloc based on BF
################################################################################


res_coloc <- coloc.wrapper(trait, out)


################################################################################
# Vanilla coloc from the summary stats
################################################################################


# # Flip alleles for exp and out if not consistent.
# # Hold on, needed? lbfs won't change...
# # Will be needed for the MR and maybe the plots but we can opt out for now
# 
# 
# # Run coloc based on lbfs computed in the format data function
# 
# trait_mat <- as.matrix(dplyr::select(trait, lbf)) # %>%
# # filter(!is.na(lbf)))
# row.names(trait_mat) <- sub("(_[ACGT]+_[ACGT]+)$", "", trait$variant[!is.na(trait$lbf)])
# trait_mat <- t(trait_mat)
# 
# 
# out_mat <- as.matrix(dplyr::select(out, lbf))
# row.names(out_mat) <- out$variant
# out_mat <- t(out_mat)
# 
# res_abf_coloc <- coloc::coloc.bf_bf(trait_mat, out_mat)
# res_abf_coloc$summary
# 
# rm(trait_mat, out_mat)


################################################################################
# Investigate causality with MR
################################################################################


res_mr <- mr.wrapper(trait, out, 780000, res_coloc)


################################################################################
# ZZ and locus plot
################################################################################


## Query LD matrix for desired region
out_name = "outcome"
trait_name = "exposure"

# Define plot window
lead_pos <- trait %>%
  filter(abs(z) == max(abs(z), na.rm = TRUE)) %>%
  pull(pos) %>%
  mean()

required_start <- lead_pos - 500000
required_end <- lead_pos + 500000

# Extract LD matrix
SNP_list <- trait$SNP[between(trait$pos, required_start, required_end)]

# LD_matrix <- get_ld_matrix_from_bim(SNP_list,
#                                     plink_loc = plink_path,
#                                     bfile_loc = bfile_path,
#                                     with_alleles = T,
#                                     temp_dir_path = temp_dir_path
# )
# 
# saveRDS(LD_matrix, paste0("data/gen_utils/ld_mat_", trait_id))

LD_matrix <- readRDS(paste0("data/gen_utils/ld_mat_", trait_id))

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
  vars_trait <- as.data.frame(res_coloc) %>% dplyr::select(!!sym(paste0("hit_", str_replace(trait_name," ","_"))))
  lead_snp <- trait$SNP[trait$variant == vars_trait[which.max(res_coloc$PP.H4.abf),1]]
  # coloc_snp = hit in hit_out_name
  vars <- as.data.frame(res_coloc) %>% dplyr::select(!!sym(paste0("hit_", str_replace(out_name," ","_"))))
  coloc_snp <- out$SNP[out$variant == vars[which.max(res_coloc$PP.H4.abf),1]]
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
  trait_name, out_name, out_type = "cc")




plots <- ggarrange(ggarrange(plot_list$zzplot, plot_list$table, ncol = 2, widths = c(1, .8)),
           ggarrange(plot_list$pvalues_at_locus,
                     ggarrange(plot_list$exposure_at_locus,
                               plot_list$outcome_at_locus,
                               nrow = 2),
                     ncol = 2, widths = c(1, .8)),
           nrow = 2)

plotfile <- paste0("results/plots_", 
                   str_replace(trait_id, " ", "_"), "_", 
                   "RenalC", ".pdf")
ggsave(plotfile, plots, dpi = 300, width = 15, height = 10)











################################################################################
################################################################################





fullres_noplot <- readRDS("~/Code/Coloc_pipeline/results/fullres_ENSG00000286106_RenalC.rds")
trait <- fullres_noplot$trait
out <- fullres_noplot$out

res_coloc <- coloc.wrapper(trait, out)


res_mr = NULL




plots <- plot.wrapper(trait, out, res_coloc, res_mr)
plotfile <- paste0("results/plots_ENSG00000286106_", 
                   "RenalC", ".pdf")
ggsave(plotfile, plots, dpi = 300, width = 15, height = 10)
