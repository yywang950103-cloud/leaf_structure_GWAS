#!/usr/bin/env Rscript

# Calculate pairwise linkage disequilibrium (LD) in a target interval and
# assign haplotypes using SNPs in strong LD with a specified lead SNP.
#
# Usage:
#   Rscript 95ci.R [vcf_file] [output_directory]
#
# Defaults:
#   vcf_file        = renamed.vcf
#   output_directory = current working directory

suppressPackageStartupMessages({
  library(vcfR)
})

# -----------------------------------------------------------------------------
# User-configurable parameters
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

vcf_file <- if (length(args) >= 1) args[[1]] else "renamed.vcf"
output_dir <- if (length(args) >= 2) args[[2]] else "."

target_chr <- "5A"
start_bp <- 587342000L
end_bp <- 591065000L

lead_snp_id <- "AVRIG27662"
vrn_snp_id <- "Vrn-A1a"

min_maf <- 0.05
max_missing <- 0.10
ld_cutoffs <- c(0.2, 0.4, 0.6, 0.8)
haplotype_r2_cutoff <- 0.8
major_haplotype_min_n <- 5L

region_label <- sprintf(
  "%s_%.3f_%.3fMb_LD",
  target_chr,
  start_bp / 1e6,
  end_bp / 1e6
)
out_prefix <- file.path(output_dir, region_label)

if (!file.exists(vcf_file)) {
  stop("VCF file not found: ", vcf_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

format_cutoff <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

make_variant_name <- function(chrom, pos, id) {
  bad_id <- is.na(id) | id == "." | id == ""
  result <- id
  result[bad_id] <- paste0(chrom[bad_id], ":", pos[bad_id])
  result
}

gt_to_dosage <- function(g) {
  g <- gsub("\\|", "/", g)
  dosage <- rep(NA_real_, length(g))
  dosage[g == "0/0"] <- 0
  dosage[g %in% c("0/1", "1/0")] <- 1
  dosage[g == "1/1"] <- 2
  dosage
}

calc_maf <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  alt_frequency <- sum(x) / (2 * length(x))
  min(alt_frequency, 1 - alt_frequency)
}

gt_to_allele <- function(gt, ref, alt) {
  if (is.na(gt) || gt == "") return("N")

  gt <- trimws(gsub("\\|", "/", gt))

  if (gt == "0/0") return(ref)
  if (gt == "1/1") return(alt)

  # The panel consists of inbred wheat accessions; heterozygous and missing
  # calls are therefore coded as N for haplotype assignment.
  "N"
}

resolve_variant_name <- function(query, variant_table, required = TRUE) {
  hits <- which(
    variant_table$ID == query |
      variant_table$SNP_name == query
  )

  if (length(hits) == 0) {
    if (required) {
      stop("Variant was not retained after filtering: ", query)
    }
    return(NA_character_)
  }

  if (length(hits) > 1) {
    stop("Variant identifier is not unique after filtering: ", query)
  }

  variant_table$SNP_name[hits]
}

# -----------------------------------------------------------------------------
# Read the VCF and extract the target interval
# -----------------------------------------------------------------------------

cat("Reading VCF: ", vcf_file, "\n", sep = "")

vcf <- read.vcfR(vcf_file, verbose = FALSE)
fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
fix_df$POS <- suppressWarnings(as.integer(fix_df$POS))
fix_df$SNP_name <- make_variant_name(fix_df$CHROM, fix_df$POS, fix_df$ID)

gt <- extract.gt(vcf, element = "GT")

idx_region <- which(
  fix_df$CHROM == target_chr &
    fix_df$POS >= start_bp &
    fix_df$POS <= end_bp
)

if (length(idx_region) == 0) {
  stop("No variants were found in the specified interval.")
}

region <- fix_df[idx_region, , drop = FALSE]
gt_region <- gt[idx_region, , drop = FALSE]

# -----------------------------------------------------------------------------
# Retain biallelic SNPs and convert genotypes to alternate-allele dosage
# -----------------------------------------------------------------------------

is_biallelic_snp <- (
  nchar(region$REF) == 1 &
    nchar(region$ALT) == 1 &
    !grepl(",", region$ALT)
)

region <- region[is_biallelic_snp, , drop = FALSE]
gt_region <- gt_region[is_biallelic_snp, , drop = FALSE]

if (nrow(region) == 0) {
  stop("No biallelic SNPs were found in the specified interval.")
}

cat("Biallelic SNPs in interval: ", nrow(region), "\n", sep = "")

dosage <- t(apply(gt_region, 1, gt_to_dosage))
colnames(dosage) <- colnames(gt_region)

maf <- apply(dosage, 1, calc_maf)
missing_rate <- apply(dosage, 1, function(x) mean(is.na(x)))

keep <- (
  !is.na(maf) &
    maf >= min_maf &
    missing_rate <= max_missing
)

region_f <- region[keep, , drop = FALSE]
gt_f <- gt_region[keep, , drop = FALSE]
dosage_f <- dosage[keep, , drop = FALSE]
maf_f <- maf[keep]
missing_f <- missing_rate[keep]

if (nrow(region_f) < 2) {
  stop("Fewer than two SNPs remained after MAF and missingness filtering.")
}

if (anyDuplicated(region_f$SNP_name)) {
  duplicate_names <- unique(region_f$SNP_name[duplicated(region_f$SNP_name)])
  stop(
    "Duplicate SNP names were found after filtering: ",
    paste(duplicate_names, collapse = ", ")
  )
}

rownames(dosage_f) <- region_f$SNP_name
rownames(gt_f) <- region_f$SNP_name

cat(
  "SNPs retained after MAF and missingness filtering: ",
  nrow(region_f),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Calculate pairwise Pearson r and r-squared
# -----------------------------------------------------------------------------

# cor() expects samples in rows and SNPs in columns.
geno_for_ld <- t(dosage_f)

r_matrix <- cor(
  geno_for_ld,
  use = "pairwise.complete.obs",
  method = "pearson"
)
r2_matrix <- r_matrix^2

write.csv(
  r2_matrix,
  paste0(out_prefix, "_r2_matrix.csv"),
  quote = FALSE
)

upper_index <- which(upper.tri(r2_matrix), arr.ind = TRUE)
snp_names <- colnames(r2_matrix)

pairwise_ld <- data.frame(
  SNP1 = snp_names[upper_index[, 1]],
  SNP2 = snp_names[upper_index[, 2]],
  r = r_matrix[upper_index],
  r2 = r2_matrix[upper_index],
  stringsAsFactors = FALSE
)

pos_map <- setNames(region_f$POS, region_f$SNP_name)
pairwise_ld$POS1 <- unname(pos_map[pairwise_ld$SNP1])
pairwise_ld$POS2 <- unname(pos_map[pairwise_ld$SNP2])
pairwise_ld$Distance_bp <- abs(pairwise_ld$POS1 - pairwise_ld$POS2)
pairwise_ld <- pairwise_ld[order(pairwise_ld$POS1, pairwise_ld$POS2), ]

write.csv(
  pairwise_ld,
  paste0(out_prefix, "_pairwise_r2.csv"),
  row.names = FALSE,
  quote = FALSE
)

snp_info_output <- data.frame(
  SNP = region_f$SNP_name,
  ID = region_f$ID,
  CHROM = region_f$CHROM,
  POS = region_f$POS,
  REF = region_f$REF,
  ALT = region_f$ALT,
  MAF = unname(maf_f),
  Missing_rate = unname(missing_f),
  stringsAsFactors = FALSE
)

write.csv(
  snp_info_output,
  paste0(out_prefix, "_SNP_info.csv"),
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------------------------------------------------------
# Extract LD between the lead SNP and all retained SNPs
# -----------------------------------------------------------------------------

lead_name <- resolve_variant_name(lead_snp_id, region_f, required = TRUE)
lead_pos <- unname(pos_map[lead_name])

lead_r2 <- data.frame(
  Lead_SNP = lead_name,
  SNP = colnames(r2_matrix),
  r = unname(r_matrix[lead_name, ]),
  r2 = unname(r2_matrix[lead_name, ]),
  stringsAsFactors = FALSE
)

lead_r2$POS <- unname(pos_map[lead_r2$SNP])
lead_r2$Distance_from_lead_bp <- abs(lead_r2$POS - lead_pos)
lead_r2 <- lead_r2[order(lead_r2$Distance_from_lead_bp), ]

write.csv(
  lead_r2,
  paste0(out_prefix, "_", lead_snp_id, "_vs_all_r2.csv"),
  row.names = FALSE,
  quote = FALSE
)

for (cutoff in ld_cutoffs) {
  linked_snps <- subset(
    lead_r2,
    SNP != lead_name &
      !is.na(r2) &
      r2 >= cutoff
  )

  cutoff_label <- format_cutoff(cutoff)

  write.csv(
    linked_snps,
    paste0(
      out_prefix,
      "_", lead_snp_id,
      "_r2_ge_", cutoff_label,
      ".csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )

  cat(
    "Lead-SNP r-squared >= ", cutoff_label, ": ",
    nrow(linked_snps), " SNPs\n",
    sep = ""
  )
}

# -----------------------------------------------------------------------------
# Report LD between the lead SNP and Vrn-A1a when both pass filtering
# -----------------------------------------------------------------------------

vrn_name <- resolve_variant_name(vrn_snp_id, region_f, required = FALSE)

if (!is.na(vrn_name)) {
  cat("\nSpecific LD comparison\n")
  cat("Lead SNP: ", lead_name, "\n", sep = "")
  cat("VRN SNP:  ", vrn_name, "\n", sep = "")
  cat(
    "Distance: ",
    abs(pos_map[lead_name] - pos_map[vrn_name]) / 1e6,
    " Mb\n",
    sep = ""
  )
  cat("r:  ", round(r_matrix[lead_name, vrn_name], 4), "\n", sep = "")
  cat("r2: ", round(r2_matrix[lead_name, vrn_name], 4), "\n", sep = "")
} else {
  cat(
    "\n",
    vrn_snp_id,
    " was not retained after filtering; lead-VRN LD was not calculated.\n",
    sep = ""
  )
}

# -----------------------------------------------------------------------------
# Assign haplotypes using SNPs above the selected lead-SNP LD threshold
# -----------------------------------------------------------------------------

selected_snps <- lead_r2$SNP[
  !is.na(lead_r2$r2) &
    lead_r2$r2 >= haplotype_r2_cutoff
]
selected_snps <- unique(c(lead_name, selected_snps))

cat("\nSNPs used for haplotype assignment:\n")
print(selected_snps)
cat("Number of haplotype SNPs: ", length(selected_snps), "\n", sep = "")

selected_index <- match(selected_snps, region_f$SNP_name)

if (anyNA(selected_index)) {
  stop(
    "Internal matching error for selected SNPs: ",
    paste(selected_snps[is.na(selected_index)], collapse = ", ")
  )
}

hap_snp_info <- region_f[selected_index, , drop = FALSE]
hap_gt <- gt_f[selected_index, , drop = FALSE]

position_order <- order(hap_snp_info$POS)
hap_snp_info <- hap_snp_info[position_order, , drop = FALSE]
hap_gt <- hap_gt[position_order, , drop = FALSE]

allele_df <- data.frame(
  SampleID = colnames(hap_gt),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(hap_snp_info))) {
  snp_name <- hap_snp_info$SNP_name[i]
  allele_df[[snp_name]] <- vapply(
    hap_gt[i, ],
    gt_to_allele,
    character(1),
    ref = hap_snp_info$REF[i],
    alt = hap_snp_info$ALT[i]
  )
}

hap_snp_columns <- hap_snp_info$SNP_name
allele_df$Haplotype_code <- apply(
  allele_df[, hap_snp_columns, drop = FALSE],
  1,
  paste0,
  collapse = ""
)

haplotype_counts <- sort(
  table(allele_df$Haplotype_code),
  decreasing = TRUE
)

haplotype_summary <- data.frame(
  Haplotype_code = names(haplotype_counts),
  n = as.integer(haplotype_counts),
  stringsAsFactors = FALSE
)
haplotype_summary$Haplotype <- paste0(
  "Hap",
  seq_len(nrow(haplotype_summary))
)
haplotype_summary$Major_haplotype <- (
  haplotype_summary$n >= major_haplotype_min_n
)

haplotype_match <- match(
  allele_df$Haplotype_code,
  haplotype_summary$Haplotype_code
)
allele_df$Haplotype <- haplotype_summary$Haplotype[haplotype_match]
allele_df$n <- haplotype_summary$n[haplotype_match]
allele_df$Major_haplotype <- haplotype_summary$Major_haplotype[haplotype_match]

allele_df <- allele_df[, c(
  "SampleID",
  hap_snp_columns,
  "Haplotype_code",
  "Haplotype",
  "n",
  "Major_haplotype"
)]

write.csv(
  allele_df,
  file.path(output_dir, "5A_AVRIG27662_genotype_haplotype_assignment.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  haplotype_summary,
  file.path(output_dir, "5A_AVRIG27662_haplotype_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  hap_snp_info[, c("SNP_name", "ID", "CHROM", "POS", "REF", "ALT")],
  file.path(output_dir, "5A_AVRIG27662_haplotype_SNPs.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nHaplotype counts:\n")
print(haplotype_summary)
cat("\nAnalysis complete. Output directory: ", output_dir, "\n", sep = "")
