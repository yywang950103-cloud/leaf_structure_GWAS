# Population genetics and GWAS analyses for the wheat panel
#
# Usage:
#   Rscript GWAS_2.R /path/to/project /path/to/output
#
# If command-line arguments are omitted, the current directory is used as the
# project directory and results are written to ./results.

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
output_dir <- if (length(args) >= 2) {
  normalizePath(args[[2]], mustWork = FALSE)
} else {
  file.path(project_dir, "results")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Select the analyses to run.
run_pca <- TRUE
run_population_structure <- TRUE
run_nj_tree <- TRUE
run_histogram <- TRUE
run_gwas <- TRUE
run_ld_decay <- TRUE

# Input files. Change only these relative paths when the project layout differs.
pca_eigenvec_file <- file.path(project_dir, "proj.pca.eigenvec")
pca_metadata_file <- file.path(project_dir, "hap", "sample information")
structure_metadata_file <- file.path(project_dir, "lines_meta.tsv")
plink_prefix <- file.path(project_dir, "step3_pruned")
gds_file <- file.path(output_dir, "step3_pruned.gds")
ld_file <- file.path(project_dir, "hap", "LD_genomewide.ld")
histogram_file <- file.path(project_dir, "trait_mean_treatment.csv")
gwas_vcf_file <- file.path(project_dir, "hap", "vcf file")
gwas_phenotype_file <- file.path(project_dir, "trait_control/heat/relative.txt")
gwas_prefix <- "trait"

assert_files <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop("Missing input file(s):\n", paste(missing_paths, collapse = "\n"), call. = FALSE)
  }
}

# -----------------------------------------------------------------------------
# 1. Principal component analysis plot
# -----------------------------------------------------------------------------
if (run_pca) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggsci)
  })

  assert_files(c(pca_eigenvec_file, pca_metadata_file))
  pca <- read.table(
    pca_eigenvec_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pca_metadata <- read.delim(
    pca_metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_metadata <- c("Genotype", "origin")
  if (!all(required_metadata %in% names(pca_metadata))) {
    stop("PCA metadata must contain: Genotype and origin.", call. = FALSE)
  }
  if (ncol(pca) < 4) {
    stop("The PCA eigenvector file must contain sample IDs, PC1, and PC2.", call. = FALSE)
  }

  pca_plot_data <- merge(pca_metadata, pca, by.x = "Genotype", by.y = "V2")
  pca_plot <- ggplot(pca_plot_data, aes(x = V3, y = V4, color = origin)) +
    geom_point(size = 4) +
    geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
    scale_color_lancet() +
    labs(x = "PC1", y = "PC2", color = NULL) +
    theme(
      legend.text = element_text(size = 16),
      legend.position = c(0.4, 0.9),
      legend.background = element_rect(),
      panel.background = element_rect(fill = "transparent", color = "grey"),
      text = element_text(family = "sans", size = 14)
    )

  ggsave(file.path(output_dir, "PCA.pdf"), pca_plot, width = 6, height = 5)
}

# -----------------------------------------------------------------------------
# 2. Population structure plot (ADMIXTURE K = 2-7)
# -----------------------------------------------------------------------------
if (run_population_structure) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(ggsci)
  })

  k_values <- 2:7
  reference_k <- 7
  fam_file <- paste0(plink_prefix, ".fam")
  q_files <- file.path(project_dir, sprintf("step3_pruned.%d.Q", k_values))
  assert_files(c(structure_metadata_file, fam_file, q_files))

  structure_metadata <- read.delim(
    structure_metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("SampleID", "origin") %in% names(structure_metadata))) {
    stop("Structure metadata must contain: SampleID and origin.", call. = FALSE)
  }

  fam <- read.table(fam_file, stringsAsFactors = FALSE)
  fam_names <- fam$V2

  read_q <- function(k) {
    q_file <- file.path(project_dir, sprintf("step3_pruned.%d.Q", k))
    q_matrix <- read.table(q_file)
    if (nrow(q_matrix) != length(fam_names)) {
      stop("The number of rows in ", q_file, " does not match the FAM file.", call. = FALSE)
    }
    colnames(q_matrix) <- paste0("Anc", seq_len(k))

    tibble(SampleID = fam_names) |>
      bind_cols(as_tibble(q_matrix)) |>
      left_join(structure_metadata, by = "SampleID") |>
      filter(!is.na(origin)) |>
      mutate(K = factor(k, levels = k_values))
  }

  structure_long <- map_dfr(k_values, read_q) |>
    pivot_longer(starts_with("Anc"), names_to = "Anc", values_to = "prop")

  # Order samples within each origin using ancestry proportions at K = 7.
  order_table <- structure_long |>
    filter(K == reference_k) |>
    pivot_wider(names_from = Anc, values_from = prop) |>
    group_by(SampleID, origin) |>
    summarise(across(starts_with("Anc"), first), .groups = "drop")

  order_by_origin <- order_table |>
    group_split(origin) |>
    map(function(group_data) {
      if (nrow(group_data) <= 1) {
        ordered_ids <- group_data$SampleID
      } else {
        ancestry_matrix <- as.matrix(select(group_data, starts_with("Anc")))
        ordered_ids <- group_data$SampleID[
          hclust(dist(ancestry_matrix), method = "average")$order
        ]
      }
      tibble(SampleID = ordered_ids, origin = first(group_data$origin))
    }) |>
    bind_rows()

  structure_long$SampleID <- factor(
    structure_long$SampleID,
    levels = order_by_origin$SampleID
  )

  origin_counts <- structure_long |>
    distinct(SampleID, origin) |>
    count(origin)
  origin_labels <- setNames(
    paste0(origin_counts$origin, " (n=", origin_counts$n, ")"),
    origin_counts$origin
  )

  structure_plot <- ggplot(
    structure_long,
    aes(x = SampleID, y = prop, fill = Anc)
  ) +
    geom_col(width = 1) +
    facet_grid(
      K ~ origin,
      scales = "free_x",
      space = "free_x",
      labeller = labeller(origin = origin_labels, K = label_both)
    ) +
    scale_fill_nejm() +
    labs(x = "Samples", y = "Ancestry proportion", fill = NULL) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = grid::unit(0.8, "lines"),
      strip.text = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)
    )

  ggsave(
    file.path(output_dir, "population_structure.pdf"),
    structure_plot,
    width = 12,
    height = 6
  )
}

# -----------------------------------------------------------------------------
# 3. Neighbor-joining tree based on IBS distance
# -----------------------------------------------------------------------------
if (run_nj_tree) {
  suppressPackageStartupMessages({
    library(SNPRelate)
    library(ape)
    library(ggplot2)
    library(ggtree)
    library(dplyr)
    library(ggsci)
  })

  maf_cutoff <- 0.05
  missing_cutoff <- 0.10
  ld_window_snps <- 500
  ld_r2_cutoff <- 0.20
  n_threads <- 4
  rebuild_gds <- FALSE

  bed_file <- paste0(plink_prefix, ".bed")
  bim_file <- paste0(plink_prefix, ".bim")
  fam_file <- paste0(plink_prefix, ".fam")
  assert_files(c(bed_file, bim_file, fam_file, structure_metadata_file))

  if (rebuild_gds && file.exists(gds_file)) {
    unlink(gds_file)
  }
  if (!file.exists(gds_file)) {
    snpgdsBED2GDS(
      bed.fn = bed_file,
      bim.fn = bim_file,
      fam.fn = fam_file,
      out.gdsfn = gds_file,
      cvt.chr = "char",
      cvt.snpid = "auto",
      verbose = TRUE
    )
  }

  local({
    genofile <- snpgdsOpen(gds_file)
    on.exit(snpgdsClose(genofile), add = TRUE)

    filtered_snps <- snpgdsSelectSNP(
      genofile,
      maf = maf_cutoff,
      missing.rate = missing_cutoff,
      autosome.only = FALSE,
      remove.monosnp = TRUE,
      verbose = TRUE
    )
    if (length(filtered_snps) == 0) {
      stop("No SNPs remained after MAF and missing-rate filtering.", call. = FALSE)
    }

    set.seed(2025)
    pruned_snps <- snpgdsLDpruning(
      genofile,
      slide.max.n = ld_window_snps,
      slide.max.bp = Inf,
      ld.threshold = ld_r2_cutoff,
      maf = maf_cutoff,
      missing.rate = missing_cutoff,
      autosome.only = FALSE,
      remove.monosnp = TRUE,
      method = "corr",
      verbose = TRUE
    )
    snps_for_tree <- unlist(pruned_snps, use.names = FALSE)
    if (length(snps_for_tree) == 0) {
      warning("LD pruning retained no SNPs; using the filtered SNP set.")
      snps_for_tree <- filtered_snps
    }

    ibs <- snpgdsIBS(
      genofile,
      snp.id = snps_for_tree,
      num.thread = n_threads,
      autosome.only = FALSE,
      verbose = TRUE
    )
    rownames(ibs$ibs) <- colnames(ibs$ibs) <- ibs$sample.id
    tree <- ape::nj(as.dist(1 - ibs$ibs))
    tree$tip.label <- make.unique(tree$tip.label)

    metadata <- read.delim(
      structure_metadata_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) |>
      distinct(SampleID, origin)
    if (!all(c("SampleID", "origin") %in% names(metadata))) {
      stop("NJ-tree metadata must contain: SampleID and origin.", call. = FALSE)
    }

    keep_ids <- intersect(tree$tip.label, metadata$SampleID)
    if (length(keep_ids) < 3) {
      stop("Fewer than three tree samples matched the metadata.", call. = FALSE)
    }
    tree_subset <- ape::keep.tip(tree, keep_ids)
    tip_metadata <- metadata |>
      filter(SampleID %in% tree_subset$tip.label) |>
      transmute(label = SampleID, origin = factor(origin))

    origin_levels <- levels(tip_metadata$origin)
    palette_values <- ggsci::pal_lancet()(max(3, length(origin_levels)))
    palette_values <- setNames(
      palette_values[seq_along(origin_levels)],
      origin_levels
    )

    nj_plot <- ggtree(tree_subset, layout = "circular", branch.length = "none", size = 0.35) %<+%
      tip_metadata +
      geom_tippoint(aes(color = origin), size = 2, alpha = 0.95) +
      geom_tiplab2(
        aes(label = label, color = origin),
        size = 1.6,
        offset = 0.004,
        linesize = 0.2,
        show.legend = FALSE
      ) +
      scale_color_manual(values = palette_values, name = NULL, na.value = "grey70") +
      ggtitle("Population structure: NJ tree (IBS distance)") +
      theme(
        legend.position = c(0.55, 0.45),
        legend.justification = c(0.5, 0.5),
        legend.background = element_rect(
          fill = scales::alpha("white", 0.85),
          color = NA
        ),
        legend.text = element_text(size = 9)
      )

    ape::write.tree(tree_subset, file = file.path(output_dir, "NJ_tree.nwk"))
    ggsave(file.path(output_dir, "NJ_tree.pdf"), nj_plot, width = 8, height = 8)
    ggsave(file.path(output_dir, "NJ_tree.png"), nj_plot, width = 8, height = 8, dpi = 300)
  })
}

# -----------------------------------------------------------------------------
# 4. Trait histogram and density curve
# -----------------------------------------------------------------------------
if (run_histogram) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })

  histogram_trait <- "aperture width"
  histogram_binwidth <- 0.15
  assert_files(histogram_file)
  histogram_data <- read.csv(
    histogram_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c("Parameter", "final_Length", "treatment")
  if (!all(required_columns %in% names(histogram_data))) {
    stop(
      "Histogram data must contain: Parameter, final_Length, and treatment.",
      call. = FALSE
    )
  }

  histogram_data <- histogram_data |>
    filter(Parameter == histogram_trait, is.finite(final_Length)) |>
    mutate(
      treatment = recode(
        tolower(as.character(treatment)),
        "control" = "Control",
        "heat" = "Heat stress",
        "treat" = "Heat stress",
        "treatment" = "Heat stress"
      ),
      treatment = factor(treatment, levels = c("Control", "Heat stress"))
    ) |>
    filter(!is.na(treatment))
  if (nrow(histogram_data) == 0) {
    stop("No observations remained for the selected histogram trait.", call. = FALSE)
  }

  binwidth <- max(histogram_binwidth, 0.10)
  histogram_plot <- ggplot(
    histogram_data,
    aes(x = final_Length, fill = treatment, color = treatment)
  ) +
    geom_histogram(binwidth = binwidth, alpha = 0.5, position = "identity") +
    geom_density(
      aes(y = after_stat(density * n * binwidth)),
      bw = binwidth,
      adjust = 2,
      fill = "transparent",
      linewidth = 1.2,
      trim = TRUE
    ) +
    scale_fill_manual(values = c("Control" = "#08479A", "Heat stress" = "#df356b")) +
    scale_color_manual(values = c("Control" = "#08479A", "Heat stress" = "#df356b")) +
    labs(
      x = expression("Aperture width (" * mu * "m)"),
      y = "Count",
      fill = NULL,
      color = NULL
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      text = element_text(family = "sans", size = 14),
      panel.background = element_rect(fill = "transparent", color = "black"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 0, size = 14, hjust = 0.5)
    )

  ggsave(
    file.path(output_dir, "aperture_width_histogram.pdf"),
    histogram_plot,
    width = 5,
    height = 5,
    dpi = 720
  )
}

# -----------------------------------------------------------------------------
# 5. Genome-wide association analysis with FarmCPU
# -----------------------------------------------------------------------------
if (run_gwas) {
  suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(rMVP)
    library(bigmemory)
  })

  assert_files(c(gwas_vcf_file, gwas_phenotype_file))

  vcf <- readVcf(gwas_vcf_file, genome = "Ta")
  phenotype_input <- read.table(
    gwas_phenotype_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  phenotype_ids <- as.character(phenotype_input[[1]])
  common_samples <- intersect(phenotype_ids, samples(header(vcf)))
  if (length(common_samples) < 3) {
    stop("Fewer than three phenotype samples matched the VCF.", call. = FALSE)
  }

  filtered_vcf_file <- file.path(output_dir, "filtered_samples.vcf")
  filtered_phenotype_file <- file.path(output_dir, "filtered_phenotype.txt")
  writeVcf(vcf[, common_samples], filename = filtered_vcf_file)
  filtered_phenotype <- phenotype_input[
    match(common_samples, phenotype_ids),
    ,
    drop = FALSE
  ]
  write.table(
    filtered_phenotype,
    filtered_phenotype_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  local({
    previous_directory <- setwd(output_dir)
    on.exit(setwd(previous_directory), add = TRUE)

    rMVP::MVP.Data(
      fileVCF = filtered_vcf_file,
      filePhe = filtered_phenotype_file,
      sep.phe = "\t",
      MAF = 0.05,
      MISS = 0.30,
      fileKin = TRUE,
      filePC = TRUE,
      out = gwas_prefix
    )

    genotype <- bigmemory::attach.big.matrix(paste0(gwas_prefix, ".geno.desc"))
    phenotype <- read.table(paste0(gwas_prefix, ".phe"), header = TRUE)
    map <- read.table(paste0(gwas_prefix, ".geno.map"), header = TRUE)
    kinship <- bigmemory::attach.big.matrix(paste0(gwas_prefix, ".kin.desc"))

    rMVP::MVP(
      phe = phenotype,
      geno = genotype,
      map = map,
      K = kinship,
      nPC.FarmCPU = 3,
      vc.method = "BRENT",
      method.bin = "static",
      method = "FarmCPU",
      threshold = 0.05,
      file.output = c("pmap", "pmap.signal", "log", "plot")
    )
  })
}

# -----------------------------------------------------------------------------
# 6. Genome-wide LD decay
# -----------------------------------------------------------------------------
if (run_ld_decay) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
  })

  # Example PLINK command used to generate the input:
  # plink --bfile step3_pruned --r2 --ld-window 99999 \
  #   --ld-window-kb 5000 --ld-window-r2 0 --out LD_genomewide

  ld_decay_plot <- function(
      ld_file,
      max_kb = 5000,
      bin_kb = 200,
      span = 0.35,
      r2_threshold = 0.20,
      out_pdf = "LD_decay_smooth.pdf") {
    assert_files(ld_file)
    ld <- data.table::fread(ld_file)
    required_columns <- c("BP_A", "BP_B", "R2")
    if (!all(required_columns %in% names(ld))) {
      stop("LD input must contain: BP_A, BP_B, and R2.", call. = FALSE)
    }

    # Retain intra-chromosomal SNP pairs when chromosome columns are available.
    if (all(c("CHR_A", "CHR_B") %in% names(ld))) {
      ld <- ld[as.character(CHR_A) == as.character(CHR_B)]
    }
    ld[, distance_kb := abs(BP_B - BP_A) / 1000]
    ld <- ld[distance_kb <= max_kb & is.finite(R2)]
    if (nrow(ld) == 0) {
      stop("No LD pairs remained after distance and R2 filtering.", call. = FALSE)
    }

    ld_binned <- ld[
      , .(mean_r2 = mean(R2, na.rm = TRUE)),
      by = .(distance_bin = floor(distance_kb / bin_kb) * bin_kb)
    ]
    data.table::setorder(ld_binned, distance_bin)
    if (nrow(ld_binned) < 3) {
      stop("At least three distance bins are required for LOESS fitting.", call. = FALSE)
    }

    loess_fit <- stats::loess(
      mean_r2 ~ distance_bin,
      data = ld_binned,
      span = span
    )
    ld_binned[, loess_r2 := pmax(
      0,
      as.numeric(stats::predict(loess_fit, distance_bin))
    )]

    distance_at_threshold <- NA_real_
    valid <- is.finite(ld_binned$loess_r2)
    x <- ld_binned$distance_bin[valid]
    y <- ld_binned$loess_r2[valid]
    crossings <- which(y[-length(y)] >= r2_threshold & y[-1] <= r2_threshold)
    if (length(crossings) > 0) {
      i <- crossings[[1]]
      distance_at_threshold <- x[i] +
        (r2_threshold - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
    }

    ld_plot <- ggplot(ld_binned, aes(x = distance_bin, y = mean_r2)) +
      geom_point(size = 1, alpha = 0.6) +
      geom_line(linewidth = 0.5, alpha = 0.5) +
      geom_line(aes(y = loess_r2), colour = "red", linewidth = 1) +
      geom_hline(
        yintercept = r2_threshold,
        linetype = "dashed",
        colour = "grey40"
      ) +
      labs(
        x = "Distance (kb)",
        y = expression(mean ~ r^2),
        title = "Genome-wide LD decay in wheat"
      ) +
      theme_classic(base_size = 12) +
      theme(
        panel.border = element_rect(fill = NA, colour = "grey50"),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    if (is.finite(distance_at_threshold)) {
      ld_plot <- ld_plot +
        geom_vline(
          xintercept = distance_at_threshold,
          linetype = "dotted",
          colour = "grey40"
        ) +
        annotate(
          "text",
          x = distance_at_threshold,
          y = r2_threshold,
          label = sprintf(
            "~%d kb at r^2 = %.2f",
            round(distance_at_threshold),
            r2_threshold
          ),
          vjust = -0.6,
          hjust = -0.1,
          size = 3.2
        )
    }

    ggsave(out_pdf, ld_plot, width = 7.5, height = 5.5, limitsize = FALSE)
    invisible(list(
      plot = ld_plot,
      summary = ld_binned,
      distance_at_threshold_kb = distance_at_threshold
    ))
  }

  ld_result <- ld_decay_plot(
    ld_file = ld_file,
    max_kb = 5000,
    bin_kb = 200,
    span = 0.35,
    r2_threshold = 0.20,
    out_pdf = file.path(output_dir, "LD_decay_smooth.pdf")
  )
}
