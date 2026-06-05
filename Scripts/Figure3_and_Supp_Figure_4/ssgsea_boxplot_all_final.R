library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)
library(effsize)  # for Cliff's delta
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(ggrepel)

###########################################################################
# PARAMETERS - Set these for your run
###########################################################################

# --- Pick a preset ---
# One of: "CM-PreTfh_vs_GC-Tfh-SAP", "Gut_CM-non-Tfh_vs_GC-Tfh", "Gut_Th17_vs_GC-Tfh"
preset <- "Gut_Th17_vs_GC-Tfh"

# --- Independent toggles (apply regardless of preset) ---
evaluate_extended  <- FALSE   # include the CoreNet-Extended set
filter_upregulated <- TRUE   # restrict gene sets to upregulated genes (LogChange > 0)

# --- Shared paths (same across all presets) ---
base_output_dir <- "outputs"
ppi_csv         <- "../../Data/PPS_1_significant_gene_list.csv"
vinuesa_csv     <- "../../Data/literature_sets/vinuesa_list.csv"
hart_csv        <- "../../Data/literature_sets/TFH_review_list.csv"
scenic_csv      <- "../../Data/scenic_all_genes.csv"
slide_csv       <- "../../Data/cm_pre_tfh_vs_gc_SLIDE_genes.csv"
extended_csv    <- "../../Data/union_propagation_genes.csv"
rna_logfc_csv   <- "../../Data/bulk_RNA/rna_earlyvsGC_log2FC.csv"

# --- Preset definitions ---
presets <- list(
  
  "CM-PreTfh_vs_GC-Tfh-SAP" = list(
    input_rds     = "/ix/djishnu/Alisa/Tfh/Tonsil_scdata/multiome/CM-PreTfh_GC-Tfh_SAP_TonsilAtlas_multiome.rds",
    filter_hiv    = FALSE,
    group_column  = "annotation_20230508",
    group_levels  = c("CM PreTfh", "GC-Tfh-SAP"),
    save_rds_path = "scGSEA_CM_PreTfhvsTfh_reference_analysis.RSD",
    plot_prefix   = "CM_PreTfh_vs_Tfh_"
  ),
  
  "Gut_CM-non-Tfh_vs_GC-Tfh" = list(
    input_rds     = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_nonTfh_Tfh_combo.rds",
    filter_hiv    = TRUE,
    group_column  = "cell_type",
    group_levels  = c("non-Tfh", "Tfh"),
    save_rds_path = "scGSEA_Gut_CM-non-Tfh_vs_GC-Tfh_analysis.RSD",
    plot_prefix   = "Gut_CMnon_vs_Tfh_"
  ),
  
  "Gut_Th17_vs_GC-Tfh" = list(
    input_rds     = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_Th17_Tfh_combo.rds",
    filter_hiv    = TRUE,
    group_column  = "cell_type",
    group_levels  = c("non-Tfh", "Tfh"),
    save_rds_path = "scGSEA_Gut_Th17_vs_GC-Tfh_analysis.RSD",
    plot_prefix   = "Gut_Th17_vs_Tfh_"
  )
)

# --- HIV filter settings (used when preset$filter_hiv is TRUE) ---
hiv_column <- "HIV"
hiv_groups <- c("PLWHIV")

###########################################################################
# Resolve preset into the variables the rest of the script uses
###########################################################################
if (!preset %in% names(presets)) {
  stop("Unknown preset: '", preset, "'. Must be one of: ",
       paste(names(presets), collapse = ", "))
}
cfg <- presets[[preset]]

input_rds     <- cfg$input_rds
filter_hiv    <- cfg$filter_hiv
group_column  <- cfg$group_column
group_levels  <- cfg$group_levels
save_rds_path <- cfg$save_rds_path
plot_prefix   <- cfg$plot_prefix
if (filter_upregulated) {
  plot_prefix <- paste0(plot_prefix, "GC_upregulated_")
}

output_dir  <- file.path(base_output_dir, preset)
summary_suffix <- if (filter_upregulated) "_GC_upregulated" else ""
summary_csv <- file.path(output_dir,
                         paste0("cliffs_delta_summary_", preset, summary_suffix, ".csv"))

###########################################################################
# END PARAMETERS - nothing below here should normally need changing
###########################################################################


# Load in data
combined <- readRDS(input_rds)

# Optional HIV filter
if (filter_hiv) {
  combined <- subset(combined, subset = !!sym(hiv_column) %in% hiv_groups)
}

# Diagnostic: show how many cells fall into each requested group
cat("\n--- Cell counts per group (after HIV filter if applied) ---\n")
group_meta <- combined[[group_column]][, 1]
cat("  Unique values in '", group_column, "': ",
    paste(unique(as.character(group_meta)), collapse = ", "), "\n", sep = "")
for (lvl in group_levels) {
  cat("  ", lvl, ": ", sum(group_meta == lvl, na.rm = TRUE), " cells\n", sep = "")
}
cat("-----------------------------------------------------------\n\n")

seurat_genes <- rownames(combined)

########################################################
# Load in the sets of interest
########################################################

# Modules
PPI_df  <- read.csv(ppi_csv)
PPI <- as.list(PPI_df$'Genes')
PPI_genes <- intersect(unlist(PPI), seurat_genes)

# Reference Sets
vinuesa_df  <- read.csv(vinuesa_csv)
vinuesa <- as.list(vinuesa_df$'Genes')
vinuesa_genes <- intersect(unlist(vinuesa), seurat_genes)

hart_df  <- read.csv(hart_csv)
hart_review <- as.list(hart_df$'Genes')
hart_genes <- intersect(unlist(hart_review), seurat_genes)

# Other analysis
scenic_tfs <- read.csv(scenic_csv)
scenic <- as.list(scenic_tfs$'Genes')
scenic <- intersect(unlist(scenic), seurat_genes)

slide_genes <- read.csv(slide_csv)
slide <- as.list(slide_genes$'names')
slide <- intersect(unlist(slide), seurat_genes)

# Optional CoreNet-Extended set
if (evaluate_extended) {
  extended_df <- read.csv(extended_csv)
  extended <- as.list(extended_df$'Genes')
  extended_genes <- intersect(unlist(extended), seurat_genes)
}

########################################################
# Optional: restrict to upregulated genes
########################################################
if (filter_upregulated) {
  rna_df <- read.csv(rna_logfc_csv)
  cat("--- filter_upregulated diagnostics ---\n")
  cat("  rna_logfc_csv columns:", paste(colnames(rna_df), collapse = ", "), "\n")
  
  # Resolve the gene-symbol column (try common names, fall back to first column)
  gene_col_candidates <- c("Genes", "Gene", "gene", "symbol", "Symbol", "GeneSymbol")
  gene_col <- intersect(gene_col_candidates, colnames(rna_df))
  if (length(gene_col) == 0) {
    gene_col <- colnames(rna_df)[1]
    warning("None of ", paste(gene_col_candidates, collapse = "/"),
            " found in ", rna_logfc_csv, "; falling back to first column: '", gene_col, "'")
  } else {
    gene_col <- gene_col[1]
  }
  
  # Resolve the logFC column (try common names)
  lfc_col_candidates <- c("LogChange", "logFC", "log2FC", "log2FoldChange",
                          "avg_log2FC", "avg_logFC")
  lfc_col <- intersect(lfc_col_candidates, colnames(rna_df))
  if (length(lfc_col) == 0) {
    stop("No logFC column found in ", rna_logfc_csv,
         ". Tried: ", paste(lfc_col_candidates, collapse = ", "),
         ". Available columns: ", paste(colnames(rna_df), collapse = ", "))
  }
  lfc_col <- lfc_col[1]
  
  cat("  Using gene column: '", gene_col, "', logFC column: '", lfc_col, "'\n", sep = "")
  
  rna_up      <- rna_df[!is.na(rna_df[[lfc_col]]) & rna_df[[lfc_col]] > 0, ]
  gc_up_genes <- as.character(rna_up[[gene_col]])
  cat("  Total genes in logFC table:        ", nrow(rna_df), "\n")
  cat("  Upregulated genes (logFC > 0):     ", length(gc_up_genes), "\n")
  
  # Sanity-check overlap *before* filtering, so we can detect ID mismatches
  # (e.g. one file has gene symbols, the other has Ensembl IDs)
  overlap_with_seurat <- length(intersect(gc_up_genes, seurat_genes))
  cat("  Overlap with seurat rownames:      ", overlap_with_seurat, "\n")
  if (overlap_with_seurat == 0) {
    stop("Zero overlap between gc_up_genes and seurat rownames - ",
         "check that both use the same gene-ID type (symbols vs Ensembl, etc).")
  }
  
  filter_one <- function(name, gs) {
    n_before <- length(gs)
    out      <- intersect(unlist(gs), gc_up_genes)
    n_after  <- length(out)
    cat("  ", format(name, width = 14), ": ", n_before, " -> ", n_after,
        " (", round(100 * n_after / max(n_before, 1), 1), "% retained)\n", sep = "")
    out
  }
  
  PPI_genes  <- filter_one("CoreNet-PPI",  PPI_genes)
  hart_genes <- filter_one("Hart",         hart_genes)
  if (evaluate_extended) {
    extended_genes <- filter_one("CoreNet-Ext", extended_genes)
  }
  
  vinuesa_genes <- list("STAT1", "STAT3","NFKB1","ASCL2", "IL21","IL21R","IL6","IL6R","IFN","IFNR","CD28",
                        "ICOSL","ICOS","ITCH","PI3K","AKT","LEF","TCF1","ERG2","ERG3","BCL6","MAF","NOTCH1",
                        "NOTCH2","IRF4","CXCR4","CXCR5","PDCD1")
  vinuesa_genes <- intersect(unlist(vinuesa_genes), seurat_genes)
  cat("--- end filter_upregulated diagnostics ---\n")
}

########################################################
# Make a list of all gene lists of interest
########################################################
pathway_names <- c("CoreNet-PPI", "Vinuesa-2016", "Hart-2022", 'SCENIC+', 'SLIDE')
gene_set <- list(
  CoreNet_PPI = PPI_genes,
  Vinuesa_2016 = vinuesa_genes,
  Hart_2022 = hart_genes,
  SCENIC = scenic,
  SLIDE = slide
)

if (evaluate_extended) {
  pathway_names <- c(pathway_names, "CoreNet-PPI-Extended")
  gene_set[["CoreNet_PPI_Extended"]] <- extended_genes
}

########################################################
# Run ssGSEA
########################################################

DefaultAssay(combined) <- "RNA"

# Print final gene-set sizes for sanity check
cat("\n--- Final gene set sizes going into runEscape ---\n")
cat("filter_upregulated =", filter_upregulated, "\n")
for (nm in names(gene_set)) {
  cat("  ", format(nm, width = 22), ": ", length(gene_set[[nm]]),
      " genes\n", sep = "")
}
cat("-------------------------------------------------\n\n")

ssg_seurat <- runEscape(combined, method = "ssGSEA", gene.sets = gene_set, min.size = 2,
                        new.assay.name = "escape.ssGSEA", normalize = FALSE)

DefaultAssay(ssg_seurat) <- "escape.ssGSEA"

# Reorder the grouping factor in the way you want it to show up on your plot
ssg_seurat[[group_column]][, 1] <- factor(ssg_seurat[[group_column]][, 1], levels = group_levels)

# Save File
if (nzchar(save_rds_path)) {
  # Append the upregulated flag so the cache doesn't get clobbered between modes
  if (filter_upregulated) {
    save_rds_path <- sub("\\.([^.]+)$", "_GC_upregulated.\\1", save_rds_path)
  }
  cat("Saving ssg_seurat to:", save_rds_path, "\n")
  saveRDS(ssg_seurat, file = save_rds_path)
}

# Gives the comparison for the p-value calculation in the plot
my_comparison <- list(group_levels)

# Convert underscores to hyphens to match ssGSEA output
# (Anything with a hyphen in the pathway gets changed by ssGSEA from _ to -)
valid_gene_sets <- names(gene_set)
converted_names <- gsub("_", "-", valid_gene_sets)

# Check which ones exist in the ssGSEA data
# Use GetAssayData (not @data) so this works for both Seurat v4 (Assay)
# and v5 (Assay5) objects -- v5 stores data in @layers$data, not @data.
existing_sets <- rownames(GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data"))
valid_gene_sets <- valid_gene_sets[converted_names %in% existing_sets]

###########################################################################
# Loop through gene sets and plot each one saving a summary.csv
###########################################################################
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results_list <- list()

# Get ssGSEA enrichment matrix once
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

# Pull the grouping vector once
group_vec <- ssg_seurat[[group_column]][, 1]

for (gene_set_name in valid_gene_sets) {
  cat("Plotting:", gene_set_name, "\n")
  
  converted_name <- gsub("_", "-", gene_set_name)
  
  # Check if gene set is present in assay
  if (!(converted_name %in% rownames(ssgsea_mat))) {
    cat("  Skipping - gene set not found\n")
    next
  }
  
  # Build data frame of enrichment scores and group labels
  df <- data.frame(
    score = ssgsea_mat[converted_name, ],
    group = group_vec
  )
  df <- df[df$group %in% group_levels, ]
  df$group <- factor(df$group, levels = group_levels)
  
  # Skip if invalid data
  if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
    cat("  Skipping - no variability\n")
    next
  }
  
  # Group sizes (diagnostic)
  g1_scores <- df[df$group == group_levels[1], "score"]
  g2_scores <- df[df$group == group_levels[2], "score"]
  cat("  Group sizes: ", group_levels[1], "=", length(g1_scores),
      ", ", group_levels[2], "=", length(g2_scores), "\n", sep = "")
  
  if (length(g1_scores) < 2 || length(g2_scores) < 2) {
    cat("  Skipping - one or both groups have <2 cells\n")
    next
  }
  
  # Run Cliff's Delta
  cliff <- tryCatch({
    cidv2(g1_scores, g2_scores, alpha = 0.05)[c("d.hat", "p.value")]
  }, error = function(e) {
    cat("  Skipping - cidv2 failed: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  })
  
  if (is.null(cliff)) next
  
  # Compute group means
  avg_g1 <- mean(df$score[df$group == group_levels[1]], na.rm = TRUE)
  avg_g2 <- mean(df$score[df$group == group_levels[2]], na.rm = TRUE)
  
  # Save result for summary
  results_list[[gene_set_name]] <- data.frame(
    gene_set = gene_set_name,
    d.hat = -cliff[["d.hat"]],
    p.value = cliff[["p.value"]],
    mean_nontfh = avg_g1,
    mean_Tfh = avg_g2
  )
  
  # Build plotting frame from full assay
  ssgsea_data <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")
  df <- as.data.frame(t(as.matrix(ssgsea_data)))
  df[[group_column]] <- group_vec
  df <- df[, c(converted_name, group_column)]
  colnames(df)[1] <- "score"
  
  df_for_points <- df %>%
    group_by(.data[[group_column]]) %>%
    sample_frac(0.7) %>%
    ungroup()
  
  delta_label <- paste0("Cliffs~Delta==", round(-cliff[["d.hat"]], 2))
  
  gs <- ggplot(df, aes(x = .data[[group_column]], y = score, fill = .data[[group_column]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_jitter(data = df_for_points, width = 0.2, size = 1, alpha = 0.1) +
    theme_classic() +
    labs(x = "Cell Type", y = paste0(gene_set_name, " Enrichment Score")) +
    stat_compare_means(
      comparisons = my_comparison,
      method = "wilcox.test",
      label = "p.signif",
      size = 6.5
    ) +
    annotate("text",
             x = 1.5,
             y = max(df$score, na.rm = TRUE),
             label = delta_label,
             parse = TRUE,
             vjust = 0.5,
             size = 6.5) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18)
    )
  
  # Save plot
  ggsave(filename = file.path(output_dir, paste0(plot_prefix, gene_set_name, "_boxplot.pdf")),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}

# Write summary CSV
if (length(results_list) > 0) {
  summary_df <- do.call(rbind, results_list)
  write.csv(summary_df, summary_csv, row.names = FALSE)
  cat("Summary CSV saved.\n")
} else {
  cat("No valid results to write to summary.\n")
}