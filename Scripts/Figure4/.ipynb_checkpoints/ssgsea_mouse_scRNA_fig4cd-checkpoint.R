library(zellkonverter)
library(Seurat)
library(SummarizedExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)
library(effsize)  # for Cliff's delta

###########################################################################
# PARAMETERS - Set these for your run
###########################################################################

# --- Input data ---
input_rds <- "/ix3/djishnu/Alisa/Tfh/mouse_scRNA/mouse_scRNA.rds"

# --- Cell group metadata ---
group_column  <- "cell_type"
group_levels  <- c("Tfh", "Resting T")     # order matters: first = left side / reference

# --- Gene set CSV paths ---
crotty_csv   <- "../../Data/literature_sets/crotty_list.csv"
vinuesa_csv  <- "../../Data/literature_sets/vinuesa_list.csv"
pps_csv      <- "../../Data/PPS_1_significant_gene_list.csv"
taiji_csv    <- "../../Data/GRN_genes_all_sets_152.csv"

# --- Human -> mouse translation table from BioMart ---
# CSV with two columns: "Mouse gene name", "Gene name" (human)
biomart_csv  <- "../../Data/mart_export.txt"

# --- Output paths ---
save_rds_path <- "outputs/scGSEA_modules_pathway_analysis_mouse_Tfh_RestingT.RSD"
output_dir    <- "outputs"
plot_prefix   <- "mouse_geyser_Tfh_vs_RestingT_"
summary_csv   <- file.path(output_dir, "cliffs_delta_summary_Tfh_vs_RestingT.csv")

###########################################################################
# END PARAMETERS
###########################################################################


# Load data
seurat_obj <- readRDS(input_rds)

# Subset to the two groups of interest
subset_seurat <- subset(seurat_obj, subset = !!sym(group_column) %in% group_levels)

# Mirror originalexp into RNA so other tools that expect "RNA" still work,
# but keep originalexp as the default for ssGSEA.
subset_seurat[["RNA"]] <- subset_seurat[["originalexp"]]
DefaultAssay(subset_seurat) <- "originalexp"

# Gene universe = rownames of the assay ssGSEA will actually use
seurat_genes <- rownames(subset_seurat[["originalexp"]])

########################################################
# Load in the sets of interest
########################################################

########################################################
# Load BioMart translation table (human -> mouse)
########################################################
biomart <- read.csv(biomart_csv, check.names = FALSE)
# Keep only rows with both a human and a mouse symbol
biomart <- biomart[
  biomart[["Mouse gene name"]] != "" & biomart[["Gene name"]] != "",
]
# Named vector: names are human symbols, values are mouse symbols
human_to_mouse <- setNames(
  biomart[["Mouse gene name"]],
  biomart[["Gene name"]]
)

translate_human_to_mouse <- function(human_genes, mapping, universe) {
  human_genes <- unique(unlist(human_genes))
  mouse_genes <- unname(mapping[human_genes])
  mouse_genes <- mouse_genes[!is.na(mouse_genes) & mouse_genes != ""]
  intersect(unique(mouse_genes), universe)
}

########################################################
# Load in the sets of interest (genes are HUMAN in the CSVs,
# translated to MOUSE here, then intersected with the assay)
########################################################
crotty_df    <- read.csv(crotty_csv)
crotty_genes <- translate_human_to_mouse(crotty_df$Genes, human_to_mouse, seurat_genes)

vinuesa_df    <- read.csv(vinuesa_csv)
vinuesa_genes <- translate_human_to_mouse(vinuesa_df$Genes, human_to_mouse, seurat_genes)

pps_df    <- read.csv(pps_csv)
pps_genes <- translate_human_to_mouse(pps_df$Genes, human_to_mouse, seurat_genes)

taiji_df    <- read.csv(taiji_csv)
taiji_genes <- translate_human_to_mouse(taiji_df$Genes, human_to_mouse, seurat_genes)

# Quick sanity check - how many genes survived translation + intersection?
cat("Gene set sizes after human->mouse translation and intersection with assay:\n")
cat("  Crotty 2020 :", length(crotty_genes), "\n")
cat("  CoreNet-PPI :", length(pps_genes), "\n")
cat("  CoreNet-GRN :", length(taiji_genes), "\n")
cat("  Vinuesa 2016:", length(vinuesa_genes), "\n")

########################################################
# Make a list of all gene lists of interest
########################################################
pathway_names <- c("Crotty 2020", "CoreNet-PPI", "CoreNet-GRN", "Vinuesa 2016")
gene_set <- list(
  Crotty_2020   = crotty_genes,
  CoreNet_PPI   = pps_genes,
  CoreNet_GRN   = taiji_genes,
  Vinuesa_2016  = vinuesa_genes
)

########################################################
# Run ssGSEA
########################################################

# IMPORTANT: ssGSEA runs on whatever the default assay is, so set it to the
# assay that holds the full expression matrix. For this mouse object that's
# `originalexp`. Switching to RNA here causes "No gene-sets meet the size
# threshold" because the gene index differs.
DefaultAssay(subset_seurat) <- "originalexp"

ssg_seurat <- runEscape(subset_seurat, method = "ssGSEA", gene.sets = gene_set, min.size = 2,
                        new.assay.name = "escape.ssGSEA", normalize = FALSE)

DefaultAssay(ssg_seurat) <- "escape.ssGSEA"

# Reorder the grouping factor for plotting
ssg_seurat[[group_column]][, 1] <- factor(ssg_seurat[[group_column]][, 1], levels = group_levels)

# Save File
if (nzchar(save_rds_path)) {
  saveRDS(ssg_seurat, file = save_rds_path)
}

# Comparison for p-value calculation in the plot
my_comparison <- list(group_levels)

# Convert underscores to hyphens to match ssGSEA output
valid_gene_sets <- names(gene_set)
converted_names <- gsub("_", "-", valid_gene_sets)

# Check which ones exist in the ssGSEA data
existing_sets <- rownames(GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data"))
valid_gene_sets <- valid_gene_sets[converted_names %in% existing_sets]

###########################################################################
# Loop through gene sets and plot each one, saving a summary CSV
###########################################################################
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results_list <- list()

# Get ssGSEA enrichment matrix once
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

# Pull the grouping vector once
group_vec <- ssg_seurat[[group_column]][, 1]

for (gene_name in valid_gene_sets) {
  cat("Plotting:", gene_name, "\n")
  
  converted_name <- gsub("_", "-", gene_name)
  
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
  
  # Run Cliff's Delta
  cliff <- tryCatch({
    cidv2(
      df[df$group == group_levels[1], "score"],
      df[df$group == group_levels[2], "score"],
      alpha = 0.05
    )[c("d.hat", "p.value")]
  }, error = function(e) {
    cat("  Skipping - cidv2 failed\n")
    return(NULL)
  })
  
  if (is.null(cliff)) next
  
  # Compute group means
  avg_g1 <- mean(df$score[df$group == group_levels[1]], na.rm = TRUE)
  avg_g2 <- mean(df$score[df$group == group_levels[2]], na.rm = TRUE)
  
  # Save result for summary
  results_list[[gene_name]] <- data.frame(
    gene_set  = gene_name,
    d.hat     = cliff[["d.hat"]],
    p.value   = cliff[["p.value"]],
    mean_g1   = avg_g1,
    mean_g2   = avg_g2
  )
  
  # Build plotting frame from full assay
  ssgsea_data <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")
  df <- as.data.frame(t(as.matrix(ssgsea_data)))
  df[[group_column]] <- group_vec
  df <- df[, c(converted_name, group_column)]
  colnames(df)[1] <- "score"
  
  df_for_points <- df %>%
    group_by(.data[[group_column]]) %>%
    sample_frac(1.0) %>%
    ungroup()
  
  delta_label <- paste0("Cliffs~Delta==", round(cliff[["d.hat"]], 2))
  
  gs <- ggplot(df, aes(x = .data[[group_column]], y = score, fill = .data[[group_column]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_jitter(data = df_for_points, width = 0.2, size = 1, alpha = 0.1) +
    theme_classic() +
    labs(x = "Cell Type", y = paste0(gene_name, " Enrichment Score")) +
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
  ggsave(filename = file.path(output_dir, paste0(plot_prefix, gene_name, "_boxplot.pdf")),
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