# =============================================================================
# ssGSEA Pathway Analysis: Gut Tfh vs non-Tfh (Th17) in PLWHIV
# =============================================================================
#
# Description:
#   Performs single-sample GSEA (ssGSEA) on a Seurat object using the escape
#   package, comparing Tfh vs non-Tfh cells in PLWHIV samples. Computes
#   Cliff's delta and Wilcoxon p-values per pathway, then plots the top 15
#   pathways by absolute Cliff's delta as a heat strip.
#
# Inputs:
#   - Seurat .rds object with `HIV` and `cell_type` metadata columns
#     (cell_type values: "Tfh", "non-Tfh"; HIV values include "PLWHIV")
#
# Outputs:
#   - Seurat object with escape.ssGSEA assay
#   - per-pathway summary statistics
# =============================================================================

library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)
library(msigdbr)
library(effsize)  # for Cliff's delta
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(ggrepel)

# -----------------------------------------------------------------------------
# Paths -- edit these for your environment
# -----------------------------------------------------------------------------
#input_rds  <- "data/labels_GUT_Th17_Tfh_combo.rds"
input_rds  <- "/../../Data/labels_GUT_Th17_Tfh_combo.rds"
output_dir <- "outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

###########################################################################
combined <- readRDS(input_rds)
hiv.groups <- c('PLWHIV')
combined <- subset(combined, subset = HIV %in% hiv.groups)
###########################################################################

#Load msigDB pathway 
###########################################################################
GS.hallmark <- getGeneSets(library = "H")
GS.pid <- getGeneSets(library = "C2", subcategory = "CP:PID")
###########################################################################

#Run analysis
###########################################################################
#ssg_seurat = runEscape(combined, method = "ssGSEA", gene.sets = GS.hallmark, min.size=2,
#                       new.assay.name = "escape.ssGSEA", normalize = FALSE)

ssg_seurat = runEscape(combined, method = "ssGSEA", gene.sets = GS.pid, min.size=2,
                       new.assay.name = "escape.ssGSEA", normalize = FALSE)

# Reorder the levels of the cell_type metadata
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("non-Tfh", "Tfh"))

#saveRDS(ssg_seurat, file.path(output_dir, "scGSEA_hallmark_pathway_analysis_Th17_PLWHIV.RSD"))
saveRDS(ssg_seurat, file.path(output_dir, "scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD"))

#Alternatively Load it in:
#ssg_seurat <- readRDS(file.path(output_dir, "scGSEA_hallmark_pathway_analysis_Th17_PLWHIV.RSD"))
#ssg_seurat <- readRDS(file.path(output_dir, "scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD"))
###########################################################################



#Loop through all pathways and make a summary statistic and saves
########################################################################################################
# Extract enrichment matrix
enrichment_matrix <- as.data.frame(t(GetAssayData(ssg_seurat, assay = "escape.ssGSEA", layer = "data")))

# Add metadata for cell type classification
enrichment_matrix$cell_type <- ssg_seurat$cell_type

# Initialize results list
results <- list()

# Loop through each pathway
for (pathway in colnames(enrichment_matrix)[!colnames(enrichment_matrix) %in% "cell_type"]) {
  
  # Subset scores by cell type
  tfh_scores <- enrichment_matrix %>% filter(cell_type == "Tfh") %>% pull(pathway)
  non_tfh_scores <- enrichment_matrix %>% filter(cell_type == "non-Tfh") %>% pull(pathway)
  
  # Compute means
  mean_tfh <- mean(tfh_scores, na.rm = TRUE)
  mean_non_tfh <- mean(non_tfh_scores, na.rm = TRUE)
  
  # Cliff's delta
  cliff <- tryCatch({
    effsize::cliff.delta(tfh_scores, non_tfh_scores)$estimate
  }, error = function(e) NA)
  
  # Wilcoxon p-value
  pval <- tryCatch({
    wilcox.test(tfh_scores, non_tfh_scores)$p.value
  }, error = function(e) NA)
  
  # Append to results
  results[[pathway]] <- data.frame(
    Pathway = pathway,
    Mean_Tfh = mean_tfh,
    mean_non_tfh = mean_non_tfh,
    Cliffs_Delta = cliff,
    P_Value = pval
  )
}

# Combine and write CSV
summary_df <- do.call(rbind, results)

print(summary_df)

#write.csv(summary_df, file.path(output_dir, "hallmark_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv"), row.names = FALSE)
write.csv(summary_df, file.path(output_dir, "pid_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv"), row.names = FALSE)
########################################################################

#Heat Plot
########################################################################
#ssg_seurat <- readRDS(file.path(output_dir, "scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD"))
#summary_df <- read.csv(file.path(output_dir, "pid_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv"))

top15 <- summary_df %>%
  arrange(desc(abs(Cliffs_Delta))) %>%
  head(15) %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))

# Add a dummy column for the x-axis (single column)
top15$column <- "Cliff's Delta"

ggplot(top15, aes(x = column, y = Pathway, fill = Cliffs_Delta)) +
  geom_tile(color = "white", width = 0.3, height = 0.8) +  # narrower width
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, name = "Cliff's Δ"
  ) +
  theme_minimal() +
  labs(
    #title = "Top 15 Pathways by Absolute Cliff's Delta",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )

########################################################################
# Per-pathway boxplots for selected gene sets
########################################################################
# Toggle which library to load by uncommenting one of the pairs below.
ssg_seurat <- readRDS(file.path(output_dir, "scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD"))
summary_df <- read.csv(file.path(output_dir, "pid_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv"))

#ssg_seurat <- readRDS(file.path(output_dir, "scGSEA_hallmark_pathway_analysis_Th17_PLWHIV.RSD"))
#summary_df <- read.csv(file.path(output_dir, "hallmark_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv"))

# Specify gene set to plot
gene_name <- "PID_IL12_2PATHWAY"
gene_name <- "PID_IL27_PATHWAY"
gene_name <- "PID_IL23_PATHWAY"
gene_name <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
gene_name <- "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
gene_name <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"
converted_name <- gsub("_", "-", gene_name)

# Get ssGSEA enrichment matrix
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", layer = "data")

# Check if gene set exists
if (!(converted_name %in% rownames(ssgsea_mat))) {
  stop("Gene set not found in ssGSEA matrix.")
}

# Create enrichment score data frame
df <- data.frame(
  score = ssgsea_mat[converted_name, ],
  group = ssg_seurat$cell_type
)
df <- df[df$group %in% c("non-Tfh", "Tfh"), ]
df$group <- factor(df$group, levels = c("non-Tfh", "Tfh"))

# Stop if no variability
if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
  stop("No variability in enrichment scores.")
}

# Calculate Cliff's Delta
cliff <- tryCatch({
  cidv2(
    df[df$group == "non-Tfh", "score"],
    df[df$group == "Tfh", "score"],
    alpha = 0.05
  )[c("d.hat", "p.value")]
}, error = function(e) {
  stop("Cliff's Delta calculation failed.")
})

# Calculate group means
avg_non_tfh <- mean(df$score[df$group == "non-Tfh"], na.rm = TRUE)
avg_tfh <- mean(df$score[df$group == "Tfh"], na.rm = TRUE)

  # Prepare plotting data; rename non-Tfh -> Th17 for display
df_plot <- as.data.frame(t(as.matrix(GetAssayData(ssg_seurat, assay = "escape.ssGSEA", layer = "data"))))
df_plot$cell_type <- ssg_seurat$cell_type
df_plot <- df_plot[, c(converted_name, "cell_type")]
colnames(df_plot)[1] <- "score"
  
df_for_points <- df_plot %>%
  group_by(cell_type) %>%
  sample_frac(0.7) %>%
  ungroup()
  
df_plot$cell_type <- factor(df_plot$cell_type, levels = c("Tfh", "non-Tfh"))
levels(df_plot$cell_type)[levels(df_plot$cell_type) == "non-Tfh"] <- "Th17"
df_for_points$cell_type <- factor(df_for_points$cell_type, levels = c("Tfh", "non-Tfh"))
levels(df_for_points$cell_type)[levels(df_for_points$cell_type) == "non-Tfh"] <- "Th17"
  
# Plot label (sign flipped so delta is reported as Tfh vs non-Tfh)
delta_label <- paste0("Cliffs~Delta==", -1 * round(cliff[["d.hat"]], 2))
  
my_comparison <- list(c("Th17", "Tfh"))
  
  # Create plot
gs <- ggplot(df_plot, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_jitter(data = df_for_points, width = 0.2, size = 1, alpha = 0.1) +
    theme_classic() +
    labs(x = "Cell Type", y = paste0(gene_name, " Enrichment Score")) +
    stat_compare_means(comparisons = my_comparison,
                       method = "wilcox.test",
                       label = "p.signif",
                       size = 6.5) +
    annotate("text", x = 1.5, y = max(df_plot$score, na.rm = TRUE),
             label = delta_label, parse = TRUE, vjust = 0.5, size = 6.5) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18)
    )
  
  # Save to file
ggsave(
    filename = file.path(output_dir,
                         paste0("gut_Th17_PLWHIV_", gene_name, "_boxplot.pdf")),
    plot = gs,
    device = "pdf",
    width = 7,
    height = 5
  )
  
print(gs)
