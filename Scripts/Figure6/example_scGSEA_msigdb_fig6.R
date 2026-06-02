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

###########################################################################
combined <- readRDS("../../Data/labels_GUT_Th17_Tfh_combo.rds")
hiv.groups <- c('PLWHIV')
combined <- subset(combined, subset = HIV %in% hiv.groups)
###########################################################################

#Load msigDB pathway 
###########################################################################
#dir.create(".cache/R/escape", recursive = TRUE, showWarnings = FALSE)
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

#saveRDS(ssg_seurat, "../outputs/scGSEA_hallmark_pathway_analysis_Th17_PLWHIV.RSD")
saveRDS(ssg_seurat, "../outputs/scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD")

#Alternatively Load it in:
#ssg_seurat <- readRDS("../outputs/scGSEA_hallmark_pathway_analysis_Th17_PLWHIV.RSD")
#ssg_seurat <- readRDS("../outputs/scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD")
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

#write.csv(summary_df, "../outputs/hallmark_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv", row.names = FALSE)
write.csv(summary_df, "../outputs/pid_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv", row.names = FALSE)
########################################################################

#Heat Plot
########################################################################
#ssg_seurat <- readRDS("../outputs/scGSEA_pid_pathway_analysis_Th17_PLWHIV.RSD")
#summary_df <- read.csv("../outputs/pid_ssGSEA_Gut_GCTfh_vs_Th17_summary_PLWHIV.csv")

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
    #title = "Top 15 Pathways by Absolute Cliff’s Delta",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )
