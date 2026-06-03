# =============================================================================
# CITE-seq ADT analysis: GC-Tfh-SAP vs CM Pre-non-Tfh
# =============================================================================
#
# Produces:
#   1. Bar plot of mean CCR3 (CD193) ADT abundance with significance stars
#      and log2FC annotation
#   2. Volcano plot of differential ADT abundance, labeling CCR3, PDCD1, BTLA
#
# Inputs (edit paths in the CONFIG section):
#   - CD4 Seurat object (CITE-seq, with ADT assay and predicted.celltype.l2)
#   - matched_full CSV mapping ADT ids to gene symbols
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tibble)

# -----------------------------------------------------------------------------
# Paths -- edit these for your environment
# -----------------------------------------------------------------------------
cd4_rds       <- "/ix/djishnu/Alisa/Tfh/testing_for_github/Data/cd4_subset_azimuth.rds"
matched_csv   <- "adt_gene_matches_full.csv"
output_dir    <- "outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
cd4_cells    <- readRDS(cd4_rds)
matched_full <- read.csv(matched_csv)

# Subset to the two cell types of interest
subset_cells <- subset(cd4_cells,
                       subset = predicted.celltype.l2 %in%
                         c("GC-Tfh-SAP", "CM Pre-non-Tfh"))
DefaultAssay(subset_cells) <- "ADT"
Idents(subset_cells) <- "predicted.celltype.l2"

# -----------------------------------------------------------------------------
# Differential ADT abundance (Wilcoxon)
# -----------------------------------------------------------------------------
adt_de <- FindMarkers(
  subset_cells,
  ident.1         = "GC-Tfh-SAP",
  ident.2         = "CM Pre-non-Tfh",
  assay           = "ADT",
  test.use        = "wilcox",
  logfc.threshold = 0,
  min.pct         = 0
) %>%
  rownames_to_column("protein") %>%
  mutate(
    p_val_adj   = p.adjust(p_val, method = "BH"),
    significant = p_val_adj < 0.05,
    direction   = case_when(
      avg_log2FC > 0 & significant ~ "Up in GC-Tfh-SAP",
      avg_log2FC < 0 & significant ~ "Up in CM Pre-non-Tfh",
      TRUE                         ~ "NS"
    )
  ) %>%
  left_join(
    matched_full %>% mutate(protein = gsub("_", "-", id)) %>%
      select(protein, matched_gene),
    by = "protein"
  )

# -----------------------------------------------------------------------------
# Plot 1: Bar plot of CCR3 (CD193) mean ADT expression with significance
# -----------------------------------------------------------------------------
ccr3_protein <- "CD193-(CCR3)"

# Extract CCR3 expression per cell
ccr3_expr <- data.frame(
  expression = as.numeric(GetAssayData(subset_cells, assay = "ADT", layer = "data")[ccr3_protein, ]),
  cell_type  = subset_cells$predicted.celltype.l2
)

# Per-group mean and SE
ccr3_summary <- ccr3_expr %>%
  group_by(cell_type) %>%
  summarise(mean = mean(expression),
            se   = sd(expression) / sqrt(n()),
            .groups = "drop")

# Stats from adt_de
ccr3_stats <- adt_de %>% filter(protein == ccr3_protein)
fc   <- round(ccr3_stats$avg_log2FC, 2)
padj <- ccr3_stats$p_val_adj
sig_label <- if (padj < 0.0001) "****" else if (padj < 0.001) "***" else if (padj < 0.01) "**" else if (padj < 0.05) "*" else "ns"

# Bracket / star y positions
y_max     <- max(ccr3_summary$mean + ccr3_summary$se) * 1.15
y_bracket <- y_max * 0.98
y_star    <- y_max * 1.02

p_bar <- ggplot(ccr3_summary, aes(x = cell_type, y = mean, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  annotate("segment", x = 1,   xend = 2, y = y_bracket, yend = y_bracket,        color = "black") +
  annotate("segment", x = 1,   xend = 1, y = y_bracket, yend = y_bracket * 0.97, color = "black") +
  annotate("segment", x = 2,   xend = 2, y = y_bracket, yend = y_bracket * 0.97, color = "black") +
  annotate("text",    x = 1.5,           y = y_star,    label = sig_label, size = 6) +
  scale_fill_manual(values = c("GC-Tfh-SAP" = "#F8766D", "CM Pre-non-Tfh" = "#00BFC4")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  theme_classic() +
  theme(axis.text.x     = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y     = element_text(size = 14),
        axis.title.y    = element_text(size = 16),
        plot.subtitle   = element_text(size = 16, color = "grey40"),
        legend.position = "none") +
  labs(x = "", y = "Mean abundance",
       subtitle = paste0("log2FC = ", fc))

ggsave(file.path(output_dir, "CCR3_abundance_barplot.pdf"),
       p_bar, width = 4, height = 5)

# -----------------------------------------------------------------------------
# Plot 2: Volcano plot of differential ADT abundance
# -----------------------------------------------------------------------------
label_proteins <- c("CCR3", "PDCD1", "BTLA")

top_label <- adt_de %>%
  filter(matched_gene %in% label_proteins) %>%
  mutate(label = matched_gene)

p_volc <- ggplot(adt_de, aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(data               = top_label,
                  aes(label          = label),
                  size               = 5,
                  max.overlaps       = 20,
                  box.padding        = 1,
                  point.padding      = 0.5,
                  min.segment.length = 0,
                  nudge_x            = 0.5,
                  show.legend        = FALSE) +
  scale_color_manual(values = c("Up in GC-Tfh-SAP"     = "#F8766D",
                                "Up in CM Pre-non-Tfh" = "#00BFC4",
                                "NS"                   = "grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5),  linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(x     = "Average log2 Fold Change\n(GC-Tfh-SAP vs CM Pre-non-Tfh)",
       y     = "-log10 adjusted p-value",
       color = "") +
  theme(legend.position = "top",
        axis.text   = element_text(size = 13),
        axis.title  = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(file.path(output_dir, "volcanoplot_ADT_top_DE.pdf"),
       p_volc, width = 6, height = 6)