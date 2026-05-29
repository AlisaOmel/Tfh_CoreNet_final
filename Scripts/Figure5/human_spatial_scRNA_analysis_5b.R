library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(ggplot2)
library(BiocParallel)
library(rogme)
library(Cairo)


#Subset to groups of interest
###########################################################################
seurat_obj <- readRDS('../../Data/sc_rds/spatial_human_tonsil_no_naive.RSD')

seurat_obj@meta.data$cell_type
seurat_obj@reductions$spatial
FeaturePlot(seurat_obj,
            #features   = "CCR3",
            reduction  = "spatial",
            group.by = 'cell_type',
            order      = TRUE,
            pt.size    = 0.5,
            label      = TRUE,
            label.size = 3,
            repel      = TRUE) +
  scale_color_gradient(low = "lightgrey", high = "darkblue") +
  ggtitle("CCR3")

DimPlot(seurat_obj,
            reduction  = "spatial",
            group.by = 'cell_type',
            order      = TRUE,
            pt.size    = 0.5,
            label      = TRUE,
            label.size = 3,
            repel      = TRUE) +
  ggtitle("CCR3")
print(seurat_obj@meta.data)
unique(seurat_obj$cell_type)

subset_seurat <- subset(seurat_obj, subset = cell_type %in% c("T_follicular_helper", "T_CD4"))

subset_seurat[["RNA"]] <- subset_seurat[["originalexp"]]
#subset_seurat[["originalexp"]] <- NULL
DefaultAssay(subset_seurat) <- "RNA"

seurat_genes <- rownames(subset_seurat)
print(rownames(subset_seurat))

DimPlot(subset_seurat, reduction = 'spatial', group.by = "cell_type", label = TRUE) +
  ggtitle("UMAP - All Cell Types")


Assays(subset_seurat)
DefaultAssay(subset_seurat) <- "originalexp"

subset_seurat[["originalexp"]]@data[1:5, 1:5]

saveRDS(subset_seurat, file = "../../Data/sc_rds/spatial_human_tonsil_no_naive_subset.RSD")
subset_seurat <- readRDS("../../Data/sc_rds/spatial_human_tonsil_no_naive_subset.RSD")

#Load in the sets of interest - I had a bunch of different .csv, .gmt, .txt etc
###########################################################################
vinuesa_df  <-read.csv("../../Data/literature_sets/vinuesa_list.csv")
vinuesa <- as.list(vinuesa_df$'Genes')
vinuesa_genes <- intersect(unlist(vinuesa), seurat_genes)

tfh_review_df  <- read.csv("../../Data/literature_sets/TFH_review_list.csv")
tfh_review <- as.list(tfh_review_df$'Genes')
tfh_review_genes <- intersect(unlist(tfh_review), seurat_genes)

CoreNet_PPI_df  <- read.csv('../../Data/PPS_1_significant_gene_list.csv')
CoreNet_PPI <- as.list(CoreNet_PPI_df$'Genes')

CoreNet_PPI_genes <- intersect(unlist(CoreNet_PPI), seurat_genes)


###########################################################################

#Make a list of all gene lists of interest
###########################################################################
pathway_names <- c("CoreNet-PPI", " Vinuesa 2016", "Hart 2022")
gene_set <- list(
  CoreNet_PPI = CoreNet_PPI_genes,
  Vinuesa_2016 = vinuesa_genes,
  Hart_2022 = tfh_review_genes
)

###########################################################################

#Run ssGSEA
###########################################################################
ssg_seurat = runEscape(subset_seurat, method = "ssGSEA", gene.sets = gene_set, min.size=2,
                       new.assay.name = "escape.ssGSEA", normalize = FALSE)

DefaultAssay(ssg_seurat)  <- "escape.ssGSEA"
row.names(ssg_seurat@assays$escape.ssGSEA)

#subset_seurat$cell_type <- factor(subset_seurat$cell_type, levels = c("T_follicular_helper", "T_CD4"))

saveRDS(ssg_seurat, file = "../../Data/sc_rds/human_tonsil_spatial_ssGSEA_analysis_Tfh.RSD")

my_comparison <- list(c("Tfh", "CD4 T-Cell"))

# Convert underscores to hyphens to match ssGSEA output (Anything with a hyphen in the pathway gets changed by ssGSEA from _ to -)
valid_gene_sets <- names(gene_set)
converted_names <- gsub("_", "-", valid_gene_sets)

# Check which ones exist in the ssGSEA data
existing_sets <- rownames(GetAssayData(ssg_seurat, assay = "escape.ssGSEA", layer = "data"))
valid_gene_sets <- valid_gene_sets[converted_names %in% existing_sets]
###########################################################################


#Plot spatial
###########################################################################
plot1 <- FeaturePlot(
  ssg_seurat,
  features = "PPI",  
  reduction = "spatial",
  cols = c("cyan2","salmon")
)


plot2 <- DimPlot(subset_seurat, reduction = 'spatial', group.by = "cell_type", label = FALSE)

plot1


ggsave(filename ="outputs/Tfh_vs_CD4T_xy_PPI_brighter_reverse.pdf",
       plot = plot1,
       device = "pdf",
       width = 4,
       height = 4)
###########################################################################

# Prepare results list
results_list <- list()

# Get ssGSEA enrichment matrix once
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

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
    group = ssg_seurat$cell_type
  )
  df <- df[df$group %in% c("T_CD4", "T_follicular_helper"), ]
  df$group <- factor(df$group, levels = c("T_CD4", "T_follicular_helper"))
  
  # Skip if invalid data
  if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
    cat("  Skipping - no variability\n")
    next
  }
  
  # Run Cliff's Delta
  cliff <- tryCatch({
    cidv2(
      df[df$group == "T_CD4", "score"],
      df[df$group == "T_follicular_helper", "score"],
      alpha = 0.05
    )[c("d.hat", "p.value")]
  }, error = function(e) {
    cat("  Skipping - cidv2 failed\n")
    return(NULL)
  })
  
  if (is.null(cliff)) next
  
  # Compute group means
  avg_nontfh <- mean(df$score[df$group == "T_CD4"], na.rm = TRUE)
  avg_tfh <- mean(df$score[df$group == "T_follicular_helper"], na.rm = TRUE)
  
  # Save result for summary
  results_list[[gene_name]] <- data.frame(
    gene_set = gene_name,
    d.hat = cliff[["d.hat"]],
    p.value = cliff[["p.value"]],
    mean_nonTfh = avg_nontfh,
    mean_Tfh = avg_tfh
  )
  
  # Plot label
  delta_label <- paste0("Cliffs~Delta==", -1*round(cliff[["d.hat"]], 2))
  
  # then in your ggplot:
  annotate("text",
           x = 1.5,
           y = max(df$score, na.rm = TRUE),
           label = delta_label,
           parse = TRUE,
           vjust = 0.5,
           size = 6.5) 
  
  df <- as.data.frame(t(as.matrix(GetAssayData(ssg_seurat, assay = "escape.ssGSEA", layer = "data"))))
  df$cell_type <- ssg_seurat$cell_type
  df <- df[, c(converted_name, "cell_type")]
  colnames(df)[1] <- "score"
  
  df_for_points <- df %>%
    group_by(cell_type) %>%
    sample_frac(1.0) %>%  # Sample 10% of each group
    ungroup()
  
  # Rename levels first
  levels(df$cell_type)[levels(df$cell_type) == "T_follicular_helper"] <- "Tfh"
  levels(df$cell_type)[levels(df$cell_type) == "T_CD4"] <- "CD4 T-Cell"
  
  # Set the desired level order
  df$cell_type <- factor(df$cell_type, levels = c("Tfh", "CD4 T-Cell"))
  
  # Repeat for df_for_points
  levels(df_for_points$cell_type)[levels(df_for_points$cell_type) == "T_follicular_helper"] <- "Tfh"
  levels(df_for_points$cell_type)[levels(df_for_points$cell_type) == "T_CD4"] <- "CD4 T-Cell"
  
  df_for_points$cell_type <- factor(df_for_points$cell_type, levels = c("Tfh", "CD4 T-Cell"))
  
  gs <- ggplot(df, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_jitter(data = df_for_points, width = 0.2, size = 1, alpha = 0.1) +
    #scale_fill_brewer(palette = "Spectral") +
    theme_classic() +
    labs(x = "Cell Type", y = paste0(gene_name, " Enrichment Score")) +
    stat_compare_means(comparisons = my_comparison, 
                       method = "wilcox.test", 
                       label = "p.signif",
                       size=6.5) +
    annotate("text", x = 1.5, y = max(df$score, na.rm = TRUE), 
             label = delta_label, vjust = 0.5, size = 6.5) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  gs
  
  # Save plot
  ggsave(filename = paste0("outputs/Tfh_vs_CD4T_spatial_", gene_name, "boxplot_smaller.pdf"),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}

# Write summary CSV
if (length(results_list) > 0) {
  summary_df <- do.call(rbind, results_list)
  write.csv(summary_df,
            "outputs/cliffs_delta_summary_Tfh_vs_CD4T_spatial.csv",
            row.names = FALSE)
  cat("Summary CSV saved.\n")
} else {
  cat("No valid results to write to summary.\n")
}

##################################################################################
