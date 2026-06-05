library(dplyr)
library(ggplot2)
library(effsize)

source("CalcCliffDelta_Helper.R")

Tfh_Tfh2_subset <- readRDS('/ix/djishnu/Alisa/Tfh/ForPaper/Resubmission/groom_comparisons/multiome_tonsil_atalas_tfh2_analysis/Tfh_Tfh2_multiome_subset_nocycT.rds')
CM_Tfh2_subset <- readRDS('/ix/djishnu/Alisa/Tfh/ForPaper/Resubmission/groom_comparisons/multiome_tonsil_atalas_tfh2_analysis/CM_Tfh2_multiome_subset_nocyclT.rds')

plot_1 <- DimPlot(CM_Tfh2_subset, reduction = "umap", group.by = "Tfh_final", label = TRUE, cols = c("red", "blue"))
plot_1

plot_2 <- DimPlot(Tfh_Tfh2_subset, reduction = "umap", group.by = "Tfh_final", label = TRUE, cols = c("red", "blue"))
plot_2

slide_boxplot_CM <- readRDS('slide_output/CM_Tfh2_SLIDECV_boxplot_data.rds')
write.csv(slide_boxplot_CM, "CM_slide_boxplot_aucs.csv", row.names = FALSE)

mean_auc_CM <- aggregate(auc ~ method, data = slide_boxplot_CM, FUN = mean)
wilcox_result_CM <- wilcox.test(auc ~ method, data = slide_boxplot_CM)
print(wilcox_result_CM)
cat("p-value:", wilcox_result_CM$p.value, "\n")

sig_genes_CM <- readRDS('slide_output/CM_Tfh2_plotSigGenes_data.RDS')
sig_genes_marginal_CM <- sig_genes_CM[sig_genes_CM$is_marginal != 'plain',]


z_CM <- as.matrix(read.csv("slide_output/CM_Tfh2_z_matrix.csv",
                           row.names = 1))
y_CM <- as.matrix(read.csv("slide_output/CM_Tfh2_y.csv",
                           row.names = 1))
SLIDE_res_CM <- readRDS("slide_output/CM_Tfh2_SLIDE_LFs.rds")


comb_CM <- GetPairwiseComb(y_CM)
sig_z_idx_CM <- c(67, 60) # 87
#sig_z_idx_CM <- SLIDE_res_CM$SLIDE_res_CM$marginal_vars # for this context, only do standalone
sig_z_CM <- z_CM[, sig_z_idx_CM]
SLIDE_cd_CM <- CalcCliffDelta(sig_z_CM, y_CM, comb_CM, sig_idx = sig_z_idx_CM)
write.csv(SLIDE_cd_CM, "CliffDelta_5d.csv")



slide_boxplot_Tfh <- readRDS('slide_output/Tfh_Tfh2_SLIDECV_boxplot_data.rds')
write.csv(slide_boxplot_Tfh, "Tfh_slide_boxplot_aucs.csv", row.names = FALSE)

mean_auc_Tfh <- aggregate(auc ~ method, data = slide_boxplot_Tfh, FUN = mean)
wilcox_result_Tfh <- wilcox.test(auc ~ method, data = slide_boxplot_Tfh)
print(wilcox_result_Tfh)
cat("p-value:", wilcox_result_Tfh$p.value, "\n")


sig_genes_Tfh <- readRDS('slide_output/Tfh_Tfh2_plotSigGenes_data.RDS')
sig_genes_marginal_Tfh <- sig_genes_Tfh[sig_genes_Tfh$is_marginal != 'plain',]

z_Tfh <- as.matrix(read.csv("slide_output/Tfh_Tfh2_z_matrix.csv",
                           row.names = 1))
y_Tfh <- as.matrix(read.csv('slide_output/Tfh_Tfh2_y.csv',
                           row.names = 1))
SLIDE_res_Tfh <- readRDS("slide_output/Tfh_Tfh2_SLIDE_LFs.rds")


comb_Tfh <- GetPairwiseComb(y_Tfh)
sig_z_idx_Tfh <- c(66,39) # 87
sig_z_Tfh <- z_Tfh[, sig_z_idx_Tfh]
SLIDE_cd_Tfh <- CalcCliffDelta(sig_z_Tfh, y_Tfh, comb_Tfh, sig_idx = sig_z_idx_Tfh)
write.csv(SLIDE_cd_Tfh, "CliffDelta_5h.csv")
