library(ggplot2)
library(ComplexUpset)

vinuesa_df  <-read.csv("../../Data/literature_sets/vinuesa_list.csv")
tfh_review_df  <- read.csv("../../Data/literature_sets/TFH_review_list.csv")

vinuesa <- as.list(vinuesa_df$'Genes')
tfh_review <- as.list(tfh_review_df$'Genes')

#Significant Genes
corenet_ext_df <- read.csv('../../Data/corenet_extended_full_genes.csv')
corenet_ext <- as.list(corenet_ext_df$'Genes')


#Random Gene Overlaps
corenet_ext_random_df <- read.csv("../../Data/corenet_extended_random_upset_human_nofd.csv")
corenet_ext_random <- as.list(corenet_ext_random_df$'Genes')


combine_gene_lists <- function(gene_lists, column_names) {
  # Ensure the number of column names matches the number of gene lists
  if (length(gene_lists) != length(column_names)) {
    stop("The number of column names must match the number of gene lists.")
  }
  
  # Get unique genes from all lists
  unique_genes <- unique(unlist(gene_lists))
  
  # Initialize a data frame with the unique genes
  df <- data.frame(Genes = unique_genes)
  
  # Add columns for each gene list indicating presence (1) or absence (0)
  for (i in seq_along(gene_lists)) {
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])
  }
  
  return(df)
}

fisher_pval <- function(df, method_col, random_col, reference_col, random_size) {
  a1 <- sum(df[[method_col]] == 1 & df[[reference_col]] == 1)
  m  <- sum(df[[method_col]] == 1)
  a2 <- sum(df[[random_col]] == 1 & df[[reference_col]] == 1)
  r  <- random_size
  mat <- matrix(c(a1, m - a1, a2, r - a2), nrow = 2, byrow = TRUE)
  fisher.test(mat, alternative = "greater")$p.value
}


gene_lists <- list(vinuesa, tfh_review, corenet_ext, corenet_ext_random)

colname <- c("Vinuesa 2016", "Hart 2022","CoreNet-Extended", 'CoreNet-Extended Random')
result <- combine_gene_lists(gene_lists, colname)
df <- subset(result, select = -Genes)


p <- upset(
  df, colname, 
  mode = 'inclusive_intersection',
  min_size = 0,
  intersections=list(
    c('CoreNet-Extended', 'Vinuesa 2016'),
    c('CoreNet-Extended Random', 'Vinuesa 2016'),
    c('CoreNet-Extended', 'Hart 2022'),
    c('CoreNet-Extended Random', 'Hart 2022')),
  #sort_sets = FALSE,
  sort_intersections = FALSE,
  base_annotations=list('Intersection size'=intersection_size( mode='inclusive_intersection')), #counts=FALSE, 
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(text = element_text(size = 20)),
      'overall_sizes' = theme(text = element_text(size = 20)),
      'set_sizes' = theme(text = element_text(size = 20)),
      'set_names' = theme(text = element_text(size = 20)),
      'intersection_size' = theme(  # ← This is key for the y-axis
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18)
      ),
      'default' = theme(text = element_text(size = 20))
    )
  ),

  set_size = FALSE,
  queries = list(
    upset_query(intersect = c('CoreNet-Extended', 'Vinuesa 2016'), color = "orange", fill = "orange"),
    upset_query(intersect = c('CoreNet-Extended', 'Hart 2022'), color = "orange", fill = "orange"),
    upset_query(intersect = c('CoreNet-Extended Random', 'Vinuesa 2016'), color = "gray", fill = "gray"),
    upset_query(intersect = c('CoreNet-Extended Random', 'Hart 2022'), color = "gray", fill = "gray")
  )
)

p

# Calculate Fisher's exact test p-values
# random_size = total number of genes in the random list (background size for comparison)
random_size <- length(corenet_ext_random)

# Fisher's p-value: CoreNet-Extended vs Vinuesa 2016 (compared to random)
p_vinuesa <- fisher_pval(
  df = result,
  method_col = "CoreNet-Extended",
  random_col = "CoreNet-Extended Random",
  reference_col = "Vinuesa 2016",
  random_size = random_size
)

# Fisher's p-value: CoreNet-Extended vs Hart 2022 (compared to random)
p_hart <- fisher_pval(
  df = result,
  method_col = "CoreNet-Extended",
  random_col = "CoreNet-Extended Random",
  reference_col = "Hart 2022",
  random_size = random_size
)

cat("Fisher's exact test p-values (one-sided, greater):\n")
cat("CoreNet-Extended vs Vinuesa 2016:", format.pval(p_vinuesa, digits = 4), "\n")
cat("CoreNet-Extended vs Hart 2022:   ", format.pval(p_hart, digits = 4), "\n")
