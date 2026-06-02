# =============================================================================
# PPI Network Permutation Analysis: MSigDB Pathway Connectivity
# =============================================================================
#
# Description:
#   For each pathway in the chosen MSigDB collection (Hallmark or PID),
#   counts the number of unique pathway genes and module genes connected
#   in the PPI network at hop distance <= 1, then compares against 1000
#   random gene sets of the same size as the module list.
#
# Usage:
#   Rscript example_permutations_for_bulk_msigdb_analysis.R <database>
#     where <database> is "hallmark" or "pid"
#
#   If no argument is supplied, defaults to "hallmark".
#
# Inputs (edit paths in the CONFIG section):
#   - PPI network edge list (3 columns: Node1, Node2, hop.distance)
#   - Module gene list (tab-delimited with a `uGenes` column)
#
# Outputs:
#   - results_patternb_<database>_<random_>.csv  (in output_dir)
# =============================================================================

library(matrixStats)
library(readxl)
library(stringr)
library(msigdbr)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# Parse command-line argument: database choice
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
database <- if (length(args) >= 1) tolower(args[1]) else "hallmark"

# Map database name to msigdbr arguments
db_config <- switch(
  database,
  "hallmark" = list(coll = "H",  sub = NULL),
  "pid"      = list(coll = "C2", sub = "CP:PID"),
  stop("database must be 'hallmark' or 'pid'")
)

message("Running with database: ", database)

# -----------------------------------------------------------------------------
# Number of random permutations per pathway
# -----------------------------------------------------------------------------
perms <- 1000

# -----------------------------------------------------------------------------
# Paths -- edit these for your environment
# -----------------------------------------------------------------------------
#ppi_file    <- "data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.txt"
ppi_file    <- "/ix/djishnu/Alisa/Tfh/Network_analysis/data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.txt"
module_file <- "../../Data/pps_taiji_unique_genes_df.txt"
output_dir  <- "outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_csv <- file.path(output_dir,
                        paste0("results_CoreNet-DI_", database, "_", perms, ".csv"))

# -----------------------------------------------------------------------------
# Load PPI network
# -----------------------------------------------------------------------------
message("[", format(Sys.time(), "%H:%M:%S"), "] Loading PPI network: ", ppi_file)
netwk <- read.delim(ppi_file, header = FALSE, sep = ' ')
hop.distance <- rep(1, nrow(netwk))
df <- cbind.data.frame(netwk, hop.distance)
names(df) <- c("Node1", "Node2", "hop.distance")
message("  -> ", nrow(df), " edges loaded")

# -----------------------------------------------------------------------------
# Load module gene list
# -----------------------------------------------------------------------------
message("[", format(Sys.time(), "%H:%M:%S"), "] Loading module gene list: ", module_file)
selectedGenes <- read.delim(module_file)
Modules <- unique(selectedGenes$uGenes)
message("  -> ", length(Modules), " unique module genes")

# -----------------------------------------------------------------------------
# Load MSigDB gene sets for the chosen database
# -----------------------------------------------------------------------------
h_gene_sets_all <- msigdbr(species = "Homo sapiens")

# Identify the collection / subcollection columns (names vary by version)
coll_col <- intersect(c("gs_collection", "gs_cat"),    colnames(h_gene_sets_all))[1]
sub_col  <- intersect(c("gs_subcollection", "gs_subcat"), colnames(h_gene_sets_all))[1]

if (is.na(coll_col)) {
  stop("Could not find collection column in msigdbr output. Columns: ",
       paste(colnames(h_gene_sets_all), collapse = ", "))
}

# Filter to the requested collection
h_gene_sets <- h_gene_sets_all[h_gene_sets_all[[coll_col]] == db_config$coll, ]

# Further filter to the requested subcollection if specified
if (!is.null(db_config$sub)) {
  h_gene_sets <- h_gene_sets[h_gene_sets[[sub_col]] == db_config$sub, ]
}

GS_NAME <- unique(h_gene_sets$gs_name)
message("  -> ", length(GS_NAME), " pathways in ", database,
        " collection (", db_config$coll,
        if (!is.null(db_config$sub)) paste0(" / ", db_config$sub) else "",
        ")")

if (length(GS_NAME) == 0) {
  stop("No pathways matched. Check db_config and the values in column ", coll_col)
}

gene_col <- if ("gene_symbol" %in% colnames(h_gene_sets)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% colnames(h_gene_sets)) {
  "human_gene_symbol"
} else {
  stop("msigdbr output has neither `gene_symbol` nor `human_gene_symbol` column.")
}

# -----------------------------------------------------------------------------
# Main analysis: per-pathway connectivity vs random gene sets
# -----------------------------------------------------------------------------
M <- matrix(0, nrow = length(GS_NAME), ncol = 6)
colnames(M) <- c('PathwayName', "Genes", "Node1", "Random1", "Node2", "Random2")

t_start <- Sys.time()
message("[", format(t_start, "%H:%M:%S"), "] Starting analysis: ",
        length(GS_NAME), " pathways x ", perms, " permutations")

for (i in 1:length(GS_NAME)) {
  
  t_pathway <- Sys.time()
  
  gs <- GS_NAME[i]
  sub_gene_set <- subset.data.frame(h_gene_sets, gs_name == gs)
  TNF <- sub_gene_set[[gene_col]]
  
  TNFSubset1 <- subset(df, ((df$Node1 %in% Modules) & (df$Node2 %in% TNF)))
  TNFSubset2 <- subset(df, ((df$Node2 %in% Modules) & (df$Node1 %in% TNF)))
  
  TNFSubset12 <- cbind.data.frame(
    c(TNFSubset1$Node1, TNFSubset2$Node2),
    c(TNFSubset1$Node2, TNFSubset2$Node1),
    c(TNFSubset1$hop.distance, TNFSubset2$hop.distance)
  )
  names(TNFSubset12) <- c("Node1", "Node2", "hop.distance")
  
  test1 <- subset(TNFSubset12, hop.distance <= 1)
  
  M[i, 3] <- length(unique(test1$Node2))
  M[i, 5] <- length(unique(test1$Node1))
  M[i, 2] <- length(sub_gene_set[[gene_col]])
  M[i, 1] <- gs
  
  Pathways <- c()
  MODs <- c()
  
  for (j in 1:perms) {
    
    Random <- sample(setdiff(unique(c(df$Node1, df$Node2)), c(Modules, TNF)),
                     size = length(Modules))
    
    randomSubset1 <- subset(df, ((df$Node1 %in% TNF) & (df$Node2 %in% Random)))
    randomSubset2 <- subset(df, ((df$Node2 %in% TNF) & (df$Node1 %in% Random)))
    
    randomSubset12 <- cbind.data.frame(
      c(randomSubset1$Node1, randomSubset2$Node2),
      c(randomSubset1$Node2, randomSubset2$Node1),
      c(randomSubset1$hop.distance, randomSubset2$hop.distance)
    )
    names(randomSubset12) <- c("Node1", "Node2", "hop.distance")
    
    testRandom1 <- subset(randomSubset12, hop.distance <= 1)
    
    Pathways[j] <- length(unique(testRandom1$Node1))
    MODs[j]     <- length(unique(testRandom1$Node2))
    
    # Heartbeat every 250 permutations so very slow pathways still show life
    if (j %% 250 == 0) {
      message("    ... pathway ", i, "/", length(GS_NAME),
              " | permutation ", j, "/", perms)
    }
  }
  
  M[i, 4] <- mean(Pathways)
  M[i, 6] <- mean(MODs)
  
  elapsed_pathway <- as.numeric(difftime(Sys.time(), t_pathway, units = "secs"))
  elapsed_total   <- as.numeric(difftime(Sys.time(), t_start,   units = "mins"))
  message(sprintf("[%s] [%d/%d] %s | %.1fs | total %.1f min",
                  format(Sys.time(), "%H:%M:%S"),
                  i, length(GS_NAME), gs, elapsed_pathway, elapsed_total))
  flush.console()
}

message("[", format(Sys.time(), "%H:%M:%S"), "] Analysis complete (",
        round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
        " min total)")

write.csv(M, output_csv)
message("Wrote: ", output_csv)
