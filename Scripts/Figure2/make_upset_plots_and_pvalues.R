# =============================================================================
# Upset plots + p-values for DI permutation pipeline outputs.
# Reads from the new folder structure: <BASE_DIR>/<PREFIX>/<PREFIX>_<TAG>_*.csv
#
# Three pipelines are visualized separately, each for both species:
#   tag = "first_DI_then_size_match"  -> Fisher's exact (real & random same size)
#   tag = "first_size_match_then_DI"  -> Binomial proportions + Empirical permutation
#   tag = "sig_only_size_match"       -> Fisher's exact (no DI step; real & random same size)
#
# Two figure variants per (pipeline, species):
#   network    : CoreNet_PPI, CoreNet_GRN, CoreNet
#   just_scores: PPS_noprop_184, logFC_rna
#
# Species:
#   human  -> references Vinuesa (2016) and Hart (2022)
#   mouse  -> reference Crotty
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexUpset)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# CONFIG -- edit paths here
# -----------------------------------------------------------------------------
BASE_DIR    <- "../../Data/upset_plot_overlaps"
OUTPUT_DIR  <- "upset_outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Reference (gold-standard) gene lists
VINUESA_CSV    <- "../../Data/literature_sets/vinuesa_list.csv"
TFH_REVIEW_CSV <- "../../Data/literature_sets/TFH_review_list.csv"
CROTTY_CSV     <- "../../Data/literature_sets/crotty_list.csv"   

# Per-pipeline tag (folder + filename middle component)
TAG_DI_FIRST  <- "first_DI_then_size_match"   # pipeline 1: Fisher's
TAG_MM_FIRST  <- "first_size_match_then_DI"   # pipeline 2: Binomial + Empirical
TAG_SIG_ONLY  <- "sig_only_size_match"        # pipeline 3: Fisher's (no DI)

# Method definitions for the two plot variants.
# Each entry has:
#   prefix       : folder + filename prefix produced by the python pipeline
#   name         : display name shown in the plot
#   color        : highlight color for non-random intersections
#   sig_input    : path to the ORIGINAL significant-gene CSV (pre-DI). Only used
#                  when tag == "sig_only_size_match", because that pipeline
#                  compares bare sig_genes (not the DI set) to randoms.
#   sig_gene_col : column name for genes in sig_input (varies per file).
NETWORK_METHODS <- list(
  list(prefix = "CoreNet",
       name   = "CoreNet-DI",
       color  = "olivedrab3",
       #union file written inline by the slurm CoreNet job
       sig_input    = file.path(BASE_DIR, "CoreNet", "CoreNet_union_genes.csv"),
       sig_gene_col = "Genes"),
  
  list(prefix = "CoreNet_GRN",
       name   = "CoreNet-GRN-DI",
       color  = "gold",
       sig_input    = "../../Data/GRN_genes_all_sets_152.csv",
       sig_gene_col = "Genes"),

  list(prefix = "CoreNet_PPI",
       name   = "CoreNet-PPI-DI",
       color  = "steelblue2",
       sig_input    = "../../Data/PPS_1_significant_gene_list.csv",
       sig_gene_col = "Genes")
)
JUST_SCORES_METHODS <- list(
  list(prefix = "PPS_noprop_184",
       name   = "PPS_noprop",
       color  = "#DDA0DD",
       sig_input    = "../../Data/EBSeq_Genes_top184.csv",
       sig_gene_col = "Gene"),
  list(prefix = "logFC",
       name   = "logFC RNA",
       color  = "#F08080",
       sig_input    = "../../Data/EBSeq_RNASeq_logFC_top184_0.01.csv",
       sig_gene_col = "gene")
)

# Reference column display names (kept identical to the original script)
# Display labels shown on the rendered plot (data columns stay short).
# Names not in this map render as-is.
REF_DISPLAY <- c("Vinuesa" = "Vinuesa 2016", "Hart" = "Hart 2022")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Build per-method file paths under <BASE_DIR>/<prefix>/.
# species = "human" or "mouse" -- selects which _random_upset_*.csv to read.
paths_for <- function(prefix, tag, species = "human") {
  base <- file.path(BASE_DIR, prefix)
  list(
    real_genes        = file.path(base, paste0(prefix, "_real_DI_unique_genes_df.csv")),
    random_upset      = file.path(base, paste0(prefix, "_", tag,
                                              "_random_upset_", species, ".csv")),
    random_overlaps   = file.path(base, paste0(prefix, "_", tag, "_random_overlaps_df.csv"))
  )
}

read_genes <- function(csv_path, gene_col = "Genes") {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (!gene_col %in% colnames(df)) {
    # try fallbacks
    for (c in c("gene", "uGenes", "Gene", "Symbol")) {
      if (c %in% colnames(df)) { gene_col <- c; break }
    }
  }
  unique(as.character(df[[gene_col]]))
}

# Parse the 'Genes' column of *_random_overlaps_df.csv into a vector of set sizes.
# Python pandas writes Python sets/frozensets like "{'A', 'B', 'C'}", and lists
# like "['A', 'B', 'C']". We just count the comma-separated tokens.
parse_random_set_sizes <- function(random_overlaps_csv) {
  df <- read.csv(random_overlaps_csv, stringsAsFactors = FALSE)
  if (!"Genes" %in% colnames(df)) {
    stop("No 'Genes' column in ", random_overlaps_csv)
  }
  sizes <- vapply(df$Genes, function(s) {
    if (is.na(s) || nchar(s) < 3) return(0L)
    # strip outer braces/brackets, split on commas, count nonempty pieces
    inner <- gsub("^[\\{\\[]|[\\}\\]]$", "", s, perl = TRUE)
    if (nchar(inner) == 0) return(0L)
    length(strsplit(inner, ",", fixed = TRUE)[[1]])
  }, integer(1))
  as.integer(sizes)
}

# Build the upset-format binary df from a named list of gene vectors
combine_gene_lists <- function(gene_lists, column_names) {
  if (length(gene_lists) != length(column_names)) {
    stop("The number of column names must match the number of gene lists.")
  }
  unique_genes <- unique(unlist(gene_lists))
  df <- data.frame(Genes = unique_genes, stringsAsFactors = FALSE)
  for (i in seq_along(gene_lists)) {
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])
  }
  df
}

# -----------------------------------------------------------------------------
# Statistical tests
# -----------------------------------------------------------------------------

# Fisher's exact, alternative='greater'.
# Used when method and random sets are the same size by construction (pipeline 1).
fisher_pval <- function(df, method_col, random_col, reference_col, random_size) {
  a1 <- sum(df[[method_col]] == 1 & df[[reference_col]] == 1)
  m  <- sum(df[[method_col]] == 1)
  a2 <- sum(df[[random_col]] == 1 & df[[reference_col]] == 1)
  r  <- random_size
  mat <- matrix(c(a1, m - a1, a2, r - a2), nrow = 2, byrow = TRUE)
  fisher.test(mat, alternative = "greater")$p.value
}

# Two-proportion z-test, alternative='greater'. Pipeline 2 historical method.
binom_prop_pval <- function(df, method_col, random_col, reference_col, random_size) {
  a1 <- sum(df[[method_col]] == 1 & df[[reference_col]] == 1)
  m  <- sum(df[[method_col]] == 1)
  a2 <- sum(df[[random_col]] == 1 & df[[reference_col]] == 1)
  r  <- random_size
  prop.test(x = c(a1, a2), n = c(m, r),
            alternative = "greater", correct = FALSE)$p.value
}

# Empirical permutation p-value. For each perm, compute overlap_random_i vs reference;
# p = (1 + #{i : overlap_i >= real_overlap}) / (num_perm + 1).
# Uses the per-perm overlaps DF directly; the reference must be one of the columns
# saved by the python pipeline ('Vinuesa', 'Hart').
empirical_pval <- function(real_overlap, random_overlap_vec) {
  num_perm <- length(random_overlap_vec)
  k <- sum(random_overlap_vec >= real_overlap)
  (1 + k) / (1 + num_perm)
}

# -----------------------------------------------------------------------------
# Load reference lists once
# -----------------------------------------------------------------------------
vinuesa    <- read_genes(VINUESA_CSV)
tfh_review <- read_genes(TFH_REVIEW_CSV)
crotty     <- read_genes(CROTTY_CSV)

cat("Vinuesa:", length(vinuesa), "genes\n")
cat("Hart:", length(tfh_review), "genes\n")
cat("Crotty (mouse):", length(crotty), "genes\n\n")

# =============================================================================
# Per-pipeline analysis
# =============================================================================

run_pipeline <- function(tag, methods_def, plot_label, species = "human") {
  cat("\n=========================================================\n")
  cat(" Pipeline: ", tag, " | variant: ", plot_label,
      " | species: ", species, "\n", sep = "")
  cat("=========================================================\n")

  # Species-dependent references and per-perm CSV column names.
  if (species == "human") {
    refs_named <- list(Vinuesa = vinuesa, Hart = tfh_review)
    ref_to_perm_col <- c("Vinuesa" = "Vinuesa", "Hart" = "Hart")
  } else if (species == "mouse") {
    refs_named <- list(Crotty = crotty)
    ref_to_perm_col <- c("Crotty" = "Crotty")
  } else {
    stop("Unknown species: ", species)
  }
  references <- names(refs_named)

  # Pipelines where real and random sets are the same size by construction
  # -> Fisher's exact applies. Pipeline 2 has random sets of varying size,
  # so we use binomial-prop + empirical permutation p instead.
  use_fisher <- tag %in% c(TAG_DI_FIRST, TAG_SIG_ONLY)

  # ---- 1. Load real and random gene lists per method ----
  real_lists   <- list()
  random_lists <- list()
  random_sizes <- list()      # avg size from per-perm overlaps df
  perm_overlaps <- list()     # per-perm raw overlap counts (for empirical p)
  real_overlaps <- list()     # real DI-set overlap with each reference
  display_names <- character()
  random_names  <- character()
  colors        <- character()

  for (m in methods_def) {
    p <- paths_for(m$prefix, tag, species = species)

    # For sig_only, the "method" is the original sig_genes (no DI step).
    # For the other two pipelines, the "method" is the DI set of sig_genes.
    if (tag == TAG_SIG_ONLY) {
      real_path     <- m$sig_input
      real_gene_col <- if (is.null(m$sig_gene_col)) "Genes" else m$sig_gene_col
    } else {
      real_path     <- p$real_genes
      real_gene_col <- "Genes"
    }

    if (!file.exists(real_path)) {
      stop("Missing: ", real_path,
           "\n  -> ", if (tag == TAG_SIG_ONLY)
                       "this is the original sig_genes input file (sig_input)"
                     else
                       "did you run with --save_di / TRUE_DI=1 ?")
    }
    if (!file.exists(p$random_upset)) {
      stop("Missing: ", p$random_upset)
    }
    if (!file.exists(p$random_overlaps)) {
      stop("Missing: ", p$random_overlaps)
    }

    real_genes <- read_genes(real_path, gene_col = real_gene_col)
    rand_genes <- read_genes(p$random_upset)

    real_lists[[m$name]]   <- real_genes
    random_lists[[paste0(m$name, " Random")]] <- rand_genes

    sizes <- parse_random_set_sizes(p$random_overlaps)
    random_sizes[[paste0(m$name, " Random")]] <- round(mean(sizes))
    cat(sprintf("  %-20s real=%d, random_avg=%d (n_perm=%d)\n",
                m$name, length(real_genes), round(mean(sizes)), length(sizes)))

    # Per-perm overlap counts for this species' references (empirical p input).
    ovr_df <- read.csv(p$random_overlaps, stringsAsFactors = FALSE)
    perm_overlaps[[m$name]] <- setNames(
      lapply(references, function(r) ovr_df[[ref_to_perm_col[[r]]]]),
      references
    )

    # Real overlap counts: real_genes intersected with each species reference.
    real_overlaps[[m$name]] <- setNames(
      lapply(references, function(r) length(intersect(real_genes, refs_named[[r]]))),
      references
    )

    display_names <- c(display_names, m$name)
    random_names  <- c(random_names,  paste0(m$name, " Random"))
    colors        <- c(colors,        m$color)
  }

  # ---- 2. Build upset binary df ----
  # Order chosen so randoms come before reals (matches original layout).
  gene_lists <- c(random_lists, refs_named, real_lists)
  colname    <- c(random_names, references,  display_names)
  df <- combine_gene_lists(gene_lists, colname)

  # ---- 3. Compute p-values ----
  # Always store all three test types; columns not applicable for this tag stay NA.
  pval_rows <- list()
  for (i in seq_along(display_names)) {
    method_col <- display_names[i]
    random_col <- random_names[i]
    rsize      <- random_sizes[[random_col]]

    for (ref in references) {
      perm_col <- ref_to_perm_col[[ref]]
      real_o   <- real_overlaps[[method_col]][[ref]]
      perm_v   <- perm_overlaps[[method_col]][[ref]]

      row <- list(
        Method      = method_col,
        Reference   = ref,
        Real_overlap     = real_o,
        Real_set_size    = sum(df[[method_col]] == 1),
        Random_avg_size  = rsize,
        Fisher_p         = NA_real_,
        Binomial_p       = NA_real_,
        Empirical_p      = NA_real_
      )

      if (use_fisher) {
        row$Fisher_p <- fisher_pval(df, method_col, random_col, ref, rsize)
      } else if (tag == TAG_MM_FIRST) {
        row$Binomial_p  <- binom_prop_pval(df, method_col, random_col, ref, rsize)
        row$Empirical_p <- empirical_pval(real_o, perm_v)
      }
      pval_rows[[length(pval_rows) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
    }
  }
  pvals_df <- do.call(rbind, pval_rows)

  pval_path <- file.path(OUTPUT_DIR,
                         sprintf("pvalues_%s_%s_%s.csv", plot_label, species, tag))
  write.csv(pvals_df, pval_path, row.names = FALSE, quote = FALSE)
  cat("\nWrote p-values: ", pval_path, "\n", sep = "")
  print(pvals_df)

  # ---- 4. Build the upset plot ----
  # Intersections: each (real or random method) x each reference, in plot order.
  # Plot top-to-bottom matches the order references are listed (reversed here so
  # the first-listed reference is closest to the top of the rendered figure).
  intersections <- list()
  for (ref in rev(references)) {
    for (i in seq_along(display_names)) {
      intersections[[length(intersections) + 1]] <- c(display_names[i], ref)
      intersections[[length(intersections) + 1]] <- c(random_names[i], ref)
    }
  }

  # Color queries: each real method gets its color; randoms get grey80.
  # Skip queries that target empty intersections -- ComplexUpset chokes on those
  # with a 'Aesthetics must be either length 1 or the same as the data' error.
  has_intersection <- function(cols) {
    sum(rowSums(df[, cols, drop = FALSE]) == length(cols)) > 0
  }

  queries <- list()
  skipped_empty <- character()
  for (i in seq_along(display_names)) {
    for (ref in references) {
      pair_real <- c(display_names[i], ref)
      pair_rand <- c(random_names[i],  ref)
      if (has_intersection(pair_real)) {
        queries[[length(queries) + 1]] <- upset_query(
          intersect = pair_real, color = colors[i], fill = colors[i]
        )
      } else {
        skipped_empty <- c(skipped_empty, paste(pair_real, collapse = " & "))
      }
      if (has_intersection(pair_rand)) {
        queries[[length(queries) + 1]] <- upset_query(
          intersect = pair_rand, color = "grey80", fill = "grey80"
        )
      } else {
        skipped_empty <- c(skipped_empty, paste(pair_rand, collapse = " & "))
      }
    }
  }
  if (length(skipped_empty)) {
    cat("  (skipped color queries for empty intersections: ",
        paste(skipped_empty, collapse = "; "), ")\n", sep = "")
  }

  # Labeller: rename reference columns on the rendered plot only.
  # Data columns stay short ('Vinuesa', 'Hart') for clean p-value CSVs.
  set_labeller <- function(x) {
    out <- ifelse(x %in% names(REF_DISPLAY), REF_DISPLAY[x], x)
    # ifelse with named-vector lookup can return NA when x is not in names();
    # the line above returns the matched display when in names(), else x itself.
    unname(out)
  }

  p <- upset(
    df, colname,
    mode = "inclusive_intersection",
    min_size = 0,
    intersections = intersections,
    sort_intersections = FALSE,
    labeller = set_labeller,
    base_annotations = list(
      "Intersection size" = intersection_size(mode = "inclusive_intersection")
    ),
    themes = upset_modify_themes(list(
      "intersections_matrix" = theme(text = element_text(size = 20)),
      "overall_sizes"        = theme(text = element_text(size = 20)),
      "set_sizes"            = theme(text = element_text(size = 20)),
      "set_names"            = theme(text = element_text(size = 20)),
      "intersection_size"    = theme(
        axis.title.y = element_text(size = 20),
        axis.text.y  = element_text(size = 18)
      ),
      "default" = theme(text = element_text(size = 20))
    )),
    set_size = FALSE,
    queries = queries
  )
  

  # Plot dimensions: scale width with number of intersections
  w_in <- max(5.5, 1.0 * length(intersections) + 1)
  pdf_path <- file.path(OUTPUT_DIR,
                        sprintf("upset_%s_%s_%s.pdf", plot_label, species, tag))
  png_path <- file.path(OUTPUT_DIR,
                        sprintf("upset_%s_%s_%s.png", plot_label, species, tag))
  pdf(pdf_path, width = w_in, height = 10); print(p); dev.off()
  png(png_path, width = w_in, height = 10, units = "in", res = 300); print(p); dev.off()
  cat("Wrote: ", pdf_path, "\n", sep = "")
  cat("Wrote: ", png_path, "\n", sep = "")

  invisible(pvals_df)
}

# =============================================================================
# Run all combinations
# =============================================================================

cat("\n=========================================================\n")
cat(" Generating all combinations: 3 pipelines x 2 variants x 2 species = 12\n")
cat("=========================================================\n")

# Network plot variants
for (sp in c("human", "mouse")) {
  run_pipeline(TAG_DI_FIRST, NETWORK_METHODS, "network", species = sp)
  run_pipeline(TAG_MM_FIRST, NETWORK_METHODS, "network", species = sp)
  run_pipeline(TAG_SIG_ONLY, NETWORK_METHODS, "network", species = sp)
}

# Just-scores plot variants
for (sp in c("human", "mouse")) {
  run_pipeline(TAG_DI_FIRST, JUST_SCORES_METHODS, "just_scores", species = sp)
  run_pipeline(TAG_MM_FIRST, JUST_SCORES_METHODS, "just_scores", species = sp)
  run_pipeline(TAG_SIG_ONLY, JUST_SCORES_METHODS, "just_scores", species = sp)
}

cat("\n=========================================================\n")
cat(" All outputs in: ", OUTPUT_DIR, "\n", sep = "")
cat("=========================================================\n")
