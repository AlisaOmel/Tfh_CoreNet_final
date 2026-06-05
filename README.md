# TFH-CoreNet

**Core Regulatory Network of T Follicular Helper Cell Differentiation**

### Table of Contents

- [Description](#description)
- [Repository Structure](#repository-structure)
- [How To Use](#how-to-use)
- [Code Documentation](#code-documentation)
- [Citation](#citation)

---

## Description

CoreNet was generated using a multi-omic computational approach for identifying and validating the core regulatory gene network underlying T follicular helper (TFH) cell differentiation. The pipeline integrates bulk RNA-seq, ATAC-seq, and CUT&RUN (H3K4me1, H3K4me3, H3K27Ac) data across four TFH differentiation states: **Naive → Early Pre-TFH → Late Pre-TFH → GC TFH**.

The dataset has 2 components:

- **CoreNet (PPI)** — Identifies direct protein–protein interactors of a significant gene set using the HINT PPI database and assesses statistical enrichment via permutation testing (HotNet2 and direct-interactor permutations).
- **CoreNet (GRN)** — Constructs a gene regulatory network using Taiji (ATAC-seq and CUT&RUN-derived transcription factor rankings) filtered by significance across TFH differentiation stages.


Validation is performed across human tonsil spatial transcriptomics, mouse scRNA-seq, CITE-seq tonsil atlas data, and MSigDB gene set enrichment analyses (ssGSEA).

---

## Repository Structure

```
tfh_corenet/
├── Data/
│   ├── biomart_corenet.csv                  # Human–mouse gene ortholog mapping for CoreNet genes
│   ├── cm_pre_tfh_vs_gc_SLIDE_genes.csv     # SLIDE significant genes (CM Pre-TFH vs GC TFH)
│   ├── corenet_extended_full_genes.csv       # Full extended CoreNet gene list
│   ├── corenet_extended_random_upset_human_nofd.csv
│   ├── CoreNet_PPI_protein_pairs_sig_genes_df.csv
│   ├── EBSeq_Genes_top184.csv               # Top 184 genes from EBSeq differential expression
│   ├── EBSeq_RNASeq_logFC_top184_0.01.csv
│   ├── GRN_genes_all_sets_152.csv           # CoreNet GRN gene list (152 genes)
│   ├── GRN_genes_all_sets_ALLGENES.csv      # All TFs tested 
│   ├── groom_tfhcore.csv
│   ├── hotnet2_heat_score_genes.csv
│   ├── merged_output.rds
│   ├── PPS_1_significant_gene_list.csv      # CoreNet-PPI propagation significant genes
│   ├── pps_taiji_unique_genes_df.txt
│   ├── scenic_all_genes.csv
│   ├── taiji_score_propagation_genes.csv
│   ├── union_propagation_genes.csv
│   ├── literature_sets/                     # Published gene sets
│   │   ├── crotty_list.csv
│   │   ├── dalit_coreTFH_genes.csv
│   │   ├── dalit_coreTFR_genes.csv
│   │   ├── TFH_review_list.csv
│   │   └── vinuesa_list.csv
│   ├── SLIDE/                                #SLIDE inputs
│   ├── Taiji_outputs/                       # Processed Taiji GeneRank and ranked outputs
│   └── upset_plot_overlaps/                 # Permutation results for CoreNet, GRN, PPI, logFC, PPS
│       ├── CoreNet/
│       ├── CoreNet_GRN/
│       ├── CoreNet_PPI/
│       ├── logFC/
│       └── PPS_noprop/
└── Scripts/
    ├── 1_bulkRNA_processing/
    │   └── process_RNA_and_logFC.ipynb       # Normalize bulk RNA-seq; compute log2FC between TFH states
    ├── 2_CoreNet_PPI_processing/
    │   ├── HotNet2.py                        # HotNet2 wrapper for network propagation
    │   ├── first_degree_interactors_pipeline.py  # Direct-interactor permutation pipeline
    │   ├── example_hotnet2_run.slurm
    │   └── run_di_perm.slurm
    ├── 3_CoreNet_GRN_processing/
    │   ├── GRN_filter_and_get_top_genes.ipynb  # Filter Taiji GeneRanks; select top TF genes
    │   ├── config_atac.yml
    │   └── mod_input_atac.yml
    ├── Figure2/
    │   └── make_upset_plots_and_pvalues.R    # UpSet plots and permutation p-values for CoreNet
    ├── Figure3_and_Supp_Figure_4/
    │   ├── run_cm_tfh_slide.R                # SLIDE analysis on CM Pre-TFH vs GC TFH SEACells
    │   ├── run_cm_tfh_slide_cv.R             # Cross-validation of SLIDE results
    │   └── ssgsea_boxplot_all_final.R        # ssGSEA enrichment boxplots
    ├── Figure4/
    │   └── ssgsea_mouse_scRNA_fig4cd.R       # ssGSEA on mouse scRNA-seq data
    ├── Figure5/
    │   └── human_spatial_scRNA_analysis_5b.R # Human tonsil spatial transcriptomics analysis
    ├── Figure6/
    │   ├── msigdb_random_permutations_fig6a.R
    │   ├── run_msigdb_random_perm.slurm
    │   └── ssGSEA_msigdb_fig6bcef_final.R    # ssGSEA against MSigDB gene sets
    ├── Supp_Figure_2/
    │   ├── module_visualization_extfig2A.R
    │   └── TFR_vs_TFH_extfig2E.ipynb
    ├── Supp_Figure_3/
    │   └── steiner_tree_extfig.py            # Steiner tree construction for network figure
    ├── Supp_Figure_5/
    │   └── SLIDE_CV_and_Cliffs_Delta_analysis.R
    ├── Supp_Figure_7/
    │   └── make_human_upset_plot_sig_corenetextended.R
    ├── Supp_Figure_8/
    │   └── cite_seq_tonsil_atlas_suppfig_bc.R
    └── Supp_Figure_9/
        └── IL12RB1_expr_fig9a.ipynb
```

---


## How To Use

### Dependencies

**Python**
- pandas
- numpy
- scipy
- seaborn
- matplotlib
- h5py
- pickle
- scikit-learn
- hotnet2

**R**
- Seurat / Signac
- SLIDE
- ggplot2
- ComplexUpset / UpSetR
- cliffsDelta
- GSVA (ssGSEA)



## Code Documentation

### `first_degree_interactors_pipeline.py`

#### `sig_set_nodot(gene_set, set_name, output_dir, save, taiji)`

Reads a gene list from a CSV file and strips Ensembl version suffixes (`.` notation).

- `gene_set` — path to CSV containing gene names
- `taiji` — if `'True'`, reads from a `'Genes'` column (Taiji format); otherwise reads column headers

Returns a deduplicated list of gene names.

---

#### `first_degree_interactors(sig_genes, hint_data, output_dir, set_name, save, save_formats)`

Queries the HINT PPI network for all direct interactors of the input significant gene set.

- `sig_genes` — list of gene names
- `hint_data` — parsed HINT HDF5 interaction data
- `save_formats` — output formats: `'csv'`, `'pkl'`, or `'txt'`

When `save='on'`, writes:
- `{set_name}_protein_pairs_df.csv` — all PPI edges involving sig genes
- `{set_name}_unique_genes_df.csv` — unique genes in the DI network

---

#### Permutation pipelines

Three permutation strategies assess whether the observed overlap between CoreNet gene sets and literature/validation sets is greater than expected by chance:

- **`pipeline_di_then_match`** — Expands sig genes to their DI set once, then draws 1000 random sets matched on size to estimate a null overlap distribution.
- **`pipeline_match_then_di`** — Draws 1000 random sets the same size as sig genes, expands each to its DI set, and computes overlaps.
- **`pipeline_sig_only`** — Computes overlap of sig genes directly without DI expansion, with 1000 size-matched random draws as the null.

---

### `GRN_filter_and_get_top_genes.ipynb`

#### `filter_GeneRanks(generank_dir, output_dir, experiment)`

Loads Taiji GeneRank and p-value TSVs for a given experiment, renames cell-state columns to standardized TFH nomenclature, and filters to genes with p < 0.01 across all states.

---

#### `rank_TFs(analysis, experiment_list)`

Rank-transforms filtered GeneRank scores within each TFH cell state. The top-ranked gene receives a score of 10⁻⁴; all others are ranked 1–N and normalized to [0, 1].

---

#### `get_top_genes(experiment_list, top_num, exclude_znf)`

Selects the top `top_num` genes by minimum rank across all experiments and cell states. Optionally excludes ZNF family genes. Returns a dataframe of −log-transformed rank values.

---

### `process_RNA_and_logFC.ipynb`

Reads normalized bulk RNA-seq counts, computes per-state means across replicates (Naive, Early Pre-TFH, Late Pre-TFH, GC TFH), and calculates log2 fold-change (GC vs Early Pre-TFH; GC vs Late Pre-TFH). Values are clipped to ±2 for a thresholded version used in downstream UpSet analyses.

---

## Citation


---

## About

Scripts and data associated with the TFH CoreNet manuscript. For questions, please open an issue or contact the corresponding author.
