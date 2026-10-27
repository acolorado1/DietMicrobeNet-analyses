# DietMicrobeNet-analyses

Analyses of (DietMicrobe[Host]Net (DM[H]N))[https://github.com/lozuponelab/DietMicrobeNet] networks in two case-study cohorts (Zimbabwe and Diet Modification), and a benchmark of DM[H]N runtime and memory across the Zimbabwe, Diet Modification and IBD datasets.

## Repository layout

```
DietMicrobeNet-analyses/
├── data/
│   ├── ZIM/                              Zimbabwe cohort
│   │   ├── Zim_KOs.tsv                   HUMAnN KO x taxonomy table
│   │   ├── FoodEquivalents.xlsx          FFQ items -> FooDB foods
│   │   ├── zim2_baseline.txt             baseline metadata and FFQ
│   │   ├── zim2_longitudinal.txt         immune markers by visit
│   │   └── participants/ZIM###/          DM[H]N inputs per participant
│   ├── DM/                               Diet Modification cohort
│   │   ├── DM_KOs_shotgun.tsv            HUMAnN KO x taxonomy table (shotgun)
│   │   ├── DM_KOs_rnaseq.txt             HUMAnN KO table (RNA-seq)
│   │   ├── host_gene_counts_tpm.csv      host gene expression
│   │   ├── metadata.xlsx                 diet arm and clinical markers by week
│   │   ├── agrarian_foodb.xlsx           FooDB foods, agrarian diet
│   │   ├── western_foodb.xlsx            FooDB foods, western diet
│   │   └── participants/DM###.1/         DM[H]N inputs per participant
│   ├── IBD/                              IBDMDB subset (benchmarking only)
│   │   ├── metadata.csv                  HMP2 metadata
│   │   ├── sub_metaG.tsv                 EC x taxonomy, metagenomes
│   │   ├── sub_metaT.tsv                 EC x taxonomy, metatranscriptomes
│   │   ├── host_tx_counts.tsv            host transcriptome counts
│   │   └── participants/<ID>_week<N>/    hg_mgx, hg_mtx, htx_mgx, htx_mtx inputs
│   ├── host_genome/                      human-genome KO list shared by all studies
│   └── benchmarking/                     benchmark results
└── scripts/
    ├── R/
    │   └── functions.R                   shared functions
    ├── preprocessing/
    │   ├── 00_prepare_inputs.R           builds the participant input folders
    │   ├── 01_run_dmhn_benchmark.sh      runs DM[H]N on IBD with runtime/memory tracking
    │   ├── run_workflow_benchmark.py     benchmarking wrapper around DM[H]N
    │   └── 02_collapse_outputs.R         collapses DM[H]N outputs, builds the benchmark table
    └── analysis/
        ├── CaseStudy_Zimbabwe.Rmd
        ├── CaseStudy_DM.Rmd
        └── Benchmarking.Rmd
```

Every script and Rmd resolves paths from the repository root with `here::here()`, so run them from inside the repository (or open it as an RStudio project). Figures are written to `figures/<analysis>/`.

## Workflow

```
00_prepare_inputs.R  ->  DM[H]N  ->  02_collapse_outputs.R  ->  analysis/*.Rmd
```

### 1. Build DM[H]N inputs: `scripts/preprocessing/00_prepare_inputs.R`

Writes one folder per participant under `data/<study>/participants/`, **overwriting** the files there. Switch studies on or off with `run_zim`, `run_dm` and `run_ibd` at the top.

| Study | Reads | Writes per participant |
|---|---|---|
| All | `data/host_genome/noquote_ko_host.txt` | `data/host_genome/host_ko_abundance.csv` (placeholder abundance of 1 per KO) |
| ZIM | `Zim_KOs.tsv`, `zim2_baseline.txt` (FFQ, visit 2), `FoodEquivalents.xlsx` | `ko_taxonomy_abundance.csv`, `noquote_ko.txt`, host files, `FFQ.csv`, `foodb_foods_dataframe.csv`* |
| DM | `DM_KOs_shotgun.tsv`, `metadata.xlsx`, `agrarian_foodb.xlsx`, `western_foodb.xlsx` | `ko_taxonomy_abundance.csv`, `noquote_ko.txt`, host files, `foodb_foods_dataframe.csv` for the assigned diet |
| IBD | `metadata.csv`, `sub_metaG.tsv`, `sub_metaT.tsv`, `host_tx_counts.tsv` | `hg_mgx`, `hg_mtx`, `htx_mgx`, `htx_mtx` folders, each with microbe and host KO files |

\* FFQ food frequencies are written into `foodb_foods_dataframe.csv` only when the participant folder contains the unweighted `foodb_foods_dataframe_og.csv`; otherwise the existing file is kept.

IBD participants are each participant's earliest timepoint with metagenome, metatranscriptome and host transcriptome data, excluding P6028. EC numbers and host gene symbols are mapped to KEGG KOs, with the mapping tables cached in `data/IBD/cache/`.

### 2. Run DM[H]N

DM[H]N is run from the DietMicrobeNet repository on each participant folder.

`scripts/preprocessing/01_run_dmhn_benchmark.sh` runs every IBD pair folder through `run_workflow_benchmark.py`, which records runtime and peak memory in `data/benchmarking/benchmark_IBD.csv`. Set `DMN_DIR` at the top to your DietMicrobeNet repository and run it from the DietMicrobeNet run-workflow directory.

### 3. Collapse outputs: `scripts/preprocessing/02_collapse_outputs.R`

DM[H]N writes many large files per participant; the analyses only need `graph/graph_results.csv`, plus `graph/network_summary.txt` for benchmarking. This script:

1. Stacks every `graph_results.csv` into one file per study, adding `run_dir` (participant folder) and `output` (DM[H]N output folder) columns:
   - `data/ZIM/graph_results_ZIM.csv.gz`
   - `data/DM/graph_results_DM.csv.gz`
   - `data/DM/graph_results_DM_ffq.csv.gz`
   - `data/IBD/graph_results_IBD.csv.gz` (when `include_ibd <- TRUE`)
2. Adds network size, node-origin counts, Shannon diversity, number of taxa and KO counts to the benchmark runtime/memory tables, written to `data/benchmarking/benchmark_all_results_with_network_stats.csv`.

Set `zim_dir`, `dm_dir` and `bench_dir` at the top to where the ZIM and DM runs were written. IBD is read from `data/IBD/`.

**Run this before deleting any DM[H]N output folders.** It warns if any `network_summary.txt` files are missing. `output_*` folders are git-ignored.

### 4. Analyses: `scripts/analysis/`

| Rmd | Reads | Figures |
|---|---|---|
| `CaseStudy_Zimbabwe.Rmd` | `data/ZIM/graph_results_ZIM.csv.gz`, `zim2_baseline.txt`, `zim2_longitudinal.txt`, `Zim_KOs.tsv` | `figures/ZIM/` |
| `CaseStudy_DM.Rmd` | `data/DM/graph_results_DM.csv.gz`, `DM_KOs_shotgun.tsv`, `metadata.xlsx`, `data/host_genome/noquote_ko_host.txt` | `figures/DM/` |
| `Benchmarking.Rmd` | `data/benchmarking/benchmark_all_results_with_network_stats.csv` | `figures/benchmarking/` |

The two case-study Rmds also need the HMDB chemical-taxonomy table from the DietMicrobeNet repository (`hmdb_path` in their setup chunk).

## Setup

1. Add the human-genome KO list (one KO per line) as `data/host_genome/noquote_ko_host.txt`. Every participant's `noquote_ko_host.txt` is a copy of this file:
   ```bash
   mkdir -p data/host_genome
   cp data/DM/participants/DM001.1/noquote_ko_host.txt data/host_genome/
   ```
2. Set `hmdb_path` in `CaseStudy_Zimbabwe.Rmd` and `CaseStudy_DM.Rmd`.
3. Set `DMN_DIR` in `01_run_dmhn_benchmark.sh` and the run folders in `02_collapse_outputs.R`.
4. Install the R packages:

```r
install.packages(c(
  "tidyverse", "readxl", "here", "vegan", "cowplot", "ggpubr", "Polychrome",
  "scales", "flextable", "broom", "knitr", "gt", "lme4", "lmerTest",
  "performance", "mediation", "sessioninfo"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("KEGGREST", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db"))
```

KEGG compound names are downloaded once and cached in `data/cache/kegg_compound_names.rds`. KEGG enrichment and the IBD ID mapping need an internet connection on first run.

## Shared functions: `scripts/R/functions.R`

| Section | Functions |
|---|---|
| Input preparation | `write_participant_ko_dirs`, `write_host_ko_abundance`, `ffq_to_foodb`, `add_ffq_frequencies`, `build_gene_to_ko_map`, `build_ec_to_ko_map`, `find_sample_col`, `participant_weeks`, `get_ext_id`, `write_microbe_files`, `write_host_tx_files`, `build_pair_dir` |
| Collapsing outputs | `collapse_graph_results`, `read_collapsed_results`, `parse_network_summary`, `count_ko_lines`, `compute_shannon`, `benchmark_paths`, `add_benchmark_stats` |
| Parsing list columns | `parse_py_list`, `parse_py_list_orgs`, `parse_py_list_safe` |
| Reference lookups | `read_hmdb`, `kegg_compound_names`, `kegg_name_lookup` |
| Enrichment | `run_ko_enrichment`, `taxonomy_class_enrichment`, `jaccard` |
| Reaction composition | `build_reaction_matrix`, `reaction_composition_test`, `plot_reaction_pcoa`, `simper_robustness` |
| Marker models | `label_zim_experiment`, `fit_lm_by_marker`, `cor_by_level`, `fit_marker_models` |
| Benchmarking | `fit_diversity_pair`, `source_boxplot` |
