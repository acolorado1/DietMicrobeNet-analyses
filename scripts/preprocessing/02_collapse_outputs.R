# =============================================================================
# Harvest everything the analyses need from the DM[H]N output folders, so the
# per-participant outputs can be deleted afterwards. Run from the repository
# root; reads the folders DM[H]N wrote to and saves the results under data/.
#   1. Collapse per-participant graph_results.csv into one file per study
#   2. Build the benchmarking table (needs network_summary.txt from the outputs)
# =============================================================================

source(here::here("scripts", "R", "functions.R"))
setwd(here::here())

# ---- where DM[H]N was run: edit for your machine ---------------------------
zim_dir   <- "~/Desktop/ZIM_dmn/participants"
dm_dir    <- "~/Desktop/DM_dmn"
bench_dir <- "~/Desktop/benchmarking"     # holds the ZIM/DM benchmark_*.csv files

# IBD is run inside the repo by 01_run_dmhn_benchmark.sh
ibd_dir       <- "data/IBD"
ibd_bench_csv <- "data/benchmarking/benchmark_IBD.csv"

include_ibd <- FALSE   # IBD collapse + benchmark rows


# =============================================================================
# 1. Collapse graph_results.csv
# =============================================================================
collapse_graph_results(zim_dir,
                       "data/ZIM/graph_results_ZIM.csv.gz",
                       outputs = c("output_dm", "output_dm_ALL", "output_dmh_ALL"))

collapse_graph_results(file.path(dm_dir, "participants"),
                       "data/DM/graph_results_DM.csv.gz",
                       outputs = "output_dmh")

collapse_graph_results(file.path(dm_dir, "dm_ffq"),
                       "data/DM/graph_results_DM_ffq.csv.gz",
                       outputs = "output_dmh_ffq")

if (include_ibd) {
  collapse_graph_results(file.path(ibd_dir, "participants"),
                         "data/IBD/graph_results_IBD.csv.gz",
                         outputs = "output_dmh")
}


# =============================================================================
# 2. Benchmarking table with network stats
# =============================================================================
hpc_dm  <- "/projects/aburkhartcolorado@xsede.org/DietMicrobeNet/dm"
hpc_zim <- "/projects/aburkhartcolorado@xsede.org/DietMicrobeNet/zim"
local_dm  <- path.expand(file.path(dm_dir, "participants"))
local_zim <- path.expand(zim_dir)

benchmark_IBD <- if (include_ibd) read_csv(ibd_bench_csv, show_col_types = FALSE) %>%
  benchmark_paths(output_dir = "output_dmh") %>%
  mutate(experiment_type     = basename(directory_clean),
         participant_id      = basename(dirname(directory_clean)),
         source              = "IBD",
         food_representation = "All Foods",
         analysis_type       = "DietMicrobeHost",
         host                = "yes")

benchmark_DMall <- read_csv(file.path(bench_dir, "benchmark_DMall_ZIMdm.csv"), show_col_types = FALSE) %>%
  filter(source == "DM") %>%
  benchmark_paths(output_dir = "output_dmh", old_prefix = hpc_dm, new_prefix = local_dm) %>%
  mutate(experiment_type     = NA_character_,
         food_representation = "All Foods",
         analysis_type       = "DietMicrobeHost",
         host                = "yes")

benchmark_DMffq <- read_csv(file.path(bench_dir, "benchmark_DMffq.csv"), show_col_types = FALSE) %>%
  benchmark_paths(output_dir = "output_dmh_ffq", old_prefix = hpc_dm, new_prefix = local_dm) %>%
  mutate(network_summary_path = file.path(dirname(local_dm), "dm_ffq", participant_id,
                                          "output_dmh_ffq", "graph", "network_summary.txt"),
         experiment_type      = NA_character_,
         food_representation  = "FFQ",
         analysis_type        = "DietMicrobeHost",
         host                 = "yes")

benchmark_ZIMdm <- read_csv(file.path(bench_dir, "benchmark_DMall_ZIMdm.csv"), show_col_types = FALSE) %>%
  filter(source == "ZIM") %>%
  benchmark_paths(output_dir = "output_dm", old_prefix = hpc_zim, new_prefix = local_zim,
                  has_host = FALSE) %>%
  mutate(experiment_type     = NA_character_,
         food_representation = "FFQ",
         analysis_type       = "DietMicrobe",
         host                = "no")

benchmark_ZIMdmh <- read_csv(file.path(bench_dir, "benchmark_ZIMdmh.csv"), show_col_types = FALSE) %>%
  benchmark_paths(output_dir = "output_dmh_ALL", old_prefix = hpc_zim, new_prefix = local_zim) %>%
  mutate(experiment_type     = NA_character_,
         food_representation = "All Foods",
         analysis_type       = "DietMicrobeHost",
         host                = "yes")

benchmark <- bind_rows(
  if (include_ibd) benchmark_IBD,
  benchmark_DMall, benchmark_DMffq, benchmark_ZIMdm, benchmark_ZIMdmh
)

missing_summaries <- sum(!file.exists(benchmark$network_summary_path))
if (missing_summaries > 0) {
  warning(missing_summaries, " network_summary.txt files not found; their stats will be NA")
}

benchmark_final <- add_benchmark_stats(benchmark)
dir.create("data/benchmarking", showWarnings = FALSE)
write_csv(benchmark_final, "data/benchmarking/benchmark_all_results_with_network_stats.csv")
