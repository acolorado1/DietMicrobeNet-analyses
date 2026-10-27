# =============================================================================
# Prepare DM[H]N inputs for each case study (run from the repository root)
#   0. Shared host (human genome) KO files
#   1. ZIM  - per-participant KO + host files, FFQ -> FooDB food frequencies
#   2. DM   - per-participant KO, host and FooDB files
#   3. IBD  - participant selection and microbe/host pair directories
#
# Writes into data/<study>/participants/, overwriting the files there.
# =============================================================================

source(here::here("scripts", "R", "functions.R"))
suppressPackageStartupMessages(library(readxl))
setwd(here::here())
set.seed(5)

run_zim <- TRUE
run_dm  <- TRUE
run_ibd <- TRUE

# human-genome KO list shared by all studies; abundance file is generated
host_ko_file    <- "data/host_genome/noquote_ko_host.txt"
host_abund_file <- "data/host_genome/host_ko_abundance.csv"


# =============================================================================
# 0. Host (human genome) KO files
# =============================================================================
if (!file.exists(host_ko_file)) {
  stop(host_ko_file, " not found. Copy the human-genome KO list there, e.g. ",
       "data/DM/participants/DM001.1/noquote_ko_host.txt")
}
write_host_ko_abundance(host_ko_file, host_abund_file)
host_files <- c(host_ko_file, host_abund_file)


# =============================================================================
# 1. ZIM
# =============================================================================
if (run_zim) {
  zim_dir <- "data/ZIM"

  # ---- per-participant KO and host files -----------------------------------
  zim_ko_table <- read_delim(file.path(zim_dir, "Zim_KOs.tsv"), delim = "\t",
                             escape_double = FALSE, trim_ws = TRUE)
  zim_ids <- str_sub(names(zim_ko_table)[-1], 1, 6)    # "ZIM089.1" -> "ZIM089"

  zim_dirs <- write_participant_ko_dirs(zim_ko_table, file.path(zim_dir, "participants"),
                                        dir_names = zim_ids, quote = TRUE)
  for (d in zim_dirs$dir) file.copy(host_files, d, overwrite = TRUE)

  # ---- FFQ -> FooDB food frequencies ---------------------------------------
  food_eq <- read_excel(file.path(zim_dir, "FoodEquivalents.xlsx")) %>%
    select(food = `Zimbabwe Food Item`, FooDB)

  # FFQ column names as recorded vs. names split into FooDB-matchable foods
  ffq_cols_raw <- c(
    "w_Butter_Margarine", "a_WholeMilk_Hodzeko_Amasi_Lacto_GoatMilk_SkimmedMilk",
    "w_Cheese_Yogurt_IceCream", "a_Grapes_Banana_Apple_Peaches_SugarCane",
    "a_Mango_Naartjie_Orange_Pineapple_PawPaw_Guava_Nyii_Maonde",
    "a_Mazhanje_Lemon_Boabab_Tsubvu_Mathowe_Masau", "a_Avocado_Plum",
    "w_Bakedbeans_Cannedmushrooms_GreenVegetables_frozen",
    "w_Fries_Friedpotatoes_Roastedpotatoes", "a_Sugarbeans_Nyemba_Peas",
    "a_Muboora_Nyevhe_Rape_Rugare_Kovo_Mowa_TsungaCabbage_Carrots_Celery_Cucumber",
    "a_SweetPotato_Cassava_Pumpkin_Yam_Gourd_Turnip_Butternut", "w_Beef_Pork",
    "w_Fried_Chicken", "w_Pilchards_Sardines", "w_Hamburger", "w_Pizza",
    "w_Peri_Peri_Chicken_Nandos", "a_Chicken_Birds_Rabbit_Mbeva_without_skin", "w_Eggs",
    "a_Matemba_Fish_Cannedtuna", "a_Caterpillar_Locust_Termites", "w_Bread_Dark_White",
    "w_WhiteRice", "w_Samp", "a_Sadza_Rapoko_Zviyo_Mhunga_Mapfunde_Mugaiwa",
    "a_Sadza_Pearlenta_RollerMeal", "w_Macaroni_Pasta_no_sauce",
    "w_Cascade_SoftDrinks_Juice", "a_Maheu_TraditionalBeer", "w_PeanutButter",
    "w_Sauces_Dressing", "w_Cake_Sweets", "w_Biscuits"
  )

  ffq_cols <- c(
    "w_Butter_Margarine", "a_WholeMilk_HodzekoAmasi_Lacto_GoatMilk_SkimmedMilk",
    "w_Cheese_Yogurt_IceCream", "a_Grapes_Banana_Apple_Peaches_SugarCane",
    "a_Mango_Naartjie_Orange_Pineapple_PawPaw_Guava_Nyii_Maonde",
    "a_Mazhanje_Lemon_Boabab_Tsubvu_Mathowe_Masau", "a_Avocado_Plum",
    "w_Bakedbeans_Cannedmushrooms_GreenVegetables",
    "w_Fries_Friedpotatoes_Roastedpotatoes", "a_Sugarbeans_Nyemba_Peas",
    "a_Muboora_Nyevhe_Rape_Rugare_Kovo_Mowa_Tsunga_Cabbage_Carrots_Celery_Cucumber",
    "a_SweetPotato_Cassava_Pumpkin_Yam_Gourd_Turnip_Butternut", "w_Beef_Pork",
    "w_FriedChicken_Chicken", "w_Pilchards_Sardines", "w_Hamburger", "w_Pizza",
    "w_PeriPeriChicken_Nandos", "a_Chicken_Birds_Rabbit_Mbeva_without_skin", "w_Eggs",
    "a_MatembaFish_Cannedtuna", "a_Caterpillar_Locust_Termites", "w_Darkbread_Whitebread",
    "w_WhiteRice", "w_Samp", "a_SadzaRapoko_ZviyoMhunga_MapfundeMugaiwa",
    "a_Sadza_Pearlenta_RollerMeal", "w_MacaroniPastanosauce",
    "w_Cascade_SoftDrinks_Juice", "a_Maheu_TraditionalBeer", "w_PeanutButter",
    "w_Sauces_Dressing", "w_Cake_Sweets", "w_Biscuits"
  )

  zim_ffq <- read_delim(file.path(zim_dir, "zim2_baseline.txt"), delim = "\t",
                        escape_double = FALSE, trim_ws = TRUE) %>%
    rename(!!!setNames(ffq_cols_raw, ffq_cols)) %>%
    filter(PID %in% zim_ids, Visit == 2)

  # FFQ.csv per participant; food frequencies are written into
  # foodb_foods_dataframe.csv when the unweighted foodb_foods_dataframe_og.csv
  # from DM[H]N is present
  for (pid in zim_ids) {
    pid_dir <- file.path(zim_dir, "participants", pid)

    ffq <- ffq_to_foodb(zim_ffq, pid = pid, food_eq = food_eq, ffq_cols = ffq_cols)
    write.csv(ffq, file = file.path(pid_dir, "FFQ.csv"), row.names = FALSE)

    # this does not need to run again since the OG df is redundant
    og_file <- file.path(pid_dir, "foodb_foods_dataframe_og.csv")
    if (file.exists(og_file)) {
      foodb_df <- add_ffq_frequencies(ffq, read_csv(og_file, show_col_types = FALSE), pid = pid)
      write.csv(foodb_df, file = file.path(pid_dir, "foodb_foods_dataframe.csv"), row.names = FALSE)
    } else {
      message(pid, ": no foodb_foods_dataframe_og.csv, kept existing foodb_foods_dataframe.csv")
    }
  }
}


# =============================================================================
# 2. DM
# =============================================================================
if (run_dm) {
  dm_dir <- "data/DM"

  dm_ko_table <- read_tsv(file.path(dm_dir, "DM_KOs_shotgun.tsv"), show_col_types = FALSE)
  dm_metadata <- read_excel(file.path(dm_dir, "metadata.xlsx"))

  dm_dirs <- write_participant_ko_dirs(dm_ko_table, file.path(dm_dir, "participants"))

  # FooDB foods for each participant's assigned diet, plus host files
  foodb_by_diet <- list(
    Agrarian = read_excel(file.path(dm_dir, "agrarian_foodb.xlsx")),
    Western  = read_excel(file.path(dm_dir, "western_foodb.xlsx"))
  )

  for (i in seq_len(nrow(dm_dirs))) {
    diet_type <- dm_metadata %>%
      filter(SampleID == dm_dirs$sample_col[i]) %>%
      pull(Diet) %>%
      dplyr::first()
    foodb <- if (diet_type == "Agrarian") foodb_by_diet$Agrarian else foodb_by_diet$Western
    write_csv(foodb, file.path(dm_dirs$dir[i], "foodb_foods_dataframe.csv"))
    file.copy(host_files, dm_dirs$dir[i], overwrite = TRUE)
  }
}


# =============================================================================
# 3. IBD
# =============================================================================
if (run_ibd) {
  ibd_dir <- "data/IBD"

  ibd_meta <- read_csv(file.path(ibd_dir, "metadata.csv"), show_col_types = FALSE)
  ibd_meta$diagnosis <- factor(ibd_meta$diagnosis, levels = c("UC", "CD", "nonIBD"))

  excluded_participants <- c("P6028")   # no microbiome
  # participant-weeks with host transcriptome + metagenome + metatranscriptome
  triple_match <- participant_weeks(ibd_meta, "host_transcriptomics") %>%
    inner_join(participant_weeks(ibd_meta, "metagenomics"),        by = c("Participant ID", "week_num")) %>%
    inner_join(participant_weeks(ibd_meta, "metatranscriptomics"), by = c("Participant ID", "week_num"))

  # earliest triple-matched timepoint per participant
  chosen <- triple_match %>%
    filter(!`Participant ID` %in% excluded_participants) %>%
    group_by(`Participant ID`) %>%
    slice_min(week_num, n = 1, with_ties = FALSE) %>%
    ungroup()

  cat("Participants:", n_distinct(chosen$`Participant ID`), "\n")

  # input tables and ID mappings
  mgx_df <- read_tsv(file.path(ibd_dir, "sub_metaG.tsv"), show_col_types = FALSE)
  mtx_df <- read_tsv(file.path(ibd_dir, "sub_metaT.tsv"), show_col_types = FALSE)
  htx_df <- read_tsv(file.path(ibd_dir, "host_tx_counts.tsv"), show_col_types = FALSE)

  ec_to_ko_map <- build_ec_to_ko_map(file.path(ibd_dir, "cache/ec_to_ko_map.rds"))
  gene_ko_map  <- build_gene_to_ko_map(htx_df[[1]],
                                       file.path(ibd_dir, "cache/gene_to_ko_map.rds"))

  host_genome_files <- c(abundance = host_abund_file, ko = host_ko_file)

  for (i in seq_len(nrow(chosen))) {
    pid <- chosen$`Participant ID`[i]
    wk  <- chosen$week_num[i]
    base_dir <- file.path(ibd_dir, "participants", paste0(pid, "_week", wk))

    mgx_id <- get_ext_id(ibd_meta, pid, wk, "metagenomics")
    mtx_id <- get_ext_id(ibd_meta, pid, wk, "metatranscriptomics")
    htx_id <- get_ext_id(ibd_meta, pid, wk, "host_transcriptomics")

    build_pair_dir(file.path(base_dir, "htx_mtx"), mtx_df, mtx_id, ec_to_ko_map,
                   host_df = htx_df, host_id = htx_id, gene_ko_map = gene_ko_map)
    build_pair_dir(file.path(base_dir, "htx_mgx"), mgx_df, mgx_id, ec_to_ko_map,
                   host_df = htx_df, host_id = htx_id, gene_ko_map = gene_ko_map)
    build_pair_dir(file.path(base_dir, "hg_mtx"),  mtx_df, mtx_id, ec_to_ko_map,
                   host_genome_files = host_genome_files)
    build_pair_dir(file.path(base_dir, "hg_mgx"),  mgx_df, mgx_id, ec_to_ko_map,
                   host_genome_files = host_genome_files)
  }
}
