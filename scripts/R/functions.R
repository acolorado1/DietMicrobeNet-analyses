# =============================================================================
# Shared functions for the DM[H]N case studies (ZIM, DM, IBD, benchmarking)
#
#   source(here::here("scripts", "R", "functions.R"))
#
# Sections
#   1. Input preparation (before DM[H]N)
#   2. Collapsing DM[H]N outputs
#   3. Parsing DM[H]N list columns
#   4. Reference lookups (HMDB, KEGG)
#   5. Enrichment and set similarity
#   6. Reaction composition (PERMANOVA / PCoA / SIMPER)
#   7. Clinical / immune marker models
#   8. Resource-use models (benchmarking)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
})


# =============================================================================
# 1. Input preparation (before DM[H]N)
# =============================================================================

# ---- ZIM and DM: per-participant KO files from a HUMAnN table ----------------

# Split a HUMAnN "KO|taxonomy" table and write, for each sample column,
#   <out_root>/<dir name>/ko_taxonomy_abundance.csv  (KO, taxonomy, Abundance_RPKs)
#   <out_root>/<dir name>/noquote_ko.txt             (unique KOs)
# UNGROUPED / UNMAPPED rows, rows without taxonomy and zero abundances are
# dropped. Returns tibble(sample_col, dir).
write_participant_ko_dirs <- function(ko_table, out_root,
                                      dir_names = names(ko_table)[-1],
                                      quote = FALSE) {
  sample_cols <- names(ko_table)[-1]
  stopifnot(length(dir_names) == length(sample_cols))

  gf <- ko_table[[1]]
  drop <- str_detect(gf, "UNGROUPED") | gf == "UNMAPPED" | startsWith(gf, "UNMAPPED|")
  ko_table <- ko_table[!drop, ]

  split_gf <- strsplit(ko_table[[1]], "|", fixed = TRUE)
  KO       <- vapply(split_gf, `[`, character(1), 1)
  taxonomy <- vapply(split_gf, function(x) if (length(x) > 1) x[2] else "", character(1))

  has_tax  <- taxonomy != ""
  KO       <- KO[has_tax]
  taxonomy <- taxonomy[has_tax]
  ko_table <- ko_table[has_tax, ]

  dirs <- file.path(out_root, dir_names)
  for (i in seq_along(sample_cols)) {
    abundance <- ko_table[[sample_cols[i]]]
    nonzero   <- abundance != 0

    sample_tbl <- data.frame(
      KO             = KO[nonzero],
      taxonomy       = taxonomy[nonzero],
      Abundance_RPKs = abundance[nonzero],
      stringsAsFactors = FALSE
    )

    dir.create(dirs[i], showWarnings = FALSE, recursive = TRUE)
    write.csv(sample_tbl, file = file.path(dirs[i], "ko_taxonomy_abundance.csv"),
              row.names = FALSE, quote = quote)
    writeLines(unique(sample_tbl$KO), con = file.path(dirs[i], "noquote_ko.txt"))
  }
  invisible(tibble(sample_col = sample_cols, dir = dirs))
}

# Host (human genome) KO list -> abundance table with a placeholder abundance
# of 1, since there are no read counts for host genes
write_host_ko_abundance <- function(ko_file, out_file) {
  read_csv(ko_file, col_names = "ko", show_col_types = FALSE) %>%
    mutate(abundance = 1) %>%
    write.csv(file = out_file, row.names = FALSE)
}

# ---- ZIM: FFQ -> FooDB -------------------------------------------------------

# Convert one participant's FFQ scores into per-food percentages joined to
# FooDB equivalents. Composite FFQ items ("a_Mango_Orange") are split and the
# score is shared equally between the foods. `ffq` should hold one row per
# participant (e.g. a single visit).
ffq_to_foodb <- function(ffq, pid, food_eq, ffq_cols) {
  ffq %>%
    filter(PID == pid) %>%
    select(all_of(ffq_cols)) %>%
    pivot_longer(cols = all_of(ffq_cols), names_to = "food", values_to = "score") %>%
    mutate(row_id = row_number(),
           food   = str_remove(food, "^(w_|a_)")) %>%
    separate_rows(food, sep = "_") %>%
    group_by(row_id) %>%
    mutate(score = score / n()) %>%
    ungroup() %>%
    select(-row_id) %>%
    mutate(raw_perc = score / sum(score) * 100,
           perc     = floor(raw_perc)) %>%
    filter(perc != 0) %>%
    left_join(food_eq, by = "food")
}

# Replace food_frequency in a FooDB foods table with one participant's FFQ
# percentages (summed over FFQ items sharing a FooDB food). Stops if the FFQ
# and FooDB food names do not match.
add_ffq_frequencies <- function(ffq, foodb_df, pid = "") {
  clean_ffq <- ffq %>%
    group_by(FooDB) %>%
    summarise(perc = sum(perc), .groups = "drop") %>%
    rename(name = FooDB, food_frequency = perc)

  foodb_df <- foodb_df %>%
    select(-food_frequency) %>%
    left_join(clean_ffq, by = "name")

  # spelling differs between the FFQ equivalents and FooDB
  # (applied after the join, as in the original runs)
  clean_ffq$name[clean_ffq$name == "Mustard Spinach"] <- "Mustard spinach"

  if (!identical(sort(clean_ffq$name), sort(foodb_df$name))) {
    stop("The names of the foods did not match for sample: ", pid)
  }
  foodb_df
}

# ---- IBD: ID mapping tables (built once, cached) -----------------------------

# Human gene symbol -> KEGG KO
build_gene_to_ko_map <- function(symbols, cache_path) {
  if (file.exists(cache_path)) return(readRDS(cache_path))
  message("Building gene -> KO mapping (first run only)...")

  sym_to_entrez <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = unique(symbols),
    keytype = "SYMBOL",
    columns = "ENTREZID"
  ) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  links <- KEGGREST::keggLink("ko", "hsa")
  entrez_to_ko <- tibble(
    ENTREZID = str_remove(names(links), "^hsa:"),
    KO       = str_remove(links, "^ko:")
  )

  gene_ko_map <- sym_to_entrez %>%
    left_join(entrez_to_ko, by = "ENTREZID") %>%
    select(SYMBOL, KO) %>%
    mutate(KO = ifelse(is.na(KO), "UNMAPPED", KO))

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(gene_ko_map, cache_path)
  gene_ko_map
}

# EC number -> KEGG KO
build_ec_to_ko_map <- function(cache_path) {
  if (file.exists(cache_path)) return(readRDS(cache_path))
  message("Building EC -> KO mapping (first run only)...")

  links <- KEGGREST::keggLink("ko", "ec")
  ec_to_ko <- tibble(
    EC = str_remove(names(links), "^ec:"),
    KO = str_remove(links, "^ko:")
  )

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ec_to_ko, cache_path)
  ec_to_ko
}

# ---- IBD: per-sample input files ----------------------------------------------

# First column whose name starts with the sample's External ID
find_sample_col <- function(df, external_id) {
  matches <- names(df)[str_starts(names(df), fixed(external_id))]
  if (length(matches) == 0) return(NULL)
  matches[1]
}

# Distinct participant-weeks that have a given data type
participant_weeks <- function(meta, dtype) {
  meta %>%
    filter(data_type == dtype) %>%
    distinct(`Participant ID`, week_num)
}

# External ID for one participant / week / data type
get_ext_id <- function(meta, pid, wk, dtype) {
  meta %>%
    filter(`Participant ID` == pid, week_num == wk, data_type == dtype) %>%
    pull(`External ID`) %>%
    dplyr::first()
}

# Microbial "EC: desc|g__Genus.s__species" table -> KO x taxonomy abundance
write_microbe_files <- function(df, external_id, ec_to_ko_map, abundance_path, ko_path) {
  if (is.null(external_id) || is.na(external_id)) return(invisible(NULL))
  feature_col <- names(df)[1]
  sample_col  <- find_sample_col(df, external_id)
  if (is.null(sample_col)) {
    warning(paste("No column match for", external_id))
    return(invisible(NULL))
  }

  df_clean <- df %>% filter(!str_starts(.data[[feature_col]], "#"))
  parts <- str_split_fixed(df_clean[[feature_col]], "\\|", 2)

  translated <- tibble(
    EC             = str_trim(str_remove(parts[, 1], ":.*$")),
    taxonomy       = ifelse(parts[, 2] == "", NA_character_, parts[, 2]),
    Abundance_RPKs = df_clean[[sample_col]]
  ) %>%
    left_join(ec_to_ko_map, by = "EC", relationship = "many-to-many") %>%
    mutate(KO = case_when(
      EC %in% c("UNMAPPED", "UNGROUPED") ~ EC,
      is.na(KO)                         ~ "UNMAPPED",
      TRUE                              ~ KO
    )) %>%
    group_by(KO, taxonomy) %>%
    summarise(Abundance_RPKs = sum(Abundance_RPKs, na.rm = TRUE), .groups = "drop") %>%
    filter(Abundance_RPKs != 0) %>%
    select(KO, taxonomy, Abundance_RPKs)

  write_csv(translated, abundance_path, na = "NA")
  writeLines(unique(translated$KO), ko_path)
}

# Host transcriptome (gene symbol counts) -> KO abundance
write_host_tx_files <- function(df, external_id, gene_ko_map, abundance_path, ko_path) {
  if (is.null(external_id) || is.na(external_id)) return(invisible(NULL))
  feature_col <- names(df)[1]
  sample_col  <- find_sample_col(df, external_id)
  if (is.null(sample_col)) {
    warning(paste("No column match for", external_id))
    return(invisible(NULL))
  }

  translated <- df %>%
    select(all_of(feature_col), all_of(sample_col)) %>%
    rename(SYMBOL = 1, Abundance = 2) %>%
    left_join(gene_ko_map, by = "SYMBOL", relationship = "many-to-many") %>%
    mutate(KO = ifelse(is.na(KO), "UNMAPPED", KO)) %>%
    group_by(KO) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(Abundance != 0)

  write_csv(translated, abundance_path, na = "NA")
  writeLines(unique(translated$KO), ko_path)
}

# One microbe + host input directory. Host is either the fixed reference genome
# (host_genome_files = c(abundance = ..., ko = ...)) or a host transcriptome
# sample (host_df + host_id + gene_ko_map).
build_pair_dir <- function(out_dir, microbe_df, microbe_id, ec_to_ko_map,
                           host_genome_files = NULL,
                           host_df = NULL, host_id = NULL, gene_ko_map = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write_microbe_files(
    df = microbe_df, external_id = microbe_id, ec_to_ko_map = ec_to_ko_map,
    abundance_path = file.path(out_dir, "ko_taxonomy_abundance.csv"),
    ko_path        = file.path(out_dir, "noquote_ko.txt")
  )

  host_abund <- file.path(out_dir, "host_ko_abundance.csv")
  host_ko    <- file.path(out_dir, "noquote_ko_host.txt")

  if (!is.null(host_genome_files)) {
    file.copy(host_genome_files[["abundance"]], host_abund, overwrite = TRUE)
    file.copy(host_genome_files[["ko"]],        host_ko,    overwrite = TRUE)
  } else if (!is.null(host_id) && !is.na(host_id)) {
    write_host_tx_files(host_df, host_id, gene_ko_map, host_abund, host_ko)
  }
  invisible(NULL)
}


# =============================================================================
# 2. Collapsing DM[H]N outputs
# =============================================================================

# Stack every <run_dir>/<output>/graph/graph_results.csv under `root` into one
# file. Adds `run_dir` (path relative to root, e.g. "ZIM089" or
# "C3006_week0/hg_mgx") and `output` (e.g. "output_dm", "output_dmh").
# `outputs` optionally restricts which output folders are included.
collapse_graph_results <- function(root, out_file, outputs = NULL) {
  rel <- list.files(root, pattern = "^graph_results\\.csv$",
                    recursive = TRUE, full.names = FALSE)
  rel <- rel[basename(dirname(rel)) == "graph"]

  output_dir <- basename(dirname(dirname(rel)))
  run_dir    <- dirname(dirname(dirname(rel)))
  if (!is.null(outputs)) {
    keep <- output_dir %in% outputs
    rel <- rel[keep]; output_dir <- output_dir[keep]; run_dir <- run_dir[keep]
  }
  if (length(rel) == 0) stop("No graph/graph_results.csv files found under ", root)

  # read as text so files with differently-typed empty columns still stack
  collapsed <- pmap_dfr(
    list(file.path(root, rel), run_dir, output_dir),
    function(f, rd, od) {
      read_csv(f, col_types = cols(.default = col_character()), progress = FALSE) %>%
        mutate(run_dir = rd, output = od, .before = 1)
    }
  ) %>%
    type_convert(col_types = cols(run_dir = col_character(),
                                  output  = col_character(),
                                  .default = col_guess()))

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  write_csv(collapsed, out_file)

  message(sprintf("%s: %d files, %d runs, %d rows -> %s",
                  root, length(rel), n_distinct(paste(run_dir, output_dir)),
                  nrow(collapsed), out_file))
  print(count(collapsed, output, name = "rows"))
  invisible(collapsed)
}

# Read a collapsed graph_results file (keeps run_dir as text, guesses types
# from every row so sparse columns are not misread)
read_collapsed_results <- function(path) {
  read_csv(path,
           col_types = cols(run_dir = col_character(),
                            output  = col_character(),
                            .default = col_guess()),
           guess_max = Inf, show_col_types = FALSE)
}

# One network_summary.txt -> one-row tibble (NA row if the file is missing)
parse_network_summary <- function(file_path) {
  origins <- c("microbe", "hostmicrobe", "none", "host", "diet",
               "all", "hostdiet", "microbediet")

  if (is.na(file_path) || !file.exists(file_path)) {
    empty <- as_tibble(setNames(
      rep(list(NA_integer_), 3 + length(origins)),
      c("n_reactions", "n_nodes", "n_edges", paste0("origin_", origins))
    ))
    return(empty)
  }

  lines <- readLines(file_path, warn = FALSE)

  extract_number <- function(pattern) {
    hit <- lines[str_detect(lines, pattern)]
    if (length(hit) == 0) return(NA_integer_)
    as.integer(str_extract(hit[1], "\\d+"))
  }
  extract_origin_count <- function(category) {
    hit <- lines[str_detect(lines, paste0("^", category, ":\\s*\\d+"))]
    if (length(hit) == 0) return(NA_integer_)
    as.integer(str_extract(hit[1], "\\d+$"))
  }

  bind_cols(
    tibble(n_reactions = extract_number("number of reactions were found"),
           n_nodes     = extract_number("nodes were created"),
           n_edges     = extract_number("edges were created")),
    as_tibble(setNames(lapply(origins, extract_origin_count), paste0("origin_", origins)))
  )
}

# Number of non-blank lines in a KO list file (NA if missing)
count_ko_lines <- function(file_path) {
  if (is.na(file_path) || !file.exists(file_path)) return(NA_integer_)
  kos <- readLines(file_path, warn = FALSE)
  length(kos[nzchar(kos)])
}

# Shannon diversity and number of taxa from a ko_taxonomy_abundance.csv,
# summing Abundance_RPKs across KOs within each taxon
compute_shannon <- function(file_path) {
  if (is.na(file_path) || !file.exists(file_path)) {
    return(tibble(shannon = NA_real_, n_taxa = NA_integer_))
  }
  taxon_profile <- read_csv(file_path, show_col_types = FALSE) %>%
    group_by(taxonomy) %>%
    summarise(total_abundance = sum(Abundance_RPKs, na.rm = TRUE), .groups = "drop")

  tibble(shannon = vegan::diversity(taxon_profile$total_abundance, index = "shannon"),
         n_taxa  = nrow(taxon_profile))
}

# Add the input/output file paths for each row of a benchmark csv.
# old_prefix/new_prefix rewrite HPC paths to local ones.
benchmark_paths <- function(bench, output_dir, old_prefix = NULL, new_prefix = NULL,
                            has_host = TRUE) {
  bench %>%
    mutate(
      directory_clean      = str_remove(directory, "/$"),
      directory_clean      = if (!is.null(old_prefix))
                               str_replace(directory_clean, fixed(old_prefix), new_prefix)
                             else directory_clean,
      participant_id       = basename(directory_clean),
      network_summary_path = file.path(directory_clean, output_dir, "graph", "network_summary.txt"),
      noquote_ko_path      = file.path(directory_clean, "noquote_ko.txt"),
      noquote_ko_host_path = if (has_host) file.path(directory_clean, "noquote_ko_host.txt")
                             else NA_character_,
      ko_tax_abund_path    = file.path(directory_clean, "ko_taxonomy_abundance.csv")
    )
}

# Add network size, origin counts, diversity and KO counts to a benchmark
# table built with benchmark_paths(); drops the path columns
add_benchmark_stats <- function(bench) {
  network_stats   <- map_dfr(bench$network_summary_path, parse_network_summary)
  diversity_stats <- map_dfr(bench$ko_tax_abund_path, compute_shannon)

  bind_cols(bench, network_stats, diversity_stats) %>%
    mutate(
      n_ko       = map_int(noquote_ko_path, count_ko_lines),
      n_ko_host  = map_int(noquote_ko_host_path, count_ko_lines),
      n_ko_total = coalesce(n_ko, 0L) + coalesce(n_ko_host, 0L)
    ) %>%
    select(-directory_clean, -network_summary_path,
           -noquote_ko_path, -noquote_ko_host_path, -ko_tax_abund_path)
}


# =============================================================================
# 3. Parsing DM[H]N list columns
# =============================================================================

# "['a', 'b']" -> list(c("a", "b")); vectorised over x
parse_py_list <- function(x) {
  str_extract_all(x, "'([^']*)'") %>%
    map(~ str_remove_all(.x, "'"))
}

# As parse_py_list, also splitting entries that hold comma-separated organisms
parse_py_list_orgs <- function(x) {
  parse_py_list(x) %>%
    map(~ unlist(str_split(.x, ",\\s*")))
}

# Single string; handles single or double quotes and escaped quotes
# (food names such as "Hen's egg")
parse_py_list_safe <- function(x) {
  m <- gregexpr("'((?:[^'\\\\]|\\\\.)*)'|\"((?:[^\"\\\\]|\\\\.)*)\"", x, perl = TRUE)
  matches <- regmatches(x, m)[[1]]
  gsub("^['\"]|['\"]$", "", matches)
}


# =============================================================================
# 4. Reference lookups
# =============================================================================

read_hmdb <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    distinct() %>%
    na.omit()
}

# All KEGG compound IDs with their first name; one API call, optionally cached
kegg_compound_names <- function(cache_path = NULL) {
  if (!is.null(cache_path) && file.exists(cache_path)) return(readRDS(cache_path))

  out <- KEGGREST::keggList("compound") %>%
    enframe(name = "kegg_id", value = "names") %>%
    mutate(kegg_id       = str_remove(kegg_id, "^cpd:"),
           compound_name = str_trim(str_extract(names, "^[^;]+"))) %>%
    select(kegg_id, compound_name)

  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_path)
  }
  out
}

# tibble(<id_col>, name) for a set of compound IDs (NA name if not in KEGG)
kegg_name_lookup <- function(ids, kegg_names, id_col = "kegg_id") {
  tibble(kegg_id = unique(ids)) %>%
    left_join(kegg_names, by = "kegg_id") %>%
    rename(name = compound_name) %>%
    rename(!!id_col := kegg_id)
}


# =============================================================================
# 5. Enrichment and set similarity
# =============================================================================

# KEGG pathway enrichment on a KO list
run_ko_enrichment <- function(ko_list,
                              universe      = NULL,
                              pvalueCutoff  = 0.05,
                              qvalueCutoff  = 0.2,
                              pAdjustMethod = "BH",
                              minGSSize     = 5,
                              maxGSSize     = 500) {
  ko_list <- unique(trimws(ko_list))
  ko_list <- ko_list[ko_list != ""]

  clusterProfiler::enrichKEGG(
    gene          = ko_list,
    organism      = "ko",
    keyType       = "kegg",
    universe      = universe,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize
  )
}

# Fisher's exact test of HMDB taxonomy-class enrichment against all of HMDB.
# class_counts: one row per taxonomy_class x group with count column `n`.
#   k = compounds of the class in the group,  n = all compounds in the group
#   K = compounds of the class in HMDB,       N = all compounds in HMDB
taxonomy_class_enrichment <- function(class_counts, hmdb, group_cols) {
  ref     <- hmdb %>% count(taxonomy_class, name = "K")
  N_total <- sum(ref$K)

  class_counts %>%
    distinct(across(all_of(c("taxonomy_class", group_cols))), n) %>%
    rename(k = n) %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(n = sum(k)) %>%
    ungroup() %>%
    left_join(ref, by = "taxonomy_class") %>%
    mutate(N = N_total) %>%
    rowwise() %>%
    mutate(
      ft = list(fisher.test(matrix(c(k, n - k, K - k, (N - n) - (K - k)),
                                   nrow = 2, byrow = TRUE))),
      p_value    = ft$p.value,
      odds_ratio = unname(ft$estimate),
      ci_low     = ft$conf.int[1],
      ci_high    = ft$conf.int[2]
    ) %>%
    ungroup() %>%
    select(-ft) %>%
    group_by(across(all_of(group_cols))) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      fold_enrichment = (k / n) / (K / N),
      log2_or         = log2(odds_ratio),
      neg_log10_p     = -log10(p_adj),
      direction = case_when(
        p_adj < 0.05 & odds_ratio > 1 ~ "Enriched",
        p_adj < 0.05 & odds_ratio < 1 ~ "Depleted",
        TRUE                          ~ "Not significant"
      )
    ) %>%
    arrange(across(all_of(group_cols)), p_adj)
}

jaccard <- function(a, b) {
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  length(intersect(a, b)) / uni
}


# =============================================================================
# 6. Reaction composition (PERMANOVA / PCoA / SIMPER)
# =============================================================================

# PID x reaction abundance matrix from diet -> microbe edges
build_reaction_matrix <- function(hits, meta_cols = c("cohort", "gender", "diet")) {
  meta <- hits %>% distinct(PID, across(all_of(meta_cols)))

  wide <- hits %>%
    distinct(PID, reaction, abundance) %>%
    pivot_wider(names_from = reaction, values_from = abundance, values_fill = 0) %>%
    left_join(meta, by = "PID") %>%
    relocate(PID, all_of(meta_cols))

  mat <- wide %>% select(-PID, -all_of(meta_cols)) %>% as.matrix()
  rownames(mat) <- wide$PID

  list(wide = wide, meta = select(wide, PID, all_of(meta_cols)), mat = mat)
}

# Bray-Curtis PERMANOVA (marginal terms), dispersion check and PCoA
reaction_composition_test <- function(rm, rhs = "cohort + gender + diet",
                                      group = "diet", permutations = 1000) {
  dist_bc <- vegan::vegdist(rm$mat, method = "bray")
  f <- as.formula(paste("dist_bc ~", rhs))

  permanova  <- vegan::adonis2(f, data = rm$meta, permutations = permutations, by = "margin")
  dispersion <- anova(vegan::betadisper(dist_bc, rm$meta[[group]]))

  pcoa <- cmdscale(dist_bc, eig = TRUE)
  eig  <- pcoa$eig
  points <- as.data.frame(pcoa$points) %>%
    setNames(c("PCo1", "PCo2")) %>%
    mutate(PID = rm$meta$PID, group = rm$meta[[group]])

  list(dist       = dist_bc,
       permanova  = permanova,
       dispersion = dispersion,
       pcoa       = points,
       pct_var    = (eig / sum(eig[eig > 0]) * 100)[1:2],
       eig_check  = c(n_negative   = sum(eig < 0),
                      sum_positive = sum(eig[eig > 0]),
                      sum_negative = sum(eig[eig < 0])))
}

plot_reaction_pcoa <- function(comp, color_lab = "diet") {
  ggplot2::ggplot(comp$pcoa, ggplot2::aes(PCo1, PCo2, color = group)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::stat_ellipse() +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(
      x = paste0("PCo1 (", round(comp$pct_var[1], 1), "%)"),
      y = paste0("PCo2 (", round(comp$pct_var[2], 1), "%)"),
      color = color_lab,
      title = "PCoA of reaction abundance profiles by diet (Bray-Curtis)"
    )
}

# SIMPER reactions driving the group difference (p <= 0.05), then PERMANOVA
# rerun without them
simper_robustness <- function(rm, group = "diet", rhs = "cohort + gender + diet",
                              simper_perm = 999, permanova_perm = 1000) {
  simper_sum <- summary(vegan::simper(rm$mat, group = rm$meta[[group]],
                                      permutations = simper_perm))

  sig <- simper_sum[[1]] %>%
    as.data.frame() %>%
    rownames_to_column("reaction") %>%
    filter(p <= 0.05) %>%
    pull(reaction)

  reduced <- rm$mat[, !(colnames(rm$mat) %in% sig), drop = FALSE]
  dist_bc <- vegan::vegdist(reduced, method = "bray")
  f <- as.formula(paste("dist_bc ~", rhs))

  list(simper        = simper_sum,
       sig_reactions = sig,
       permanova_reduced = vegan::adonis2(f, data = rm$meta,
                                          permutations = permanova_perm, by = "margin"))
}


# =============================================================================
# 7. Clinical / immune marker models
# =============================================================================

# ZIM experiment names for figure labels
label_zim_experiment <- function(x) {
  recode(x, DietMicrobe_AllFoods = "DietMicrobe AllFoods",
            DietMicrobe_FFQ      = "DietMicrobe FFQ")
}

# Fit the same lm separately for each marker; BH across all returned terms
fit_lm_by_marker <- function(data, formula, marker_col, markers) {
  map_dfr(markers, function(i) {
    lm(formula, data = filter(data, .data[[marker_col]] == i)) %>%
      broom::tidy() %>%
      mutate(marker = i)
  }) %>%
    mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
    arrange(p_adj)
}

# ---- DM: host KOs vs clinical marker change ----------------------------------

# Spearman correlation between a patient's number of distinct KOs in each
# group (category / subcategory / pathway) and the marker change
cor_by_level <- function(level, sig_ko_host, deltas_long, patients, min_pairs) {
  sig_ko_host %>%
    distinct(RecordID, group = .data[[level]], ko) %>%
    filter(!is.na(group)) %>%
    count(RecordID, group, name = "n_ko") %>%
    complete(RecordID = patients, group, fill = list(n_ko = 0)) %>%
    inner_join(deltas_long, by = "RecordID", relationship = "many-to-many") %>%
    group_by(group, marker, change) %>%
    filter(n() >= min_pairs, n_distinct(n_ko) > 1, n_distinct(delta) > 1) %>%
    summarise(n   = n(),
              rho = cor(n_ko, delta, method = "spearman"),
              p   = suppressWarnings(cor.test(n_ko, delta, method = "spearman")$p.value),
              .groups = "drop") %>%
    group_by(marker, change) %>%
    mutate(q = p.adjust(p, "BH")) %>%
    ungroup() %>%
    mutate(q_all = p.adjust(p, "BH"))
}

# Rank models for one pathway x marker x comparison:
#   main:        rank(delta) ~ rank(frac) + Diet + rank(wk1)
#   interaction: rank(delta) ~ rank(frac) * Diet + rank(wk1)
fit_marker_models <- possibly(function(d) {
  m1 <- summary(lm(rank(delta) ~ rank(frac) + Diet + rank(wk1), data = d))$coefficients
  m2 <- summary(lm(rank(delta) ~ rank(frac) * Diet + rank(wk1), data = d))$coefficients
  tibble(beta_main     = m1["rank(frac)", "Estimate"],
         p_main        = m1["rank(frac)", "Pr(>|t|)"],
         beta_agrarian = m2["rank(frac)", "Estimate"],
         p_agrarian    = m2["rank(frac)", "Pr(>|t|)"],
         beta_diff     = m2["rank(frac):DietWestern", "Estimate"],
         p_int         = m2["rank(frac):DietWestern", "Pr(>|t|)"])
}, otherwise = tibble(beta_main = NA_real_, p_main = NA_real_,
                      beta_agrarian = NA_real_, p_agrarian = NA_real_,
                      beta_diff = NA_real_, p_int = NA_real_))


# =============================================================================
# 8. Resource-use models (benchmarking)
# =============================================================================

# Fit the same log-scale mixed model with Shannon vs number of taxa as the
# diversity term, for comparison with performance::compare_performance()
fit_diversity_pair <- function(response, data, covariates) {
  rhs <- paste(c("log(n_ko)", covariates, "(1 | participant_id)"), collapse = " + ")
  list(
    shan = lmer(as.formula(paste0("log(", response, ") ~ shannon + ", rhs)),
                data = data, REML = FALSE),
    taxa = lmer(as.formula(paste0("log(", response, ") ~ log(n_taxa) + ", rhs)),
                data = data, REML = FALSE)
  )
}

# Boxplot + jitter of one metric per benchmark source
source_boxplot <- function(data, y, ylab, legend = "none") {
  ggplot2::ggplot(data, ggplot2::aes(x = source, y = .data[[y]], fill = source)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.5) +
    ggplot2::labs(x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = legend,
                   axis.text.x  = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank())
}
