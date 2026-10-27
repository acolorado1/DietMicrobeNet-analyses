# Zimbabwe Cohort Files 

These are all files references in the scripts associated with the Zimbabwe cohort anaysis, from data gathering/formatting, to `DM[H]N` usage, and following analyses. 

## General files 

- `Zim_KOs.tsv`: Full humann3 output, per sample organism and ko abundances in RPKs
- `zim2_baseline.txt`: Baseline metadata (see publication)
- `zim2_longitudinal.txt`: Metadata containing paired longitudinal data (see publication)
- `FoodEquivalents.xlsx`: Original Zimbabwe FFQ foods items and their FooDB equivalents 
- `graph_results_ZIM.csv.gz`: aggregated graph_results.csv (DM[H]N output) across all samples

## Per participant files 

- `FFQ.csv`: FFQ with foods, scores, calculated percentages, and FooDB equivalents
- `foodb_foods_dataframe.csv`: Web application output
- `host_ko_abundance.csv`: Host KOs (KEGG reference genomes) and abundance 
- `ko_taxonomy_abundance.csv`: Microbiome organims and read abundance information 
- `noquote_ko_host`: List of host KOs 
- `noquote_ko`: List of microbe KOs

Every file here besides `FFQ.csv` is used as an input in DM[H]N.  
