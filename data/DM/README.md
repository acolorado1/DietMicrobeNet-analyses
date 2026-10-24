# Diet Modification Cohort Files 

These are all files references in the scripts associated with the Diet Modification cohort anaysis, from data gathering/formatting, to `DM[H]N` usage, and following analyses. 

## General files 

- `DM_KOs_rnaseq.txt`: Microbiome mmetatranscriptome
- `DM_KOs_shotgun.tsv`: Microbiome metagenome
- `host_gene_counts_tmp.csv`: Host transcriptome
- `agrarian_foodb.xlsx`: Publication has sample agrarian meals, this contains foodb equivalents 
- `western_foodb.xlsx`: Publication has sample western meals, this contains foodb equivalents 
- `metadata.xlsx`: Cohort metadata

## Per participant files 

- `foodb_foods_dataframe.csv`: Web application output (this will either be the sample agrarian or sample western meal)
- `host_ko_abundance.csv`: Host KOs (KEGG reference genomes) and abundance 
- `ko_taxonomy_abundance.csv`: Microbiome organims and read abundance information 
- `noquote_ko_host`: List of host KOs 
- `noquote_ko`: List of microbe KOs

Every file here is used as an input in DM[H]N.  
