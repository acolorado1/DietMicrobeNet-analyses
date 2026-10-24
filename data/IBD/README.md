# Diet Modification Cohort Files 

These are all files references in the scripts associated with the Inflammatory Bowel Disease cohort anaysis, from data gathering/formatting, to `DM[H]N` usage, and following analyses. 

## General files 

These files were all subset from the publicly available [site](https://ibdmdb.org/) associated with this study. 

- `sub_metaG.tsv`: Microbiome metagenome, subset to relevant participants 
- `sub_metaT.tsv`: Microbiome metatranscriptome, subset to relevant participants
- `metadata.csv`: Cohort metadata, subset to only relevant participants 
- `host_tx_counts.tsv`: human transcriptomes 

## Per participant files 

- `hg_mgx`: human reference genome + microbiome metagenome 
- `hg_mtx`: human reference genome + microbiome metatranscriptome
- `htx_mgx`: human transcriptome + microbiome metagenome 
- `htx_mtx`: human transcriptome + microbiome metatranscriptome

Within each of these files there is: 

- `host_ko_abundance.csv`: Host KOs (KEGG reference genomes) and abundance 
- `ko_taxonomy_abundance.csv`: Microbiome organims and read abundance information 
- `noquote_ko_host`: List of host KOs 
- `noquote_ko`: List of microbe KOs

Every file here is used as an input in DM[H]N.  