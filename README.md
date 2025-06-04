# ISLAND2022-WATER COLUMN DYNAMICS
Bacterial and eukaryotic microbiota revealed through 16S and 18S rRNA metabarcoding of water column at Lastovo Island, South Adriatic Sea.

## PAPER SUBMITTED TO Marine Environmental Research
## Influence of Internal Island-Trapped Waves on Plankton Structure and Trophic Networks in Stratified Oligotrophic Coastal Waters  
Maja Mucko(1*), Luca Russo(2), Antonija Matek(1), Filip Grgurević(1), Branka Pestorić(3), Eric P. Achterberg(4), Domenico D’Alelio(2), Zrinka Ljubešić(1)

(1)University of Zagreb, Faculty of Science, Biology Department
(2)Department of Integrative Marine Ecology, Stazione Zoologica Anton Dohrn, Naples, Italy
(3)University of Montenegro, Institute of Marine Biology
(4)4GEOMAR Helmholtz Centre for Ocean Research
*corresponding author: maja.mucko@biol.pmf.hr


## ABSTRACT: 
During periods of water column stratification in marine oligotrophic ecosystems, physical forcings such as internal island-trapped waves (ITWs) can facilitate nutrient fluxes to surface waters and determine fine-scale changes in microbial communities. During a two-week in situ experiment at Lastovo Island, South Adriatic Sea, we conducted water-column community diversity assessments with environmental metabarcoding of plankton during and after ITWs events. Bacteria and eukaryotes communities showed significant dissimilarities between size-fraction and depths, with ITW event significantly contributing to clustering of eukaryotic communities. Major bacterial contributors were Gammaproteobacteria, Cyanobacteria, Bacteroidota and Verrucomicrobiota following ITW event. Bacterial functional profiling indicated that ureolysis, aerobic ammonia oxidation, nitrification and ectoparasitic or predatory roles were directly linked to shifting environmental parameters in the water column. Metazoan sequences (mainly Arthropoda, class Copepoda) dominated the micro size fraction, while various dinoflagellates (with high contribution of parasitic Syndiniales) dominated nano and pico size fractions. Primary producers were Mamiellophyceae, Prymnesiophyceae and Bacillariophaceae, with highest relative abundances in the deep chlorophyll maximum (DCM). Assessment of genus-level networks in surface waters and the DCM revealed that about 34% and 39% co-occurrences, respectively, were attributable to putative trophic interactions with protists dominating over zooplankton taxa in both communities. Important network taxa connectors were mainly identified among autotrophic protists within DCM, while in surface network connectors number generally decreased.

## KEYWORDS:
oligotrophic stratified ecosystem, bacteria, eukaryotes, trophic associations, internal island-trapped waves

## CONTENT IN GITHUB REPOSITORY:
### Raw files 
(24 samples, paired-end reads) submmited to European Nucleotide Archive (ENA: https://www.ebi.ac.uk/ena/browser/search) under project number PRJEB63220; samples W1-W24
### Qiime2 code
DATA FILE: Mucko_etal_2025_WORKFLOW_QIIME2.nb

Whole bioinformatic pipeline composed of steps:
1. Environment activation
2. Data import
3. Data denoise
4. Feature table generation
5. Taxonomy assignment
6. Phylogenetic tree reconstruction
8. Alpha & Beta Diversity
### Analysis output visualizations
.qzv artifact visualizations of all bioinformatic pipeline steps available for viewing interactively via view.qiime2.org

BACTERIAL DATASET

16S_demux-paired-end.qzv - interactive plots of imported raw reads

WATER-denoising-stats_final.qzv - denoising statistics after DADA2 tool application on raw reads

WATER_16s_table_final.qzv - feature table of identified ASVs

WATER_16s_repseqs_final.qzv - representative sequences of identified ASVs

WATER_16S_taxonomy_gg2.qzv - identified taxonomy of ASVs

taxa-barplots_16S_final.qzv - taxonomy barplots on sample level

taxa-barplots_16S_depth.qzv - taxonomy barplots on grouped depth category

taxa-barplots_16S_fraction.qzv - taxonomy barplots on grouped fraction category


EUCARYOTIC DATASET

demux-paired-end_18S.qzv - interactive plots of imported raw reads

WATER_18s_denoising-stats_trunc220.qzv - denoising statistics after DADA2 tool application on raw reads

WATER_18s_table_trunc220.qzv - feature table of identified ASVs

WATER_18s_repseqs_trunc220.qzv - representative sequences of identified ASVs

WATER_18s_taxonomy_pr2.qzv - identified taxonomy of ASVs

taxa-barplot_18S.qzv - taxonomy barplots on sample level

taxa-barplots_18S_depth.qzv - taxonomy barplots on grouped depth category

taxa-barplots_18S_fraction.qzv - taxonomy barplots on grouped fraction category

### Metadata table
Full metadata table describing each sample, variable for their discriminations and sample groups
metadata_water_ISLAND2022_envs.tsv

### R analysis scripts: 
16S_18S_cluster_tanglegram_and_heatmap.R

16S_downstream_analysis.R

16S_FunctionalGroup_Spearman_heatmap.R

18S_downstream_analysis.R

Genus_NETWORKs.R

map.R

modularity analyses.R

