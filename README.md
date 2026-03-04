# experimental-gingivitis-metagenomics
The complete R code for data analysis used in published paper: Metagenomic characterization fo the supragingival and salivary microbiota during initiation and resolution of experimental gingivitis

## Repository purpose 

This repository contains the workflow for taxonomic and functional profiling of microbiome data used in the study: 

* Importing **MetaPhlAn** and **HUMAnN** profiling datasets
* Alpha-diversity analysis
* Beta-diversity analysis
* Differential abundance analysis using **MaAsLin3**

## Short description of R-code (in order of execution) 

**1. Data import and data preprocessing**

_File:_ `R/metaphlan_to_phyloseq.Rmd`

_File:_ `R/humann_to_phyloseq.Rmd`
* Load MetaPhlAn4 taxonomic profiles
* Load HUMAnN functional profiles
* Import and merge metadata

**2. Alpha-diversity analysis**

_File:_ `R/alpha_v2.Rmd` 
* Calculate diversity metrics (Shannon index)
* Statistical comparison across sample type and timepoints 
* Visualization

**3. Beta-diversity analysis**

_File:_ `R/beta_v5.Rmd`
* Compute distance matrices (robust Aitchison distance)
* Ordination analysis and visualization (PCoA)
* Statistical comparisons across sample type and timepoints

**4. Differential abundance analysis** 
_File:_ `R/DA_v2.Rmd`
* Model setup using MaAsLin3
* Identification of significantly associated species and pathways 
