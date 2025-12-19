# Analysis of ILCs in Human Psoriasis (GSE228421)

This repository contains the R analysis pipeline used for the single-cell RNA sequencing (scRNA-seq) reanalysis in:

PD-L1hiSca-1+ innate lymphoid cells restrain psoriatic inflammation through IL-10-dependent regulation of CCR2+ γδT17 trafficking  
(Jo et al., Journal of Investigative Dermatology, submitted)

## Overview

We reanalyzed the public dataset GSE228421 to quantify IL-10+ and PD-L1+ innate lymphoid cells (ILCs) in human psoriatic skin.  
The pipeline performs:

- Loading of 10x Genomics count matrices
- Basic cell filtering (gene/UMI thresholds)
- Lineage-based ILC gating
- Detection of PD-L1+ and IL-10+ ILC subsets
- Sample- and group-level summaries for downstream statistics and plotting

## Software requirements

- R (tested with 4.5.1)
- R packages:
  - Matrix (e.g. 1.7.4)
  - readr
  - dplyr
  - stringr
  - tibble

All packages can be installed from CRAN using `install.packages()`.

## Input data

Raw data are available from NCBI GEO:

- Accession: GSE228421
- Expected input: standard 10x Genomics output for each sample  
  - `matrix.mtx` or `matrix.mtx.gz`  
  - `features.tsv` or `features.tsv.gz` (or `genes.tsv`)  
  - `barcodes.tsv` or `barcodes.tsv.gz`

You do **not** need to modify the script for individual sample names, as long as the 10x files follow this pattern.

## Usage

1. Download the raw data files from GEO into a local directory.
2. Place `GSE228421_ILC_pipeline.R` in the same or another working directory.
3. Open R or RStudio and set the working directory to the location of the script.
4. Run the script, for example:
   ```r
   source("GSE228421_ILC_pipeline.R")
