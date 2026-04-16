# Analysis of ILCs in Human Psoriasis (GSE228421)

This repository contains the R analysis pipeline used for reanalysis of the public human psoriasis skin scRNA-seq dataset GSE228421.

## Overview

We reanalyzed GSE228421 to quantify IL-10+ and PD-L1+ innate lymphoid cells (ILCs) in human psoriatic skin.

This repository includes code for:

- loading 10x Genomics count matrices
- basic cell filtering using gene/UMI thresholds
- lineage-based ILC gating
- detection of PD-L1+ and IL-10+ ILC subsets
- sample- and group-level summaries
- additional subtype annotation of ILC2-, ILC3-, and ILC1-like populations
- exploratory γδ T-cell analysis in the same dataset

## Software requirements

- R (tested with 4.5.1)

Required R packages:

- Matrix
- readr
- dplyr
- stringr
- tibble
- tidyr

Packages can be installed from CRAN using:

```r
install.packages(c("Matrix", "readr", "dplyr", "stringr", "tibble", "tidyr"))
```

## Input data

Raw data are publicly available from NCBI GEO:

- Accession: GSE228421

Expected input for each sample is standard 10x Genomics output:

- `matrix.mtx` or `matrix.mtx.gz`
- `features.tsv` or `features.tsv.gz` (or `genes.tsv`)
- `barcodes.tsv` or `barcodes.tsv.gz`

## Main pipeline

### `GSE228421_ILC_pipeline.R`

This script performs the core reanalysis of GSE228421:
- loading raw 10x matrices
- basic barcode filtering
- merging filtered matrices across samples
- lineage-based ILC gating
- quantification of total ILC, IL-10+ ILC, PD-L1+ ILC, and DP ILC (PD-L1+IL-10+ ILC)
- generation of summary tables

### Usage

Run:

```r
source("GSE228421_ILC_pipeline.R")
```

## Additional revision analysis

### `02_GSE228421_ILC_subtype_annotation.R`

This script performs additional subtype annotation of human skin ILCs in GSE228421.

Subtype definitions used:

- ILC2-like: GATA3, IL1RL1 (ST2), KLRG1, HPGDS, PTGDR2, IL7R
- ILC3-like: RORC
- ILC1-like: TBX21

This script:
- loads the merged filtered count matrix
- reproduces ILC gating
- annotates ILC2-, ILC3-, and ILC1-like subsets
- summarizes subtype composition in:
  - All ILC
  - PD-L1+ ILC
  - DP ILC (PD-L1+IL-10+ ILC)

### Usage

Run after the main pipeline:

```r
source("02_GSE228421_ILC_subtype_annotation.R")
```

Processed output files:
- `GSE228421_ILC2_ILC3_ILC1_in_PDL1_and_DP.csv`
- `GSE228421_ILC_subtype_by_group4_long.csv`
- `GSE228421_ILC_subtype_by_sample_long.csv`
- `GSE228421_ILC_per_cell_with_subtype.csv`

### `03_GSE228421_gdT_exploratory_gating.R`

This script performs an exploratory γδ T-cell analysis in GSE228421 using a strict transcript-based gate.

Strict γδ T-cell definition used:
- TRDC > 0
- TRGC1/2 > 0
- TRAC == 0

This script:
- loads the merged filtered count matrix
- applies a strict γδ T-cell gating strategy
- summarizes γδ T-cell counts and frequencies at the sample level with group labels

### Usage

Run after the main pipeline:

```r
source("03_GSE228421_gdT_exploratory_gating.R")
```

Processed output files:
- `gdT_by_sample_with_group.csv`

## Notes

- `DP ILC` in this repository refers to PD-L1+IL-10+ ILC.
- This repository contains analysis code and processed summary tables generated from the public dataset GSE228421.
- The original raw scRNA-seq data are hosted by GEO and are not duplicated here.
- The γδ T-cell analysis is exploratory and highly conservative, and may underestimate γδ T-cell abundance in sparse droplet-based scRNA-seq data.
