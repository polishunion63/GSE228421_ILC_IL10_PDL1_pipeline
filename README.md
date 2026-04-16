# Analysis of ILCs in Human Psoriasis (GSE228421)

This repository contains the R analysis pipeline used for the single-cell RNA sequencing (scRNA-seq) reanalysis in:

**PD-L1hiSca-1+ innate lymphoid cells restrain psoriatic inflammation through IL-10-dependent regulation of CCR2+ γδT17 trafficking**  
(Jo et al., *Journal of Investigative Dermatology*, submitted)

## Overview

We reanalyzed the public dataset **GSE228421** to quantify **IL-10+** and **PD-L1+** innate lymphoid cells (ILCs) in human psoriatic skin.

The repository includes code for:

- loading 10x Genomics count matrices
- basic cell filtering (gene/UMI thresholds)
- lineage-based ILC gating
- detection of PD-L1+ and IL-10+ ILC subsets
- sample- and group-level summaries for downstream statistics and plotting
- additional subtype annotation of ILC2-, ILC3-, and ILC1-like populations for revised Figure 3 / supplementary analyses

## Software requirements

- R (tested with **4.5.1**)

Required R packages:

- Matrix
- readr
- dplyr
- stringr
- tibble
- tidyr

All packages can be installed from CRAN using:

```r
install.packages(c("Matrix", "readr", "dplyr", "stringr", "tibble", "tidyr"))
