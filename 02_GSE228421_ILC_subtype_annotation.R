############################################################
## 02_GSE228421_ILC_subtype_annotation.R
## Purpose:
##   - Load merged filtered count matrix (RDS)
##   - Reproduce ILC gating used in GSE228421 analysis
##   - Annotate ILC2/ILC3/ILC1-like subsets
##   - Summarize subtype composition in:
##       1) All ILC
##       2) PD-L1+ ILC
##       3) PD-L1+IL-10+ ILC (DP_ILC)
##   - Save processed summary tables for manuscript/Supplement
############################################################

rm(list = ls()); invisible(gc())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

OUTDIR <- file.path(getwd(), "gse228421_out")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

logi <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%H:%M:%S"),
              paste(..., collapse = " ")))
}

## ------------------------------------------------------------
## 0. Load merged filtered matrix
## ------------------------------------------------------------
rds_path <- file.path(OUTDIR, "GSE228421_merged_counts_filtered.rds")
if (!file.exists(rds_path)) {
  stop("Cannot find merged filtered RDS: ", rds_path)
}

M <- readRDS(rds_path)
if (!inherits(M, "dgCMatrix")) {
  stop("Loaded object is not a dgCMatrix.")
}

logi(sprintf("Loaded matrix: %d genes x %d cells", nrow(M), ncol(M)))

## ------------------------------------------------------------
## 1. Basic metadata from cell barcodes
## ------------------------------------------------------------
cn <- colnames(M)

sample_id <- sub("_[^_]+$", "", cn, perl = TRUE)

site <- dplyr::if_else(stringr::str_detect(sample_id, "_NL$"), "NL", "L")

visit <- dplyr::case_when(
  stringr::str_detect(sample_id, "_V1_") ~ "V1",
  stringr::str_detect(sample_id, "_V2_") ~ "V2",
  stringr::str_detect(sample_id, "_V3_") ~ "V3",
  TRUE ~ NA_character_
)

group4 <- paste0(site, "_", visit)

## ------------------------------------------------------------
## 2. Gene helper functions
## ------------------------------------------------------------
uniq_pos <- function(x) unique(x[!is.na(x) & x > 0])

find_any <- function(symbol, ensg = NULL) {
  rn <- rownames(M)
  hits <- integer(0)

  ## exact symbol match
  hits <- c(hits, which(rn == symbol))

  ## composite rownames like SYMBOL|ENSG / SYMBOL_ENSG / SYMBOL.ENSG
  patt <- paste0("(^|[|;,_\\s\\.])", symbol, "($|[|;,_\\s\\.])")
  hits <- c(hits, grep(patt, rn))

  ## optional Ensembl match
  if (!is.null(ensg)) {
    patt2 <- paste0("^", ensg, "(\\.|$)")
    hits <- c(hits, grep(patt2, rn))
  }

  uniq_pos(hits)
}

get_any <- function(sym, ensg = NULL) {
  hit <- find_any(sym, ensg)
  if (!length(hit)) return(rep(0, ncol(M)))
  Matrix::colSums(M[hit, , drop = FALSE])
}

gt0 <- function(v) v > 0
ge2 <- function(v) v >= 2
ge3 <- function(v) v >= 3

v <- function(k) {
  switch(
    k,
    IL10   = get_any("IL10",   "ENSG00000136634"),
    CD274  = get_any("CD274",  "ENSG00000120217"),
    IL7R   = get_any("IL7R",   "ENSG00000140992"),
    ID2    = get_any("ID2",    "ENSG00000134644"),
    GATA3  = get_any("GATA3",  "ENSG00000107485"),
    RORC   = get_any("RORC",   "ENSG00000143365"),
    TBX21  = get_any("TBX21",  "ENSG00000073861"),
    KIT    = get_any("KIT",    "ENSG00000157404"),
    IL1RL1 = get_any("IL1RL1", "ENSG00000115602"),
    KLRG1  = get_any("KLRG1",  "ENSG00000139193"),
    HPGDS  = get_any("HPGDS",  "ENSG00000170209"),
    PTGDR2 = get_any("PTGDR2", "ENSG00000163659"),
    TRAC   = get_any("TRAC",   "ENSG00000277734"),
    TRBC1  = get_any("TRBC1",  "ENSG00000109541"),
    TRBC2  = get_any("TRBC2",  "ENSG00000211677"),
    TRDC   = get_any("TRDC",   "ENSG00000198835"),
    TRGC1  = get_any("TRGC1",  "ENSG00000211711"),
    TRGC2  = get_any("TRGC2",  "ENSG00000211680"),
    CD3D   = get_any("CD3D",   "ENSG00000167286"),
    EOMES  = get_any("EOMES",  "ENSG00000163508"),
    PRF1   = get_any("PRF1",   "ENSG00000180644"),
    NKG7   = get_any("NKG7",   "ENSG00000105374"),
    FCER1A = get_any("FCER1A", "ENSG00000179639"),
    IL3RA  = get_any("IL3RA",  "ENSG00000185291"),
    LST1   = get_any("LST1",   "ENSG00000154165"),
    FCGR3A = get_any("FCGR3A", "ENSG00000203747"),
    S100A8 = get_any("S100A8", "ENSG00000143556"),
    S100A9 = get_any("S100A9", "ENSG00000163220"),
    MS4A1  = get_any("MS4A1",  "ENSG00000156738"),
    CD79A  = get_any("CD79A",  "ENSG00000105369"),
    CD79B  = get_any("CD79B",  "ENSG00000007312"),
    TPSAB1 = get_any("TPSAB1", "ENSG00000172236"),
    CPA3   = get_any("CPA3",   "ENSG00000163751"),
    HDC    = get_any("HDC",    NULL)
  )
}

MS4A2_v  <- get_any("MS4A2",  NULL)
FCER1B_v <- get_any("FCER1B", NULL)

## ------------------------------------------------------------
## 3. ILC gating
## ------------------------------------------------------------
## lineage exclusion
is_T <- gt0(v("TRAC")) | gt0(v("TRBC1")) | gt0(v("TRBC2")) |
        gt0(v("TRDC")) | gt0(v("TRGC1")) | gt0(v("TRGC2")) |
        gt0(v("CD3D"))

is_NK <- gt0(v("EOMES")) | (gt0(v("PRF1")) & gt0(v("NKG7")))

is_pDC <- gt0(v("FCER1A")) & gt0(v("IL3RA"))

is_my_pair <- (gt0(v("LST1")) & gt0(v("FCGR3A"))) |
              (gt0(v("S100A8")) & gt0(v("S100A9")))
my_high3   <- ge3(v("LST1")) | ge3(v("FCGR3A")) |
              ge3(v("S100A8")) | ge3(v("S100A9"))
is_my      <- is_my_pair | my_high3

B_pair <- (gt0(v("CD79A")) & gt0(v("CD79B"))) |
          (gt0(v("MS4A1")) & (gt0(v("CD79A")) | gt0(v("CD79B"))))
B_ge2  <- ge2(v("MS4A1")) | ge2(v("CD79A")) | ge2(v("CD79B"))
is_B   <- B_pair | B_ge2

MAST_pair  <- gt0(v("TPSAB1")) & gt0(v("CPA3"))
MAST_high  <- ge2(v("TPSAB1")) | ge2(v("CPA3"))
MAST_extra <- ge2(v("HDC")) | ge2(MS4A2_v) | ge2(FCER1B_v)
is_Mast    <- MAST_pair | MAST_high | MAST_extra

## anchor
core_any <- gt0(v("GATA3")) | gt0(v("RORC")) |
            gt0(v("TBX21")) | gt0(v("KIT"))  |
            gt0(v("ID2"))

ilc2_any <- gt0(v("IL1RL1")) | gt0(v("KLRG1")) |
            gt0(v("HPGDS"))  | gt0(v("PTGDR2"))

il7r_any <- gt0(v("IL7R"))

id2_only <- gt0(v("ID2")) &
            !(gt0(v("GATA3")) | gt0(v("RORC")) |
              gt0(v("TBX21")) | gt0(v("KIT")) | ilc2_any)

anchor <- (core_any | ilc2_any) & !(id2_only & !il7r_any)

exclude  <- is_T | is_NK | is_pDC | is_my | is_B | is_Mast
ILC_flag <- anchor & !exclude

logi("ILC-flagged cells:", sum(ILC_flag))

## ------------------------------------------------------------
## 4. PD-L1 / IL-10 / DP flags
## ------------------------------------------------------------
IL10_ge1 <- v("IL10")  >= 1
PDL1_ge1 <- v("CD274") >= 1
DP       <- IL10_ge1 & PDL1_ge1

## ------------------------------------------------------------
## 5. ILC subtype annotation
##    Marker set used in revision:
##      ILC2: GATA3, IL1RL1(ST2), KLRG1, HPGDS, PTGDR2, IL7R
##      ILC3: RORC
##      ILC1: TBX21
##
##    Mutually exclusive rule:
##      - ILC3_like: RORC+
##      - ILC1_like: TBX21+ and not ILC2-marker+
##      - ILC2_like: any ILC2 marker+, not RORC+
##      - else: other_or_mixed_ILC
## ------------------------------------------------------------
ILC2_markers_pos <- gt0(v("GATA3"))  |
                    gt0(v("IL1RL1")) |
                    gt0(v("KLRG1"))  |
                    gt0(v("HPGDS"))  |
                    gt0(v("PTGDR2")) |
                    gt0(v("IL7R"))

ILC3_markers_pos <- gt0(v("RORC"))
ILC1_markers_pos <- gt0(v("TBX21"))

subtype <- dplyr::case_when(
  !ILC_flag ~ NA_character_,
  ILC3_markers_pos ~ "ILC3_like",
  !ILC3_markers_pos & ILC1_markers_pos & !ILC2_markers_pos ~ "ILC1_like",
  !ILC3_markers_pos & ILC2_markers_pos ~ "ILC2_like",
  TRUE ~ "other_or_mixed_ILC"
)

cell_df <- tibble(
  cell_id    = colnames(M),
  sample_id  = sample_id,
  site       = site,
  visit      = visit,
  group4     = group4,
  ILC        = as.integer(ILC_flag),
  IL10_ge1   = as.integer(IL10_ge1),
  PDL1_ge1   = as.integer(PDL1_ge1),
  DP         = as.integer(DP),
  subtype    = subtype,
  GATA3      = v("GATA3"),
  IL1RL1     = v("IL1RL1"),
  KLRG1      = v("KLRG1"),
  HPGDS      = v("HPGDS"),
  PTGDR2     = v("PTGDR2"),
  IL7R       = v("IL7R"),
  RORC       = v("RORC"),
  TBX21      = v("TBX21")
)

## keep only ILCs for subtype summaries
ilc_df <- cell_df %>%
  dplyr::filter(ILC == 1)

logi("Subtype counts among all ILC:")
print(table(ilc_df$subtype, useNA = "ifany"))

## ------------------------------------------------------------
## 6. Helper to summarize subtype composition in a given set
## ------------------------------------------------------------
summarise_subtype_set <- function(df, label) {
  n_total <- nrow(df)

  n_ilc2  <- sum(df$subtype == "ILC2_like", na.rm = TRUE)
  n_ilc3  <- sum(df$subtype == "ILC3_like", na.rm = TRUE)
  n_ilc1  <- sum(df$subtype == "ILC1_like", na.rm = TRUE)
  n_other <- sum(df$subtype == "other_or_mixed_ILC", na.rm = TRUE)

  tibble(
    set        = label,
    n_cells    = n_total,
    n_ILC2     = n_ilc2,
    n_ILC3     = n_ilc3,
    n_ILC1     = n_ilc1,
    n_other    = n_other,
    pct_ILC2   = ifelse(n_total > 0, 100 * n_ilc2  / n_total, NA_real_),
    pct_ILC3   = ifelse(n_total > 0, 100 * n_ilc3  / n_total, NA_real_),
    pct_ILC1   = ifelse(n_total > 0, 100 * n_ilc1  / n_total, NA_real_),
    pct_other  = ifelse(n_total > 0, 100 * n_other / n_total, NA_real_)
  )
}

## manuscript-facing summary table
summary_main <- dplyr::bind_rows(
  summarise_subtype_set(ilc_df, "All_ILC"),
  summarise_subtype_set(ilc_df %>% dplyr::filter(PDL1_ge1 == 1), "PDL1pos_ILC"),
  summarise_subtype_set(ilc_df %>% dplyr::filter(DP == 1), "DP_ILC")
)

print(summary_main)

## ------------------------------------------------------------
## 7. Optional: group4-level subtype composition
## ------------------------------------------------------------
summary_by_group4 <- ilc_df %>%
  dplyr::mutate(
    set = dplyr::case_when(
      DP == 1 ~ "DP_ILC",
      PDL1_ge1 == 1 ~ "PDL1pos_ILC",
      TRUE ~ "All_ILC"
    )
  ) %>%
  dplyr::group_by(group4, set, subtype) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(total = sum(n), pct = 100 * n / total) %>%
  dplyr::ungroup()

## per-sample version if needed
summary_by_sample <- ilc_df %>%
  dplyr::mutate(
    set = dplyr::case_when(
      DP == 1 ~ "DP_ILC",
      PDL1_ge1 == 1 ~ "PDL1pos_ILC",
      TRUE ~ "All_ILC"
    )
  ) %>%
  dplyr::group_by(sample_id, group4, set, subtype) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(total = sum(n), pct = 100 * n / total) %>%
  dplyr::ungroup()

## ------------------------------------------------------------
## 8. Save outputs
## ------------------------------------------------------------
readr::write_csv(
  summary_main,
  file.path(OUTDIR, "GSE228421_ILC2_ILC3_ILC1_in_PDL1_and_DP.csv")
)

readr::write_csv(
  summary_by_group4,
  file.path(OUTDIR, "GSE228421_ILC_subtype_by_group4_long.csv")
)

readr::write_csv(
  summary_by_sample,
  file.path(OUTDIR, "GSE228421_ILC_subtype_by_sample_long.csv")
)

readr::write_csv(
  ilc_df,
  file.path(OUTDIR, "GSE228421_ILC_per_cell_with_subtype.csv")
)

logi("Saved subtype annotation outputs:")
logi(" - GSE228421_ILC2_ILC3_ILC1_in_PDL1_and_DP.csv")
logi(" - GSE228421_ILC_subtype_by_group4_long.csv")
logi(" - GSE228421_ILC_subtype_by_sample_long.csv")
logi(" - GSE228421_ILC_per_cell_with_subtype.csv")
