############################################################
## 03_GSE228421_gdT_exploratory_gating.R
## Purpose:
##   - Load merged filtered count matrix (RDS)
##   - Apply a strict γδ T-cell gating strategy
##   - Save per-cell γδ T flags
##   - Save per-sample summary with group labels
##
## Strict γδ T definition used here:
##   TRDC > 0 AND TRGC1/2 > 0 AND TRAC == 0
##
## Note:
##   This is a highly conservative exploratory gate and may
##   underestimate γδ T-cell abundance in sparse scRNA-seq data.
############################################################

rm(list = ls()); invisible(gc())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
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
## 2. Helper functions: find genes and extract expression
## ------------------------------------------------------------
find_gene <- function(mat, candidates) {
  hit <- intersect(candidates, rownames(mat))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_expr <- function(mat, gene_symbol) {
  if (is.na(gene_symbol)) return(rep(0, ncol(mat)))
  mat[gene_symbol, , drop = TRUE]
}

## ------------------------------------------------------------
## 3. γδ T-related genes
## ------------------------------------------------------------
trdc_sym  <- find_gene(M, c("TRDC",  "Trdc"))
trgc1_sym <- find_gene(M, c("TRGC1", "Trgc1"))
trgc2_sym <- find_gene(M, c("TRGC2", "Trgc2"))
trac_sym  <- find_gene(M, c("TRAC",  "Trac"))
trbc1_sym <- find_gene(M, c("TRBC1", "Trbc1"))
trbc2_sym <- find_gene(M, c("TRBC2", "Trbc2"))

cat("TRDC :",  trdc_sym,  "\n")
cat("TRGC1:",  trgc1_sym, "\n")
cat("TRGC2:",  trgc2_sym, "\n")
cat("TRAC :",  trac_sym,  "\n")
cat("TRBC1:",  trbc1_sym, "\n")
cat("TRBC2:",  trbc2_sym, "\n")

## ------------------------------------------------------------
## 4. Extract expression vectors
## ------------------------------------------------------------
expr_trdc  <- get_expr(M, trdc_sym)
expr_trgc1 <- get_expr(M, trgc1_sym)
expr_trgc2 <- get_expr(M, trgc2_sym)
expr_trac  <- get_expr(M, trac_sym)
expr_trbc1 <- get_expr(M, trbc1_sym)
expr_trbc2 <- get_expr(M, trbc2_sym)

## gamma-chain combined, beta-chain combined
expr_trgc <- pmax(expr_trgc1, expr_trgc2)
expr_trbc <- pmax(expr_trbc1, expr_trbc2)

## ------------------------------------------------------------
## 5. Strict γδ T-cell gating
## ------------------------------------------------------------
## γδ T candidate: TRDC+ and (TRGC1 or TRGC2)+
is_gd_candidate <- (expr_trdc > 0) & (expr_trgc > 0)

## αβ-like exclusion: TRAC > 0
is_ab_like <- expr_trac > 0

## final strict γδ T-cell gate
is_gdT <- is_gd_candidate & !is_ab_like

cat("Total cells:", ncol(M), "\n")
cat("γδ T candidates:", sum(is_gd_candidate), "\n")
cat("Final strict γδ T cells:", sum(is_gdT), "\n")

## ------------------------------------------------------------
## 6. γδ T subset object
## ------------------------------------------------------------
M_gdT <- M[, is_gdT, drop = FALSE]
gd_cell_ids <- colnames(M)[is_gdT]

logi("Number of strict γδ T cells:", length(gd_cell_ids))
logi("First few γδ T cell IDs:")
print(utils::head(gd_cell_ids))

## ------------------------------------------------------------
## 7. Per-cell flag table
## ------------------------------------------------------------
gdT_flag_per_cell <- tibble(
  cell_id   = colnames(M),
  sample_id = sample_id,
  site      = site,
  visit     = visit,
  group4    = group4,
  is_gdT    = as.integer(is_gdT)
)

## ------------------------------------------------------------
## 8. Per-sample summary with group labels
## ------------------------------------------------------------
gdT_by_sample_with_group <- gdT_flag_per_cell %>%
  dplyr::group_by(sample_id, site, visit, group4) %>%
  dplyr::summarise(
    n_cells        = dplyr::n(),
    n_gdT          = sum(is_gdT),
    pooled_pct_gdT = 100 * n_gdT / n_cells,
    .groups = "drop"
  ) %>%
  dplyr::arrange(group4, sample_id)

print(gdT_by_sample_with_group)

## ------------------------------------------------------------
## 9. Save outputs
## ------------------------------------------------------------
readr::write_csv(
  gdT_flag_per_cell,
  file.path(OUTDIR, "gdT_flag_per_cell.csv")
)

readr::write_csv(
  gdT_by_sample_with_group,
  file.path(OUTDIR, "gdT_by_sample_with_group.csv")
)

logi("Saved outputs:")
logi(" - gdT_flag_per_cell.csv")
logi(" - gdT_by_sample_with_group.csv")