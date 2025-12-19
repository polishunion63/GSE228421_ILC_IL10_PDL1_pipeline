############################################################
## GSE228421 — Human psoriasis skin scRNA-seq
## End-to-end pipeline (submission version)
##  1) raw 10x .tsv/.mtx(.gz) → cell-filtered merged counts (dgCMatrix)
##  2) ILC gating + PD-L1/IL-10 markers, nILC >= 100 QC
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

## ------------------------------
## Global parameters
## ------------------------------
MIN_GENES_PER_BC   <- 200   # cell-calling: minimum detected genes per barcode
MIN_UMI_PER_BC     <- 500   # cell-calling: minimum UMI per barcode
MIN_ILC_PER_SAMPLE <- 100   # minimum ILC count to include a sample in summaries

############################################################
## Step 1 — raw 10x .tsv/.mtx(.gz) → filtered merged counts
############################################################

read_tsv_maybe_gz <- function(path) {
  if (is.na(path) || !file.exists(path)) stop("Missing TSV: ", path)
  x <- try(
    readr::read_tsv(path, col_names = FALSE,
                    show_col_types = FALSE, progress = FALSE),
    silent = TRUE
  )
  if (!inherits(x, "try-error")) return(x)
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  readr::read_tsv(con, col_names = FALSE,
                  show_col_types = FALSE, progress = FALSE)
}

readMM_maybe_gz <- function(path) {
  if (is.na(path) || !file.exists(path)) stop("Missing MTX: ", path)
  x <- try(Matrix::readMM(path), silent = TRUE)
  if (!inherits(x, "try-error")) return(x)
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  Matrix::readMM(con)
}

## Example file names:
##   GSM7120449_P1_V1_L.barcodes.tsv
##   GSM7120449_P1_V1_L_barcodes.tsv
##   GSM7120449_P1_V1_L.matrix.mtx.gz
## Stub becomes: <dir>/GSM7120449_P1_V1_L
stub_from_path <- function(p) {
  d <- dirname(p)
  b <- basename(p)
  b0 <- sub("([._](barcodes|features|genes|matrix))\\.(tsv|mtx)(\\.gz)?$",
            "", b, perl = TRUE, ignore.case = TRUE)
  file.path(d, b0)
}

resolve_triplets <- function(paths) {
  stubs <- unique(vapply(paths, stub_from_path, FUN.VALUE = character(1)))

  pick_first <- function(stub, kind) {
    d    <- dirname(stub)
    base <- basename(stub)
    cands <- c(
      file.path(d, paste0(base, "_", kind, ".tsv.gz")),
      file.path(d, paste0(base, "_", kind, ".tsv")),
      file.path(d, paste0(base, ".",  kind, ".tsv.gz")),
      file.path(d, paste0(base, ".",  kind, ".tsv")),
      file.path(d, paste0(base, "_", kind, ".mtx.gz")),
      file.path(d, paste0(base, "_", kind, ".mtx")),
      file.path(d, paste0(base, ".",  kind, ".mtx.gz")),
      file.path(d, paste0(base, ".",  kind, ".mtx"))
    )
    existing <- cands[file.exists(cands)]
    if (length(existing) == 0) return(NA_character_)
    existing[1]
  }

  rows <- lapply(stubs, function(s) {
    list(
      stub     = s,
      barcodes = pick_first(s, "barcodes"),
      features = {
        x <- pick_first(s, "features")
        if (is.na(x)) pick_first(s, "genes") else x
      },
      matrix   = pick_first(s, "matrix")
    )
  })

  dplyr::bind_rows(rows)
}

## Load 10x triplet for a single stub and apply simple cell-calling
read_filter_10x_from_row <- function(row) {
  logi("Loading 10x sample:", basename(row$stub))

  bc   <- read_tsv_maybe_gz(row$barcodes)[[1]]
  feat <- read_tsv_maybe_gz(row$features)
  gene_ids <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]

  mm <- readMM_maybe_gz(row$matrix)
  dimnames(mm) <- list(make.unique(gene_ids), make.unique(bc))

  genes_per_bc <- Matrix::colSums(mm > 0)
  umi_per_bc   <- Matrix::colSums(mm)

  keep <- (genes_per_bc >= MIN_GENES_PER_BC) &
          (umi_per_bc   >= MIN_UMI_PER_BC)

  n_total <- ncol(mm)
  n_keep  <- sum(keep)

  logi(sprintf("  total barcodes = %d, kept (genes >= %d & UMI >= %d) = %d (%.3f%%)",
               n_total, MIN_GENES_PER_BC, MIN_UMI_PER_BC,
               n_keep, ifelse(n_total > 0, 100 * n_keep / n_total, 0)))

  if (n_keep == 0L) {
    warning("No barcodes passed filtering in sample: ", basename(row$stub))
    mm_f <- mm[, 0, drop = FALSE]
  } else {
    mm_f <- mm[, keep, drop = FALSE]
  }

  token <- sub("^GSM\\d+[_-]", "", basename(row$stub), perl = TRUE)
  colnames(mm_f) <- paste0(token, "_", colnames(mm_f))

  list(
    M_filt          = as(mm_f, "dgCMatrix"),
    genes           = rownames(mm_f),
    sample_id       = token,
    n_genes_all     = nrow(mm),
    n_bc_all        = n_total,
    n_bc_keep       = n_keep,
    tot_UMI_all     = sum(mm),
    tot_UMI_keep    = sum(mm_f),
    mean_genes_keep = if (n_keep > 0) mean(genes_per_bc[keep]) else NA_real_,
    mean_UMI_keep   = if (n_keep > 0) mean(umi_per_bc[keep])   else NA_real_
  )
}

pad_to_union <- function(mats, genes_list) {
  allg <- unique(unlist(genes_list))
  logi("Union genes (after filtering):", length(allg))
  out <- vector("list", length(mats))
  pb2 <- utils::txtProgressBar(min = 0, max = length(mats), style = 3)
  for (i in seq_along(mats)) {
    mm <- mats[[i]]
    if (ncol(mm) == 0L) {
      pad <- Matrix::Matrix(0, nrow = length(allg), ncol = 0, sparse = TRUE)
      rownames(pad) <- allg
    } else {
      idx <- match(rownames(mm), allg)
      pad <- Matrix::Matrix(0, nrow = length(allg), ncol = ncol(mm), sparse = TRUE)
      pad[idx, ] <- mm
      rownames(pad) <- allg
      colnames(pad) <- colnames(mm)
    }
    out[[i]] <- pad
    utils::setTxtProgressBar(pb2, i)
  }
  close(pb2)
  do.call(cbind, out)
}

pick_files <- function() {
  if (tolower(Sys.info()[["sysname"]]) == "windows") {
    utils::choose.files(
      caption = "Select GSE228421 barcodes/features/matrix files"
    )
  } else {
    message("Enter space-separated paths:")
    scan(what = "character")
  }
}

## Step 1A — collect triplets
paths <- pick_files()
paths <- paths[nzchar(paths)]
stopifnot(length(paths) >= 3)

logi("Number of selected files:", length(paths))

logi("Resolving barcodes/features/matrix triplets ...")
trip <- resolve_triplets(paths)

miss <- trip %>%
  dplyr::mutate(
    miss_barcodes = is.na(barcodes),
    miss_features = is.na(features),
    miss_matrix   = is.na(matrix)
  ) %>%
  dplyr::filter(miss_barcodes | miss_features | miss_matrix)

if (nrow(miss) > 0) {
  print(
    miss %>%
      dplyr::transmute(
        stub     = basename(stub),
        barcodes = basename(barcodes),
        features = basename(features),
        matrix   = basename(matrix)
      )
  )
  stop("Some samples are missing barcodes/features/matrix files.")
}

logi("Number of complete samples:", nrow(trip))
print(
  trip %>%
    dplyr::transmute(
      stub     = basename(stub),
      barcodes = basename(barcodes),
      features = basename(features),
      matrix   = basename(matrix)
    )
)

## Step 1B — per-sample loading and filtering
logi("Loading and filtering 10x samples ...")

mats_filt  <- list()
genes_list <- list()
qc_list    <- list()

pb <- utils::txtProgressBar(min = 0, max = nrow(trip), style = 3)
for (i in seq_len(nrow(trip))) {
  r <- read_filter_10x_from_row(trip[i, ])

  mats_filt[[i]]  <- r$M_filt
  genes_list[[i]] <- r$genes

  qc_list[[i]] <- tibble::tibble(
    sample_token      = r$sample_id,
    stub              = basename(trip$stub[i]),
    n_genes_all       = r$n_genes_all,
    n_bc_all          = r$n_bc_all,
    n_bc_keep         = r$n_bc_keep,
    total_UMI_all     = r$tot_UMI_all,
    total_UMI_keep    = r$tot_UMI_keep,
    mean_genes_keep   = r$mean_genes_keep,
    mean_UMI_keep     = r$mean_UMI_keep
  )

  utils::setTxtProgressBar(pb, i)
}
close(pb)

sample_qc <- dplyr::bind_rows(qc_list) %>%
  dplyr::arrange(sample_token)

logi("Per-sample pre/post filtering QC:")
print(sample_qc)

## Step 1C — union padding and merge
logi("Union-padded merge of filtered matrices ...")
M_filt <- pad_to_union(mats_filt, genes_list)
stopifnot(inherits(M_filt, "dgCMatrix"))

logi(sprintf("Merged filtered counts: %d genes × %d cells",
             nrow(M_filt), ncol(M_filt)))

## Step 1D — save outputs
rds_path <- file.path(OUTDIR, "GSE228421_merged_counts_filtered.rds")
saveRDS(M_filt, rds_path)
logi("Saved filtered count matrix RDS to:", rds_path)

qc_path <- file.path(OUTDIR, "GSE228421_sample_qc_pre_post.csv")
readr::write_csv(sample_qc, qc_path)
logi("Saved per-sample QC table to:", qc_path)

## Use M for downstream
M <- M_filt

############################################################
## Step 2 — ILC gating + PD-L1/IL-10 summary
############################################################

logi("Starting ILC gating and summary ...")

## 2A — basic metadata from column names
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

## 2B — gene helpers
uniq_pos  <- function(x) unique(x[!is.na(x) & x > 0])

find_any <- function(symbol, ensg = NULL) {
  rn <- rownames(M)
  hits <- integer(0)
  hits <- c(hits, which(rn == symbol))
  patt <- paste0("(^|[|;,_\\s\\.])", symbol, "($|[|;,_\\s\\.])")
  hits <- c(hits, grep(patt, rn))
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

ci_pct_safe <- function(x, n, which = 1) {
  if (is.na(n) || is.na(x) || n < 1) return(NA_real_)
  if (x > n) {
    warning(sprintf("x > n in binom.test (x = %s, n = %s); clamped", x, n))
    x <- n
  }
  100 * stats::binom.test(as.integer(x), as.integer(n))$conf.int[which]
}

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

## 2C — lineage exclusion
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

## 2D — ILC anchor
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

## 2E — marker flags
IL10_ge1 <- v("IL10")  >= 1
PDL1_ge1 <- v("CD274") >= 1
DP_lite  <- IL10_ge1 & PDL1_ge1

## 2F — cell-level table
cell_df <- tibble::tibble(
  sample_id = sample_id,
  site      = site,
  visit     = visit,
  group4    = group4,
  ILC       = as.integer(ILC_flag),
  IL10_ge1  = as.integer(IL10_ge1),
  PDL1_ge1  = as.integer(PDL1_ge1),
  DP_lite   = as.integer(DP_lite)
)

## 2G — by-sample counts (all samples)
by_sample_core_all <- cell_df %>%
  dplyr::group_by(sample_id, site, visit, group4) %>%
  dplyr::summarise(
    nILC        = sum(ILC),
    n_IL10_ge1  = sum(IL10_ge1[ILC == 1]),
    n_PDL1_ge1  = sum(PDL1_ge1[ILC == 1]),
    n_DP_lite   = sum(DP_lite[ILC == 1]),
    .groups = "drop"
  ) %>%
  dplyr::arrange(group4, sample_id)

## 2H — ILC-based sample QC (nILC >= 100)
ILC_QC <- by_sample_core_all %>%
  dplyr::mutate(
    pass_nILC100 = nILC >= MIN_ILC_PER_SAMPLE
  ) %>%
  dplyr::arrange(group4, sample_id)

## 2I — restrict to samples with nILC >= 100 for percentages
by_sample_core <- ILC_QC %>%
  dplyr::filter(pass_nILC100) %>%
  dplyr::select(-pass_nILC100)

by_sample_pct <- by_sample_core %>%
  dplyr::mutate(
    IL10_ge1_pct = 100 * n_IL10_ge1 / nILC,
    PDL1_ge1_pct = 100 * n_PDL1_ge1 / nILC,
    DP_lite_pct  = 100 * n_DP_lite   / nILC
  ) %>%
  dplyr::arrange(group4, sample_id)

by_sample_presence <- by_sample_core %>%
  dplyr::transmute(
    sample_id, site, visit, group4, nILC,
    has_IL10_ge1 = n_IL10_ge1 > 0,
    has_PDL1_ge1 = n_PDL1_ge1 > 0,
    has_DP_lite  = n_DP_lite  > 0
  ) %>%
  dplyr::arrange(group4, sample_id)

## 2J — pooled summaries (group4)
pooled <- by_sample_core %>%
  dplyr::group_by(group4) %>%
  dplyr::summarise(
    n_samples_included = dplyr::n(),
    total_ILC          = sum(nILC),
    total_IL10_ge1     = sum(n_IL10_ge1),
    total_PDL1_ge1     = sum(n_PDL1_ge1),
    total_DP_lite      = sum(n_DP_lite),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    IL10_ge1_pct = 100 * total_IL10_ge1 / total_ILC,
    IL10_ge1_lo  = mapply(ci_pct_safe, total_IL10_ge1, total_ILC,
                          MoreArgs = list(which = 1)),
    IL10_ge1_hi  = mapply(ci_pct_safe, total_IL10_ge1, total_ILC,
                          MoreArgs = list(which = 2)),
    PDL1_ge1_pct = 100 * total_PDL1_ge1 / total_ILC,
    PDL1_ge1_lo  = mapply(ci_pct_safe, total_PDL1_ge1, total_ILC,
                          MoreArgs = list(which = 1)),
    PDL1_ge1_hi  = mapply(ci_pct_safe, total_PDL1_ge1, total_ILC,
                          MoreArgs = list(which = 2)),
    DP_lite_pct  = 100 * total_DP_lite  / total_ILC,
    DP_lite_lo   = mapply(ci_pct_safe, total_DP_lite,  total_ILC,
                          MoreArgs = list(which = 1)),
    DP_lite_hi   = mapply(ci_pct_safe, total_DP_lite,  total_ILC,
                          MoreArgs = list(which = 2))
  ) %>%
  dplyr::arrange(group4)

## 2K — gating audit (for documentation)
kept <- ILC_flag
pct  <- function(n, d) ifelse(d > 0, 100 * n / d, NA_real_)

mast_single1 <- (gt0(v("TPSAB1")) | gt0(v("CPA3"))) &
                !(MAST_pair | MAST_high | MAST_extra)
pdc_single   <- (gt0(v("FCER1A")) | gt0(v("IL3RA"))) & !is_pDC
my_single_hi <- (!is_my_pair) & my_high3

gating_audit <- tibble::tibble(
  kept_ILC = sum(kept),
  drop_T   = sum(anchor & is_T),
  drop_NK  = sum(anchor & is_NK),
  drop_pDC = sum(anchor & is_pDC),
  drop_my_pair  = sum(anchor & is_my_pair & !my_high3),
  drop_my_high3 = sum(anchor & my_high3),
  drop_B   = sum(anchor & is_B),
  drop_Mast_pair  = sum(anchor & MAST_pair),
  drop_Mast_high  = sum(anchor & MAST_high),
  drop_Mast_extra = sum(anchor & MAST_extra),
  pct_ILC_mast_single1 = pct(sum(kept & mast_single1), sum(kept)),
  pct_ILC_pdc_single   = pct(sum(kept & pdc_single),   sum(kept)),
  pct_ILC_my_single_hi = pct(sum(kept & my_single_hi), sum(kept))
)

## 2L — write outputs
readr::write_csv(by_sample_core_all,
                 file.path(OUTDIR, "GSE228421_by_sample_core_allSamples.csv"))
readr::write_csv(ILC_QC,
                 file.path(OUTDIR, "GSE228421_ILC_QC_nILC_ge100.csv"))
readr::write_csv(by_sample_core,
                 file.path(OUTDIR, "GSE228421_by_sample_core_nILC_ge100.csv"))
readr::write_csv(by_sample_pct,
                 file.path(OUTDIR, "GSE228421_by_sample_pct_nILC_ge100.csv"))
readr::write_csv(by_sample_presence,
                 file.path(OUTDIR, "GSE228421_by_sample_presence_nILC_ge100.csv"))
readr::write_csv(pooled,
                 file.path(OUTDIR, "GSE228421_pooled_summary_nILC_ge100.csv"))
readr::write_csv(gating_audit,
                 file.path(OUTDIR, "GSE228421_gating_audit_ILC.csv"))

logi("Finished: raw 10x -> filtered counts -> single ILC gate with nILC >= 100 QC.")
