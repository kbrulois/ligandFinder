#!/usr/bin/env Rscript
## ---- the FULL step-1 scan as a CSV --------------------------------------------
## 10_9b's table has 41,870 rows because it is the peptide-end set's CANDIDATE
## windows -- dibasic d-1, amidation d-2, precursor termini, knowns. That split is
## confusingly named `all` in the data, but it is not the scan.
##
## This writes the actual exhaustive scan: every mature residue of every
## precursor, step 1, 2,980,852 rows for C over 5,616 precursors.
##
## Deliberately LEAN. At the candidate table's 45-column width the same rows come
## to ~965 MB, and most of that is wasted: win_type, motif, pep_name and the db_*
## block exist only for candidate anchors, so ~2.94M rows would carry NA there.
## Kept here:
##   * the window and its scores, with the per-seed summaries
##   * is_peak / rank_peak, so the de-duplicated set is `is_peak` away
##   * the EEC expression block, which is per GENE and so applies to every row
##   * candidate metadata and the db model's scores where the position IS a
##     candidate, NA elsewhere -- that is information, not a gap
##
## A .csv.gz is written alongside; it is ~5x smaller and fread/pandas read it
## directly.
##
##   Rscript inst/scripts/10_9d_pepend_scan_csv.R --term C
##   ... --min-score 0.05        # drop the background, far smaller file
##   ... --peaks-only            # just the local maxima
##
## Output: ~/AF2_analysis/lf_pepend_scan_windows_<term>.csv(.gz)
## ------------------------------------------------------------------------------

suppressMessages({ library(data.table); library(dplyr) })
setDTthreads(0L)

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
.flag <- function(flag) flag %in% .args

term     <- toupper(.opt("--term", "C"))
scan_p   <- path.expand(.opt("--scan", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
cand_csv <- path.expand(.opt("--candidates", sprintf("~/AF2_analysis/lf_pepend_all_windows_%s_eec.csv", term)))
eec_p    <- path.expand(.opt("--eec", "~/AF2_analysis/eec_gene_classification.tsv"))
cache_p  <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
pa_p     <- path.expand(.opt("--proteinatlas", "~/AF2_analysis/hpa_cache/proteinatlas.tsv"))
out_csv  <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_scan_windows_%s.csv", term)))
min_sc   <- as.numeric(.opt("--min-score", "-1"))
no_gz    <- .flag("--no-gzip")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "eec_expression.R"))
source(file.path(ROOT, "R", "pepend_windows.R"))

## ---- 1. the scan --------------------------------------------------------------
s  <- readRDS(scan_p)
dt <- as.data.table(s$scan)
message(sprintf("scan: %s windows over %d precursors; columns %s",
                format(nrow(dt), big.mark = ","), uniqueN(dt$accession),
                paste(names(dt), collapse = ", ")))

## the local maxima, as a flag rather than a separate file
pk <- as.data.table(s$peaks)[, .(accession, anchor, rank_peak)]
dt <- merge(dt, pk, by = c("accession", "anchor"), all.x = TRUE, sort = FALSE)
dt[, is_peak := !is.na(rank_peak)]
message(sprintf("  local peaks (+/-%s): %s", s$peak_w,
                format(sum(dt$is_peak), big.mark = ",")))

if (.flag("--peaks-only")) {
  dt <- dt[is_peak == TRUE]
  out_csv <- sub("\\.csv$", "_peaks.csv", out_csv)
  message("  --peaks-only: ", format(nrow(dt), big.mark = ","), " rows kept")
}
if (min_sc > 0) {
  n0 <- nrow(dt); dt <- dt[score >= min_sc]
  message(sprintf("  --min-score %.3g: %s of %s rows kept", min_sc,
                  format(nrow(dt), big.mark = ","), format(n0, big.mark = ",")))
}

## ---- 2. candidate metadata + the db model, where the position is a candidate ---
if (file.exists(cand_csv)) {
  keep <- c("accession", "anchor", "peps", "win_type", "motif", "at_candidate",
            "known", "pep_name", "pep_start", "pep_end", "len", "db_pos",
            grep("^db_(mean|n_seeds|sd)", names(fread(cand_csv, nrows = 0)), value = TRUE))
  cand <- fread(cand_csv, select = keep)
  setnames(cand, "peps", "candidate_peps")
  ## the candidate table's `motif` exists for the knowns only; keep it aside to
  ## check the one computed for every window below against it
  setnames(cand, "motif", "motif_candidate")
  before <- nrow(dt)
  dt <- merge(dt, cand, by = c("accession", "anchor"), all.x = TRUE, sort = FALSE)
  stopifnot(nrow(dt) == before)
  dt[, is_candidate := !is.na(win_type)]
  message(sprintf("\ncandidate metadata on %s of %s rows (%.2f%%); db scores on %s",
                  format(sum(dt$is_candidate), big.mark = ","),
                  format(nrow(dt), big.mark = ","), 100 * mean(dt$is_candidate),
                  format(sum(!is.na(dt$db_mean_allseeds)), big.mark = ",")))
} else {
  message("\n", cand_csv, " not found -- candidate/db columns omitted")
  dt[, is_candidate := FALSE]
}

## ---- 2b. motif: what sits just outside EVERY anchor ---------------------------
## The candidate table carries `motif` for the knowns only. Here it is computed
## for every scanned position, with lf_pepend_motif_vec -- validated against the
## scalar lf_pepend_motif over 157,833 anchors per terminus, and re-checked below
## against whatever the candidate table already had.
fc    <- readRDS(cache_p)
seqs  <- setNames(fc$prec$seq,    fc$prec$accession)
nprot <- setNames(fc$prec$n_prot, fc$prec$accession)
miss  <- setdiff(unique(dt$accession), names(seqs))
if (length(miss)) stop(length(miss), " accession(s) absent from ", cache_p, call. = FALSE)

dt[, motif := lf_pepend_motif_vec(seqs[[.BY$accession]], term, anchor,
                                  nprot[[.BY$accession]]), by = accession]
## the four classes asked for: the amidation variant is a dibasic pair with a
## glycine in front, so it folds into `dibasic` rather than becoming a fifth
dt[, motif4 := fifelse(motif == "G + dibasic", "dibasic", motif)]

if ("motif_candidate" %in% names(dt)) {
  ## fread gives an absent character field "" rather than NA, and the candidate
  ## table only has a motif for the knowns -- test for both or the comparison
  ## silently covers every candidate row instead of the 86 that carry one
  chk <- dt[!is.na(motif_candidate) & nzchar(motif_candidate)]
  if (nrow(chk) && !all(chk$motif == chk$motif_candidate))
    stop(sprintf("computed motif disagrees with the candidate table on %d of %d rows",
                 sum(chk$motif != chk$motif_candidate), nrow(chk)), call. = FALSE)
  message(sprintf("\nmotif: computed for all %s rows; matches the candidate table on the %s that had one",
                  format(nrow(dt), big.mark = ","), format(nrow(chk), big.mark = ",")))
  dt[, motif_candidate := NULL]
}
message("  motif (5 classes, as elsewhere in the codebase):")
print(as.data.frame(dt[, .(rows = .N, peaks = sum(is_peak),
                           median_score = round(median(score), 4),
                           max_score = round(max(score), 4)), by = motif][order(-rows)]),
      row.names = FALSE)
message("  motif4 (the four asked for):")
print(as.data.frame(dt[, .(rows = .N, pct = round(100 * .N / nrow(dt), 2),
                           peaks = sum(is_peak),
                           median_score = round(median(score), 4)), by = motif4][order(-rows)]),
      row.names = FALSE)

## ---- 3. EEC expression: per GENE, so it lands on every row ---------------------
## Join on GENE SYMBOL -- the HPA table has no UniProt accession or entry name,
## and id_mapping.rds has no Ensembl, so the symbol is the only shared key
## (see R/eec_expression.R).
eec <- fread(eec_p)
dt  <- as.data.table(lf_eec_attach(dt, eec, pa_p))
message(sprintf("EEC attached: %s rows with a class other than no_expression_data",
                format(sum(dt$eec_class != "no_expression_data"), big.mark = ",")))

## ---- 4. write -----------------------------------------------------------------
front <- intersect(c("gene", "accession", "anchor", "score", "sd",
                     grep("^score_top", names(dt), value = TRUE),
                     grep("^n_seeds_gt_", names(dt), value = TRUE),
                     "rank_all", "is_peak", "rank_peak",
                     "motif4", "motif", "is_candidate", "candidate_peps", "win_type",
                     "known", "pep_name"), names(dt))
setcolorder(dt, c(front, setdiff(names(dt), front)))
setorderv(dt, "rank_all")

fwrite(dt, out_csv)
sz <- file.info(out_csv)$size
message(sprintf("\n%s rows x %d columns -> %s  (%.0f MB)",
                format(nrow(dt), big.mark = ","), ncol(dt), out_csv, sz / 1e6))
if (!no_gz) {
  gz <- paste0(out_csv, ".gz")
  fwrite(dt, gz, compress = "gzip")
  message(sprintf("  gzipped -> %s  (%.0f MB, %.1fx smaller)", gz,
                  file.info(gz)$size / 1e6, sz / file.info(gz)$size))
}

## ---- 5. what is in it ---------------------------------------------------------
ngt <- grep("^n_seeds_gt_", names(dt), value = TRUE)[1]
message("\nrows by EEC class:")
dt[, .(rows = .N, genes = uniqueN(gene), peaks = sum(is_peak),
       median_score = round(median(score), 4), max_score = round(max(score), 4),
       all_seeds_over = sum(get(ngt) == max(get(ngt), na.rm = TRUE), na.rm = TRUE)),
   by = eec_class][order(-rows)] %>% as.data.frame() %>% print(row.names = FALSE)

message("\ntop 20 local peaks in EEC-specific genes:")
dt[eec_class == "specifically_expressed" & is_peak][order(-score)][1:20,
   .(gene, anchor, score = round(score, 4),
     top_k = round(get(grep("^score_top", names(dt), value = TRUE)[1]), 4),
     n_gt = get(ngt), rank_peak, win_type,
     eec_ntpm = round(eec_ntpm, 1), eec_tier,
     known = ifelse(is.na(pep_name), "", pep_name))] %>%
  as.data.frame() %>% print(row.names = FALSE)
