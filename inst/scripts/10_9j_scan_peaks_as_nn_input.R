#!/usr/bin/env Rscript
## ---- scan peaks -> an nn_input_comb the protein page can render ----------------
## make_protein_plot_win draws every window as a segment on a stacked layer, so
## it cannot take the step-1 scan whole: NPY's 69 overlapping windows would stack
## 69 layers. This writes the LOCAL MAXIMA only (3 for NPY) in the shape
## 10_4_plot_proteins.R expects, so the full annotated page renders with the
## peptide-end model's windows in place of the production ones.
##
##   Rscript inst/scripts/10_9j_scan_peaks_as_nn_input.R --gene NPY
##   LF_NN_INPUT_COMB=~/AF2_analysis/lf_scan_peaks_nn_input_NPY.rds \
##     LF_PLOT_DIR=~/AF2_analysis/scan_peak_plots \
##     Rscript inst/scripts/plot_all.R NPY
##
## TWO GEOMETRY CAVEATS, since the page was written for the dibasic-anchored set:
##
##  1. The page derives each window's anchor mark as `end - 6`, i.e. window
##     positions 30-31, because that is where 9.2 puts the dibasic pair. A
##     peptide-end window anchors the peptide's last residue at position 28, so
##     positions 30-31 are residues anchor+2 and anchor+3. For an AMIDATED site
##     that is exactly the dibasic pair (anchor+1 is the glycine), so NPY's mark
##     lands correctly; for a plain dibasic anchor it sits one residue late.
##  2. `per_index` carries the peptide-end model's six classes. DB and gap do not
##     exist in that vocabulary, so those two lines are simply absent from the
##     detail panels -- correctly, not as a rendering failure.
##
## Output: ~/AF2_analysis/lf_scan_peaks_nn_input_<gene>.rds
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tibble); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
gene   <- .opt("--gene", "NPY")
term   <- toupper(.opt("--term", "C"))
scan_p <- path.expand(.opt("--scan", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
cache_p<- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
run_p  <- path.expand(.opt("--run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
out_p  <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_scan_peaks_nn_input_%s.rds", gene)))
top_n  <- as.integer(.opt("--top", "0"))        # 0 = local maxima only
max_pk <- as.integer(.opt("--max-peaks", "0"))  # cap the peaks, best first

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "dcnn_bridge.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)

ANCHOR <- c(N = 8L, C = 28L)[[term]]

## ---- 1. which windows -----------------------------------------------------------
fc <- readRDS(cache_p)
p  <- fc$prec %>% filter(gene == !!gene) %>% slice(1)
if (!nrow(p)) stop("no precursor called ", gene, call. = FALSE)
feat <- fc$feats[[match(p$accession, fc$prec$accession)]]
sc <- readRDS(scan_p)
pk <- sc$peaks %>% filter(accession == p$accession) %>% arrange(desc(score))
if (top_n > 0) {
  pk <- sc$scan %>% filter(accession == p$accession) %>% slice_max(score, n = top_n)
  message("using the top ", nrow(pk), " windows by score")
} else {
  message("using the ", nrow(pk), " local maxima")
  ## a long precursor has many local maxima -- KNG1 has 24 -- and every window is
  ## its own overlap layer, so cap them rather than stack two dozen rows
  if (max_pk > 0 && nrow(pk) > max_pk) {
    message("  capping to the top ", max_pk, " by score (dropped ",
            nrow(pk) - max_pk, ", best dropped scores ",
            sprintf("%.3f", pk$score[max_pk + 1L]), ")")
    pk <- pk %>% slice_head(n = max_pk)
  }
}
stopifnot(nrow(pk) > 0)
kn <- sc$knowns %>% filter(accession == p$accession)

## ---- 2. build and score them ----------------------------------------------------
cfg <- mod$Config$from_json(file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"),
                                      "in", "config.json"))
models <- mod$scan$load_member_models(
  file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out"), cfg, term)
pi_names <- as.character(cfg$pi_names)
w_start <- pk$anchor - ANCHOR + 1L
X <- array(0, c(nrow(pk), 36L, ncol(feat)))
for (i in seq_len(nrow(pk)))
  X[i, , ] <- as.matrix(lf_pepend_slice(feat, w_start[i], p$n_prot, p$c_prot)[, fc$all_params3])
xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))

g <- matrix(0, nrow(pk), length(models))
PI <- array(0, c(length(models), nrow(pk), 36L, length(pi_names)))
for (m in seq_along(models)) {
  out <- mod$pipeline$predict_all(models[[m]], xp)
  g[, m] <- as.numeric(reticulate::py_to_r(out[["global"]]))
  PI[m, , , ] <- reticulate::py_to_r(out[["per_index_cat"]])
}
## the scan's stored score is the same ensemble mean -- check, don't assume
stopifnot(max(abs(rowMeans(g) - pk$score)) < 1e-5)
message(sprintf("scored %d window(s) with %d members; max|mean - scan score| %.2g",
                nrow(pk), length(models), max(abs(rowMeans(g) - pk$score))))

aa <- strsplit(p$seq, "")[[1]]
pi_tbl <- function(m) as_tibble(setNames(as.data.frame(m), pi_names)) %>%
  mutate(index = row_number(), .before = 1)

## ---- 3. the nn_input_comb shape --------------------------------------------------
## `peps` records the REAL residue range, excluding padding, as 9.2 writes it --
## the page re-derives start/end from this string.
real <- lapply(seq_len(nrow(pk)), function(i) {
  r <- w_start[i] + 0:35
  r[r >= p$n_prot & r <= p$c_prot]
})
nn <- tibble(
  peps     = vapply(seq_len(nrow(pk)), function(i)
               sprintf("%s_w%d-%d", gene, min(real[[i]]), max(real[[i]])), ""),
  win_type = "scan_peak",
  gene     = gene,
  target   = term,
  pep_id   = NA_character_,
  known    = as.integer(pk$anchor %in% kn$anchor),
  ## `index` here is the PRECURSOR coordinate, not the window position 1..36:
  ## make_detail_panel binds it as `index_og` and uses it as the x axis of the
  ## click-through panel, which is the protein's residue axis. Numbering it
  ## 1..36 drew every panel over the first 36 residues of the protein instead of
  ## over its own window. (per_index keeps a window-local `index`, but that one
  ## is dropped by the panel -- only meta_data's is load-bearing.)
  meta_data = lapply(seq_len(nrow(pk)), function(i) {
    r <- w_start[i] + 0:35
    tibble(index = r,
           AA = ifelse(r >= p$n_prot & r <= p$c_prot,
                       aa[pmin(pmax(r, 1L), length(aa))], NA_character_),
           win_pos = seq_len(36L))
  }),
  pred_raw = rowMeans(g),
  pred_sd  = apply(g, 1, sd),
  pred     = rowMeans(g),
  per_index    = lapply(seq_len(nrow(pk)), function(i) pi_tbl(apply(PI[, i, , ], c(2, 3), mean))),
  per_index_sd = lapply(seq_len(nrow(pk)), function(i) pi_tbl(apply(PI[, i, , ], c(2, 3), sd))),
  nn_closest_peptide = NA_character_,
  nn_closest_sim     = NA_real_,
  model    = term,
  category = paste0(term, "_scan_peak"),
  rank     = pk$rank_all,
  rank_cat = pk$rank_peak,
  end_type = NA_character_)

saveRDS(nn, out_p)
message("\nnn_input_comb -> ", out_p)
print(as.data.frame(nn %>% mutate(anchor = pk$anchor) %>%
        select(peps, anchor, known, pred, pred_sd, rank_peak = rank_cat, rank_all = rank)),
      row.names = FALSE, digits = 4)
message("\nnow run:\n  LF_NN_INPUT_COMB=", out_p,
        " \\\n    LF_PLOT_DIR=~/AF2_analysis/scan_peak_plots \\\n",
        "    Rscript inst/scripts/plot_all.R ", gene)
