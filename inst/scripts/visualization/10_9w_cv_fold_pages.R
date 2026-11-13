#!/usr/bin/env Rscript
## ---- protein pages drawn by the 2-fold CV trunks ---------------------------------
## Every track on these pages comes from a model trained on half the knowns, so a
## page can show what a model that never saw the peptide says. Two modes:
##
##   held out (default) -- one page per gene in the fold map, each terminus taken
##     from the fold that held THAT terminus out, and the residue models from the
##     lf_resid fold that held the gene out. Those splits are independent (36 of
##     66 genes land on opposite sides), so the window fold letter must never pick
##     the residue fold.
##
##       Rscript inst/scripts/10_9w_cv_fold_pages.R
##
##   fixed fold -- every track from one fold, for genes that are not knowns (e.g.
##     CXCL14, CXCL17), where both folds are out-of-sample and both are worth a
##     page. Refuses a gene either fold trained on, window or residue side.
##
##       Rscript inst/scripts/10_9w_cv_fold_pages.R --genes CXCL14,CXCL17 --fold A
##
##     The fold xgb npz files (lf_winxgb_cv_*) only cover the knowns, so score the
##     xgb head for new genes first, into a suffixed file that leaves them alone:
##
##       cd inst/python
##       python -m lf_winxgb.scan_genes --genes CXCL14,CXCL17 --term C \
##         --run ~/AF2_analysis/lf_pepend_run_C_foldA.rds \
##         --groups ~/AF2_analysis/lf_pepend_window_groups_C_foldA.npz \
##         --pm-rounds 100 --skip-stacked \
##         --out ~/AF2_analysis/lf_winxgb_cv_C_foldA_cxcl.npz
##
##     (x N/C x A/B) and pass --xgb-suffix _cxcl. --pm-rounds 100 is what the
##     held-out set used; the heads reproduce its val metrics exactly.
##
## Output: <out>/<GENE>.html (held out) or <out>/fold<F>/<GENE>.html (fixed)
## ------------------------------------------------------------------------------
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
A2      <- "~/AF2_analysis"
fold    <- .opt("--fold", NA_character_)
genes   <- .opt("--genes", NA_character_)
xsuf    <- .opt("--xgb-suffix", "")
meta_d  <- path.expand(.opt("--meta", file.path(A2, "cv_meta")))
out_d   <- .opt("--out", file.path(A2, if (is.na(fold)) "cv_heldout" else "cv_fold_pages"))
if (!is.na(genes)) genes <- strsplit(genes, ",", fixed = TRUE)[[1]]
if (!is.na(fold) && !fold %in% c("A", "B")) stop("--fold must be A or B", call. = FALSE)
if (!is.na(fold) && anyNA(genes)) stop("--fold needs --genes", call. = FALSE)

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()

## ---- which fold draws which track, per gene -------------------------------------
if (is.na(fold)) {
  asg <- readRDS(file.path(meta_d, "gene_folds.rds")) %>%
    left_join(readRDS(file.path(meta_d, "resid_fold_map.rds")) %>%
                select(gene, resid_fold = resid_fold_heldout), by = "gene")
  stopifnot(!any(is.na(asg$resid_fold)))
  if (!anyNA(genes)) asg <- asg %>% filter(gene %in% genes)
  asg <- asg %>% mutate(sub = "")
} else {
  ## a fixed fold is only honest for a gene that fold never trained on.
  ## The bridge goes up BEFORE numpy: reticulate binds whichever python it meets
  ## first, and a bare numpy import binds one without tensorflow -- the class
  ## tracks then drop off every page without an error.
  source(file.path(ROOT, "R", "dcnn_bridge.R"))
  lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
  np <- reticulate::import("numpy", convert = FALSE)
  fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
  acc <- fc$prec$accession[match(genes, fc$prec$gene)]
  win_known <- unlist(lapply(c("N", "C"), function(t) {
    x <- readRDS(path.expand(file.path(A2, sprintf("lf_pepend_nn_input_fold%s.rds", fold))))$nn_input[[t]]
    x$train$accession[x$train$known == 1]
  }))
  z  <- np$load(path.expand(file.path(A2, sprintf("lf_resid_preds_fold%s.npz", fold))), allow_pickle = TRUE)
  ra <- as.character(reticulate::py_to_r(z[["accession"]]))
  rp <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
  rt <- as.integer(reticulate::py_to_r(z[["in_train"]]))
  res_known <- unique(ra[rp[rt > 0] + 1L])
  bad <- genes[acc %in% c(win_known, res_known)]
  if (length(bad)) stop("fold ", fold, " trained on ", paste(bad, collapse = ", "),
                        " -- use the held-out mode for these", call. = FALSE)
  asg <- tibble(gene = genes, fold_C = fold, fold_N = fold, resid_fold = fold,
                clean = TRUE, n_held = NA_integer_, n_ends = NA_integer_,
                sub = paste0("fold", fold))
}
## the bundle is ~1 GB; load it only once the request is known to be honest
.b <- readRDS(file.path(A2, "plot_bundle.rds"))
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
source(file.path(ROOT, "R", "pepend_tracks.R"))
gs <- purrr::map_chr(the_input, \(x) x$gene)

asg <- asg %>% filter(gene %in% gs)
message(sprintf("rendering %d gene(s) -> %s", nrow(asg), out_d))

for (r in seq_len(nrow(asg))) {
  g <- asg$gene[r]; fC <- asg$fold_C[r]; fN <- asg$fold_N[r]
  plot_dir <- file.path(out_d, asg$sub[r])
  dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)
  options(
    lf.track_run_C    = sprintf("%s/lf_pepend_run_C_fold%s.rds",  A2, fC),
    lf.track_run_N    = sprintf("%s/lf_pepend_run_N_fold%s.rds",  A2, fN),
    lf.track_scan_C   = sprintf("%s/lf_pepend_scan_C_fold%s.rds", A2, fC),
    lf.track_scan_N   = sprintf("%s/lf_pepend_scan_N_fold%s.rds", A2, fN),
    lf.track_winxgb_C = sprintf("%s/lf_winxgb_cv_C_fold%s%s.npz", A2, fC, xsuf),
    lf.track_winxgb_N = sprintf("%s/lf_winxgb_cv_N_fold%s%s.npz", A2, fN, xsuf),
    lf.track_resid    = sprintf("%s/lf_resid_preds_fold%s.npz",   A2, asg$resid_fold[r]))
  ## a fresh cache each gene: the holders key on file basename, but the env
  ## persists, and a stale array from the previous gene's fold would be silent
  assign(".lf_tracks_env", new.env(parent = emptyenv()), envir = globalenv())
  i <- which(gs == g)[1]
  tm <- tryCatch(system.time(
    make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]]))[["elapsed"]],
    error = function(e) { message("!!! ", g, ": ", conditionMessage(e)); NA_real_ })
  message(sprintf("=== %-10s C=fold%s N=fold%s resid=fold%s  %s  %5.1f s", g, fC, fN,
                  asg$resid_fold[r],
                  if (asg$clean[r]) "all held out" else
                    sprintf("%d/%d held out", asg$n_held[r], asg$n_ends[r]), tm))
}
