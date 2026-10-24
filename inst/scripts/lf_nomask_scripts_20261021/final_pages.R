## Protein pages from the FLATTEN window head only -- no attention anywhere
## (Kevin, 2026-10-21: the x1 flatten head is the model going forward, and the
## CNN attention tracks come off the pages).
##
##   MODE=2fold  the 2-fold CV flatten runs (foldA/foldB x 5 seeds), pooled 10
##   MODE=cv5    the 5-fold CV flatten runs (cv5_f1..f5 x 10 seeds), pooled 50
##
## Every model-derived track comes from those runs:
##   class-probability tracks + the page's own popups  <- the pooled ensemble
##                                                        (lf.track_run_*)
##   CNN linear head, xgb window head                 <- one thin series per fold,
##                                                        the pooled ensemble heavy
## The page's own scan-backed window series are switched OFF (no scan rds / npz):
## those were the attention trunk's. The window-free MLP/XGB residue tracks stay.
##
##   MODE=2fold NPZ_SFX=_brinp3 Rscript final_pages.R BRINP3
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W  <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7"
WP <- file.path(W, "inst", "python")
MODE <- Sys.getenv("MODE", "cv5"); stopifnot(MODE %in% c("2fold", "cv5"))
SFX  <- Sys.getenv("NPZ_SFX", "")        # which per-gene xgb npz set (see arm_post.sh)
A2 <- path.expand("~/AF2_analysis")
plot_dir <- Sys.getenv("LF_PLOT_DIR", if (MODE == "cv5") "~/AF2_analysis/pepend_pages_flatten_cv5"
                                      else "~/AF2_analysis/pepend_pages_flatten")
.b <- readRDS("~/AF2_analysis/plot_bundle.rds")
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source(file.path(W, "R/plot_proteins_win_new.R")); source(file.path(W, "R/pepend_tracks.R"))
## The peptide-annotation helpers live in R/plot_proteins.R, which the page does
## NOT source -- it uses the INSTALLED package's copy. Take this worktree's
## definitions of the four it needs (the "never seen" line), and nothing else.
local({
  e <- parse(file.path(W, "R/plot_proteins.R"))
  keep <- c("anno_feats", "anno_feat_colors", "pep_nudges", "assemble_anno_feats")
  hit <- vapply(e, function(x) is.call(x) && length(x) >= 3 &&
                  tryCatch(as.character(x[[1]])[1], error = function(...) "") %in% c("<-", "=") &&
                  tryCatch(as.character(x[[2]])[1], error = function(...) "") %in% keep, logical(1))
  for (i in which(hit)) eval(e[[i]], envir = globalenv())
  stopifnot("never seen" %in% anno_feats$feature)
  message("annotation helpers from the worktree: ", sum(hit), " definitions")
})
## The never-seen peptides: the evaluation set the window metrics use
## (cv5_unseen_ends.csv -- ends whose gene no fold trained a window on), one row
## per END; rebuilt here into one span per peptide.
never_seen <- local({
  u <- utils::read.csv(path.expand("~/AF2_analysis/lf_nomask_scripts_20261021/cv5_unseen_ends.csv"),
                       stringsAsFactors = FALSE)
  u$start <- ifelse(u$term == "N", u$anchor, u$anchor - u$len + 1L)
  u$end   <- ifelse(u$term == "N", u$anchor + u$len - 1L, u$anchor)
  u <- unique(u[, c("gene", "accession", "name", "start", "end", "len", "set")])
  u$description <- sprintf("never seen: %s %d-%d (%d aa%s; no CV fold trained on this gene)",
                           u$name, u$start, u$end, u$len,
                           ifelse(grepl("^large", u$set), ", large cleavage product", ""))
  u
})
## the ChimeraX helpers, definitions only (sourcing the file whole runs a batch)
local({
  e <- parse(file.path(W, "inst/scripts/generate_cm_sct.R"))
  isfn <- vapply(e, function(x) {
    if (!is.call(x) || length(x) < 3) return(FALSE)
    op <- tryCatch(as.character(x[[1]])[1], error = function(...) "")
    if (!op %in% c("<-", "=")) return(FALSE)
    is.call(x[[3]]) && identical(tryCatch(as.character(x[[3]][[1]])[1],
                                          error = function(...) ""), "function")
  }, logical(1))
  for (i in which(isfn)) eval(e[[i]], envir = globalenv())
  message("chimerax helpers loaded: ", sum(isfn), " functions")
})

if (MODE == "cv5") {
  folds <- as.character(1:5)
  fold_run <- function(tm, f) sprintf("lf_pepend_run_%s_cv5_f%s_nomask_flatten", tm, f)
  fold_npz <- function(tm, f) sprintf("lf_winxgb_scan_%s_cv5_f%s_nomask_flatten%s.npz", tm, f, SFX)
  ens_run  <- function(tm) sprintf("lf_pepend_run_%s_cv5all50_nomask_flatten", tm)
  ens_npz  <- function(tm) sprintf("lf_winxgb_scan_%s_cv5all50_nomask_flatten%s.npz", tm, SFX)
  ens_name <- "all 50"
} else {
  folds <- c("A", "B")
  fold_run <- function(tm, f) sprintf("lf_pepend_run_%s_fold%s_nomask_flatten", tm, f)
  fold_npz <- function(tm, f) sprintf("lf_winxgb_scan_%s_fold%s_nomask_flatten%s.npz", tm, f, SFX)
  ens_run  <- function(tm) sprintf("lf_pepend_run_%s_2fold_nomask_flatten", tm)
  ens_npz  <- function(tm) sprintf("lf_winxgb_scan_%s_2fold_nomask_flatten%s.npz", tm, SFX)
  ens_name <- "both folds"
}
## one thin series per fold, and the pooled ensemble as the track's heavy line.
## The ensemble is named "<ens>, N terminus" and carries no popup of its own: its
## surface IS the page's own N series (same members, lf.track_run_N), and a click
## resolves by finding "N terminus" in the tooltip -- so it opens that, once.
specs <- unlist(lapply(c("C", "N"), function(tm) c(
  lapply(folds, function(f)
    list(name = sprintf("%s fold %s", tm, f), term = tm, group = "per-fold", popup = TRUE,
         cnn_label = "CNN linear head",
         run = file.path(A2, paste0(fold_run(tm, f), ".rds")),
         winxgb = file.path(A2, fold_npz(tm, f)))),
  list(list(name = sprintf("%s, %s terminus", ens_name, tm), term = tm, group = "per-fold",
            popup = FALSE, group_ensemble = TRUE, cnn_label = "CNN linear head",
            run = file.path(A2, paste0(ens_run(tm), ".rds")),
            winxgb = file.path(A2, ens_npz(tm)))))), recursive = FALSE)
for (s in specs) stopifnot(dir.exists(sub("\\.rds$", "_isolated", s$run)), file.exists(s$winxgb))

## held-out marks: the fold whose val split holds a known end is the honest
## series on that precursor (triangles)
held <- if (MODE == "cv5") NULL else do.call(rbind, lapply(folds, function(f) {
  k <- as.data.frame(readRDS(file.path(A2, sprintf("lf_pepend_nn_input_fold%s.rds", f)))$knowns)
  k <- k[k$split == "val", , drop = FALSE]
  data.frame(term = sprintf("%s fold %s", k$term, f), terminus = k$term,
             accession = k$accession, residue = k$anchor)
}))
options(lf.pepend_python_path = WP,
        lf.class_tracks_by_term = TRUE,
        lf.pepend_per_index = TRUE,
        lf.pepend_popup_max = as.integer(Sys.getenv("LF_POPUP_MAX", "2000")),
        lf.win_overlay_all50 = TRUE,
        lf.track_cnn_label = "CNN linear head",
        ## the page's own trunk = the pooled flatten ensemble: class tracks + popups
        lf.track_run_C = file.path(A2, paste0(ens_run("C"), ".rds")),
        lf.track_run_N = file.path(A2, paste0(ens_run("N"), ".rds")),
        ## no scan-backed own window series: those were the attention trunk's
        lf.track_scan_C = "none", lf.track_scan_N = "none",
        lf.track_winxgb_C = "none", lf.track_winxgb_N = "none",
        lf.cv_knowns = if (MODE == "cv5") list(N = file.path(A2, "lf_pepend_cv5_knowns_N.csv"),
                                               C = file.path(A2, "lf_pepend_cv5_knowns_C.csv")),
        lf.heldout_extra = held,
        lf.track_resid = file.path(A2, "lf_resid_preds_v2.npz"),
        lf.track_extra_win = specs)
dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)
want <- commandArgs(trailingOnly = TRUE)
gs <- purrr::map_chr(the_input, \(x) x$gene)
stopifnot(all(want %in% gs))
for (i in which(gs %in% want)) {
  message("--- ", gs[i], " ---")
  ns <- never_seen[never_seen$accession == the_input[[i]]$accession, , drop = FALSE]
  if (nrow(ns)) {
    the_input[[i]]$features[[1]] <- dplyr::bind_rows(
      the_input[[i]]$features[[1]],
      tibble::tibble(source = "never_seen", type = "peptide", start = as.integer(ns$start),
                     end = as.integer(ns$end), description = ns$description))
    message(sprintf("    %d never-seen peptide(s): %s", nrow(ns),
                    paste(sprintf("%d-%d", ns$start, ns$end), collapse = ", ")))
  }
  tm <- system.time(make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]]))[["elapsed"]]
  message(sprintf("    %s: %.1f s", gs[i], tm))
}
## provenance on up to 3 genes: the page scores every series on the fly from its
## members; arm_post.sh / scan_genes scored the same members separately. They
## must agree, or the page drew some other model.
np <- reticulate::import("numpy", convert = FALSE)
fc <- .lf_tracks_env[["feat_cache"]]
for (g in head(want, as.integer(Sys.getenv("PROV_N", "3")))) for (s in specs) {
  i <- match(g, fc$prec$gene)
  ex <- lf_track_extra_head(modifyList(s, list(popup = FALSE)), fc$prec$accession[i],
                            fc$prec$n_prot[i], fc$prec$c_prot[i])
  z <- np$load(s$winxgb)
  a <- reticulate::py_to_r(z[[paste0(g, "__anchor")]]); v <- reticulate::py_to_r(z[[paste0(g, "__cnn")]])
  m <- match(ex$win$residue, a)
  message(sprintf("provenance %-22s %-7s: %d/%d anchors, max |page - scan| CNN = %.2e",
                  s$name, g, sum(!is.na(m)), nrow(ex$win), max(abs(ex$win$mean - v[m]), na.rm = TRUE)))
}
