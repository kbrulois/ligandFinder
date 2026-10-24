## Protein pages: the cv5 pages (pepend_pages_cv5, from cv5_pages.R) plus one
## more track group, "no mask" -- the trunk retrained WITHOUT the class-position
## mask (10_8 --no-mask), 2-fold CV x 5 seeds, both folds pooled into one
## 10-member ensemble per terminus. It gets a CNN track and an xgb track like
## the per-fold group, and the masked all-50 ensemble is overlaid into it, so
## the two models read on one pair of axes.
##
## Page code: the modest-satoshi worktree's R files, read-only, exactly as the
## cv5 pages used them (its track changes are uncommitted there).
## Python: THIS worktree's inst/python. The no-mask members' config.json says
## class_position_mask = false; the older lf_dcnn would drop that field as
## unknown and rebuild those members WITH the mask -- silently wrong scores.
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W  <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
WP <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7/inst/python"
plot_dir <- Sys.getenv("LF_PLOT_DIR", "~/AF2_analysis/pepend_pages_nomask")
.b <- readRDS("~/AF2_analysis/plot_bundle.rds")
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source(file.path(W, "R/plot_proteins_win_new.R")); source(file.path(W, "R/pepend_tracks.R"))
## the ChimeraX helpers, definitions only (see cv5_pages.R: sourcing the file
## whole runs a batch that errors on this bundle)
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
A2 <- path.expand("~/AF2_analysis")
## the cv5 per-fold series, as on the cv5 pages -- their xgb from the four-gene
## npz files (cv5_npz_4genes.R), which add CXCL17 and are otherwise identical
perfold <- unlist(lapply(1:5, function(k) lapply(c("C","N"), function(tm)
  list(name = sprintf("%s fold %d", tm, k), term = tm, group = "per-fold", popup = TRUE,
       run = file.path(A2, sprintf("lf_pepend_run_%s_cv5_f%d_slack8.rds", tm, k)),
       winxgb = file.path(A2, sprintf("lf_winxgb_scan_%s_cv5_f%d_nomaskpages.npz", tm, k))))),
  recursive = FALSE)
## the new group: the run path only names the `_isolated` dir holding the 10
## pooled members (nomask_post.sh); there is no rds behind it
nomask <- lapply(c("C","N"), function(tm)
  list(name = sprintf("%s no mask", tm), term = tm, group = "no mask", popup = TRUE,
       run = file.path(A2, sprintf("lf_pepend_run_%s_2fold_nomask.rds", tm)),
       winxgb = file.path(A2, sprintf("lf_winxgb_scan_%s_2fold_nomask.npz", tm))))
for (s in nomask) stopifnot(dir.exists(sub("\\.rds$", "_isolated", s$run)), file.exists(s$winxgb))
options(lf.pepend_python_path = WP,
        lf.class_tracks_by_term = TRUE,
        lf.pepend_per_index = TRUE,
        lf.pepend_popup_max = as.integer(Sys.getenv("LF_POPUP_MAX", "2000")),
        lf.win_overlay_all50 = TRUE,
        lf.cv_knowns = list(N = file.path(A2, "lf_pepend_cv5_knowns_N.csv"),
                            C = file.path(A2, "lf_pepend_cv5_knowns_C.csv")),
        lf.track_run_C = file.path(A2, "lf_pepend_run_C_cv5all50_slack8.rds"),
        lf.track_run_N = file.path(A2, "lf_pepend_run_N_cv5all50_slack8.rds"),
        lf.track_scan_C   = file.path(A2, "lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N   = file.path(A2, "lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C = file.path(A2, "lf_winxgb_scan_C_cv5all50_nomaskpages.npz"),
        lf.track_winxgb_N = file.path(A2, "lf_winxgb_scan_N_cv5all50_nomaskpages.npz"),
        lf.track_resid = file.path(A2, "lf_resid_preds_v2.npz"),
        lf.track_extra_win = c(perfold, nomask))
dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)
want <- commandArgs(trailingOnly = TRUE)
gs <- purrr::map_chr(the_input, \(x) x$gene)
stopifnot(all(want %in% gs))
for (i in which(gs %in% want)) {
  message("--- ", gs[i], " ---")
  tm <- system.time(make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]]))[["elapsed"]]
  message(sprintf("    %s: %.1f s", gs[i], tm))
}
## provenance: the page scores the no-mask CNN series on the fly from the pooled
## members; nomask_post.sh scored the same members through lf_winxgb.scan_genes.
## Two routes, one model -- they must agree, or the page drew something else.
np <- reticulate::import("numpy", convert = FALSE)
for (s in nomask) for (g in want) {
  fc <- .lf_tracks_env[["feat_cache"]]
  i <- match(g, fc$prec$gene)
  ex <- lf_track_extra_head(modifyList(s, list(popup = FALSE)), fc$prec$accession[i],
                            fc$prec$n_prot[i], fc$prec$c_prot[i])
  z <- np$load(s$winxgb)
  a <- reticulate::py_to_r(z[[paste0(g, "__anchor")]]); v <- reticulate::py_to_r(z[[paste0(g, "__cnn")]])
  m <- match(ex$win$residue, a)
  message(sprintf("provenance %s %-6s: %d/%d anchors matched, max |page - scan| CNN = %.2e",
                  s$name, g, sum(!is.na(m)), nrow(ex$win), max(abs(ex$win$mean - v[m]), na.rm = TRUE)))
}
