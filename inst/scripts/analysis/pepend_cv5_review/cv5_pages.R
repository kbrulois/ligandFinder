## Protein pages from the 5-fold CV tracks.
## Main window tracks  = the all-50 ensemble (CNN from the scan rds, xgb from
##                       the per-gene npz) -- both from the SAME 50 models.
## `per-fold` group    = one series per fold, from that fold's own run and npz.
## The per-index class tracks, popups and ChimeraX script are ON, and the
## per-index surface IS all 50 models: `..._cv5all50_slack8_isolated` holds
## symlinks to every fold's member weights under one members/ dir, which is all
## `load_member_models` globs for. So the class tracks, the click-to-open window
## popups and the window tracks all come from the same 50 models -- no mixing,
## and nothing is lost.
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W  <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
plot_dir <- Sys.getenv("LF_PLOT_DIR", "~/AF2_analysis/pepend_pages_cv5")
.b <- readRDS("~/AF2_analysis/plot_bundle.rds")
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source(file.path(W, "R/plot_proteins_win_new.R")); source(file.path(W, "R/pepend_tracks.R"))
## make_cm_script_text() lives in generate_cm_sct.R; without it the page
## silently drops the ChimeraX script ("unavailable (source generate_cm_sct.R)").
## That file is NOT a pure function file -- sourcing it runs a batch over
## `the_input[c(1:5435, 5441:5614)]`, which errors on a 5,180-entry bundle. So
## evaluate only its function definitions.
local({
  e <- parse(file.path(W, "inst/scripts/visualization/generate_cm_sct.R"))
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
extra <- unlist(lapply(1:5, function(k) lapply(c("C","N"), function(tm)
  ## popup = TRUE builds that fold's per-index surface too, so a dot on the
  ## per-fold track opens its window like the main track's dots do. Without it
  ## lf_track_extra_head returns pix = NULL, the series is missing from the
  ## popup payload, and seriesOf() finds nothing -- the dot is inert.
  list(name = sprintf("%s fold %d", tm, k), term = tm, group = "per-fold",
       popup = TRUE,
       run = file.path(A2, sprintf("lf_pepend_run_%s_cv5_f%d_slack8.rds", tm, k)),
       winxgb = file.path(A2, sprintf("lf_winxgb_scan_%s_cv5_f%d.npz", tm, k))))),
  recursive = FALSE)
options(lf.pepend_python_path = file.path(W, "inst", "python"),
        lf.class_tracks_by_term = TRUE,
        lf.pepend_per_index = TRUE,
        ## default 600: on a precursor with more anchors than that, only the
        ## top 600 BY SCORE carry window data and every other dot is a dead
        ## click. ANO8 has 1,232 anchors and BRINP3 766, so both were half
        ## inert. Raised to cover the largest gene here; it is a page-size
        ## tradeoff, not a correctness one, so lower it again for bulk runs.
        ## LF_POPUP_MAX overrides this per run. 2000 keeps window data for every
        ## anchor on the largest gene here; the five biggest precursors (C5 1658
        ## anchors, C3, EGF, COPA, ANO8) then produce 80-120 MB pages, so they are
        ## rebuilt at 600 -- the top 600 anchors BY SCORE stay clickable and the
        ## rest become inert dots. A page-size tradeoff, not a correctness one.
        lf.pepend_popup_max = as.integer(Sys.getenv("LF_POPUP_MAX", "2000")),
        ## the all-50 ensemble stops being its own track and is drawn into each
        ## per-fold track as a heavier line with square points
        lf.win_overlay_all50 = TRUE,
        ## one row per known end with the fold that held it out: that fold's
        ## point gets a downward triangle on the per-fold track
        lf.cv_knowns = list(N = file.path(A2, "lf_pepend_cv5_knowns_N.csv"),
                            C = file.path(A2, "lf_pepend_cv5_knowns_C.csv")),
        lf.track_run_C = file.path(A2, "lf_pepend_run_C_cv5all50_slack8.rds"),
        lf.track_run_N = file.path(A2, "lf_pepend_run_N_cv5all50_slack8.rds"),
        lf.track_scan_C   = file.path(A2, "lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N   = file.path(A2, "lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C = file.path(A2, "lf_winxgb_scan_C_cv5all50.npz"),
        lf.track_winxgb_N = file.path(A2, "lf_winxgb_scan_N_cv5all50.npz"),
        ## the per-residue models RETRAINED in this session: inner stop set, the
        ## regularised xgb preset, context kept, all four targets. Omitting this
        ## silently falls back to LF_TRACK_PATHS$resid -- the Oct 5 file -- and
        ## puts stale per-residue tracks beside fresh CV window tracks.
        lf.track_resid = file.path(A2, "lf_resid_preds_v2.npz"),
        lf.track_extra_win = if (nzchar(Sys.getenv("LF_NO_FOLDS"))) list() else extra)
dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)
want <- commandArgs(trailingOnly = TRUE)
gs <- purrr::map_chr(the_input, \(x) x$gene)
for (i in which(gs %in% want)) {
  message("--- ", gs[i], " ---")
  tm <- system.time(make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]]))[["elapsed"]]
  message(sprintf("    %s: %.1f s", gs[i], tm))
}
