## plot_only.R, with ONE difference: lf.track_extra_win is never set, so the
## CT-/NT-padded arms contribute no tracks. Every other input keeps its default.
## Writes to its own directory; ~/AF2_analysis/pepend_pages is not touched.
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W  <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
plot_dir <- Sys.getenv("LF_PLOT_DIR", "~/AF2_analysis/pepend_pages_nopad")
.b <- readRDS("~/AF2_analysis/plot_bundle.rds")
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source(file.path(W, "R/plot_proteins_win_new.R"))
source(file.path(W, "R/pepend_tracks.R"))      # the optimised copy (ec5506f)
options(lf.pepend_python_path = file.path(W, "inst", "python"),
        ## one class track per TERMINUS, not their average -- the existing
        ## pepend_pages set this (their axis reads "N/C terminus windows"), and
        ## without it the N and C class tracks collapse into one averaged track
        lf.class_tracks_by_term = TRUE)
stopifnot(length(getOption("lf.track_extra_win", list())) == 0)   # no padded arms
if (nzchar(Sys.getenv("LF_REBUILD"))) {
  ## point every model input at the rebuilt (all-residue, 29-channel) arm.
  ## LF_TRACK_PATHS is a plain object in this env, so the cache is overridden by
  ## assignment rather than by editing pepend_tracks.R -- the members are
  ## 29-channel and MUST be fed the 29-channel cache.
  options(lf.track_run_C    = "~/AF2_analysis/lf_pepend_run_C_allres_slack8.rds",
          lf.track_run_N    = "~/AF2_analysis/lf_pepend_run_N_allres_slack8.rds",
          lf.track_scan_C   = "~/AF2_analysis/lf_pepend_scan_C_allres.rds",
          lf.track_scan_N   = "~/AF2_analysis/lf_pepend_scan_N_allres.rds",
          lf.track_winxgb_C = "~/AF2_analysis/lf_winxgb_scan_C_allres_review.npz",
          lf.track_winxgb_N = "~/AF2_analysis/lf_winxgb_scan_N_allres_review.npz",
          lf.track_resid    = "~/AF2_analysis/lf_resid_preds_allres.npz")
  LF_TRACK_PATHS$cache <- "~/AF2_analysis/lf_pepend_residue_cache_allres.rds"
  for (k in c("run_C","run_N","scan_C","scan_N","winxgb_C","winxgb_N","resid")) {
    f <- path.expand(getOption(paste0("lf.track_", k)))
    message(sprintf("  %-10s %s %s", k, if (file.exists(f)) "OK " else "MISSING", basename(f)))
  }
  message(sprintf("  %-10s %s %s", "cache",
    if (file.exists(path.expand(LF_TRACK_PATHS$cache))) "OK " else "MISSING",
    basename(LF_TRACK_PATHS$cache)))
}
dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)
want <- commandArgs(trailingOnly = TRUE)
gs   <- purrr::map_chr(the_input, \(x) x$gene)
idx  <- if (length(want)) which(gs %in% want) else seq_along(the_input)
if (!length(idx)) stop("none of those genes are in the bundle")
for (i in idx) {
  message("--- ", gs[i], " ---")
  tm <- system.time(
    make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]])
  )[["elapsed"]]
  message(sprintf("    %s: %.1f s", gs[i], tm))
}
