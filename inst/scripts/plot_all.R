## =============================================================================
## Render every protein page from files on disk. No RStudio session, no keras,
## no models, no retraining -- the 5-seed ensemble predictions are already saved.
##
## Run from the repo root:
##   /usr/local/bin/Rscript inst/scripts/plot_all.R
##   /usr/local/bin/Rscript inst/scripts/plot_all.R ANO8 NPY PENK   # just these
##
## Output dir and inputs can be overridden with env vars:
##   LF_PLOT_DIR, LF_NN_INPUT_COMB, LF_SECRETOME, LF_SECRETOME_AA
##
## This is the "cold start" entry point. `plot_only.R` is the faster one, but it
## needs plot_bundle.rds -- which this script rewrites at the end, covering every
## gene it rendered. So: run this once, then use plot_only.R for touch-ups.
## =============================================================================
if (!file.exists("inst/scripts/10_4_plot_proteins.R"))
  stop("run this from the repo root (", getwd(), " is not it)", call. = FALSE)

.env <- function(k, default) { v <- Sys.getenv(k); if (nzchar(v)) v else default }

nn_input_comb_path <- .env("LF_NN_INPUT_COMB", "~/AF2_analysis/nn_input_comb_ensemble.rds")
secretome_path     <- .env("LF_SECRETOME",     "~/AF2_analysis/secretome.rds")
secretome_aa_path  <- .env("LF_SECRETOME_AA",
                           "~/peptide_alg/build_residue_db/processed/secretome_aa_1.rds")
plot_dir           <- .env("LF_PLOT_DIR",      "~/AF2_analysis/ensemble_plots_full")

## Check every input up front: these files take minutes to load, so a missing one
## should not surface twenty minutes in.
local({
  missing <- Filter(\(p) !file.exists(path.expand(p)),
                    c(nn_input_comb_path, secretome_path, secretome_aa_path))
  if (length(missing))
    stop("missing input(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  pkgs <- Filter(\(p) !requireNamespace(p, quietly = TRUE),
                 c("tidyverse", "ggiraph", "patchwork", "ggnewscale", "ligandFinder"))
  if (length(pkgs))
    stop("missing package(s): ", paste(pkgs, collapse = ", "),
         "\n  R being used: ", file.path(R.home("bin"), "Rscript"),
         "\n  if that is the pixi R, run /usr/local/bin/Rscript explicitly", call. = FALSE)
})

suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))

.load <- function(path, what) {
  t0 <- Sys.time()
  x  <- readRDS(path.expand(path))
  message(sprintf("loaded %-14s %s  (%.0f s, %.1f GB on disk)", what, path,
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  file.size(path.expand(path)) / 2^30))
  x
}

nn_input_comb <- .load(nn_input_comb_path, "nn_input_comb")
secretome     <- .load(secretome_path,     "secretome")
secretome_aa  <- .load(secretome_aa_path,  "secretome_aa")
peps_tp       <- .load("~/AF2_analysis/peps_tp.rds", "peps_tp")

## ---- AlphaFold-predicted disulfides ----------------------------------------
## secretome.rds carries UniProt's disulfides only (type "disulfide bond", which
## assemble_anno_feats renames to unidsb). The AlphaFold-derived ones live in a
## side table written by get_afdsb.R; merge them in here rather than rewriting
## the 1.8 GB secretome, so removing the file reverts the plots.
##
## Prepended, matching get_dsb.R. A pair predicted by AlphaFold AND annotated by
## UniProt collapses to feature "both" inside assemble_anno_feats, which is what
## puts all three keys in the disulfide legend.
local({
  f <- path.expand(.env("LF_AFDSB", "~/AF2_analysis/afdsb.rds"))
  if (!file.exists(f)) {
    message("no afdsb.rds -- disulfides will be UniProt-only. ",
            "Run: Rscript inst/scripts/get_afdsb.R")
    return(invisible())
  }
  af <- readRDS(f)
  ## idempotent: never add a second copy to a secretome that already has them
  if (any(purrr::map_lgl(utils::head(secretome[["features"]], 200),
                         \(x) !is.null(x) && "type" %in% names(x) && any(x$type == "afdsb")))) {
    message("secretome already carries afdsb -- not merging again")
    return(invisible())
  }
  by_acc <- split(dplyr::select(af, -accession), af$accession)
  secretome[["features"]] <<- purrr::map2(
    secretome[["features"]], secretome[["accession"]],
    \(x, a) if (is.null(by_acc[[a]])) x else dplyr::bind_rows(by_acc[[a]], x))
  message(sprintf("merged %d AlphaFold disulfides across %d proteins",
                  nrow(af), dplyr::n_distinct(af$accession)))
})

if (!"per_index_sd" %in% names(nn_input_comb))
  warning("nn_input_comb has no per_index_sd -- this is a single-model run, so ",
          "the pages will carry no ensemble spread (no ribbons, no +/- sd)")

.want <- commandArgs(trailingOnly = TRUE)
if (length(.want)) plot_genes <- .want

## 10_4 rewrites plot_bundle.rds from whatever `the_input` holds, and defaults to
## doing so. On a subset run that would replace a full 5,180-gene bundle with a
## two-gene one -- so only write it when rendering everything.
write_bundle <- !length(.want)
if (!write_bundle)
  message("subset run: leaving plot_bundle.rds alone")

message("rendering ", if (length(.want)) length(.want) else dplyr::n_distinct(nn_input_comb$gene),
        " gene(s) -> ", plot_dir)
source("inst/scripts/10_4_plot_proteins.R")
