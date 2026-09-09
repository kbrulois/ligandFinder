## =============================================================================
## Plot proteins in a FRESH R session. No training, no secretome_aa, no keras.
## Needs only the bundle written by save_plot_bundle.R.
##
##   Rscript inst/scripts/plot_only.R            # plots every gene in the bundle
##   Rscript inst/scripts/plot_only.R ANO8 NPY   # or just these
## =============================================================================
## The Rscript first on PATH may be a pixi install with no packages; say so
## plainly instead of dying on `library(tidyverse)`.
local({
  missing <- Filter(\(p) !requireNamespace(p, quietly = TRUE),
                    c("tidyverse", "ggiraph", "patchwork", "ligandFinder"))
  if (length(missing))
    stop("missing package(s): ", paste(missing, collapse = ", "),
         "\n  R being used: ", file.path(R.home("bin"), "Rscript"),
         "\n  if that is not the R you install packages into, run it explicitly:",
         "\n    /usr/local/bin/Rscript inst/scripts/plot_only.R", call. = FALSE)
})
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))

bundle_path <- "~/AF2_analysis/plot_bundle.rds"
plot_dir    <- "~/AF2_analysis/new_meth_plot2"

if (!file.exists(path.expand(bundle_path)))
  stop(bundle_path, " does not exist.\n",
       "  It is written by 10_4_plot_proteins.R, which needs a session with\n",
       "  nn_input_comb, secretome, secretome_aa and peps_tp:\n",
       "    source(\"inst/scripts/10_4_plot_proteins.R\")\n",
       "  Run that once; afterwards this script re-renders with no keras.",
       call. = FALSE)

.b <- readRDS(bundle_path)
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source("~/R_projects/ligandFinder/R/plot_proteins_win_new.R")
dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)

want <- commandArgs(trailingOnly = TRUE)
gs   <- purrr::map_chr(the_input, \(x) x$gene)
idx  <- if (length(want)) which(gs %in% want) else seq_along(the_input)
if (!length(idx)) stop("none of those genes are in the bundle: ", paste(gs, collapse = ", "))

for (i in idx) {
  message("--- ", gs[i], " ---")
  tm <- system.time(
    make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]])
  )[["elapsed"]]
  message(sprintf("    %s: %.1f s", gs[i], tm))
}
