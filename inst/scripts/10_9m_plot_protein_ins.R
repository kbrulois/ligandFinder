#!/usr/bin/env Rscript
## ---- the protein page, with the insertion head's call on each window ----------
## plot_only.R renders the page from the plot bundle. This does the same, with
## two differences:
##
##   1. it joins the peptide-end model's 3-way insertion predictions onto
##      `pred_to_plot` by `peps`, which make_protein_plot_win then draws as a
##      stacked probability bar under each window and adds to the tooltip;
##   2. it sources THIS checkout's R/plot_proteins_win_new.R. plot_only.R
##      hardcodes ~/R_projects/ligandFinder/R/..., i.e. the main checkout, so a
##      worktree's plotting changes are invisible to it (see the path-shadowing
##      note in the repo's own comments).
##
## The join is left-only and additive: a window with no insertion prediction
## keeps NA and renders exactly as before, and
## options(lf.nn_show_ins = FALSE) turns the whole thing off.
##
##   Rscript inst/scripts/10_9m_plot_protein_ins.R NPY
##   ... --ins <csv> --out-dir <dir> --bundle <rds>
##
## The insertion csv is written by the step that reads the 20 members of
## lf_pepend_run_C_ins: one row per candidate window, keyed on `peps`.
##
## Output: <out-dir>/<GENE>.html
## ------------------------------------------------------------------------------
suppressMessages({ library(tidyverse); library(ggiraph); library(patchwork) })
suppressWarnings(suppressMessages(library(ligandFinder)))

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
bundle_p <- path.expand(.opt("--bundle", "~/AF2_analysis/plot_bundle.rds"))
ins_p    <- path.expand(.opt("--ins", "~/AF2_analysis/lf_pepend_ins_by_window_C.csv"))
out_dir  <- path.expand(.opt("--out-dir", "~/AF2_analysis/new_meth_plot_ins"))
## --windows REPLACES the bundle's window set with one 10_9j built from the scan,
## which already carries its own ins_* columns -- so no csv join happens at all.
win_p    <- .opt("--windows", NA_character_)
## positional args are gene symbols; everything after a -- flag is its value
genes <- local({
  a <- .args; drop <- integer(0)
  for (f in c("--bundle", "--ins", "--out-dir", "--windows")) {
    i <- match(f, a); if (!is.na(i)) drop <- c(drop, i, i + 1L)
  }
  if (length(drop)) a <- a[-drop]
  if (length(a)) a else "NPY"
})

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()

if (!file.exists(bundle_p)) stop(bundle_p, " does not exist; see plot_only.R", call. = FALSE)
if (!file.exists(ins_p))    stop(ins_p, " does not exist", call. = FALSE)

message("loading bundle ", bundle_p, " ...")
.b <- readRDS(bundle_p)
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
stopifnot(exists("the_input"), exists("pred_to_plot"), exists("pep_input"))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))

## ---- the windows to draw -----------------------------------------------------
if (!is.na(win_p)) {
  ## Scan windows from 10_9j: they already carry ins_* (when scored against an
  ## --ins-head run), so they go straight in.
  win_p <- path.expand(win_p)
  if (!file.exists(win_p)) stop(win_p, " does not exist", call. = FALSE)
  pred_to_plot <- readRDS(win_p)
  message(sprintf("window set <- %s (%d window(s), %s)", win_p, nrow(pred_to_plot),
                  if ("ins_class" %in% names(pred_to_plot))
                    sprintf("insertion head on %d", sum(!is.na(pred_to_plot$ins_class)))
                  else "no insertion columns"))
} else {

## ---- join the insertion head's call onto the windows -------------------------
ins <- readr::read_csv(ins_p, show_col_types = FALSE)
stopifnot(all(c("peps", "ins_class", "ins_p_inserting", "ins_p_loop",
                "ins_p_non_inserting") %in% names(ins)))
## `peps` is the window id in BOTH tables (GENE_w<start>-<end>), and it is unique
## in the insertion table -- assert rather than silently fan the page's rows out.
stopifnot(!anyDuplicated(ins$peps))

## drop any previous join before re-joining, so sourcing this twice in one
## session cannot produce ins_class.x / ins_class.y and lose the columns
pred_to_plot <- pred_to_plot %>%
  dplyr::select(-dplyr::any_of(c(names(ins)[names(ins) != "peps"],
                                 paste0(names(ins), ".x"), paste0(names(ins), ".y")))) %>%
  dplyr::left_join(ins, by = "peps")

.n_hit <- sum(!is.na(pred_to_plot$ins_class))
message(sprintf("insertion predictions joined: %d of %d windows (%.1f%%)",
                .n_hit, nrow(pred_to_plot), 100 * .n_hit / nrow(pred_to_plot)))
if (!.n_hit) stop("no window matched on `peps` -- are these the same window set?",
                  call. = FALSE)
}

source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
message("sourced ", file.path(ROOT, "R", "plot_proteins_win_new.R"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gs  <- purrr::map_chr(the_input, \(x) x$gene)
idx <- which(gs %in% genes)
if (!length(idx)) stop("none of those genes are in the bundle: ",
                       paste(utils::head(gs, 20), collapse = ", "), call. = FALSE)

for (i in idx) {
  message("--- ", gs[i], " ---")
  tm <- system.time(
    make_protein_plot_win(the_input[[i]], pred_to_plot, out_dir, pep_input[[i]])
  )[["elapsed"]]
  message(sprintf("    %s: %.1f s -> %s", gs[i], tm, file.path(out_dir, paste0(gs[i], ".html"))))
}
