#!/usr/bin/env Rscript
## ---- the G + dibasic C ends, before and after --amid-g-anchor -----------------
## Kevin, 2026-11-12: "new plot of aligned known peptide windows showing new
## positioning of G dibasic". 10_7b draws every window of one build; this puts
## the moved C ends of a G-anchored build beside their old-convention windows:
## the old one (mature end at 28, G at 29, basic pair at 30-31) next to the new
## one (G at 28, pair at 29-30), and a position-frequency strip of G and K/R
## around the anchor for the moved ends (old and new anchor) and the rest.
##
## The two builds are paired by accession + mature end, and lf_gshift_pairs()
## stops unless they differ by exactly the shift and the merge it implies --
## so the counts on the figure are checked, not just reported.
##
##   Rscript inst/scripts/10_7c_plot_gshift_windows.R \
##     --new ~/AF2_analysis/lf_pepend_nn_input_cv5_f1_padg.rds \
##     --old ~/AF2_analysis/lf_pepend_nn_input_cv5_f1.rds --out-dir ~/AF2_analysis/lf_viz_20261112
##
## Any fold will do: a build's five folds hold the same knowns, only `split`
## differs, and this figure does not draw the split.
##
## Output: <out-dir>/pepend_gshift_windows_<stamp>.{png,svg}
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
new_path <- path.expand(.opt("--new", "~/AF2_analysis/lf_pepend_nn_input_cv5_f1_padg.rds"))
old_path <- path.expand(.opt("--old", "~/AF2_analysis/lf_pepend_nn_input_cv5_f1.rds"))
out_dir  <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
stamp    <- format(Sys.time(), "%Y%m%d_%H%M%S")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
## nn_class_cols / nn_class_labels: the shared class palette and figure names
source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
source(file.path(ROOT, "R", "viz_theme.R"))
source(file.path(ROOT, "R", "pepend_gshift_plots.R"))

## the knowns only: a fold input is ~0.8 GB in memory, so one at a time
knowns_of <- function(f) {
  x <- readRDS(f); out <- list(knowns = x$knowns, params = x$params)
  rm(x); invisible(gc()); out
}
nw <- knowns_of(new_path)
if (!isTRUE(nw$params$amid_g)) stop(basename(new_path), " was not built with --amid-g-anchor")
od <- knowns_of(old_path)
if (isTRUE(od$params$amid_g)) stop(basename(old_path), " is G-anchored too; --old needs the old convention")

gp <- lf_gshift_pairs(nw$knowns, od$knowns)
cn <- gp$counts
message(sprintf(paste0("C ends: %d old -> %d new | %d moved onto the G (%d peptides), %d unchanged, ",
                       "%d merged (%s) | N windows: %d of %d g_amid relabelled at the G"),
                cn$old_C, cn$new_C, cn$moved, cn$peptides, cn$unchanged, cn$merged,
                paste(sprintf("%s -> %s", gp$merged$id, gp$merged$into), collapse = ", "),
                cn$N_relabelled, cn$N_gamid))
if (nrow(gp$regrouped))
  message("10_7b insertion group changed: ",
          paste(sprintf("%s %s %s %s -> %s", gp$regrouped$term, gp$regrouped$gene, gp$regrouped$name,
                        gp$regrouped$old, gp$regrouped$new), collapse = "; "))

p <- lf_gshift_figure(gp, old_name = basename(old_path), new_name = basename(new_path))
f <- file.path(out_dir, sprintf("pepend_gshift_windows_%s", stamp))
## (the fixed part is the title block, the strip and the legend; it grew with
## the regrouping line and the caption's strip note)
h <- 5.1 + 0.2 * (nrow(gp$pairs) + nrow(gp$merged))
ggsave(paste0(f, ".png"), p, width = 15, height = h, dpi = 200, bg = LF_VIZ$surface)
ggsave(paste0(f, ".svg"), p, width = 15, height = h, bg = LF_VIZ$surface)
message("wrote ", f, ".{png,svg}")
