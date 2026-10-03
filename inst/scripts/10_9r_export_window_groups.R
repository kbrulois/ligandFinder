#!/usr/bin/env Rscript
## ---- window identities for the exported arrays --------------------------------
## `<run>_isolated/in/arrays.npz` holds x and y only, so anything that needs to
## know WHICH window a row is -- an inner split that must not put two windows of
## the same precursor on opposite sides, say -- cannot get it from there. This
## writes the accession and peps for each split, in the array row order.
##
##   Rscript inst/scripts/10_9r_export_window_groups.R --term C
## Output: ~/AF2_analysis/lf_pepend_window_groups_<term>.npz
## ------------------------------------------------------------------------------
suppressMessages(library(dplyr))
.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
term   <- toupper(.opt("--term", "C"))
in_p   <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
out_p  <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_window_groups_%s.npz", term)))
np <- reticulate::import("numpy", convert = FALSE)

x <- readRDS(in_p)$nn_input[[term]]
args <- list(out_p)
for (s in c("train", "val", "all")) {
  d <- x[[s]]
  args[[paste0(term, "__", s, "__accession")]] <- np$asarray(as.character(d$accession))
  args[[paste0(term, "__", s, "__peps")]]      <- np$asarray(as.character(d$peps))
  args[[paste0(term, "__", s, "__known")]]     <- np$asarray(as.integer(d$known))
  ## the anchor locates the window's defining residue in a per-residue table
  ## such as the PLM embeddings, which are keyed (accession, index)
  args[[paste0(term, "__", s, "__anchor")]]    <- np$asarray(as.integer(d$anchor))
  message(sprintf("%-5s %6d rows, %d precursors, %d positives",
                  s, nrow(d), dplyr::n_distinct(d$accession), sum(d$known == 1)))
}
do.call(np$savez_compressed, args)
message("\ngroups -> ", out_p)
