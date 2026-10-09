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
## A run trained with 10_8 --positives has FEWER rows than the full nn_input
## (inserting-only drops train 457 -> 429, val 429 -> 414), so its arrays.npz and
## this table only line up if the same filter is applied here. Mismatched lengths
## are an error downstream rather than a silent misalignment, but only because
## something checks -- so pass the same --positives you trained with.
pos_cls <- .opt("--positives", NA_character_)
out_p  <- path.expand(.opt("--out",
  sprintf("~/AF2_analysis/lf_pepend_window_groups_%s%s.npz", term,
          if (!is.na(pos_cls)) paste0("_pos", pos_cls) else "")))
np <- reticulate::import("numpy", convert = FALSE)

x <- readRDS(in_p)$nn_input[[term]]
if (!is.na(pos_cls)) {
  source(file.path(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), "..", "..",
    "R", "pepend_windows.R"))
  for (s_ in c("train", "val")) {
    d <- x[[s_]]; pos <- which(d$known == 1L)
    cls <- lf_pepend_insertion_class(d$known_idx[pos], d$term[pos], d$len[pos])
    drop <- pos[cls != pos_cls]
    if (length(drop)) x[[s_]] <- d[-drop, ]
  }
  message("filtered positives to '", pos_cls, "'")
}
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
