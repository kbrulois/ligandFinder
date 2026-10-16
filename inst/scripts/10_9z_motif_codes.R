#!/usr/bin/env Rscript
## ---- the cleavage motif at every scanned anchor ---------------------------------
## Classifies what sits just OUTSIDE each anchor, using the package's own
## lf_pepend_motif_vec, so the labels match the candidate tables exactly:
##
##   C: dibasic = K/R at +1,+2 | G + dibasic = G at +1 then K/R at +2,+3 (the
##      amidation motif) | monobasic = K/R at +1, or G at +1 with K/R at +2 |
##      terminus = the anchor IS the protein's last residue | other
##   N: the mirror, read at -1,-2,-3.
##
## Walks the scan's own prot_idx/anchor arrays, so the codes come back in the
## scan's row order and can be used as a mask over it with no further joining.
##
##   Rscript inst/scripts/10_9z_motif_codes.R --term C
##
## Writes ~/AF2_analysis/lf_pepend_motif_<term>.npz (motif, levels) and, for the
## python half, <tables>/prec_order.csv.
## ------------------------------------------------------------------------------
suppressMessages(library(dplyr))

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
term   <- .opt("--term", "C")
A2     <- path.expand(.opt("--a2", "~/AF2_analysis"))
scan_f <- .opt("--scan", file.path(A2, sprintf("lf_pepend_scanmem_%s_near.npz", term)))
tables <- .opt("--tables", NA_character_)
if (!term %in% c("N", "C")) stop("--term must be N or C", call. = FALSE)

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))

fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))$prec
np <- reticulate::import("numpy", convert = FALSE)
z    <- np$load(path.expand(scan_f))
prot <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
anc  <- as.integer(reticulate::py_to_r(z[["anchor"]]))

lv  <- c("dibasic", "G + dibasic", "monobasic", "terminus", "other")
out <- integer(length(anc))
## the scan emits each precursor's anchors as one contiguous ascending block
starts <- c(0L, which(diff(prot) != 0L)) + 1L
stops  <- c(which(diff(prot) != 0L), length(prot))
for (b in seq_along(starts)) {
  lo <- starts[b]; hi <- stops[b]; i <- prot[lo] + 1L
  m <- lf_pepend_motif_vec(fc$seq[i], term, anc[lo:hi], fc$n_prot[i])
  out[lo:hi] <- match(m, lv)
  if (b == 1L) stopifnot(identical(anc[lo:hi], seq.int(fc$n_prot[i], fc$c_prot[i])))
}
stopifnot(!any(is.na(out)))

np$savez_compressed(file.path(A2, sprintf("lf_pepend_motif_%s.npz", term)),
                    motif = np$asarray(as.integer(out)), levels = np$asarray(lv))
if (!is.na(tables))
  utils::write.csv(fc[, c("accession", "gene", "n_prot", "c_prot")],
                   file.path(path.expand(tables), "prec_order.csv"), row.names = FALSE)
message(sprintf("%s: %s anchors", term, format(length(out), big.mark = ",")))
print(setNames(as.integer(table(factor(lv[out], levels = lv))), lv))
