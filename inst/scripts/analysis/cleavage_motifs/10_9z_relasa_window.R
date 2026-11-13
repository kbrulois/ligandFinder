#!/usr/bin/env Rscript
## ---- solvent exposure around every scanned anchor -------------------------------
## Mean relASA over a window centred on the anchor (default +/-3, so 7 residues),
## clipped at the protein's ends. Emitted in the scan's own row order, so it can
## be used as a mask over the scan with no further joining -- the same contract
## 10_9z_motif_codes.R follows.
##
##   Rscript inst/scripts/10_9z_relasa_window.R --term C --half 3
##
## Writes ~/AF2_analysis/lf_pepend_relasa_<term>.npz (asa_w, half).
## ------------------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
term   <- .opt("--term", "C")
half   <- as.integer(.opt("--half", "3"))
A2     <- path.expand(.opt("--a2", "~/AF2_analysis"))
scan_f <- .opt("--scan", file.path(A2, sprintf("lf_pepend_scanmem_%s_near.npz", term)))
if (!term %in% c("N", "C")) stop("--term must be N or C", call. = FALSE)

fc   <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
prec <- fc$prec
np   <- reticulate::import("numpy", convert = FALSE)
z    <- np$load(path.expand(scan_f))
prot <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
anc  <- as.integer(reticulate::py_to_r(z[["anchor"]]))

out    <- numeric(length(anc))
starts <- c(0L, which(diff(prot) != 0L)) + 1L
stops  <- c(which(diff(prot) != 0L), length(prot))
for (b in seq_along(starts)) {
  lo <- starts[b]; hi <- stops[b]; i <- prot[lo] + 1L
  v  <- fc$feats[[prec$accession[i]]][, "relASA"]
  L  <- length(v)
  ## cumulative sum so the window mean is O(1) per anchor; clip at the protein ends
  cs <- c(0, cumsum(v))
  a  <- anc[lo:hi]
  w0 <- pmax(1L, a - half); w1 <- pmin(L, a + half)
  out[lo:hi] <- (cs[w1 + 1L] - cs[w0]) / (w1 - w0 + 1L)
}
stopifnot(!any(is.na(out)))

np$savez_compressed(file.path(A2, sprintf("lf_pepend_relasa_%s.npz", term)),
                    asa_w = np$asarray(out), half = np$asarray(half))
message(sprintf("%s: %s anchors, mean relASA over +/-%d", term,
                format(length(out), big.mark = ","), half))
print(round(quantile(out, c(0, .05, .25, .5, .75, .9, .95, .99, 1)), 3))
for (t in c(0.2, 0.4, 0.6, 0.8))
  message(sprintf("  >= %.1f : %9s anchors (%.1f%%)", t,
                  format(sum(out >= t), big.mark = ","), 100 * mean(out >= t)))
