## The never-seen set for arms trained on the _padg windows (Kevin, 2026-11-10).
##
## cv5_unseen_ends.csv holds the annotated ends of genes no fold trained on
## whose own 36-mer no TRAINING-NEGATIVE window overlaps (cv5_sets.R's rule).
## The _padg inputs add negatives -- near precursor termini and at monobasic
## sites -- and redraw some C negatives, so an end that was clean may not be.
## This keeps the old set's rule and drops every end overlapped by a negative
## of any of the five _padg inputs (train + val, both termini), testing the
## end's window at its mature anchor AND at its G anchor (arm_anchor), since a
## padg arm is scored at the latter.
##
## The result is a SUBSET of the old set, so it is the set to score EVERY arm
## on when a padg arm is in the comparison: UNSEEN_SFX=_padg.
##
##   Rscript cv5_unseen_padg.R            -> cv5_unseen_ends_padg.csv
suppressMessages({library(data.table)})
S  <- path.expand("~/AF2_analysis/lf_nomask_scripts_20261021")   # the end tables (data)
A2 <- path.expand("~/AF2_analysis")
SFX <- Sys.getenv("INPUT_SFX", "_padg")
## arm_conv.R is code, beside this script (LF_SCRIPTS overrides the directory)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(Sys.getenv("LF_SCRIPTS", .lf_here), "arm_conv.R"))
W  <- 36L
start_of <- function(term, anchor) ifelse(term == "N", anchor - 7L, anchor - 27L)   # LF_PEPEND anchors 8 / 28

negs <- rbindlist(lapply(1:5, function(k) {
  x <- readRDS(file.path(A2, sprintf("lf_pepend_nn_input_cv5_f%d%s.rds", k, SFX)))
  stopifnot(isTRUE(x$params$amid_g))
  rbindlist(lapply(c("N", "C"), function(tm) rbindlist(lapply(c("train", "val"), function(sp) {
    d <- x$nn_input[[tm]][[sp]]; i <- d$known == 0
    data.table(accession = d$accession[i], lo = start_of(tm, d$anchor[i]))
  }))))
}))
negs <- unique(negs)[, hi := lo + W - 1L]; setkey(negs, accession)
cat(sprintf("negative windows over the five %s inputs: %d (%d precursors)\n",
            SFX, nrow(negs), uniqueN(negs$accession)))

un <- fread(file.path(S, "cv5_unseen_ends.csv"))
hit <- function(acc, term, anchor) {
  n <- negs[.(acc), nomatch = 0L]
  if (!nrow(n)) return(FALSE)
  lo <- start_of(term, anchor); any(n$lo <= lo + W - 1L & lo <= n$hi)
}
g <- arm_anchor("curated_AAswap_padg", un$accession, un$term, un$anchor)
un[, dirty_padg := mapply(hit, accession, term, anchor) | mapply(hit, accession, term, g)]
cat(sprintf("old set %d ends; %d now overlapped -> %d kept\n", nrow(un), sum(un$dirty_padg), sum(!un$dirty_padg)))
print(un[, .(old = .N, kept = sum(!dirty_padg)), by = .(set, term)])
out <- un[dirty_padg == FALSE][, dirty_padg := NULL]
fwrite(out, file.path(S, sprintf("cv5_unseen_ends%s.csv", SFX)))
cat("wrote", file.path(S, sprintf("cv5_unseen_ends%s.csv", SFX)), "\n")
