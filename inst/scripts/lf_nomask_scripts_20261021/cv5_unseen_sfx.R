## The never-seen set for arms trained on their own windows (Kevin, 2026-11-10;
## any window set since 2026-11-12, when this was cv5_unseen_padg.R).
##
## cv5_unseen_ends.csv holds the annotated ends of genes no fold trained on
## whose own 36-mer no TRAINING-NEGATIVE window overlaps (cv5_sets.R's rule).
## A window set built with other negatives (INPUT_SFX: lf_pepend_nn_input_cv5_
## f<k><sfx>.rds, e.g. the --pad-match / --mono-match top-ups near precursor
## termini and at monobasic sites) may overlap an end that was clean. This keeps
## the old set's rule and drops every end whose window, at the anchor the end
## table carries, overlaps a negative of any of the five inputs (train + val,
## both termini).
##
## The result is a SUBSET of the old set, so it is the set to score EVERY arm
## on when an arm on those windows is in the comparison: UNSEEN_SFX=<sfx>.
##
## Kevin, 2026-11-12: it also checks the set's POSITIVES. Every evaluator scores
## an arm at the known ends, anchors and held-out folds of the shared tables
## (lf_pepend_cv5_knowns_<T>.csv); none re-reads an arm's own folds any more,
## and the rule above assumes the set trains on the same genes. So the ends the
## five inputs hold out, and in which fold, must BE that table, or this stops:
## an arm on a set that differs would be graded on ends its folds trained on.
## KNOWNS_DIFFER_OK=1 writes the set anyway -- the _padg set, whose 30 G +
## dibasic C ends are anchored one past the mature end and whose two GLP-1 C
## ends are one, and whose arms are scored at the shared anchors since 2026-11-12.
##
##   cd <scripts dir> && INPUT_SFX=_msa LC_ALL=en_US.UTF-8 Rscript cv5_unseen_sfx.R
##     -> $S/cv5_unseen_ends_msa.csv   (OUT overrides the output path)
suppressMessages({library(data.table)})
S  <- path.expand("~/AF2_analysis/lf_nomask_scripts_20261021")   # the end tables (data)
A2 <- path.expand("~/AF2_analysis")
## required: "" would read the shared inputs and write over cv5_unseen_ends.csv
SFX <- Sys.getenv("INPUT_SFX", "")
if (!nzchar(SFX)) stop("set INPUT_SFX, the window set (e.g. INPUT_SFX=_msa)", call. = FALSE)
OUT <- path.expand(Sys.getenv("OUT", file.path(S, sprintf("cv5_unseen_ends%s.csv", SFX))))
W  <- 36L
start_of <- function(term, anchor) ifelse(term == "N", anchor - 7L, anchor - 27L)   # LF_PEPEND anchors 8 / 28

ins <- lapply(1:5, function(k) {
  x <- readRDS(file.path(A2, sprintf("lf_pepend_nn_input_cv5_f%d%s.rds", k, SFX)))
  ## five folds of ONE build: a stray input from another set would drop the wrong ends
  if (!identical(as.integer(x$params$n_folds), 5L) || !identical(as.integer(x$params$fold), k))
    stop(sprintf("lf_pepend_nn_input_cv5_f%d%s.rds is fold %s of %s, not fold %d of 5", k, SFX,
                 format(x$params$fold), format(x$params$n_folds), k), call. = FALSE)
  held <- as.data.table(x$knowns)[split == "val", .(term, accession, anchor = as.integer(anchor))]
  held[, fold := k]
  list(held = held,
       negs = rbindlist(lapply(c("N", "C"), function(tm) rbindlist(lapply(c("train", "val"), function(sp) {
         d <- x$nn_input[[tm]][[sp]]; i <- d$known == 0
         data.table(accession = d$accession[i], lo = start_of(tm, d$anchor[i]))
       })))))
})
held <- rbindlist(lapply(ins, `[[`, "held"))
negs <- rbindlist(lapply(ins, `[[`, "negs")); rm(ins)
negs <- unique(negs)[, hi := lo + W - 1L]; setkey(negs, accession)
cat(sprintf("negative windows over the five %s inputs: %d (%d precursors)\n",
            SFX, nrow(negs), uniqueN(negs$accession)))

## the set's positives: the ends its five inputs hold out, and in which fold,
## must be the shared tables' (see the header)
ref <- rbindlist(lapply(c("N", "C"), function(tm)
  fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))[
    , .(term = tm, accession, anchor = as.integer(anchor), fold = as.integer(fold))]))
if (anyDuplicated(held, by = c("term", "accession", "anchor")))
  stop("a known end is held out by more than one of the five ", SFX, " inputs", call. = FALSE)
cmp <- merge(ref, held, by = c("term", "accession", "anchor"), all = TRUE, suffixes = c("", "_sfx"))
bad <- cmp[is.na(fold) | is.na(fold_sfx) | fold != fold_sfx]
cat(sprintf("known ends held out: %d in lf_pepend_cv5_knowns_{N,C}.csv, %d over the five %s inputs; %d differ\n",
            nrow(ref), nrow(held), SFX, nrow(bad)))
if (nrow(bad)) {
  print(bad, nrows = 40L)
  why <- sprintf(paste0("%d known ends of the %s inputs are not the shared tables' (%d only in ",
                        "lf_pepend_cv5_knowns_<T>.csv, %d only in the inputs, %d in another fold), and ",
                        "every evaluator scores an arm at the tables' ends and folds"),
                 nrow(bad), SFX, sum(is.na(bad$fold_sfx)), sum(is.na(bad$fold)),
                 sum(!is.na(bad$fold) & !is.na(bad$fold_sfx)))
  if (!identical(Sys.getenv("KNOWNS_DIFFER_OK"), "1"))
    stop(why, "; KNOWNS_DIFFER_OK=1 writes the never-seen set anyway", call. = FALSE)
  message("KNOWNS_DIFFER_OK=1, writing anyway: ", why)
}

un <- fread(file.path(S, "cv5_unseen_ends.csv"))
hit <- function(acc, term, anchor) {
  n <- negs[.(acc), nomatch = 0L]
  if (!nrow(n)) return(FALSE)
  lo <- start_of(term, anchor); any(n$lo <= lo + W - 1L & lo <= n$hi)
}
un[, dirty_sfx := mapply(hit, accession, term, anchor)]
cat(sprintf("old set %d ends; %d now overlapped -> %d kept\n", nrow(un), sum(un$dirty_sfx), sum(!un$dirty_sfx)))
print(un[, .(old = .N, kept = sum(!dirty_sfx)), by = .(set, term)])
out <- un[dirty_sfx == FALSE][, dirty_sfx := NULL]
fwrite(out, OUT)
cat("wrote", OUT, "\n")
