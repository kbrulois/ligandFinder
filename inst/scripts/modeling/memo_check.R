## Memorisation check (Kevin, 2026-10-22: "is that because the model overfits").
## Every seen end was trained on by 4 of the 5 folds and held out by one. Score
## and call it on each kind of profile, one fold (10 models) at a time:
##   in-sample   the 4 folds that trained on it
##   held out    the fold that did not
## and the never-seen ends on every fold (all out-of-sample), so all three are
## single-fold numbers. Same call rule as end_calls.R (T = 0.2, +/-2).
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
TOL <- 2L; thr <- 0.2
segments <- function(b, s, thr) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  d <- data.table(sid, accession = b$accession, anchor = b$anchor, s = s, up = up)
  d[, .(up = up[1], centre = anchor[which.max(s)]), by = .(sid, accession)][up == TRUE]
}
caught <- function(pk, ev) vapply(seq_len(nrow(ev)), function(r)
  any(pk$accession == ev$accession[r] & abs(pk$centre - ev$anchor[r]) <= TOL), logical(1))

m <- fread(file.path(S, "end_windows_meta.csv"))
out <- list()
for (tm in c("C", "N")) {
  b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm)))
  setDT(b); setorder(b, accession, anchor)
  b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
  ev <- unique(m[term == tm, .(set, accession, anchor, fold)])
  sc <- b[ev, on = .(accession, anchor), mult = "first"]
  for (k in 1:5) {
    pk <- segments(b, b[[sprintf("cnn__f%d", k)]], thr)
    hit <- caught(pk, ev)
    out[[length(out) + 1]] <- data.table(term = tm, ev[, .(set, accession, anchor, fold)], k = k,
                                         score = sc[[sprintf("cnn__f%d", k)]], tp = hit)
  }
  rm(b); invisible(gc())
}
R <- rbindlist(out)
R[, kind := fifelse(set != "seen", "never seen (every fold out-of-sample)",
                    fifelse(fold == k, "seen, held-out fold", "seen, in-sample folds"))]
R[set != "seen", kind := paste0(set, " (every fold)")]
res <- R[, .(ends = uniqueN(paste(accession, anchor)), median_score = median(score, na.rm = TRUE),
             tp_rate = mean(tp)), by = .(term, kind)][order(term, kind)]
fwrite(res, file.path(A2, "lf_headcmp", "memo_check.csv"))
options(width = 160); print(res, row.names = FALSE)
