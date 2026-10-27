## Per-END calls at one threshold (Kevin, 2026-10-22: "annotate TP and FN
## unseens at a threshold of 0.2"). Same rule as segments_roc.R v2:
##   TP  a run of anchors >= T whose maximum is within +/-2 of the end
##   FN  otherwise -- "below T" when the end's own run is < T, "off-centre"
##       when the end sits inside a run >= T whose maximum lands > 2 away
## Never-seen ends: the all-50 profile (what the page's heavy squares show).
## Seen (cross-validated) ends: the profile of the fold that held each one out.
## Both heads. Writes lf_headcmp/end_calls_T<thr>.csv.
##   THRESH=0.2 Rscript end_calls.R
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
## ARM selects which arm's tracks to call on and suffixes the output; TOL is
## no longer fixed at 2. TOL=2 with ARM unset reproduces the original exactly.
TOL <- as.integer(Sys.getenv("TOL", "4"))
ARM <- Sys.getenv("ARM", "")
thr <- as.numeric(Sys.getenv("THRESH", "0.2"))

segments <- function(b, s, thr) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  d <- data.table(sid, accession = b$accession, anchor = b$anchor, s = s, up = up)
  seg <- d[, .(start = anchor[1], end = anchor[.N], up = up[1], height = max(s),
               centre = anchor[which.max(s)]), by = .(sid, accession)]
  list(seg = seg, sid = sid)
}
## one row per end: call, reason, the end's score, the nearest peak
call_ends <- function(b, s, ev) {
  sg <- segments(b, s, thr); seg <- sg$seg; pk <- seg[up == TRUE]
  key <- data.table(accession = b$accession, anchor = b$anchor, sid = sg$sid, s = s)
  e <- key[ev, on = .(accession, anchor)]
  out <- rbindlist(lapply(seq_len(nrow(e)), function(r) {
    p <- pk[accession == e$accession[r]]
    p <- p[order(abs(centre - e$anchor[r]))]
    hit <- nrow(p) && abs(p$centre[1] - e$anchor[r]) <= TOL
    own <- seg[sid == e$sid[r]]
    data.table(call = if (hit) "TP" else "FN",
               reason = if (hit) "peak on the end" else if (is.na(e$sid[r])) "unscored" else
                 if (isTRUE(own$up)) "off-centre peak" else "below T",
               score = e$s[r],
               peak_centre = if (nrow(p)) p$centre[1] else NA_integer_,
               peak_height = if (nrow(p)) p$height[1] else NA_real_)
  }))
  cbind(ev, out)
}

un <- fread(file.path(S, "cv5_unseen_ends.csv"))
un[, set := ifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]
out <- list()
for (tm in c("C", "N")) {
  b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, ARM)))
  setDT(b); setorder(b, accession, anchor)
  b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  kn[, gene := sub("_w-?[0-9]+--?[0-9]+$", "", peps)]
  for (hd in c("cnn", "xgb_pm")) {
    ev <- unique(un[term == tm, .(set, gene, accession, name, len, term, anchor)])
    out[[length(out) + 1]] <- cbind(head = hd, call_ends(b, b[[paste0(hd, "__all50")]], ev))
    for (k in 1:5) {
      ev <- unique(kn[fold == k, .(set = "cross-validated", gene, accession, name = peps, len = NA_integer_,
                                    term = tm, anchor)])
      out[[length(out) + 1]] <- cbind(head = hd, fold = k, call_ends(b, b[[sprintf("%s__f%d", hd, k)]], ev))
    }
  }
  rm(b); invisible(gc())
}
R <- rbindlist(out, fill = TRUE)
R[, thr := thr]
SUF <- paste0(ARM, if (TOL == 2L) "" else paste0("_tol", TOL))
fwrite(R, file.path(A2, "lf_headcmp", sprintf("end_calls_T%s%s.csv", thr, SUF)))
print(R[, .(ends = .N, TP = sum(call == "TP"), below_T = sum(reason == "below T"),
            off_centre = sum(reason == "off-centre peak"), unscored = sum(reason == "unscored")),
        by = .(set, term, head)][order(set, term, head)], row.names = FALSE)
