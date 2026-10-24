## Threshold as a classifier over SEGMENTS (Kevin, 2026-10-22).
##
## At threshold T a precursor's anchor profile splits into alternating runs:
##   positive segment  a run of consecutive anchors with score >= T
##                     TP if its centre (the run's maximum) is within +/-TOL of
##                     an evaluation end, FP otherwise
##   negative segment  a run with score < T; TN if it holds no real end
## Ends are counted per END (v2, 2026-10-22), so TP + FN = the ends:
##   TP  an end with a positive segment's centre within TOL
##   FN  every other end -- in a negative segment, inside a positive segment
##       whose maximum lands more than TOL away (that segment is also an FP, as
##       a mislocalised detection is in object detection), or unscored
## A negative segment holding an evaluation end is neither TN nor FN (its end
## is already counted, as a TP or an FN). Segments touching an annotated end
## from OUTSIDE the evaluation set (other folds' train ends, ligand_list,
## UniProt peptide/propeptide, the other never-seen set) are held out of every
## count, as in peaks_final.R.
## v1 (segment-level FN; an end inside an off-centre peak was not an FN) is in
## lf_headcmp/segment_*_v1.csv.
##
## Each T gives one (FPR, TPR) point. Per T: the single-point ROC AUC
## (1 + TPR - FPR) / 2 -- the area under the curve through (0,0), that point
## and (1,1), i.e. balanced accuracy. Over T: the curve traced by the grid.
## At low T peaks merge into a few long segments, so TPR falls again as T -> 0
## and the path bends back; it is reported along the path (closed flat to
## FPR = 1) and as its upper envelope (best TPR reachable at each FPR).
## Cross-validated: fold k's ends on fold k's own profile, counts summed over
## folds; never seen: the all-50 profile.
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7"
TOL <- 2L
GRID <- sort(unique(round(c(seq(0, 1, by = 0.025), seq(0.1, 0.5, by = 0.05)), 3)))

## segments of one profile at threshold thr; `b` is sorted by accession, anchor
## and carries `brk0` (TRUE where a new precursor or an anchor gap starts)
segments <- function(b, s, thr) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  d <- data.table(sid, accession = b$accession, anchor = b$anchor, s = s, up = up)
  seg <- d[, .(start = anchor[1], end = anchor[.N], up = up[1],
               centre = anchor[which.max(s)]), by = .(sid, accession)]
  list(seg = seg, sid = sid)
}
## segment id of each end (ends outside the scanned range get NA)
end_sid <- function(b, sid, ends) {
  key <- data.table(accession = b$accession, anchor = b$anchor, sid = sid)
  key[ends, on = .(accession, anchor), sid]
}
## TRUE where a positive segment's centre is within TOL of any end in `ends`
centre_near <- function(seg, ends) {
  if (!nrow(seg) || !nrow(ends)) return(rep(FALSE, nrow(seg)))
  e <- unique(ends[, .(accession, anchor)]); setkey(e, accession, anchor)
  q <- data.table(accession = seg$accession, anchor = seg$centre, i = seq_len(nrow(seg)))
  j <- e[q, on = .(accession, anchor), roll = "nearest", .(i, hit = x.anchor)]
  out <- abs(j$hit - q$anchor) <= TOL; out[is.na(out)] <- FALSE
  out[order(j$i)]
}
confusion <- function(b, s, thr, ev, oth) {
  sg <- segments(b, s, thr); seg <- sg$seg
  pos <- seg$up; neg <- !pos
  tp <- pos & centre_near(seg, ev)
  ## which evaluation ends a TP segment caught (centre within TOL)
  tps <- seg[tp]
  caught <- if (nrow(tps)) vapply(seq_len(nrow(ev)), function(r)
    any(tps$accession == ev$accession[r] & abs(tps$centre - ev$anchor[r]) <= TOL), logical(1))
    else rep(FALSE, nrow(ev))
  ev_sid  <- end_sid(b, sg$sid, ev)
  oth_sid <- end_sid(b, sg$sid, oth)
  ev_up   <- seg$up[match(ev_sid, seg$sid)]          # NA: the end has no scored anchor
  miss    <- !caught
  held_pos <- pos & !tp & centre_near(seg, oth)
  held_neg <- neg & seg$sid %in% c(ev_sid, oth_sid)
  c(TP = sum(caught), FN = sum(miss), FP = sum(pos & !tp & !held_pos), TN = sum(neg & !held_neg),
    TP_seg = sum(tp),
    FN_neg = sum(miss & !is.na(ev_up) & !ev_up),      # missed: in a negative segment
    FN_peak = sum(miss & !is.na(ev_up) & ev_up),      # missed: inside an off-centre positive segment
    FN_unscored = sum(miss & is.na(ev_up)),
    ends = nrow(ev))
}

ll <- as.data.table(readRDS(file.path(W, "inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
un <- fread(file.path(S, "cv5_unseen_ends.csv"))
un[, set := ifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]

out <- list()
for (tm in c("C", "N")) {
  b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm)))
  setDT(b); setorder(b, accession, anchor)
  b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))[, .(accession, anchor, fold)]
  ann <- unique(rbind(ll[, .(accession, anchor = if (tm == "N") start else end)],
                      up[, .(accession, anchor = if (tm == "N") start else end)]))
  for (hd in c("cnn", "xgb_pm")) {
    for (thr in GRID) {
      ## cross-validated: per fold, summed
      cv <- Reduce(`+`, lapply(1:5, function(k) {
        ev  <- kn[fold == k, .(accession, anchor)]
        oth <- unique(rbind(ann, kn[fold != k, .(accession, anchor)], un[, .(accession, anchor)]))
        confusion(b, b[[sprintf("%s__f%d", hd, k)]], thr, ev, oth)
      }))
      out[[length(out) + 1]] <- data.table(term = tm, set = "cross-validated", head = hd, thr = thr, t(cv))
      for (st in c("never seen: peptides", "never seen: large products")) {
        ev  <- unique(un[set == st & term == tm, .(accession, anchor)])
        oth <- unique(rbind(ann, kn[, .(accession, anchor)], un[set != st, .(accession, anchor)]))
        cf <- confusion(b, b[[paste0(hd, "__all50")]], thr, ev, oth)
        out[[length(out) + 1]] <- data.table(term = tm, set = st, head = hd, thr = thr, t(cf))
      }
    }
    message(sprintf("%s %-6s done (%d thresholds)", tm, hd, length(GRID)))
  }
  rm(b); invisible(gc())
}
R <- rbindlist(out)
stopifnot(R[, all(TP + FN == ends)])
R[, `:=`(TPR = TP / (TP + FN), FPR = FP / (FP + TN), precision = TP_seg / (TP_seg + FP))]
R[, auc_T := (1 + TPR - FPR) / 2]
## AUC over T. Path: threshold descending from (0,0) at T = 1, closed flat to
## FPR = 1 at the last point's TPR. Envelope: the best TPR reachable at each FPR.
trap <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
auc <- R[order(-thr), {
  ok <- is.finite(FPR) & is.finite(TPR); x <- FPR[ok]; y <- TPR[ok]
  px <- c(0, x, 1); py <- c(0, y, tail(y, 1))
  o <- order(x, -y); ex <- c(0, x[o], 1); ey <- cummax(c(0, y[o], 0)); ey[length(ey)] <- max(y)
  .(roc_auc_path = trap(px, py), roc_auc_env = trap(ex, ey),
    monotone = all(diff(x) >= -1e-12) && all(diff(y) >= -1e-12),
    tpr_max = max(y), thr_tpr_max = thr[ok][which.max(y)],
    auc_T_max = max(auc_T[ok]), thr_auc_T_max = thr[ok][which.max(auc_T[ok])])
}, by = .(term, set, head)]
fwrite(R, file.path(A2, "lf_headcmp", "segment_confusion.csv"))
fwrite(auc, file.path(A2, "lf_headcmp", "segment_roc_auc.csv"))
options(width = 220)
print(auc[order(set, term, head)], row.names = FALSE)
print(R[thr %in% c(0, 0.025, 0.05, 0.1, 0.3, 0.5) & set == "cross-validated",
        .(term, head, thr, TP, FN, FN_neg, FN_peak, FN_unscored, FP, TN, TPR = round(TPR, 3), FPR = round(FPR, 4),
          precision = round(precision, 4), auc_T = round(auc_T, 3))],
      row.names = FALSE)
