## Peaks, not anchors: score the final window model over PEAKS of its profile.
##
## A peak = a run of consecutive anchors on one precursor with score >= T; its
## height is the run's maximum and its CENTRE the anchor where that maximum sits
## (HCRT's C peak over 64-68 has centre 66). T = 0.3 and 0.5.
## A peak is POSITIVE if its centre is within +/-TOL (2) of an end in the
## evaluation set, NEGATIVE otherwise -- except a peak near a real end from
## OUTSIDE the set (other folds' train ends, ligand_list, UniProt peptide /
## propeptide, the other never-seen set), which is held out: it is neither a
## clean positive nor an honest negative (10_9z_locus_pr's rule).
## Sets and scores as in metrics_final.R: cross-validated = the 5-fold knowns,
## fold k's ends scored on fold k's OWN profile; never seen = all 50.
##
## Reported per set x terminus x head x T:
##   peaks         the new search space (non-held-out peaks over the proteome)
##   ceiling       ends that have ANY peak within TOL (an end below T is missed)
##   ROC / AP / norm PR   over peaks, ranked by height
##   peaks for 50% / 80% of ends  each end's rank = 1 + negative peaks above its
##                 best matching peak (pooled over folds) -- the budget curve,
##                 which 10_9z found comparable across peak definitions; AP is
##                 not (it moves with how many peaks a definition emits)
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7"
TOL <- 2L
## thresholds: THRESH="0.1,0.15,..." in the environment, else 0.3 and 0.5
THRESH <- as.numeric(strsplit(Sys.getenv("THRESH", "0.3,0.5"), ",")[[1]])

fast_auc <- function(score, y) {
  r <- frank(score, ties.method = "average"); n1 <- sum(y); n0 <- length(y) - n1
  if (!n1 || !n0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fast_ap <- function(s, y) { if (!sum(y)) return(NA_real_); o <- order(-s); yy <- y[o]
  sum((cumsum(yy) / seq_along(yy))[yy == 1]) / sum(y) }
fast_ap_at <- function(s, y, target = 0.01) {
  if (!sum(y)) return(NA_real_)
  o <- order(-s); yy <- y[o]; n1 <- sum(y); n0 <- length(y) - n1
  w <- n1 * (1 - target) / (target * n0); tp <- cumsum(yy); fp <- cumsum(1 - yy)
  sum((tp / (tp + w * fp))[yy == 1]) / n1
}

## peaks of one profile: runs of consecutive anchors with s >= thr
call_peaks <- function(b, s, thr) {
  up <- !is.na(s) & s >= thr
  i <- which(up)
  if (!length(i)) return(data.table())
  brk <- c(TRUE, diff(i) != 1L | b$accession[i[-1]] != b$accession[i[-length(i)]] |
                   diff(b$anchor[i]) != 1L)
  run <- cumsum(brk)
  d <- data.table(run, accession = b$accession[i], anchor = b$anchor[i], s = s[i])
  d[, .(start = min(anchor), end = max(anchor), height = max(s), centre = anchor[which.max(s)]),
    by = .(run, accession)]
}
## distance from each peak centre to the nearest end of `ends` on its precursor
near <- function(pk, ends) {
  if (!nrow(pk) || !nrow(ends)) return(rep(FALSE, nrow(pk)))
  e <- unique(ends[, .(accession, anchor)]); setkey(e, accession, anchor)
  q <- data.table(accession = pk$accession, anchor = pk$centre, i = seq_len(nrow(pk)))
  j <- e[q, on = .(accession, anchor), roll = "nearest", .(i, hit = x.anchor)]
  out <- abs(j$hit - q$anchor) <= TOL; out[is.na(out)] <- FALSE
  out[order(j$i)]
}
## for each end: the height of its best peak within TOL (NA if none)
best_peak <- function(pk, ends) {
  if (!nrow(ends)) return(numeric(0))
  vapply(seq_len(nrow(ends)), function(r) {
    h <- pk$height[pk$accession == ends$accession[r] & abs(pk$centre - ends$anchor[r]) <= TOL]
    if (length(h)) max(h) else NA_real_
  }, numeric(1))
}

ll <- as.data.table(readRDS(file.path(W, "inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
un <- fread(file.path(S, "cv5_unseen_ends.csv"))
un[, set := ifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]

res <- list(); curves <- list(); hcrt <- list()
for (tm in c("C", "N")) {
  b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm)))
  setDT(b); setorder(b, accession, anchor)
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))[, .(accession, anchor, fold)]
  ann <- unique(rbind(ll[, .(accession, anchor = if (tm == "N") start else end)],
                      up[, .(accession, anchor = if (tm == "N") start else end)]))
  for (hd in c("cnn", "xgb_pm")) for (thr in THRESH) {
    ## ---- cross-validated: fold k's ends on fold k's own profile -----------
    fm <- list(); ranks <- numeric(0); nfound <- 0L; npk <- 0
    for (k in 1:5) {
      pk <- call_peaks(b, b[[sprintf("%s__f%d", hd, k)]], thr)
      ev <- kn[fold == k, .(accession, anchor)]
      oth <- unique(rbind(ann, kn[fold != k, .(accession, anchor)], un[, .(accession, anchor)]))
      pos <- near(pk, ev); held <- !pos & near(pk, oth)
      P <- pk[!held]; y <- as.integer(pos[!held])
      fm[[k]] <- c(auc = fast_auc(P$height, y), ap = fast_ap(P$height, y), ap1 = fast_ap_at(P$height, y))
      npk <- npk + nrow(P) / 5
      bh <- best_peak(pk[pos], ev); nfound <- nfound + sum(!is.na(bh))
      negh <- sort(P$height[y == 0], decreasing = TRUE)
      ranks <- c(ranks, ifelse(is.na(bh), Inf, findInterval(-bh, -negh) + 1))
    }
    fm <- do.call(rbind, fm); n_end <- nrow(kn)
    res[[length(res) + 1]] <- data.table(term = tm, set = "cross-validated", head = hd, thr = thr,
      n_ends = n_end, peaks = round(npk), ends_with_peak = nfound,
      roc_auc = mean(fm[, "auc"], na.rm = TRUE), ap = mean(fm[, "ap"], na.rm = TRUE),
      norm_pr_auc = mean(fm[, "ap1"], na.rm = TRUE),
      p50 = if (mean(is.finite(ranks)) >= 0.5) as.numeric(quantile(ranks, 0.5, type = 1)) else NA_real_,
      p80 = if (mean(is.finite(ranks)) >= 0.8) as.numeric(quantile(ranks, 0.8, type = 1)) else NA_real_)
    curves[[length(curves) + 1]] <- data.table(term = tm, set = "cross-validated", head = hd, thr = thr, rank = ranks)
    ## ---- never seen: the all-50 profile ----------------------------------
    pk <- call_peaks(b, b[[paste0(hd, "__all50")]], thr)
    if (tm == "C" && hd == "cnn")
      hcrt[[length(hcrt) + 1]] <- pk[accession == "O43612"][order(-height)][1:3][, thr := thr]
    for (st in c("never seen: peptides", "never seen: large products")) {
      ev <- unique(un[set == st & term == tm, .(accession, anchor)])
      oth <- unique(rbind(ann, kn[, .(accession, anchor)], un[set != st, .(accession, anchor)]))
      pos <- near(pk, ev); held <- !pos & near(pk, oth)
      P <- pk[!held]; y <- as.integer(pos[!held])
      bh <- best_peak(pk[pos], ev)
      negh <- sort(P$height[y == 0], decreasing = TRUE)
      ranks <- ifelse(is.na(bh), Inf, findInterval(-bh, -negh) + 1)
      res[[length(res) + 1]] <- data.table(term = tm, set = st, head = hd, thr = thr,
        n_ends = nrow(ev), peaks = nrow(P), ends_with_peak = sum(!is.na(bh)),
        roc_auc = fast_auc(P$height, y), ap = fast_ap(P$height, y), norm_pr_auc = fast_ap_at(P$height, y),
        p50 = if (mean(is.finite(ranks)) >= 0.5) as.numeric(quantile(ranks, 0.5, type = 1)) else NA_real_,
        p80 = if (mean(is.finite(ranks)) >= 0.8) as.numeric(quantile(ranks, 0.8, type = 1)) else NA_real_)
      curves[[length(curves) + 1]] <- data.table(term = tm, set = st, head = hd, thr = thr, rank = ranks)
    }
    message(sprintf("%s %-6s T=%.2f done", tm, hd, thr))
  }
  rm(b); invisible(gc())
}
R <- rbindlist(res); Cv <- rbindlist(curves)
cat("\nHCRT (O43612) C, linear head, all-50 -- top peaks:\n"); print(rbindlist(hcrt), row.names = FALSE)
options(width = 220)
print(R[, .(term, set, head, thr, n_ends, peaks, ends_with_peak, roc_auc = round(roc_auc, 4), ap = round(ap, 4),
            norm_pr_auc = round(norm_pr_auc, 3), p50, p80)][order(set, term, head, thr)], row.names = FALSE)
fwrite(R, file.path(A2, "lf_headcmp", "peak_metrics.csv"))
fwrite(Cv, file.path(A2, "lf_headcmp", "peak_ranks.csv"))
