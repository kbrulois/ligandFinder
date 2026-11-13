## Held-out comparison of arms from the trunks' OWN summaries, before any scan
## (2026-11-10, the first read of the two _padg arms; superseded the next day by
## the full-proteome metrics in lf_headcmp/unseen_padg, which are the numbers to
## quote -- this is the quick look that needs nothing but training).
##
## Per arm and terminus, over the five folds (10 seeds each):
##   * ensemble val AP / ROC AUC: each fold's `val` split (that fold's fifth of
##     the knowns + its ~400-830 negatives), scored by the ensemble at the
##     weights EarlyStopping kept (summary$val, as 10_8 saved it), not the
##     per-epoch histories, which misreport the shipped model
##   * the held-out knowns' rank among the arm's candidate pool (rank_cand,
##     ~40,000 candidate windows per terminus), every known end ranked once, by
##     the fold that held it out
##   * PAIR: the same ends, arm by arm, matched on gene + anchor -- so two arms
##     on different C conventions pair only the C ends that did not move
##   * the C ends anchored on the G (g_amid, arms built with --amid-g-anchor)
##     against the rest
##
## READ WITH CARE. Val AP is NOT comparable across window sets: the _padg val
## splits carry about twice the negatives (the top-ups), so curated_AAswap's N
## 0.552 against the _padg arms' 0.450 / 0.475 is a harder set, not a worse
## model. rank_cand is over the same candidate pool for every arm (N identical;
## C the same ~41,800 windows but for ~30 the G convention moves by a residue).
## On C much of the _padg gain is the G re-anchoring: the 30 moved ends have
## median rank 17 / 22, the 55 unmoved 313 / 248 -- a G next to a dibasic pair
## is intrinsically easier to find. No N end moved, so N's gain is the windows'.
## 2026-11-10 (C / N): val AP 0.616 / 0.552 (old windows), 0.693 / 0.450 (padg),
## 0.689 / 0.475 (padg + MSA cons); median rank 337 / 409, 103 / 178.5, 77 / 191.5.
##
##   LC_ALL=en_US.UTF-8 Rscript trunk_val_compare.R
##   ARMS=curated_AAswap,esmonly36 PAIR= HEADCMP=/tmp/tv LC_ALL=en_US.UTF-8 Rscript trunk_val_compare.R
##
## Reads ~/AF2_analysis/lf_pepend_run_<T>_cv5_f<k>_nomask_flatten[_<arm>]_summary.rds
## (stage 1's per-trunk summaries; "curated" = the unsuffixed curated baseline).
## Writes to HEADCMP (default ~/AF2_analysis/lf_headcmp):
##   trunk_val_metrics.csv   arm x term x fold: n_val, n_pos, val_ap, val_auc
##   trunk_val_ranks.csv     one row per held-out known end: rank_cand, pred
suppressMessages(library(data.table))
A2 <- path.expand("~/AF2_analysis")
ARMS <- strsplit(Sys.getenv("ARMS", "curated_AAswap,curated_AAswap_padg,curated_AAswap_msacons_padg"), ",")[[1]]
PAIR <- strsplit(Sys.getenv("PAIR", "curated_AAswap_padg,curated_AAswap_msacons_padg"), ",")[[1]]
HC <- path.expand(Sys.getenv("HEADCMP", file.path(A2, "lf_headcmp")))
dir.create(HC, showWarnings = FALSE, recursive = TRUE)
## the 2026-11-10 table's names for the three arms; any other arm keeps its own
SHORT <- c(curated_AAswap = "AAswap_old", curated_AAswap_padg = "AAswap_padg",
           curated_AAswap_msacons_padg = "msacons_padg")
short <- function(a) unname(ifelse(a %in% names(SHORT), SHORT[a], a))

ap  <- function(s, y) { o <- order(-s); yy <- y[o]; sum((cumsum(yy) / seq_along(yy))[yy == 1]) / sum(y) }
auc <- function(s, y) {
  r <- frank(s, ties.method = "average"); n1 <- sum(y)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1))
}

out <- list(); rk <- list()
for (a in ARMS) for (tm in c("C", "N")) for (k in 1:5) {
  sfx <- if (a %in% c("", "curated", "flatten")) "" else paste0("_", a)
  f <- file.path(A2, sprintf("lf_pepend_run_%s_cv5_f%d_nomask_flatten%s_summary.rds", tm, k, sfx))
  if (!file.exists(f)) { message("missing ", basename(f)); next }
  s <- readRDS(f); v <- as.data.table(s$val); y <- as.integer(as.character(v$truth))
  out[[length(out) + 1]] <- data.table(arm = short(a), term = tm, fold = k, n_val = nrow(v), n_pos = sum(y),
                                       val_ap = ap(v$score, y), val_auc = auc(v$score, y))
  kr <- as.data.table(s$knowns_ranked)[split == "val"]
  rk[[length(rk) + 1]] <- data.table(arm = short(a), term = tm, fold = k, gene = kr$gene, anchor = kr$anchor,
                                     g_amid = if ("g_amid" %in% names(kr)) kr$g_amid else NA,
                                     rank = as.double(kr$rank_cand), pred = kr$pred)
}
R <- rbindlist(out); K <- rbindlist(rk, fill = TRUE)
lev <- short(ARMS)
R[, arm := factor(arm, lev)]; K[, arm := factor(arm, lev)]

cat("=== ensemble val metrics, mean over the 5 folds (10 seeds each) ===\n")
print(R[, .(n_val = mean(n_val), n_pos = sum(n_pos), val_AP = round(mean(val_ap), 4), val_AUC = round(mean(val_auc), 4)),
        keyby = .(term, arm)])
cat("\n=== held-out knowns: rank among the arm's candidate pool (lower is better) ===\n")
print(K[, .(n = .N, median_rank = as.double(median(rank)), mean_rank = round(mean(rank), 1), top10 = round(mean(rank <= 10), 3),
            top100 = round(mean(rank <= 100), 3), median_pred = round(median(pred), 3)), keyby = .(term, arm)])

if (length(PAIR) == 2 && all(short(PAIR) %in% K$arm)) {
  p0 <- short(PAIR[1]); p1 <- short(PAIR[2])
  cat(sprintf("\n=== paired, same ends: %s -> %s ===\n", p0, p1))
  w <- dcast(K[arm %in% c(p0, p1)][, arm := as.character(arm)], term + gene + anchor ~ arm, value.var = "rank")
  setnames(w, c(p0, p1), c("r0", "r1"))
  print(w[!is.na(r0) & !is.na(r1), .(n = .N, better = sum(r1 < r0), worse = sum(r1 > r0), tied = sum(r1 == r0),
                                     median_delta = median(r1 - r0)), keyby = term])
}
if (any(!is.na(K$g_amid))) {
  cat("\n=== C ends anchored on the G (arms built with --amid-g-anchor) ===\n")
  print(K[!is.na(g_amid) & term == "C", .(n = .N, median_rank = as.double(median(rank)), top10 = round(mean(rank <= 10), 3)),
          keyby = .(arm, g_amid)])
}
fwrite(R, file.path(HC, "trunk_val_metrics.csv"))
fwrite(K, file.path(HC, "trunk_val_ranks.csv"))
cat("\nwrote", file.path(HC, c("trunk_val_metrics.csv", "trunk_val_ranks.csv")), sep = "\n  ")
