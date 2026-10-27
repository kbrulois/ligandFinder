## Window-model metrics, attention trunk (current pages) vs flatten trunk (final),
## each with its CNN head and its per-member xgb head -- four models, scored over
## EVERY mature anchor of the proteome (the full scans), with the definitions of
## win_stats.R / rr_curves.R so the numbers line up with the deck:
##
##   ROC AUC, AP   pool = all anchors minus every OTHER annotated end (ligand_list
##                 + the CV knowns + the never-seen ends): those are neither clean
##                 positives nor honest negatives
##   norm. PR-AUC  AP re-derived at a common 1% prevalence by reweighting the
##                 negatives along the same ranking (fast_ap_at)
##   recall@rank   each positive's rank among the anchors that are not an annotated
##                 end anywhere; reported as anchors needed for 50% / 80% recall
##
## Sets: "cross-validated" = the 5-fold knowns, each scored ONLY by the fold that
## held it out (metrics per fold, then averaged; ranks pooled); "never seen" =
## ends whose gene no fold trained on (cv5_unseen_ends.csv), scored by all 50.
suppressMessages({library(data.table); library(dplyr); library(reticulate)})
S <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
np <- import("numpy", convert = TRUE)
ll <- as.data.frame(readRDS(lf_src("inst/extdata/ligand_list.rds")))
un <- fread(file.path(S, "cv5_unseen_ends.csv"))
prec <- fread(file.path(A2, "pepend_locus", "prec_order.csv"))

fast_auc <- function(score, y) {
  r <- frank(score, ties.method = "average"); n1 <- sum(y); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fast_ap <- function(s, y) { o <- order(-s); yy <- y[o]; sum((cumsum(yy) / seq_along(yy))[yy == 1]) / sum(y) }
fast_ap_at <- function(s, y, target = 0.01) {
  o <- order(-s); yy <- y[o]; n1 <- sum(y); n0 <- length(y) - n1
  w <- n1 * (1 - target) / (target * n0); tp <- cumsum(yy); fp <- cumsum(1 - yy)
  sum((tp / (tp + w * fp))[yy == 1]) / n1
}
rank_in <- function(sp, neg) findInterval(-sp, -sort(neg, decreasing = TRUE)) + 1L

## the flatten trunk's per-anchor tracks, from its ten full scans
flat_tracks <- function(tm) {
  p <- file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm))
  if (file.exists(p)) return(readRDS(p))
  fs <- sprintf("%s/lf_winxgb_fullscan_%s_cv5_f%d_nomask_flatten.npz", A2, tm, 1:5)
  stopifnot(all(file.exists(fs)))
  z1 <- np$load(fs[1], allow_pickle = TRUE)
  b <- data.table(accession = prec$accession[z1[["prot_idx"]] + 1L],
                  anchor = as.integer(z1[["anchor"]]))
  for (k in 1:5) {
    z <- np$load(fs[k], allow_pickle = TRUE)
    stopifnot(identical(as.integer(z[["anchor"]]), b$anchor))
    b[[sprintf("cnn__f%d", k)]] <- as.numeric(z[["cnn"]])
    b[[sprintf("xgb_pm__f%d", k)]] <- as.numeric(z[["xgb_pm"]])
  }
  b[, cnn__all50 := rowMeans(.SD), .SDcols = sprintf("cnn__f%d", 1:5)]
  b[, xgb_pm__all50 := rowMeans(.SD), .SDcols = sprintf("xgb_pm__f%d", 1:5)]
  saveRDS(b, p); b
}

ARMS <- strsplit(Sys.getenv("ARMS", "flatten,esmonly36"), ",")[[1]]   # attention is retired
NGRID <- unique(round(10^seq(0, 6.48, length.out = 150)))
res <- list(); curves <- list()
for (tm in c("C", "N")) {
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  kn[, key := paste(accession, anchor)]
  ann <- ll %>% filter(!is.na(start), !is.na(end)) %>%
    transmute(key = paste(accession, if (tm == "N") start else end)) %>% pull(key) %>% unique()
  u <- un[term == tm]; u[, key := paste(accession, anchor)]
  all_ends <- unique(c(ann, kn$key, u$key))
  ## ARMS is comma-separated so a new arm is an env var, not an edit. The ESM-C
  ## arm reads the tracks rds tracks_esmonly36.R already built from its full
  ## scans -- nothing is rescanned here.
  for (trunk in ARMS) {
    cols <- c("accession", "anchor", sprintf("cnn__f%d", 1:5), sprintf("xgb_pm__f%d", 1:5),
              "cnn__all50", "xgb_pm__all50")
    b <- if (trunk == "attention")
      as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))[, cols, with = FALSE]
    else if (trunk == "flatten") flat_tracks(tm)
    else as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten_%s.rds",
                                                    tm, trunk))))[, cols, with = FALSE]
    b[, key := paste(accession, anchor)]
    not_end <- !(b$key %in% all_ends)
    for (hd in c("cnn", "xgb_pm")) {
      ## ---- cross-validated: fold k's own ends, fold k's own column ----------
      fm <- list(); rk <- integer(0)
      for (k in 1:5) {
        col <- sprintf("%s__f%d", hd, k); pk <- kn[fold == k]$key
        keep <- !(b$key %in% setdiff(all_ends, pk)) & !is.na(b[[col]])
        y <- as.integer(b$key[keep] %in% pk); s <- b[[col]][keep]
        fm[[k]] <- c(auc = fast_auc(s, y), ap = fast_ap(s, y), ap1 = fast_ap_at(s, y),
                     npos = sum(y), pool = length(y))
        neg <- b[[col]][not_end]; neg <- neg[!is.na(neg)]
        rk <- c(rk, rank_in(b[[col]][b$key %in% pk], neg))
      }
      fm <- do.call(rbind, fm)
      res[[length(res) + 1]] <- data.table(term = tm, set = "cross-validated", trunk = trunk, head = hd,
        n = sum(fm[, "npos"]), roc_auc = mean(fm[, "auc"]), ap = mean(fm[, "ap"]), norm_pr_auc = mean(fm[, "ap1"]),
        prev_per_M = 1e6 * mean(fm[, "npos"] / fm[, "pool"]),
        n50 = NGRID[which(sapply(NGRID, \(x) mean(rk <= x)) >= 0.5)[1]],
        n80 = NGRID[which(sapply(NGRID, \(x) mean(rk <= x)) >= 0.8)[1]])
      curves[[length(curves) + 1]] <- data.table(term = tm, set = "cross-validated", trunk = trunk, head = hd,
        N = NGRID, recall = sapply(NGRID, \(x) mean(rk <= x)))
      ## ---- never seen: all 50, which is honest for genes no fold trained on --
      for (st in unique(u$set)) {
        col <- paste0(hd, "__all50"); pk <- u[set == st]$key
        keep <- !(b$key %in% setdiff(all_ends, pk)) & !is.na(b[[col]])
        y <- as.integer(b$key[keep] %in% pk); s <- b[[col]][keep]
        neg <- b[[col]][not_end]; neg <- neg[!is.na(neg)]
        rk <- rank_in(b[[col]][b$key %in% pk], neg)
        lab <- if (grepl("^never seen", st)) "never seen: peptides" else "never seen: large products"
        res[[length(res) + 1]] <- data.table(term = tm, set = lab, trunk = trunk, head = hd,
          n = sum(y), roc_auc = fast_auc(s, y), ap = fast_ap(s, y), norm_pr_auc = fast_ap_at(s, y),
          prev_per_M = 1e6 * mean(y),
          n50 = NGRID[which(sapply(NGRID, \(x) mean(rk <= x)) >= 0.5)[1]],
          n80 = NGRID[which(sapply(NGRID, \(x) mean(rk <= x)) >= 0.8)[1]])
        curves[[length(curves) + 1]] <- data.table(term = tm, set = lab, trunk = trunk, head = hd,
          N = NGRID, recall = sapply(NGRID, \(x) mean(rk <= x)))
      }
    }
    rm(b); invisible(gc())
    message(tm, " ", trunk, " done")
  }
}
R <- rbindlist(res); Cv <- rbindlist(curves)
fwrite(R, file.path(A2, "lf_headcmp", "final_metrics.csv"))
fwrite(Cv, file.path(A2, "lf_headcmp", "final_recall_rank.csv"))
options(width = 200)
print(R[, .(term, set, trunk, head, n, roc_auc = round(roc_auc, 4), ap = round(ap, 4),
            norm_pr_auc = round(norm_pr_auc, 3), prev_per_M = round(prev_per_M, 1), n50, n80)][
  order(set, term, head, trunk)], row.names = FALSE)
