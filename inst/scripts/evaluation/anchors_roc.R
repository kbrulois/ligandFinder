## Threshold as a classifier over ANCHORS (Kevin, 2026-11-02) -- the per-anchor
## counterpart of segments_roc.R, with its end sets, held-out ends, profiles and
## threshold grid, so the two sweeps compare like for like.
##
## At threshold T every scored anchor with score >= T is a CALL:
##   TP_anc  a call within +/-TOL of an evaluation end
##   FP      a call within TOL of no annotated end
##   TN      an anchor below T within TOL of no annotated end
## A call within TOL of an annotated end from OUTSIDE the evaluation set (other
## folds' ends, ligand_list, UniProt peptide/propeptide, the other never-seen
## set) is held out of every count, as segments_roc.R holds out such segments;
## an uncalled anchor near an evaluation end is neither TN nor FN.
## Ends are counted per END, as in segments_roc.R, so TP + FN = the ends:
##   TP  an end with at least one call within TOL
##   FN  every other end (FN_unscored: no scored anchor within TOL at all)
## precision = TP_anc / (TP_anc + FP), over calls; TPR = TP / ends; FPR = FP / (FP + TN).
## TOL = 0 is metrics_final.R's own per-anchor classification (one anchor per end).
## Cross-validated: fold k's ends on fold k's own profile, counts summed over
## folds; never seen: the all-50 profile.
##
##   TOL=4 ARM=_esmonly36 Rscript anchors_roc.R   -> lf_headcmp/anchor_confusion_esmonly36_tol4.csv
##   TOL=4 Rscript anchors_roc.R                  -> lf_headcmp/anchor_confusion_tol4.csv (curated)
## segment_T_sweep_plot.R LEVEL=anchor draws them.
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- path.expand(Sys.getenv("LF_RUN", file.path(A2, "lf_run")))
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "..", "lf_paths.R"))
TOL  <- as.integer(Sys.getenv("TOL", "4"))
ARM  <- Sys.getenv("ARM", "")              # which tracks to read ("" = curated)
GRID <- sort(unique(round(c(seq(0, 1, by = 0.025), seq(0.1, 0.5, by = 0.05)), 3)))   # segments_roc.R's

## the same inputs as segments_roc.R
ll <- as.data.table(readRDS(lf_src("inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
un <- fread(file.path(S, "cv5_unseen_ends.csv"))
un[, set := ifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]

## count of x >= t for every t in GRID (x sorted ascending, no NA)
n_ge <- function(x, t) length(x) - findInterval(t, x, left.open = TRUE)

out <- list()
for (tm in c("C", "N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, ARM))))
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))[, .(accession, anchor, fold)]
  ann <- unique(rbind(ll[, .(accession, anchor = if (tm == "N") start else end)],
                      up[, .(accession, anchor = if (tm == "N") start else end)]))
  ## integer keys: accession index * 1e6 + anchor
  ACC <- sort(unique(c(b$accession, kn$accession, ann$accession, un$accession)))
  key <- function(a, x) match(a, ACC) * 1e6 + x
  bk <- key(b$accession, b$anchor)
  ## rows of b within TOL of any end in `e`
  near <- function(e) {
    if (!nrow(e)) return(rep(FALSE, nrow(b)))
    k <- unique(as.vector(outer(key(e$accession, e$anchor), -TOL:TOL, `+`)))
    bk %in% k
  }
  ## each end's best score within TOL (-Inf: no scored anchor there)
  best_near <- function(e, s) {
    k <- outer(key(e$accession, e$anchor), -TOL:TOL, `+`)
    v <- matrix(s[match(k, bk)], nrow = nrow(e))
    m <- suppressWarnings(apply(v, 1, max, na.rm = TRUE))
    m[!is.finite(m)] <- -Inf
    m
  }
  counts <- function(s, ev, oth) {
    nev <- near(ev); free <- !nev & !near(oth); ok <- !is.na(s)
    s_ev <- sort(s[nev & ok]); s_free <- sort(s[free & ok])
    bn <- best_near(ev, s)
    data.table(thr = GRID, TP = vapply(GRID, function(t) sum(bn >= t), 0), ends = nrow(ev),
               FN_unscored = sum(!is.finite(bn)),
               TP_anc = n_ge(s_ev, GRID), FP = n_ge(s_free, GRID),
               TN = length(s_free) - n_ge(s_free, GRID))
  }
  for (hd in c("cnn", "xgb_pm")) {
    ## cross-validated: per fold, summed
    cv <- rbindlist(lapply(1:5, function(k) {
      ev  <- unique(kn[fold == k, .(accession, anchor)])
      oth <- unique(rbind(ann, kn[fold != k, .(accession, anchor)], un[, .(accession, anchor)]))
      counts(b[[sprintf("%s__f%d", hd, k)]], ev, oth)
    }))[, lapply(.SD, sum), by = thr]
    out[[length(out) + 1]] <- data.table(term = tm, set = "cross-validated", head = hd, cv)
    for (st in c("never seen: peptides", "never seen: large products")) {
      ev  <- unique(un[set == st & term == tm, .(accession, anchor)])
      oth <- unique(rbind(ann, kn[, .(accession, anchor)], un[set != st, .(accession, anchor)]))
      out[[length(out) + 1]] <- data.table(term = tm, set = st, head = hd,
                                           counts(b[[paste0(hd, "__all50")]], ev, oth))
    }
    message(sprintf("%s %-6s done (%d thresholds)", tm, hd, length(GRID)))
  }
  rm(b); invisible(gc())
}
R <- rbindlist(out)
R[, FN := ends - TP]
R[, `:=`(TPR = TP / ends, FPR = FP / (FP + TN), precision = TP_anc / (TP_anc + FP))]
R[, auc_T := (1 + TPR - FPR) / 2]
R[, F1 := fifelse(is.na(precision) | precision + TPR == 0, 0, 2 * precision * TPR / (precision + TPR))]
setcolorder(R, c("term", "set", "head", "thr", "TP", "FN", "FN_unscored", "ends", "TP_anc", "FP", "TN"))
f <- file.path(A2, "lf_headcmp", sprintf("anchor_confusion%s_tol%d.csv", ARM, TOL))
fwrite(R, f)
cat("wrote", f, "\n")
options(width = 200)
print(R[, .SD[which.max(F1)], by = .(term, set, head)][
  , .(term, set, head, thr, F1 = round(F1, 3), precision = round(precision, 3), TPR = round(TPR, 3),
      TP, ends, TP_anc, FP)][order(set, term, head)], row.names = FALSE)
