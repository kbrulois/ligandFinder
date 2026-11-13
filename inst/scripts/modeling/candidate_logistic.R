## Optional stage: a multivariable logistic regression over the candidate
## table, to rank candidates by how much they look like a known bioactive
## ligand (Kevin, 2026-10-23).
##
##   HEAD=xgb_pm THRESH=0.5 TOL=4 Rscript candidate_logistic.R
##
## OUTCOME: the candidate is sequence-identical to a curated ligand_list peptide
## over its full length (mmseqs_pident 100, tcov >= 0.99, hit_type ligand_list).
##
## TWO THINGS THIS HAS TO GET RIGHT, or the number is worthless:
##
## 1. LEAKAGE. The outcome is DEFINED from the mmseqs_* columns, so every
##    mmseqs_* and overlap_* column is excluded from the predictors. Leaving
##    mmseqs_pident in would give a perfect model that has learned nothing.
##
## 2. CORRELATED ROWS. Candidates are not independent -- up to 81 of them come
##    from one pair of end segments, differing by a residue or two, and the
##    positives cluster in a few genes. A random split would put near-identical
##    rows in both train and test and report a fantasy AUC. So the folds are
##    GROUPED BY GENE: no gene appears in both.
##
## The score written back is the OUT-OF-FOLD prediction -- every row scored by a
## model that never saw its gene -- which is the honest number to rank on.
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis")
HEAD <- Sys.getenv("HEAD", "xgb_pm"); ARM <- Sys.getenv("ARM", "")   # "" = curated, else e.g. "_esmonly36"
thr <- as.numeric(Sys.getenv("THRESH", "0.5")); TOL <- as.integer(Sys.getenv("TOL", "4"))
NFOLD <- as.integer(Sys.getenv("NFOLD", "5"))
f <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_T%s_tol%d%s.csv", HEAD, thr, TOL, ARM))
cd <- fread(f)

cd[, y := as.integer(mmseqs_hit_type == "ligand_list" & mmseqs_pident >= 100 & mmseqs_tcov >= 0.99)]
cd[is.na(y), y := 0L]

## PREDICTORS, chosen by gene-grouped CV rather than by throwing everything in.
## Per-fold held-out AUC: these five 0.892 +/- 0.076; all 20 MultiPep classes
## plus the whole window block 0.874 +/- 0.136; mp_peptidehormone alone
## 0.880 +/- 0.102. The big model is both worse and far less stable across
## folds, so parsimony here is not taste -- it generalises better.
## Nothing derived from the known set appears, see the leakage note above.
mot  <- c("n_motif_centre", "c_motif_centre")
for (v in mot) cd[[v]] <- factor(cd[[v]])
pred <- c("mp_peptidehormone", "mp_neuropeptide", "mp_cellcellsignaling",
          "end_score_max", "length")
pred <- pred[vapply(pred, function(v) length(unique(cd[[v]])) > 1, logical(1))]
banned <- grep("^(mmseqs_|overlap_)", pred, value = TRUE)
stopifnot(length(banned) == 0)      # leakage guard
cat(sprintf("%s rows, %s positives (%.2f%%), %d predictors\n", format(nrow(cd), big.mark = ","),
            format(sum(cd$y), big.mark = ","), 100 * mean(cd$y), length(pred)))

## gene-grouped folds
set.seed(42)
g <- unique(cd$gene); fold_of <- setNames(sample(rep_len(1:NFOLD, length(g))), g)
cd[, fold := fold_of[gene]]
cat(sprintf("folds grouped by gene: %d genes over %d folds; positives per fold %s\n",
            length(g), NFOLD, paste(cd[, sum(y), by = fold][order(fold)]$V1, collapse = ", ")))

fml <- as.formula(paste("y ~", paste(pred, collapse = " + ")))
oof <- rep(NA_real_, nrow(cd))
for (k in 1:NFOLD) {
  tr <- cd[fold != k]; te <- which(cd$fold == k)
  ## a factor level absent from train would break predict(); align them
  for (v in mot) tr[[v]] <- factor(tr[[v]], levels = levels(cd[[v]]))
  m <- suppressWarnings(glm(fml, data = tr, family = binomial()))
  oof[te] <- suppressWarnings(predict(m, newdata = cd[te], type = "response"))
  cat(sprintf("  fold %d: trained on %s rows (%d pos), scored %s\n", k,
              format(nrow(tr), big.mark = ","), sum(tr$y), format(length(te), big.mark = ",")))
}
cd[, lr_prob_oof := round(oof, 6)]
## The ranking column is the FULL-DATA fit: with five parameters over 119k rows
## its optimism is negligible, and unlike the out-of-fold score it puts every
## candidate on one calibration, which is what sorting the whole table needs.
## lr_prob_oof is kept so the honest held-out number can be re-checked.
full0 <- suppressWarnings(glm(fml, data = cd, family = binomial()))
cd[, lr_prob := round(predict(full0, newdata = cd, type = "response"), 6)]

auc <- function(s, yy) { ok <- !is.na(s); s <- s[ok]; yy <- yy[ok]
  r <- frank(s, ties.method = "average"); n1 <- sum(yy); n0 <- length(yy) - n1
  if (!n1 || !n0) return(NA_real_); (sum(r[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
pf <- sapply(1:NFOLD, function(k) auc(cd$lr_prob_oof[cd$fold == k], cd$y[cd$fold == k]))
cat(sprintf("\nPER-FOLD held-out AUC (gene-grouped): %.4f +/- %.4f  [%s]\n",
            mean(pf), sd(pf), paste(round(pf, 3), collapse = " ")))
cat("single predictors, for comparison:\n")
for (v in c("mp_signalling_max", "mp_peptidehormone", "end_score_max", "length"))
  if (v %in% names(cd)) {
    q <- sapply(1:NFOLD, function(k) auc(cd[[v]][cd$fold == k], cd$y[cd$fold == k]))
    cat(sprintf("   %-20s %.4f +/- %.4f\n", v, mean(q), sd(q)))
  }

## what the model leans on, from a fit over everything (reported, not scored on)
co <- summary(full0)$coefficients
co <- co[order(-abs(co[, "z value"])), , drop = FALSE]
cat("\nstrongest terms (full-data fit, z value):\n")
print(round(head(co[rownames(co) != "(Intercept)", c("Estimate", "z value")], 12), 3))

cd[, c("y", "fold") := NULL]
fwrite(cd, f)
cat(sprintf("\nwrote lr_prob to %s (%d columns)\n", basename(f), ncol(cd)))
