## Optional stage: attach MultiPep (Scheele lab) class probabilities to the
## candidate table. MultiPep is multi-label over 20 bioactivity classes; the
## three that match our biology are peptidehormone / neuropeptide /
## cellcellsignaling, combined into mp_signalling_max.
##
## MultiPep keys its output by SEQUENCE (first column), so the merge is on
## sequence and survives a pipeline re-run at another THRESH/TOL.
##   HEAD=xgb_pm THRESH=0.5 TOL=4 Rscript candidate_multipep_merge.R
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis")
HEAD <- Sys.getenv("HEAD", "xgb_pm")
thr <- as.numeric(Sys.getenv("THRESH", "0.5")); TOL <- as.integer(Sys.getenv("TOL", "4"))
f <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_T%s_tol%d.csv", HEAD, thr, TOL))
cd <- fread(f)
mp <- fread(file.path(A2, "multipep_ext", "multipep_scores.csv"))
setnames(mp, 1, "sequence")
cls <- setdiff(names(mp), "sequence")
setnames(mp, cls, paste0("mp_", gsub("[^A-Za-z0-9]+", "_", cls)))
mpc <- setdiff(names(mp), "sequence")
if (length(intersect(names(cd), c(mpc, "mp_signalling_max"))))
  cd[, (intersect(names(cd), c(mpc, "mp_signalling_max"))) := NULL]
cd[mp, on = "sequence", (mpc) := mget(paste0("i.", mpc))]
sig <- c("mp_peptidehormone", "mp_neuropeptide", "mp_cellcellsignaling")
cd[, mp_signalling_max := round(pmax(get(sig[1]), get(sig[2]), get(sig[3])), 4)]
for (c_ in mpc) set(cd, j = c_, value = round(cd[[c_]], 4))
fwrite(cd, f)
cat(sprintf("merged %d MultiPep classes onto %s rows (%s unscored), %d columns total\n",
            length(mpc), format(nrow(cd), big.mark = ","),
            format(sum(is.na(cd$mp_peptidehormone)), big.mark = ","), ncol(cd)))
