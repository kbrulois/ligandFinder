## Optional stage: attach PDeepPP scores to the candidate table from the cache
## that score_candidates.py writes (keyed by SEQUENCE, so it survives a re-run
## of the pipeline at a different THRESH/TOL without re-scoring).
##   TASK=neuro HEAD=xgb_pm THRESH=0.5 TOL=4 Rscript candidate_pdeeppp_merge.R
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis")
TASK <- Sys.getenv("TASK", "neuro"); HEAD <- Sys.getenv("HEAD", "xgb_pm")
thr <- as.numeric(Sys.getenv("THRESH", "0.5")); TOL <- as.integer(Sys.getenv("TOL", "4"))
f <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_T%s_tol%d.csv", HEAD, thr, TOL))
cd <- fread(f)
sc <- fread(file.path(A2, "pdeeppp_run", sprintf("pdeeppp_%s_scores.csv", TASK)))
lc <- paste0("pdeeppp_", TASK, "_logit"); pc <- paste0("pdeeppp_", TASK, "_prob")
if (length(intersect(names(cd), c(lc, pc)))) cd[, (intersect(names(cd), c(lc, pc))) := NULL]
cd[sc, on = "sequence", (c(lc, pc)) := .(round(i.logit, 3), round(i.prob, 6))]
fwrite(cd, f)
cat(sprintf("merged %s scores onto %s rows (%s unscored)\n", TASK,
            format(nrow(cd), big.mark = ","), format(sum(is.na(cd[[lc]])), big.mark = ",")))
