## Stage 7 of the candidate pipeline: up to FOUR representative cuts per group
## (Kevin, 2026-10-27).
##
##   HEAD=xgb_pm THRESH=0.5 TOL=4 ARM=_arms3 Rscript candidate_best_per_group.R
##
## A group is one (N es, C es) pair; every member is the same two segments cut
## at a different offset, so the group needs a "which cut" answer, and four
## pieces of evidence give four different ones. Per group, the cut with the
## highest
##   xgb       xgb_score       the classifier's global score (out of fold)
##   end       end_score_max   the stronger of its two end scores (all-50 max)
##   multipep  mp_best_score   its best MultiPep class score
##   mmseqs    mmseqs_sim      its similarity to the closest known 8-60 aa peptide
##                             (ties: higher mmseqs_bits_per_res)
## Remaining ties go to xgb_score, then group_rank, then the id. A cut that wins
## several criteria is ONE row, `picked_by` naming them all ("xgb+end"), so a
## group has 1-4 rows. Groups are ordered by their best xgb_score, and within a
## group by criterion in the order above.
##
## Reads the stage-6 scored CSV (and the table, for group_size); writes
## candidate_peptides_<...>_best_per_group.csv. Cheap: no model, no MMseqs.
suppressMessages({library(data.table)})
A2   <- path.expand("~/AF2_analysis")
HEAD <- Sys.getenv("HEAD", "xgb_pm"); ARM <- Sys.getenv("ARM", "")
thr  <- as.numeric(Sys.getenv("THRESH", "0.5")); TOL <- as.integer(Sys.getenv("TOL", "4"))
key  <- sprintf("%s_T%s_tol%d%s", HEAD, thr, TOL, ARM)
f_cd <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s.csv", key))
f_sc <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_scored.csv", key))
f_bg <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_best_per_group.csv", key))

sc <- fread(f_sc)
sc[fread(f_cd, select = c("peptide_id", "group_size")), on = "peptide_id", group_size := i.group_size]

## criterion -> the columns it sorts on (all descending), before the shared
## tie-breaks
CRIT <- list(xgb = "xgb_score", end = "end_score_max", multipep = "mp_best_score",
             mmseqs = c("mmseqs_sim", "mmseqs_bits_per_res"))
miss <- setdiff(unlist(CRIT), names(sc))
if (length(miss)) stop("scored table lacks: ", paste(miss, collapse = ", "))
picks <- rbindlist(lapply(names(CRIT), function(nm) {
  by <- unique(c(CRIT[[nm]], "xgb_score", "group_rank", "peptide_id"))
  o <- sc[, c("group_id", "peptide_id", by[!by %in% c("group_id", "peptide_id")]), with = FALSE]
  setorderv(o, by, ifelse(by %in% c("group_rank", "peptide_id"), 1L, -1L), na.last = TRUE)
  o[!duplicated(group_id), .(group_id, peptide_id, crit = nm, crit_i = match(nm, names(CRIT)))]
}))
win <- picks[, .(picked_by = paste(crit[order(crit_i)], collapse = "+"), n_criteria = .N,
                 first_crit = min(crit_i)), by = .(group_id, peptide_id)]

bg <- sc[win, on = .(group_id, peptide_id)]
bg[, group_best_xgb := max(xgb_score), by = group_id]
setorder(bg, -group_best_xgb, group_id, first_crit)
bg[, c("group_best_xgb", "first_crit") := NULL]
front <- names(sc)[seq_len(match("group_rank", names(sc)))]
setcolorder(bg, c(front, "group_size", "picked_by", "n_criteria"))
fwrite(bg, f_bg)

cat(sprintf("%s groups -> %s rows (%.2f per group) in %s\n", format(uniqueN(bg$group_id), big.mark = ","),
            format(nrow(bg), big.mark = ","), nrow(bg) / uniqueN(bg$group_id), basename(f_bg)))
cat("rows per group:\n"); print(bg[, .N, by = group_id][, .(groups = .N), keyby = .(rows = N)], row.names = FALSE)
cat("how often each pair of criteria picks the same cut (of all groups):\n")
G <- uniqueN(bg$group_id)
for (a in names(CRIT)) for (b in names(CRIT)) if (match(a, names(CRIT)) < match(b, names(CRIT)))
  cat(sprintf("   %-9s = %-9s %5.1f%%\n", a, b,
              100 * bg[grepl(a, picked_by, fixed = TRUE) & grepl(b, picked_by, fixed = TRUE), .N] / G))
