## Drop the precursors with NO ESM C embedding from the TRAINING windows
## (Kevin, 2026-10-24), so the ESM-C-only arm can be retrained without the rows
## that carry an all-zero 36-channel block and therefore no information at all.
##
##   Rscript nn_input_esmcov.R
##
## Writes ~/AF2_analysis/lf_pepend_nn_input_cv5_f<k>_esmcov.rds for k = 1..5,
## which train_esmonly36x.sh passes to 10_8 as --input. 10_8 is NOT modified:
## the cut is made in its input, which is also the only place it can be made
## without changing a script that is the record of the existing numbers.
##
## WHAT IS CUT, AND WHAT IS NOT
##   train, val   filtered to the 5,182 accessions present in esm_c_pca36.parquet
##   all          LEFT WHOLE, deliberately
##
## `all` is left whole for the same reason 10_8's --positives leaves it whole:
## it is the scoring set, not a fitting set. lf_dcnn's pipeline.py touches
## splits["all"] only at lines 284 and 295 -- predict_all and embed, both AFTER
## the fit -- so it cannot influence a weight. Filtering it would change the
## per-run summary rds and nothing else, at the cost of making this arm's
## summary non-comparable to the existing one. So the cut is exactly the fitting
## set and nothing more.
##
## A consequence worth stating: 10_8 builds its zero blocks from `need`, which
## spans every split, so it will still print "433 get ZERO channels". That is
## correct and expected -- the blocks exist for the `all` split. The proof the
## cut landed is the train/val window counts 10_8 prints next, which must match
## the "after" column below.
suppressMessages({library(data.table); library(arrow); library(dplyr)})
A2 <- path.expand("~/AF2_analysis")

cov <- unique(as.character(as.data.frame(
  open_dataset(file.path(A2, "lf_plm/esm_c_pca36.parquet")) |>
    select(accession) |> collect())$accession))
cat(sprintf("esm_c_pca36.parquet covers %d accessions\n\n", length(cov)))

rep <- list()
for (k in 1:5) {
  ip <- file.path(A2, sprintf("lf_pepend_nn_input_cv5_f%d.rds", k))
  op <- file.path(A2, sprintf("lf_pepend_nn_input_cv5_f%d_esmcov.rds", k))
  if (file.exists(op)) { cat(sprintf("fold %d: %s exists, skipping\n", k, basename(op))); next }
  x <- readRDS(ip)
  for (tm in names(x$nn_input)) for (sp in c("train", "val")) {
    d <- x$nn_input[[tm]][[sp]]
    keep <- as.character(d$accession) %in% cov
    rep[[length(rep)+1]] <- data.table(
      fold = k, term = tm, split = sp, before = nrow(d), after = sum(keep),
      dropped = sum(!keep),
      pos_before = sum(d$known == 1L), pos_after = sum(d$known[keep] == 1L),
      prec_dropped = uniqueN(as.character(d$accession)[!keep]))
    x$nn_input[[tm]][[sp]] <- d[keep, ]
  }
  ## every split must still have both classes or lf_dcnn's sampler raises
  for (tm in names(x$nn_input)) for (sp in c("train", "val")) {
    d <- x$nn_input[[tm]][[sp]]
    if (sum(d$known == 1L) == 0L || sum(d$known == 0L) == 0L)
      stop(sprintf("fold %d %s %s lost a class: %d pos, %d neg", k, tm, sp,
                   sum(d$known == 1L), sum(d$known == 0L)), call. = FALSE)
  }
  saveRDS(x, op)
  cat(sprintf("fold %d -> %s\n", k, basename(op)))
  rm(x); invisible(gc(FALSE))
}

r <- rbindlist(rep)
cat("\n=== windows removed from the FITTING splits ===\n")
print(r, row.names = FALSE, nrows = 100)
cat("\n=== totals ===\n")
print(r[, .(before = sum(before), after = sum(after), dropped = sum(dropped),
            pct = round(100*sum(dropped)/sum(before), 2),
            pos_before = sum(pos_before), pos_after = sum(pos_after)),
        by = .(term, split)], row.names = FALSE)
fwrite(r, file.path(A2, "lf_nomask_scripts_20261021", "nn_input_esmcov_report.csv"))
