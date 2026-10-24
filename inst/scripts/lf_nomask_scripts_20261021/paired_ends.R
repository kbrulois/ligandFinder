## Is the ESM-C-only never-seen advantage bigger than ~45 ends can resolve?
## (Kevin, 2026-10-24)
##
## segments_roc.R reports TP and FN as independent counts per arm, which cannot
## answer that: the two arms are scored on the SAME ends, so the right test is
## paired. This script recomputes, end by end, whether each arm's profile caught
## it at T = 0.5 / TOL = 4 -- using segments_roc.R's own segment logic, copied
## verbatim below so the per-end verdicts agree with the published counts -- and
## then runs McNemar's exact test on the discordant pairs.
##
##   ARMS="_esmonly36" Rscript paired_ends.R
##   ARMS="_esmonly36 _esmonly36x" Rscript paired_ends.R
##
## Writes lf_headcmp/paired_ends_<arm>_T0.5_tol4.csv (one row per end per head)
## and prints the contingency tables.
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7"
TOL <- as.integer(Sys.getenv("TOL", "4")); THR <- as.numeric(Sys.getenv("THR", "0.5"))
ARMS <- strsplit(trimws(Sys.getenv("ARMS", "_esmonly36")), "\\s+")[[1]]

## ---- segments_roc.R's logic, verbatim -------------------------------------
segments <- function(b, s, thr) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  d <- data.table(sid, accession = b$accession, anchor = b$anchor, s = s, up = up)
  seg <- d[, .(start = anchor[1], end = anchor[.N], up = up[1],
               centre = anchor[which.max(s)]), by = .(sid, accession)]
  list(seg = seg, sid = sid)
}
centre_near <- function(seg, ends) {
  if (!nrow(seg) || !nrow(ends)) return(rep(FALSE, nrow(seg)))
  e <- unique(ends[, .(accession, anchor)]); setkey(e, accession, anchor)
  q <- data.table(accession = seg$accession, anchor = seg$centre, i = seq_len(nrow(seg)))
  j <- e[q, on = .(accession, anchor), roll = "nearest", .(i, hit = x.anchor)]
  out <- abs(j$hit - q$anchor) <= TOL; out[is.na(out)] <- FALSE
  out[order(j$i)]
}
## per-END verdict, the `caught` vector inside segments_roc.R's confusion()
caught_ends <- function(b, s, thr, ev) {
  sg <- segments(b, s, thr); seg <- sg$seg
  tp  <- seg$up & centre_near(seg, ev)
  tps <- seg[tp]
  if (!nrow(tps)) return(rep(FALSE, nrow(ev)))
  vapply(seq_len(nrow(ev)), function(r)
    any(tps$accession == ev$accession[r] & abs(tps$centre - ev$anchor[r]) <= TOL), logical(1))
}

zc <- fread(file.path(S, "zero_esm_precursors.csv"))$accession
un <- fread(file.path(S, "cv5_unseen_ends.csv"))
un[, set := ifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]

res <- list()
for (arm in c("", ARMS)) {
  for (tm in c("C", "N")) {
    b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, arm)))
    setDT(b); setorder(b, accession, anchor)
    b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
    for (hd in c("cnn", "xgb_pm")) {
      for (st in c("never seen: peptides", "never seen: large products")) {
        ev <- unique(un[set == st & term == tm, .(accession, anchor)])
        res[[length(res)+1]] <- data.table(
          arm = if (arm == "") "curated" else sub("^_", "", arm),
          term = tm, head = hd, set = st,
          accession = ev$accession, anchor = ev$anchor,
          caught = caught_ends(b, b[[paste0(hd, "__all50")]], THR, ev))
      }
    }
    rm(b); invisible(gc(FALSE))
    cat(sprintf("  %s %s done\n", if (arm=="") "curated" else arm, tm))
  }
}
R <- rbindlist(res)
R[, zero_esm := accession %in% zc]
fwrite(R, file.path(A2, "lf_headcmp", sprintf("paired_ends_T%s_tol%d.csv", THR, TOL)))

## ---- paired comparison -----------------------------------------------------
cat(sprintf("\n==== paired per-end comparison, T = %s, TOL = %d ====\n", THR, TOL))
cat("b = caught by the ESM arm only, c = caught by curated only; McNemar exact (binom)\n")
for (a in ARMS) {
  an <- sub("^_", "", a)
  for (restrict in c(FALSE, TRUE)) {
    cat(sprintf("\n--- %s vs curated%s ---\n", an,
                if (restrict) "  [restricted to ESM-covered precursors]" else ""))
    out <- list()
    for (st in unique(R$set)) for (tm in c("C","N")) for (hd in c("cnn","xgb_pm")) {
      x <- R[set == st & term == tm & head == hd]
      if (restrict) x <- x[zero_esm == FALSE]
      cur <- x[arm == "curated"][order(accession, anchor)]
      new <- x[arm == an][order(accession, anchor)]
      stopifnot(identical(cur$accession, new$accession), identical(cur$anchor, new$anchor))
      n <- nrow(cur)
      bb <- sum(new$caught & !cur$caught); cc <- sum(!new$caught & cur$caught)
      p  <- if (bb + cc == 0) 1 else binom.test(bb, bb + cc, 0.5)$p.value
      out[[length(out)+1]] <- data.table(
        set = st, term = tm, head = hd, ends = n,
        cur_TP = sum(cur$caught), new_TP = sum(new$caught),
        both = sum(cur$caught & new$caught), neither = sum(!cur$caught & !new$caught),
        b_new_only = bb, c_cur_only = cc, delta_TP = bb - cc,
        p_mcnemar = round(p, 4))
    }
    o <- rbindlist(out); setorder(o, set, term, head)
    print(o, row.names = FALSE, nrows = 100)
    cat(sprintf("  pooled over the 8 cells: new-only %d, curated-only %d, exact p = %.4g\n",
                sum(o$b_new_only), sum(o$c_cur_only),
                binom.test(sum(o$b_new_only), sum(o$b_new_only)+sum(o$c_cur_only), 0.5)$p.value))
    cat(sprintf("  pooled over the 4 never-seen-PEPTIDE cells: new-only %d, curated-only %d, exact p = %.4g\n",
                sum(o[grepl("peptides", set)]$b_new_only), sum(o[grepl("peptides", set)]$c_cur_only),
                binom.test(sum(o[grepl("peptides", set)]$b_new_only),
                           sum(o[grepl("peptides", set)]$b_new_only)+sum(o[grepl("peptides", set)]$c_cur_only), 0.5)$p.value))
  }
}
