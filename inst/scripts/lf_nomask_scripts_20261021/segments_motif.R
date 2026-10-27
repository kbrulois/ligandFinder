## Segment confusion at ONE threshold, split by cleavage motif (Kevin,
## 2026-10-22: "segment rocauc ... at a threshold of 0.3, subdivide the segments
## into overlapping a dibasic, monobasic, g-dibasic, or terminus").
##
## Motif per anchor is lf_pepend_motif_vec -- the rule 10_7 labels the knowns
## with: terminus beats dibasic beats G + dibasic beats monobasic, else other.
## C looks forward from the anchor (the cut is after it), N looks back.
##
## For a motif class M, at threshold T:
##   TP_M / FN_M  the EVALUATION ENDS whose own anchor has motif M; TP when a
##                positive segment's maximum lands within +/-TOL of the end,
##                FN otherwise (so TP + FN = the ends of that motif).
##   FP_M / TN_M  the SEGMENTS that overlap at least one anchor of motif M and
##                hold no real end: FP when the run is >= T and its maximum is
##                near no end, TN when the run is < T. A segment spanning two
##                motif classes is counted in both rows -- that is what
##                "overlapping" means, and a long negative run usually overlaps
##                several.
## Segments touching an annotated end from OUTSIDE the evaluation set (other
## folds' train ends, ligand_list, UniProt peptide/propeptide, the other
## never-seen set) are held out, as in segments_roc.R.
## ROC AUC at T = (1 + TPR - FPR) / 2: the area under the curve through (0,0),
## this threshold's one (FPR, TPR) point and (1,1).
## Several thresholds share one pass: loading the tracks and labelling 3M
## anchors is the slow part and does not depend on T.
## TOL is the half-width of the window in which a positive segment's maximum
## counts as landing on a known end; several values share one pass too.
##   THRESH=0.1,0.3,0.5,0.7 TOL=2 Rscript segments_motif.R
##   THRESH=0.5 TOL=4,6     Rscript segments_motif.R
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
source(lf_src("R/pepend_windows.R"))     # lf_pepend_motif_vec
TOLS <- as.integer(strsplit(Sys.getenv("TOL", "4"), ",")[[1]])
THRESH <- as.numeric(strsplit(Sys.getenv("THRESH", "0.5"), ",")[[1]])
## ARM picks which arm's tracks to read AND suffixes the output, exactly as in
## segments_roc.R. The default "" reads the curated tracks and keeps the output
## names the already-drawn tables use (TOL=2 THRESH=0.3 reproduces those).
##   THRESH=0.5 TOL=4 ARM=_esmonly36 Rscript segments_motif.R
ARM <- Sys.getenv("ARM", "")
CLASSES <- c("terminus", "dibasic", "G + dibasic", "monobasic", "other")

fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds")); prec <- as.data.table(fc$prec); rm(fc)

## TRUE where a segment's centre is within TOL of any end in `ends`
centre_near <- function(seg, ends, tol) {
  if (!nrow(seg) || !nrow(ends)) return(rep(FALSE, nrow(seg)))
  e <- unique(ends[, .(accession, anchor)]); setkey(e, accession, anchor)
  q <- data.table(accession = seg$accession, anchor = seg$centre, i = seq_len(nrow(seg)))
  j <- e[q, on = .(accession, anchor), roll = "nearest", .(i, hit = x.anchor)]
  out <- abs(j$hit - q$anchor) <= tol; out[is.na(out)] <- FALSE
  out[order(j$i)]
}
## one row per (motif class) for one profile
confusion <- function(b, s, thr, tol, ev, oth) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  d <- data.table(sid, accession = b$accession, anchor = b$anchor, s = s, up = up, motif = b$motif)
  seg <- d[, .(up = up[1], centre = anchor[which.max(s)]), by = .(sid, accession)]
  pos <- seg$up; neg <- !pos
  tp_seg <- pos & centre_near(seg, ev, tol)
  tps <- seg[tp_seg]
  caught <- if (nrow(tps)) vapply(seq_len(nrow(ev)), function(r)
    any(tps$accession == ev$accession[r] & abs(tps$centre - ev$anchor[r]) <= tol), logical(1))
    else rep(FALSE, nrow(ev))
  key <- data.table(accession = b$accession, anchor = b$anchor, sid = sid)
  ev_sid  <- key[ev,  on = .(accession, anchor), sid]
  oth_sid <- key[oth, on = .(accession, anchor), sid]
  ## a segment is out of the negative pool when it holds a real end (its end is
  ## already counted at end level), or when a positive one peaks on another end
  held <- (pos & !tp_seg & centre_near(seg, oth, tol)) |
          (neg & seg$sid %in% c(ev_sid, oth_sid))
  seg[, `:=`(fp = pos & !tp_seg & !held, tn = neg & !held)]
  ## which motif classes each segment overlaps
  pairs <- unique(d[, .(sid, motif)])
  segm  <- pairs[seg[, .(sid, fp, tn)], on = "sid"]
  by_seg <- segm[, .(FP = sum(fp), TN = sum(tn)), by = motif]
  ## the ends, by their own anchor's motif
  evm <- data.table(motif = ev$motif, caught = caught)[, .(TP = sum(caught), FN = sum(!caught)), by = motif]
  out <- merge(by_seg, evm, by = "motif", all = TRUE)
  out[is.na(out)] <- 0L
  out
}

ll <- as.data.table(readRDS(lf_src("inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
un <- fread(file.path(S, "cv5_unseen_ends.csv"))
un[, set := ifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]

res <- list()
for (tm in c("C", "N")) {
  b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, ARM)))
  setDT(b); setorder(b, accession, anchor)
  b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
  ## motif of every anchor, one precursor at a time
  b[prec, on = "accession", `:=`(seq = i.seq, n_prot = i.n_prot)]
  b[, motif := lf_pepend_motif_vec(seq[1], tm, anchor, n_prot[1]), by = accession]
  b[, c("seq", "n_prot") := NULL]
  message(sprintf("%s: %s anchors; motif %s", tm, format(nrow(b), big.mark = ","),
                  paste(sprintf("%s %s", b[, .N, by = motif]$motif,
                                format(b[, .N, by = motif]$N, big.mark = ",")), collapse = ", ")))
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))[, .(accession, anchor, fold)]
  ## the knowns' own motif comes from the same rule, so a page and this table agree
  kn[b, on = .(accession, anchor), motif := i.motif]
  unm <- unique(un[term == tm, .(set, accession, anchor)])
  unm[b, on = .(accession, anchor), motif := i.motif]
  ann <- unique(rbind(ll[, .(accession, anchor = if (tm == "N") start else end)],
                      up[, .(accession, anchor = if (tm == "N") start else end)]))
  for (hd in c("cnn", "xgb_pm")) {
    for (thr in THRESH) for (tol in TOLS) {
      cv <- rbindlist(lapply(1:5, function(k) {
        ev  <- kn[fold == k, .(accession, anchor, motif)]
        oth <- unique(rbind(ann, kn[fold != k, .(accession, anchor)], unm[, .(accession, anchor)]))
        confusion(b, b[[sprintf("%s__f%d", hd, k)]], thr, tol, ev, oth)
      }))[, lapply(.SD, sum), by = motif]
      res[[length(res) + 1]] <- data.table(term = tm, set = "cross-validated", head = hd,
                                           thr = thr, tol = tol, cv)
      for (st in c("never seen: peptides", "never seen: large products")) {
        ev  <- unm[set == st, .(accession, anchor, motif)]
        oth <- unique(rbind(ann, kn[, .(accession, anchor)], unm[set != st, .(accession, anchor)]))
        res[[length(res) + 1]] <- data.table(term = tm, set = st, head = hd, thr = thr, tol = tol,
                                             confusion(b, b[[paste0(hd, "__all50")]], thr, tol, ev, oth))
      }
      message(sprintf("  %s %-6s T=%.2f tol=%d done", tm, hd, thr, tol))
    }
  }
  rm(b); invisible(gc())
}
R <- rbindlist(res)
R <- R[motif %in% CLASSES]
R[, motif := factor(motif, levels = CLASSES)]
R[, `:=`(ends = TP + FN, TPR = TP / (TP + FN), FPR = FP / (FP + TN))]
R[, auc_T := (1 + TPR - FPR) / 2]
setorder(R, thr, tol, term, set, head, motif)
## tol 2 keeps the plain name, so the tables already drawn stay where they are
for (t in THRESH) for (v in TOLS)
  fwrite(R[thr == t & tol == v], file.path(A2, "lf_headcmp",
         sprintf("segment_motif_T%s%s%s.csv", t, ARM, if (v == 2L) "" else paste0("_tol", v))))
options(width = 220)
print(dcast(R[set == "cross-validated"], term + head + motif ~ thr + tol, value.var = "auc_T")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
