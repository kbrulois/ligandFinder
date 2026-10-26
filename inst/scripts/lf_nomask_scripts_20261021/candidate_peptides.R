## Candidate peptides, by joining an N end segment to a C end segment
## (Kevin, 2026-10-23).
##
## An END SEGMENT ("es") is what segments_roc.R scores: a run of consecutive
## anchors whose score is >= T, with its CENTRE at the run's maximum. An N es
## proposes where a peptide starts, a C es where it ends -- the N anchor IS the
## peptide's first residue and the C anchor its last, so a pair spans
## n_pos .. c_pos and the peptide is c_pos - n_pos + 1 residues long.
##
## TOL is deliberately not a filter here but a FAN-OUT: each es contributes
## every anchor within +/-TOL of its centre as a possible end, so one (N, C)
## pair yields up to (2*TOL+1)^2 overlapping candidates and the exact cut is
## left open rather than guessed. Candidates outside [LEN_MIN, LEN_MAX] drop.
##
## Scores are the ensemble mean at that anchor (`<head>__all50`, the mean over
## all 50 members -- the heavy squares on a protein page), taken from the
## terminus' own profile: an N position is scored by the N model, a C position
## by the C model.
##
##   HEAD=xgb_pm THRESH=0.5 TOL=4 LEN_MIN=8 LEN_MAX=60 Rscript candidate_peptides.R
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
source(lf_src("R/pepend_windows.R"))     # lf_pepend_motif_vec
HEAD    <- Sys.getenv("HEAD", "xgb_pm")
thr     <- as.numeric(Sys.getenv("THRESH", "0.5"))
TOL     <- as.integer(Sys.getenv("TOL", "4"))
LEN_MIN <- as.integer(Sys.getenv("LEN_MIN", "8"))
LEN_MAX <- as.integer(Sys.getenv("LEN_MAX", "60"))
COL     <- paste0(HEAD, "__all50")

## every anchor of every es, plus which es it belongs to and that es's centre
es_anchors <- function(b, s) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  d <- data.table(sid, accession = b$accession, anchor = b$anchor, s = s, up = up)[up == TRUE]
  if (!nrow(d)) return(d[0, .(accession, es, centre, centre_score)])
  cen <- d[, .(centre = anchor[which.max(s)], centre_score = max(s), n_anchor = .N), by = .(sid, accession)]
  setnames(cen, "sid", "es")
  cen
}
## the anchors a given es offers as a possible end: centre +/- TOL, scored
fan <- function(cen, b, s) {
  ## replicate whole ROWS, then add the offset -- building `pos` with rep(..., each)
  ## inside .() would leave the sibling columns recycled with `times` and every
  ## offset but the first attached to the wrong es
  g <- cen[rep(seq_len(nrow(cen)), each = 2L * TOL + 1L)]
  g[, pos := centre + rep(-TOL:TOL, times = nrow(cen))]
  stopifnot(all(abs(g$pos - g$centre) <= TOL))
  sc <- data.table(accession = b$accession, anchor = b$anchor, score = s, motif = b$motif)
  g[sc, on = .(accession, pos = anchor), `:=`(score = i.score, motif_end = i.motif)]
  g[sc, on = .(accession, centre = anchor), motif_centre := i.motif]
  g[!is.na(score)]                       # a position off the mature range has none
}

fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
prec <- as.data.table(fc$prec)[, .(accession, gene, seq, n_prot)]; rm(fc)

side <- list()
for (tm in c("N", "C")) {
  b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm)))
  setDT(b); setorder(b, accession, anchor)
  b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
  ## cleavage motif of every anchor, the rule the knowns are labelled with
  b[prec, on = "accession", `:=`(seq = i.seq, n_prot = i.n_prot)]
  b[, motif := lf_pepend_motif_vec(seq[1], tm, anchor, n_prot[1]), by = accession]
  b[, c("seq", "n_prot") := NULL]
  s <- b[[COL]]
  cen <- es_anchors(b, s)
  side[[tm]] <- fan(cen, b, s)
  message(sprintf("%s: %s end segments at T=%.2f (%s anchors offered within +/-%d)",
                  tm, format(nrow(cen), big.mark = ","), thr,
                  format(nrow(side[[tm]]), big.mark = ","), TOL))
  rm(b, s, cen); invisible(gc())
}

## pair them up on one precursor: N start before C end, length in range
setnames(side$N, c("es", "centre", "centre_score", "pos", "score", "n_anchor", "motif_end", "motif_centre"),
         c("n_es", "n_centre", "n_centre_score", "n_pos", "n_pos_score", "n_es_anchors", "n_motif_end", "n_motif_centre"))
setnames(side$C, c("es", "centre", "centre_score", "pos", "score", "n_anchor", "motif_end", "motif_centre"),
         c("c_es", "c_centre", "c_centre_score", "c_pos", "c_pos_score", "c_es_anchors", "c_motif_end", "c_motif_centre"))
P <- side$N[side$C, on = "accession", allow.cartesian = TRUE, nomatch = 0]
P[, len := c_pos - n_pos + 1L]
P <- P[len >= LEN_MIN & len <= LEN_MAX]
## two es pairs can propose the same cut; keep it once
setorder(P, accession, n_pos, c_pos, -n_pos_score, -c_pos_score)
P <- unique(P, by = c("accession", "n_pos", "c_pos"))
message(sprintf("%s candidate peptides, %s precursors",
                format(nrow(P), big.mark = ","), format(uniqueN(P$accession), big.mark = ",")))

P[prec, on = "accession", `:=`(gene = i.gene, seq = i.seq)]
P[, peptide_id := sprintf("%sx%dx%d", gene, n_pos, c_pos)]
## ---- groups of overlapping candidates -------------------------------------
## One (N es, C es) pair is one group: every candidate in it is the same pair of
## segments cut at a different offset, so they all overlap. The group is named
## by its two CENTRES -- i.e. the id of the centre-to-centre peptide, which is
## the group's own `peptide_id` when that cut is in the length range.
## `pair_score` treats the two ends as independent evidence, so `group_rank` 1
## is the cut in that group where BOTH ends are strongest; taking rank 1 per
## group collapses the table to one peptide per pair of segments.
P[, group_id := sprintf("%sx%dx%d", gene, n_centre, c_centre)]
P[, pair_score := n_pos_score * c_pos_score]
P[, group_size := .N, by = group_id]
P[, group_rank := frank(-pair_score, ties.method = "first"), by = group_id]
P[, sequence := substr(seq, n_pos, c_pos)]
P[, `:=`(n_dist = n_pos - n_centre, c_dist = c_pos - c_centre)]
## the stronger of the two ends. pair_score (their product) rewards cuts where
## BOTH ends are good; this surfaces a candidate carried by one confident end.
P[, end_score_max := round(pmax(n_pos_score, c_pos_score), 4)]
## how far this cut sits from the two es centres in total. 0 means both ends
## are exactly on their segment maxima; it rises to 2*TOL at the corners of the
## fan-out, so it says how far a candidate has been pushed from what the model
## actually pointed at.
P[, total_offset := abs(n_pos - n_centre) + abs(c_pos - c_centre)]

## two motif columns per terminus: the es's own type (its centre, which is what
## names the segment) and the type at THIS candidate's cut, which is what you
## would actually filter a candidate on
out <- P[, .(peptide_id, group_id, group_size, group_rank, pair_score = round(pair_score, 4),
             gene, accession, start = n_pos, end = c_pos, length = len, sequence,
             end_score_max, total_offset,
             n_dist, n_score_end = round(n_pos_score, 4), n_score_centre = round(n_centre_score, 4),
             n_motif_centre, n_motif_end,
             c_dist, c_score_end = round(c_pos_score, 4), c_score_centre = round(c_centre_score, 4),
             c_motif_centre, c_motif_end,
             n_centre, c_centre, n_es_anchors, c_es_anchors)]
setorder(out, -n_score_centre, -c_score_centre, group_id, group_rank)
f <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_T%s_tol%d.csv", HEAD, thr, TOL))
fwrite(out, f)
message(sprintf("wrote %s", f))
options(width = 230)
cat("\ntop 15 by the two es centre scores:\n")
print(out[1:12, .(peptide_id, length, n_dist, n_score_end, n_score_centre, n_motif_centre,
                  c_dist, c_score_end, c_score_centre, c_motif_centre)], row.names = FALSE)
cat("\nes centre motif, N x C (candidates):\n")
print(dcast(out[, .N, by = .(n_motif_centre, c_motif_centre)], n_motif_centre ~ c_motif_centre,
            value.var = "N", fill = 0), row.names = FALSE)
cat(sprintf("\nlength: median %d, range %d-%d | per precursor: median %d, max %d\n",
            as.integer(median(out$length)), min(out$length), max(out$length),
            as.integer(median(out[, .N, by = accession]$N)), max(out[, .N, by = accession]$N)))
g <- out[, .N, by = group_id]
cat(sprintf("groups (one per N es x C es pair): %s | size median %d, max %d | rank-1 rows = %s\n",
            format(nrow(g), big.mark = ","), as.integer(median(g$N)), max(g$N),
            format(nrow(out[group_rank == 1]), big.mark = ",")))
cat("\nthe EDN2 group, every cut it offers:\n")
print(out[group_id == out[gene == "EDN2"][1]$group_id][order(group_rank),
          .(peptide_id, group_rank, length, n_dist, n_score_end, c_dist, c_score_end, pair_score)][1:9],
      row.names = FALSE)
