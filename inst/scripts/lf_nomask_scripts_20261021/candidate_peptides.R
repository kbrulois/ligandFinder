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
##
## ARM picks which arm's tracks to segment, exactly as in segments_roc.R: ""
## (the default) is the curated arm, "_esmonly36" reads
## lf_pepend_cv5_tracks_<N|C>_flatten_esmonly36.rds. It also suffixes both
## output names, so arms never overwrite each other. An arm's threshold is its
## own (arms/<ARM>.sh THRESH_DEFAULT: 0.65 for esmonly36) -- pass it as THRESH.
##   ARM=_esmonly36 THRESH=0.65 Rscript candidate_peptides.R
##
## ARMS_MAX segments SEVERAL arms at once (Kevin, 2026-10-27): the profile is,
## per anchor, the MAX over the listed arms of each arm's all-50 mean, and every
## es records which arms it rests on -- an arm SUPPORTS an es when its own score
## is at or above THRESH somewhere in the es's run -- as `n_arms` (1, 2, 3...)
## and `arms` ("curated+esmonly36"), plus `centre_arm`, the arm holding the max
## at the centre. Candidates carry the same per end (n_es_arms, c_es_arms, ...).
## "curated" names the curated arm. One THRESH for the max (an arm's own
## threshold, e.g. esmonly36's 0.65, does not apply here). ARM is required with
## it, as the output suffix, so this mode can never overwrite a one-arm table.
##   ARMS_MAX=curated,esmonly36,curated_AAswap ARM=_max3 Rscript candidate_peptides.R
##
## Also writes the END SEGMENTS themselves, one row per es with its run of
## anchors (start..end) and centre, as end_segments_<HEAD>_T<THRESH>.csv -- the
## protein page draws them, and reading them back from here (rather than
## re-deriving them on the page) guarantees the page shows the segments the
## candidates were actually built from. No TOL in that name: an es does not
## depend on it. Segments with no partner in length range are kept.
##
## OUT_DIR (default ~/AF2_analysis/lf_headcmp) redirects both files. Stages 3-5
## add their columns to the candidate CSV IN PLACE, so re-running this stage
## alone over an annotated table would drop them -- point OUT_DIR elsewhere to
## regenerate segments for a table that already exists.
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
ARM     <- Sys.getenv("ARM", "")
ARMS_MAX <- Filter(nzchar, trimws(strsplit(Sys.getenv("ARMS_MAX", ""), ",")[[1]]))
if (length(ARMS_MAX) && !nzchar(ARM))
  stop("ARMS_MAX needs an ARM suffix of its own (e.g. ARM=_max3): without one it would ",
       "overwrite the curated arm's table")
arm_sfx <- function(a) if (a == "curated") "" else paste0("_", a)
OUT_DIR <- path.expand(Sys.getenv("OUT_DIR", file.path(A2, "lf_headcmp")))

## every es: which run it is, its first and last anchor, and its centre. With M
## (one column of scores per arm, rows aligned with b), also which arms support
## it: an arm whose own score reaches thr anywhere in the run.
es_anchors <- function(b, s, M = NULL) {
  up <- !is.na(s) & s >= thr
  sid <- cumsum(b$brk0 | c(TRUE, up[-1] != up[-length(up)]))
  keep <- which(up)
  d <- data.table(sid = sid[keep], accession = b$accession[keep], anchor = b$anchor[keep], s = s[keep])
  if (!nrow(d)) return(d[0, .(accession, es = sid, centre = anchor, centre_score = s,
                              n_anchor = integer(0), seg_start = anchor, seg_end = anchor)])
  cen <- d[, .(centre = anchor[which.max(s)], centre_score = max(s), n_anchor = .N,
               seg_start = anchor[1], seg_end = anchor[.N]), by = .(sid, accession)]
  if (!is.null(M)) {
    hit <- M[keep, , drop = FALSE] >= thr
    hit[is.na(hit)] <- FALSE
    sup <- rowsum(hit * 1L, d$sid) > 0                  # es x arm
    sup <- sup[match(cen$sid, as.integer(rownames(sup))), , drop = FALSE]
    stopifnot(all(rowSums(sup) >= 1))                   # the max reached thr, so some arm did
    ci <- keep[d[, .I[which.max(s)], by = sid]$V1]      # the centre's row in b
    cen[, `:=`(es_arms = as.integer(rowSums(sup)),
               es_arm_list = apply(sup, 1, function(r) paste(colnames(M)[r], collapse = "+")),
               centre_arm = colnames(M)[max.col(M[ci, , drop = FALSE], ties.method = "first")])]
  }
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

side <- list(); segs <- list()
for (tm in c("N", "C")) {
  M <- NULL
  if (length(ARMS_MAX)) {
    ## every arm's tracks, which must agree row for row on (accession, anchor)
    bs <- lapply(ARMS_MAX, function(a) {
      x <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds",
                                                       tm, arm_sfx(a)))))[, c("accession", "anchor", COL), with = FALSE]
      setorder(x, accession, anchor); x
    })
    for (x in bs[-1]) stopifnot(identical(x$accession, bs[[1]]$accession),
                                identical(x$anchor, bs[[1]]$anchor))
    M <- vapply(bs, function(x) x[[COL]], numeric(nrow(bs[[1]])))
    colnames(M) <- ARMS_MAX
    b <- bs[[1]][, .(accession, anchor)]
    b[, (COL) := do.call(pmax, c(lapply(seq_len(ncol(M)), function(j) M[, j]), na.rm = TRUE))]
    rm(bs); invisible(gc())
  } else {
    b <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, ARM)))
    setDT(b); setorder(b, accession, anchor)
  }
  b[, brk0 := c(TRUE, accession[-1] != accession[-.N] | diff(anchor) != 1L)]
  ## cleavage motif of every anchor, the rule the knowns are labelled with
  b[prec, on = "accession", `:=`(seq = i.seq, n_prot = i.n_prot)]
  b[, motif := lf_pepend_motif_vec(seq[1], tm, anchor, n_prot[1]), by = accession]
  b[, c("seq", "n_prot") := NULL]
  s <- b[[COL]]
  cen <- es_anchors(b, s, M)
  segs[[tm]] <- cen[, .(term = tm, accession, es, start = seg_start, end = seg_end,
                        centre, centre_score, n_anchors = n_anchor)]
  if (!is.null(M))
    segs[[tm]][, `:=`(n_arms = cen$es_arms, arms = cen$es_arm_list, centre_arm = cen$centre_arm)]
  segs[[tm]][b, on = .(accession, centre = anchor), motif_centre := i.motif]
  cen[, c("seg_start", "seg_end") := NULL]
  side[[tm]] <- fan(cen, b, s)
  message(sprintf("%s: %s end segments at T=%.2f (%s anchors offered within +/-%d)",
                  tm, format(nrow(cen), big.mark = ","), thr,
                  format(nrow(side[[tm]]), big.mark = ","), TOL))
  if (!is.null(M)) print(cen[, .(terminus = tm, .N), by = .(n_arms = es_arms)][order(n_arms)], row.names = FALSE)
  rm(b, s, cen, M); invisible(gc())
}

## pair them up on one precursor: N start before C end, length in range
## (the arm columns exist only under ARMS_MAX, hence skip_absent)
setnames(side$N, c("es", "centre", "centre_score", "pos", "score", "n_anchor", "motif_end", "motif_centre",
                   "es_arms", "es_arm_list", "centre_arm"),
         c("n_es", "n_centre", "n_centre_score", "n_pos", "n_pos_score", "n_es_anchors", "n_motif_end",
           "n_motif_centre", "n_es_arms", "n_es_arm_list", "n_centre_arm"), skip_absent = TRUE)
setnames(side$C, c("es", "centre", "centre_score", "pos", "score", "n_anchor", "motif_end", "motif_centre",
                   "es_arms", "es_arm_list", "centre_arm"),
         c("c_es", "c_centre", "c_centre_score", "c_pos", "c_pos_score", "c_es_anchors", "c_motif_end",
           "c_motif_centre", "c_es_arms", "c_es_arm_list", "c_centre_arm"), skip_absent = TRUE)
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
## under ARMS_MAX: how many arms (and which) each end's segment rests on --
## row-aligned with P, so assigned before the sort
if (length(ARMS_MAX))
  out[, `:=`(n_es_arms = P$n_es_arms, n_es_arm_list = P$n_es_arm_list, n_centre_arm = P$n_centre_arm,
             c_es_arms = P$c_es_arms, c_es_arm_list = P$c_es_arm_list, c_centre_arm = P$c_centre_arm)]
setorder(out, -n_score_centre, -c_score_centre, group_id, group_rank)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUT_DIR, sprintf("candidate_peptides_%s_T%s_tol%d%s.csv", HEAD, thr, TOL, ARM))
fwrite(out, f)
message(sprintf("wrote %s", f))
## the segments, with how many candidates each one ended up in
seg <- rbindlist(segs)
seg[prec, on = "accession", gene := i.gene]
## an es is identified by (term, accession, centre): runs are disjoint, so no
## two share a centre
nc <- rbind(out[, .(term = "N", accession, centre = n_centre)],
            out[, .(term = "C", accession, centre = c_centre)])[, .(n_candidates = .N),
                                                                by = .(term, accession, centre)]
seg[nc, on = .(term, accession, centre), n_candidates := i.n_candidates]
seg[is.na(n_candidates), n_candidates := 0L]
setcolorder(seg, c("gene", "accession", "term", "es", "start", "end", "centre", "centre_score",
                   "n_anchors", "motif_centre", "n_candidates"))
setorder(seg, accession, term, start)
fs <- file.path(OUT_DIR, sprintf("end_segments_%s_T%s%s.csv", HEAD, thr, ARM))
fwrite(seg[, centre_score := round(centre_score, 4)], fs)
message(sprintf("wrote %s (%s segments, %s in no candidate)", fs, format(nrow(seg), big.mark = ","),
                format(seg[n_candidates == 0, .N], big.mark = ",")))
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
if (length(ARMS_MAX)) {
  cat("\ncandidates by how many arms support each end (rows N, columns C):\n")
  print(dcast(out[, .N, by = .(n_es_arms, c_es_arms)], n_es_arms ~ c_es_arms, value.var = "N", fill = 0),
        row.names = FALSE)
}
cat("\nthe EDN2 group, every cut it offers:\n")
print(out[group_id == out[gene == "EDN2"][1]$group_id][order(group_rank),
          .(peptide_id, group_rank, length, n_dist, n_score_end, c_dist, c_score_end, pair_score)][1:9],
      row.names = FALSE)
