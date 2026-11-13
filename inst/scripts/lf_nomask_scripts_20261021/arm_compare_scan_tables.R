## One pass over each arm's full-proteome tracks, written as the small tables the
## msa-arm figures are drawn from (Kevin, 2026-11-13: are msa / msa_gap better
## than the site's curated_AAswap, and why).
##
## Two ways of ranking an end, both among the ~3M scanned anchors that are no
## annotated end:
##
##   rank_std  metrics_final.R's: a cross-validated end by its fold's own column
##             against that column over every non-end anchor; a never-seen end
##             by all 50 members against all 50 over every non-end anchor.
##   rank_ho   the same ends against HELD-OUT negatives only: anchors no
##             training window reaches (no train-split window, positive or
##             negative, within SPAN residues in the same precursor, so the two
##             windows share no residue) of ANY compared arm, and no UniProt
##             peptide / propeptide end either (scan-pr-hold-out-other-ends).
##             A cross-validated end against its own fold's column over the
##             anchors that fold never trained near; a never-seen end against
##             all 50 over the anchors no fold trained near. One pool per fold,
##             common to every arm. Why: the msa arms train on ~4x the
##             negatives of the original windows, chosen AT motif sites, and a
##             negative a member trained on is a memorised low score -- under
##             rank_std those anchors drop out of the ranking above the ends,
##             a gain for having seen them, not for generalising.
##
## Writes to OUT (default ~/AF2_analysis/lf_viz_20261113/msa_arms_data):
##   end_ranks.csv       term, set, arm, head, accession, gene, anchor, fold,
##                       motif, score, rank_std, rank_ho
##   fp_top.csv          the top non-end anchors per arm, head and scoring:
##                       held_out = per fold, its top TOPK/5 among the anchors it
##                       never trained near (fold, rank within the fold); all50 =
##                       the top TOPK of all 50 over every non-end anchor (as
##                       deployed). accession, anchor, score, motif, near_end
##                       (within 3 residues of an annotated end)
##   pool_motif.csv      the motif make-up of the held-out pools (mean over folds)
##                       and of every non-end anchor
##   tracks_examples.csv every anchor of the GENES precursors: all50 per arm and head
##   pool_summary.csv    anchors per term: scanned, non-end, in the held-out pool
##
##   LF_SRC=<checkout> LC_ALL=en_US.UTF-8 RETICULATE_PYTHON=~/.virtualenvs/lf-plm/bin/python \
##     Rscript arm_compare_scan_tables.R
## env: ARMS         curated,curated_AAswap,curated_AAswap_msacons_padg,msa,msa_gap
##                   ("curated" = the original flatten trunk, tracks without suffix)
##      UNSEEN_SETS  _msa[,_padg]: the never-seen set is the INTERSECTION of
##                   cv5_unseen_ends<sfx>.csv over these (default _msa, plus _padg
##                   when a _padg arm is compared), so every arm is scored on the
##                   same ends and none of them overlaps any arm's negatives;
##                   "" = the original set
##      GROUPS       arm:sfx[,...] window suffixes of arms the table below lacks
##      GENES        POMC,GCG,KNG1,AGT,NPY,CARTPT,PENK,TAC1    TOPK  10000    SPAN  35
## Memory: one arm's tracks at a time (~1 GB); the run reads ~10 tracks files.
suppressMessages({library(data.table); library(reticulate)})
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
source(lf_src("R/pepend_windows.R"))
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
ARMS <- strsplit(Sys.getenv("ARMS", "curated,curated_AAswap,curated_AAswap_msacons_padg,msa,msa_gap"), ",")[[1]]
## each arm's training windows (its window-groups npz suffix, = its INPUT_SFX);
## GROUPS=arm:sfx[,arm:sfx] adds or overrides, an arm not named anywhere is taken
## to be on the original windows ("")
GSFX <- c(curated = "", flatten = "", curated_AAswap = "", curated_wAAid = "", curated_AAswap_padg = "_padg",
          curated_AAswap_msacons_padg = "_padg", msa = "_msa", msa_gap = "_msagap")
for (kv in strsplit(strsplit(Sys.getenv("GROUPS", ""), ",")[[1]], ":")) GSFX[kv[1]] <- if (length(kv) > 1) kv[2] else ""
if (length(new <- setdiff(ARMS, names(GSFX)))) {
  message("no window suffix known for ", paste(new, collapse = ", "), ": taken to be the original windows")
  GSFX[new] <- ""
}
## UNSEEN_SETS="" = the original never-seen set (cv5_unseen_ends.csv)
us_env <- Sys.getenv("UNSEEN_SETS", paste(unique(c("_msa", if (any(GSFX[ARMS] == "_padg")) "_padg")), collapse = ","))
US <- if (nzchar(us_env)) strsplit(us_env, ",")[[1]] else ""
GENES <- strsplit(Sys.getenv("GENES", "POMC,GCG,KNG1,AGT,NPY,CARTPT,PENK,TAC1"), ",")[[1]]
TOPK <- as.integer(Sys.getenv("TOPK", "10000")); SPAN <- as.integer(Sys.getenv("SPAN", "35"))
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113/msa_arms_data")))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
np <- import("numpy", convert = TRUE)
heads <- c("cnn", "xgb_pm")

prec <- fread(file.path(A2, "pepend_locus", "prec_order.csv"))
sq_p <- file.path(OUT, "prec_seqs.rds")             # the sequences, cut from the v2 feature cache once
if (!file.exists(sq_p)) {
  fc <- readRDS(file.path(A2, "lf_pepend_residue_cache_v2.rds"))
  saveRDS(setNames(fc$prec$seq, fc$prec$accession), sq_p); rm(fc); invisible(gc())
}
seqs <- readRDS(sq_p)
ll <- as.data.table(readRDS(lf_src("inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
un_all <- fread(file.path(S, "cv5_unseen_ends.csv"))     # every never-seen end, of any set
## the never-seen csvs repeat an end once per isoform name (CXCL12's N end x4), and
## the intersection join would square that (x16): one row per end before joining
un1 <- function(s) unique(fread(file.path(S, paste0("cv5_unseen_ends", s, ".csv"))), by = c("term", "accession", "anchor"))
un <- Reduce(function(a, b) a[b[, .(term, accession, anchor)], on = .(term, accession, anchor), nomatch = 0L],
             lapply(US, un1))
un[, set := fifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]
message("never-seen set = intersection of ", paste0("cv5_unseen_ends", US, ".csv", collapse = ", "), ": ",
        paste(un[, .N, by = .(term, set)][, sprintf("%s %s %d", term, set, N)], collapse = "; "))

rank_in <- function(sp, neg) findInterval(-sp, -sort(neg, decreasing = TRUE)) + 1L
tracks_path <- function(tm, a) file.path(A2, if (a %in% c("curated", "flatten")) sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm)
                                         else sprintf("lf_pepend_cv5_tracks_%s_flatten_%s.rds", tm, a))
## train-split window anchors of one arm's fold
train_windows <- function(tm, a, k) {
  g <- np$load(sprintf("%s/lf_pepend_window_groups_%s_cv5_f%d%s.npz", A2, tm, k, GSFX[[a]]), allow_pickle = TRUE)
  data.table(accession = as.character(g[[sprintf("%s__train__accession", tm)]]),
             anchor = as.integer(g[[sprintf("%s__train__anchor", tm)]]))
}

ER <- list(); FP <- list(); PM <- list(); TX <- list(); PS <- list()
for (tm in c("C", "N")) {
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)), select = c("accession", "anchor", "gene.x", "motif", "fold"))
  setnames(kn, "gene.x", "gene")
  u <- un[term == tm]
  key_of <- function(acc, anc) paste(acc, anc)
  ends_std <- unique(c(key_of(ll$accession, if (tm == "N") ll$start else ll$end), key_of(kn$accession, kn$anchor),
                       un_all[term == tm, key_of(accession, anchor)]))
  ends_ho <- unique(c(ends_std, key_of(up$accession, if (tm == "N") up$start else up$end)))
  read_tracks <- function(a) {
    x <- as.data.table(readRDS(tracks_path(tm, a)))[, c("accession", "anchor", sprintf("cnn__f%d", 1:5),
                                                         sprintf("xgb_pm__f%d", 1:5), "cnn__all50", "xgb_pm__all50"), with = FALSE]
    setkey(x, accession, anchor); x
  }
  ## ---- the anchors, their motif and flags, and which folds trained near them --------
  x <- read_tracks(ARMS[1])
  b <- x[, .(accession, anchor)]
  b[, key := key_of(accession, anchor)]
  b[, `:=`(end_std = key %in% ends_std, end_ho = key %in% ends_ho)]
  np_ <- setNames(prec$n_prot, prec$accession)
  miss <- setdiff(unique(b$accession), intersect(names(seqs), names(np_)))
  if (length(miss)) message(length(miss), " scanned precursors have no sequence or n_prot: motif NA")
  b[!accession %in% miss, motif := lf_pepend_motif_vec(seqs[[.BY$accession]], tm, anchor, np_[[.BY$accession]]),
    by = accession]
  near <- unique(rbindlist(lapply(c(-3:-1, 1:3), function(d)
    data.table(accession = c(ll$accession, kn$accession, up$accession, un_all[term == tm]$accession),
               anchor = c(if (tm == "N") ll$start else ll$end, kn$anchor, if (tm == "N") up$start else up$end,
                          un_all[term == tm]$anchor) + d))))
  b[, near_end := key %in% near[, key_of(accession, anchor)] & !end_ho]
  b[, key := NULL]
  seen <- matrix(FALSE, nrow(b), 5)                    # fold k of ANY arm trained near the anchor
  for (a in ARMS) for (k in 1:5) {
    iv <- unique(train_windows(tm, a, k)[, .(accession, lo = anchor - SPAN, hi = anchor + SPAN)])
    seen[unique(b[iv, on = .(accession, anchor >= lo, anchor <= hi), which = TRUE, nomatch = 0L]), k] <- TRUE
  }
  pool_k <- lapply(1:5, function(k) !b$end_ho & !seen[, k])
  pool_all <- !b$end_ho & rowSums(seen) == 0
  PS[[tm]] <- data.table(term = tm, scanned = nrow(b), non_end = sum(!b$end_std),
                         pool_fold_mean = round(mean(vapply(pool_k, sum, 0))), pool_all = sum(pool_all),
                         seen_by_some_fold = sum(rowSums(seen) > 0))
  u[, motif := b[.(u$accession, u$anchor), on = .(accession, anchor)]$motif]
  u[, seen_any := rowSums(seen[b[.(u$accession, u$anchor), on = .(accession, anchor), which = TRUE], , drop = FALSE]) > 0]
  if (any(u$seen_any, na.rm = TRUE)) message(tm, ": ", sum(u$seen_any, na.rm = TRUE),
                                             " never-seen ends lie within a training window's reach of some arm")
  pm_ho <- rbindlist(lapply(1:5, function(k) b[pool_k[[k]], .N, by = .(motif, near_end)][, share := N / sum(N)]))[
    , .(share = mean(share)), by = .(motif, near_end)]
  PM[[tm]] <- rbind(pm_ho[, pool := "held_out"], b[!b$end_ho, .N, by = .(motif, near_end)][, .(motif, near_end, share = N / sum(N), pool = "all50")])[
    , term := tm]
  acc_ex <- unique(c(kn[gene %in% GENES]$accession, prec[gene %in% GENES]$accession))
  ex <- b$accession %in% acc_ex
  TX[[tm]] <- b[ex, .(term = tm, accession, anchor, motif)]

  for (a in ARMS) {
    t0 <- Sys.time()
    if (a != ARMS[1]) x <- read_tracks(a)
    stopifnot(identical(x$anchor, b$anchor), identical(x$accession, b$accession))
    for (hd in heads) {
      fk <- sprintf("%s__f%d", hd, 1:5); fa <- paste0(hd, "__all50")
      ## cross-validated ends: the fold's own column; never-seen: all 50
      kk <- kn[, .(term = tm, set = "cross-validated", accession, anchor, gene, motif, fold)]
      kk[, score := as.matrix(x[.(kk$accession, kk$anchor), on = .(accession, anchor), fk, with = FALSE])[
        cbind(seq_len(nrow(kk)), kk$fold)]]
      nn <- u[, .(term = tm, set, accession, anchor, gene, motif, fold = NA_integer_)]
      nn[, score := x[.(nn$accession, nn$anchor), on = .(accession, anchor)][[fa]]]
      e <- rbind(kk, nn, use.names = TRUE)
      e[, `:=`(arm = a, head = hd, rank_std = NA_integer_, rank_ho = NA_integer_)]
      for (k in 1:5) {
        s_k <- x[[fk[k]]]
        e[set == "cross-validated" & fold == k, `:=`(
          rank_std = rank_in(score, s_k[!b$end_std & !is.na(s_k)]),
          rank_ho  = rank_in(score, s_k[pool_k[[k]] & !is.na(s_k)]))]
        ## the top anchors this fold never trained near
        p <- which(pool_k[[k]] & !is.na(s_k)); o <- p[order(-s_k[p])[seq_len(min(TOPK %/% 5L, length(p)))]]
        FP[[length(FP) + 1]] <- b[o, .(term = tm, arm = a, head = hd, scoring = "held_out", fold = k, rank = seq_along(o),
                                       accession, anchor, score = s_k[o], motif, near_end)]
      }
      s_a <- x[[fa]]
      e[set != "cross-validated", `:=`(rank_std = rank_in(score, s_a[!b$end_std & !is.na(s_a)]),
                                       rank_ho  = rank_in(score, s_a[pool_all & !is.na(s_a)]))]
      e[is.na(score), `:=`(rank_std = NA_integer_, rank_ho = NA_integer_)]
      ER[[length(ER) + 1]] <- e
      p <- which(!b$end_ho & !is.na(s_a)); o <- p[order(-s_a[p])[seq_len(min(TOPK, length(p)))]]
      FP[[length(FP) + 1]] <- b[o, .(term = tm, arm = a, head = hd, scoring = "all50", fold = NA_integer_, rank = seq_along(o),
                                     accession, anchor, score = s_a[o], motif, near_end)]
      ## the example precursors: all 50, and its rank among the non-end anchors
      ## (the arms are calibrated differently; ranks are comparable)
      TX[[tm]][, (sprintf("%s__%s__all50", a, hd)) := s_a[ex]]
      TX[[tm]][, (sprintf("%s__%s__all50_rank", a, hd)) := rank_in(s_a[ex], s_a[!b$end_std & !is.na(s_a)])]
    }
    message(sprintf("%s %s: %.0f s", tm, a, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    rm(x); invisible(gc())
  }
  rm(b, seen, pool_k, pool_all); invisible(gc())
}
ER <- rbindlist(ER, use.names = TRUE, fill = TRUE)
setcolorder(ER, c("term", "set", "arm", "head", "accession", "gene", "anchor", "fold", "motif", "score", "rank_std", "rank_ho"))
fwrite(ER, file.path(OUT, "end_ranks.csv"))
fwrite(rbindlist(FP), file.path(OUT, "fp_top.csv"))
fwrite(rbindlist(PM), file.path(OUT, "pool_motif.csv"))
TX <- rbindlist(TX, use.names = TRUE, fill = TRUE)
TX[, gene := prec$gene[match(accession, prec$accession)]]
fwrite(TX, file.path(OUT, "tracks_examples.csv"))
fwrite(rbindlist(PS), file.path(OUT, "pool_summary.csv"))
options(width = 200)
print(rbindlist(PS))
print(ER[, .(n = .N, med_std = as.double(median(rank_std, na.rm = TRUE)), med_ho = as.double(median(rank_ho, na.rm = TRUE)),
             n_ho = sum(!is.na(rank_ho))), keyby = .(term, set, head, arm)])
message("wrote ", OUT)
