## Stage 6 of the candidate pipeline: an xgboost classifier over end, MultiPep
## and MMseqs evidence, and the scored candidate CSV (Kevin, 2026-10-27).
##
##   HEAD=xgb_pm THRESH=0.5 TOL=4 ARM=_arms3 \
##   ARMS_MAX=curated:0.5,esmonly36:0.65,curated_AAswap:0.5 Rscript candidate_xgb.R
##
## OUTCOME: the candidate IS a known peptide on its own precursor -- a curated
## ligand (ligand_list) or a UniProt `peptide` feature whose start AND end are
## each within LABEL_TOL (2) residues of the candidate's. A candidate matching
## only a PROPEPTIDE (a real cleavage product, not a signalling peptide) is left
## out of training -- neither positive nor a fair negative -- but still scored.
##
## THREE LEAKS THIS STAGE IS BUILT AROUND
## 1. MMseqs. A candidate on a known peptide's precursor always finds that very
##    peptide, so the self-inclusive hit columns ARE the label. Only the
##    mmseqs_ns_* columns (best hit on ANOTHER precursor) and *_ns_bits_per_res
##    enter; an assertion below refuses anything else from MMseqs, and every
##    overlap_* column (defined from the same annotations as the label).
## 2. End scores. The all-50 mean includes the fold that TRAINED on each known
##    end, so it is inflated exactly where the label is 1. The features use each
##    arm's HELD-OUT fold instead: on a precursor with CV knowns, the fold that
##    held out its nearest known end (lf_pepend_cv5_knowns_<N|C>.csv); elsewhere
##    no fold trained on a positive there, and the all-50 mean is used. That
##    excludes everything built from the all-50 profile -- the end scores, the
##    centre and distance-to-centre, the segment width and its arm count --
##    replaced by held-out analogues: each arm's score at the cut, how many arms
##    clear their own threshold there, and how far the cut stands above its
##    +/-2 neighbours (a held-out "is this the peak").
## 3. Correlated rows. Up to 81 near-identical candidates share a segment pair,
##    so folds are grouped BY GENE; every candidate's xgb_score comes from the
##    one model that never saw its gene (out of fold). No early stopping on the
##    held-out fold: a fixed number of rounds.
##
## Writes, beside the table:
##   xgb_score / xgb_fold / label / label_match onto the candidate table (in place)
##   xgb_nomp_score                        the same classifier, folds and rounds with every
##                                         MultiPep (mp_*) feature removed: end, motif and
##                                         MMseqs evidence only
##   curated_relation (exact | overlap | uniprot | none) and its curated_* details: the
##                                         curated known (ligand_list) of 8-60 aa on the
##                                         same precursor the candidate matches or overlaps;
##                                         "uniprot" = no curated one, but a UniProt peptide
##                                         of 8-60 aa there (so "none" = no known peptide)
##   uniprot_relation (exact | overlap | none) and its uniprot_* details: the same
##                                         against UniProt `peptide` features of 8-60 aa
##   candidate_peptides_<...>_scored.csv  the tidy table: ids, end, MultiPep and
##                                         MMseqs metrics, and xgb_score, best first
##   candidate_xgb_<...>.model             the full-data model, for new candidates
##   candidate_xgb_nomp_<...>.model        the same without MultiPep
## (stage 7, candidate_best_per_group.R, cuts the scored table to the best cuts per group)
suppressMessages({library(data.table); library(xgboost)})
A2 <- path.expand("~/AF2_analysis")
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
HEAD <- Sys.getenv("HEAD", "xgb_pm"); ARM <- Sys.getenv("ARM", "")
thr  <- as.numeric(Sys.getenv("THRESH", "0.5")); TOL <- as.integer(Sys.getenv("TOL", "4"))
LABEL_TOL <- as.integer(Sys.getenv("LABEL_TOL", "2"))
NFOLD     <- as.integer(Sys.getenv("NFOLD", "5"))
NROUNDS   <- as.integer(Sys.getenv("XGB_ROUNDS", "300"))
SEED      <- as.integer(Sys.getenv("XGB_SEED", "1"))
key  <- sprintf("%s_T%s_tol%d%s", HEAD, thr, TOL, ARM)
f    <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s.csv", key))
f_sc <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_scored.csv", key))
f_md <- file.path(A2, "lf_headcmp", sprintf("candidate_xgb_%s.model", key))
f_md_nomp <- file.path(A2, "lf_headcmp", sprintf("candidate_xgb_nomp_%s.model", key))

## the arms whose scores are features: ARMS_MAX's (with their own thresholds),
## else the one arm this table was cut from
am <- Filter(nzchar, trimws(strsplit(Sys.getenv("ARMS_MAX", ""), ",")[[1]]))
if (!length(am)) am <- if (nzchar(ARM)) sub("^_", "", ARM) else "curated"
am <- strsplit(am, ":", fixed = TRUE)
ARMS <- vapply(am, `[`, "", 1)
ARM_THR <- stats::setNames(vapply(am, function(x) if (length(x) > 1) as.numeric(x[2]) else thr, 0), ARMS)
arm_sfx <- function(a) if (a == "curated") "" else paste0("_", a)

cd <- fread(f)
rel_cols <- c("relation", "match", "match_id", "match_len", "shared_res",
              "frac_of_known", "frac_of_cand", "start_diff", "end_diff")
out_cols <- c("xgb_score", "xgb_nomp_score", "xgb_fold", "label", "label_match", "label_type",
              paste0("curated_", rel_cols), paste0("uniprot_", rel_cols))
if (length(intersect(names(cd), out_cols))) cd[, (intersect(names(cd), out_cols)) := NULL]
cd[, row := .I]
cat(sprintf("%s candidates on %s precursors; arms: %s\n", format(nrow(cd), big.mark = ","),
            format(uniqueN(cd$accession), big.mark = ","),
            paste(sprintf("%s>=%s", ARMS, ARM_THR), collapse = ", ")))

## ---- the label ----------------------------------------------------------------
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds")); prec <- as.data.table(fc$prec)[, .(accession, gene)]; rm(fc)
ll <- as.data.table(readRDS(lf_src("inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
ll[, nm := fifelse(!is.na(final_name) & nzchar(final_name), final_name, uniprot_name)]
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
kn <- rbind(ll[, .(accession, kstart = as.integer(start), kend = as.integer(end), kname = nm, ktype = "ligand_list")],
            up[, .(accession, kstart = as.integer(start), kend = as.integer(end),
                   kname = fifelse(!is.na(description) & nzchar(description), description, type), ktype = type)])
kn <- unique(kn[!is.na(kstart) & !is.na(kend)])
j <- kn[cd[, .(row, accession, start, end)], on = "accession", allow.cartesian = TRUE, nomatch = NULL][
  abs(kstart - start) <= LABEL_TOL & abs(kend - end) <= LABEL_TOL]
j[, sig := ktype %in% c("ligand_list", "peptide")]
lab <- j[order(-sig, abs(kstart - start) + abs(kend - end)), .(
  sig = any(sig), label_match = kname[1], label_type = ktype[1]), by = row]
cd[lab, on = "row", `:=`(label_match = i.label_match, label_type = i.label_type,
                         label = fifelse(i.sig, 1L, NA_integer_))]
cd[is.na(label) & is.na(label_type), label := 0L]
cat(sprintf("label (both ends within +/-%d of a known peptide on the same precursor): %s positive, %s propeptide-only (left out of training), %s negative\n",
            LABEL_TOL, format(cd[label %in% 1L, .N], big.mark = ","),
            format(cd[is.na(label), .N], big.mark = ","), format(cd[label %in% 0L, .N], big.mark = ",")))

## ---- curated knowns 8-60 aa on the same precursor: exact or overlapping ----------
## Reported, never a feature (it is the label's own annotation). Only curated
## ligands (ligand_list) of peptide length (KNOWN_LEN_MIN..MAX, 8-60): exact =
## same start and end; overlap = any shared residue. Where several overlap,
## the exact one, else the one sharing most residues, is reported.
CUR_MIN <- as.integer(Sys.getenv("KNOWN_LEN_MIN", "8"))
CUR_MAX <- as.integer(Sys.getenv("KNOWN_LEN_MAX", "60"))
cur <- unique(ll[, .(accession, kstart = as.integer(start), kend = as.integer(end), kname = nm)])[
  kend - kstart + 1L >= CUR_MIN & kend - kstart + 1L <= CUR_MAX]
cur[prec, on = "accession", kgene := i.gene]
cur[, kid := sprintf("%sx%dx%d", fifelse(is.na(kgene), accession, kgene), kstart, kend)]
jo <- cur[cd[, .(row, accession, start, end, length)], on = "accession", allow.cartesian = TRUE,
          nomatch = NULL][pmin(end, kend) >= pmax(start, kstart)]
jo[, `:=`(shared = pmin(end, kend) - pmax(start, kstart) + 1L, exact = start == kstart & end == kend)]
setorder(jo, row, -exact, -shared, kstart)
jb <- jo[jo[, .I[1], by = row]$V1]
cd[, curated_relation := "none"]
cd[jb, on = "row", `:=`(curated_relation = fifelse(i.exact, "exact", "overlap"),
                        curated_match = i.kname, curated_match_id = i.kid,
                        curated_match_len = i.kend - i.kstart + 1L, curated_shared_res = i.shared,
                        curated_frac_of_known = round(i.shared / (i.kend - i.kstart + 1L), 3),
                        curated_frac_of_cand = round(i.shared / length, 3),
                        curated_start_diff = start - i.kstart, curated_end_diff = end - i.kend)]
cat(sprintf("curated knowns %d-%d aa (%d on %d precursors): %s candidates match one exactly, %s overlap one, %s neither\n",
            CUR_MIN, CUR_MAX, nrow(cur), uniqueN(cur$accession),
            format(cd[curated_relation == "exact", .N], big.mark = ","),
            format(cd[curated_relation == "overlap", .N], big.mark = ","),
            format(cd[curated_relation == "none", .N], big.mark = ",")))

## ---- UniProt peptides 8-60 aa on the same precursor: exact or overlapping ----------
## The same comparison against UniProt `peptide` features (not propeptides, not
## chains), into uniprot_*, whatever the curated relation. And a candidate with
## no curated relation that matches or overlaps one moves out of "none" into a
## curated_relation category of its own, "uniprot" -- so "none" means no known
## peptide, curated or UniProt, at any of its residues (Kevin, 2026-10-28: "none"
## let 222 UniProt peptides into the first novel-candidate gallery). Reported,
## never a feature: it is the label's own annotation. Length window
## UNIPROT_LEN_MIN..MAX, default the curated one.
UP_MIN <- as.integer(Sys.getenv("UNIPROT_LEN_MIN", CUR_MIN))
UP_MAX <- as.integer(Sys.getenv("UNIPROT_LEN_MAX", CUR_MAX))
upk <- unique(up[type == "peptide", .(accession, kstart = as.integer(start), kend = as.integer(end),
                                      kname = fifelse(!is.na(description) & nzchar(description),
                                                      description, "peptide"))])[
  !is.na(kstart) & !is.na(kend) & kend - kstart + 1L >= UP_MIN & kend - kstart + 1L <= UP_MAX]
upk[prec, on = "accession", kgene := i.gene]
upk[, kid := sprintf("%sx%dx%d", fifelse(is.na(kgene), accession, kgene), kstart, kend)]
ju <- upk[cd[, .(row, accession, start, end, length)], on = "accession", allow.cartesian = TRUE,
          nomatch = NULL][pmin(end, kend) >= pmax(start, kstart)]
ju[, `:=`(shared = pmin(end, kend) - pmax(start, kstart) + 1L, exact = start == kstart & end == kend)]
setorder(ju, row, -exact, -shared, kstart)
jub <- ju[ju[, .I[1], by = row]$V1]
cd[, uniprot_relation := "none"]
cd[jub, on = "row", `:=`(uniprot_relation = fifelse(i.exact, "exact", "overlap"),
                         uniprot_match = i.kname, uniprot_match_id = i.kid,
                         uniprot_match_len = i.kend - i.kstart + 1L, uniprot_shared_res = i.shared,
                         uniprot_frac_of_known = round(i.shared / (i.kend - i.kstart + 1L), 3),
                         uniprot_frac_of_cand = round(i.shared / length, 3),
                         uniprot_start_diff = start - i.kstart, uniprot_end_diff = end - i.kend)]
cd[curated_relation == "none" & uniprot_relation != "none", curated_relation := "uniprot"]
cat(sprintf("UniProt peptides %d-%d aa (%d on %d precursors): %s candidates match one exactly, %s overlap one; %s with no curated relation -> curated_relation \"uniprot\", %s left \"none\"\n",
            UP_MIN, UP_MAX, nrow(upk), uniqueN(upk$accession),
            format(cd[uniprot_relation == "exact", .N], big.mark = ","),
            format(cd[uniprot_relation == "overlap", .N], big.mark = ","),
            format(cd[curated_relation == "uniprot", .N], big.mark = ","),
            format(cd[curated_relation == "none", .N], big.mark = ",")))

## ---- held-out end scores, per arm ---------------------------------------------
## For every anchor: the fold that held out its precursor's NEAREST known end
## (precursors with CV knowns), else none (all-50). One fold column per anchor.
ho_profile <- function(tm) {
  kf <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))[!is.na(fold), .(accession, anchor, fold)]
  setkey(kf, accession, anchor)
  P <- NULL
  for (a in ARMS) {
    b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, arm_sfx(a)))))
    fcols <- grep(sprintf("^%s__f[0-9]+$", HEAD), names(b), value = TRUE)
    b <- b[, c("accession", "anchor", paste0(HEAD, "__all50"), fcols), with = FALSE]
    setkey(b, accession, anchor)
    if (is.null(P)) {
      ## nearest known end's fold, only on precursors that have one
      P <- b[, .(accession, anchor)]
      P[, fold := kf[P, on = .(accession, anchor), roll = "nearest", x.fold]]
      P[!accession %in% kf$accession, fold := NA_integer_]
    } else stopifnot(identical(P$accession, b$accession), identical(P$anchor, b$anchor))
    v <- b[[paste0(HEAD, "__all50")]]
    for (k in sort(unique(stats::na.omit(P$fold)))) {
      w <- which(P$fold == k)
      v[w] <- b[[sprintf("%s__f%d", HEAD, k)]][w]
    }
    set(P, j = paste0("ho_", a), value = v)
    rm(b); invisible(gc())
  }
  hc <- paste0("ho_", ARMS)
  P[, ho_max := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = hc]
  hit <- sweep(as.matrix(P[, hc, with = FALSE]), 2, ARM_THR[ARMS], ">=")
  hit[is.na(hit)] <- FALSE
  P[, ho_arms := as.integer(rowSums(hit))]
  ## how far the cut stands above its +/-2 neighbours on the held-out max
  ## profile (within a precursor; positive = a local peak)
  P[, nb := pmax(shift(ho_max, 1), shift(ho_max, 2), shift(ho_max, -1), shift(ho_max, -2), na.rm = TRUE),
    by = accession]
  P[, ho_peak := ho_max - nb][, nb := NULL]
  P[, fold := NULL]
  P
}
for (tm in c("N", "C")) {
  P <- ho_profile(tm)
  pre <- tolower(tm)
  pos <- if (tm == "N") "start" else "end"
  nc <- setdiff(names(P), c("accession", "anchor"))
  cd[P, on = c("accession", setNames("anchor", pos)), (paste0(pre, "_", nc)) := mget(paste0("i.", nc))]
  rm(P); invisible(gc())
}

## ---- features ---------------------------------------------------------------------
mp <- setdiff(grep("^mp_", names(cd), value = TRUE), "mp_signalling_max")
if (length(mp)) {
  M <- as.matrix(cd[, mp, with = FALSE]); k <- max.col(replace(M, is.na(M), -Inf), ties.method = "first")
  cd[, `:=`(mp_best_class = sub("^mp_", "", mp[k]), mp_best_score = M[cbind(seq_len(nrow(M)), k)])]
}
cd[, mmseqs_ns_log10_evalue := pmin(log10(pmax(mmseqs_ns_evalue, 1e-30)), 4)]
for (tt in c("ligand_list", "peptide", "propeptide", "none"))
  set(cd, j = paste0("ns_type_", tt), value = as.integer(cd$mmseqs_ns_hit_type == tt))
motifs <- sort(unique(c(cd$n_motif_end, cd$c_motif_end)))
for (tm in c("n", "c")) for (mo in motifs)
  set(cd, j = sprintf("%s_motif_%s", tm, gsub("[^A-Za-z0-9]+", "_", mo)), value = as.integer(cd[[paste0(tm, "_motif_end")]] == mo))
feat <- c("length",
          grep("^[nc]_ho_", names(cd), value = TRUE),
          grep("^[nc]_motif_(?!end$|centre$)", names(cd), value = TRUE, perl = TRUE),
          mp, "mp_signalling_max",
          "mmseqs_ns_bits_per_res", "mmseqs_ns_sim", "mmseqs_ns_pident", "mmseqs_ns_tcov", "mmseqs_ns_qcov",
          "mmseqs_ns_log10_evalue", grep("^ns_type_", names(cd), value = TRUE),
          grep("^(ll|pep|prop)_ns_bits_per_res$", names(cd), value = TRUE))
feat <- intersect(feat, names(cd))
## LEAKAGE GUARD: nothing self-inclusive from MMseqs, nothing from the overlap
## annotation, no other ranker, no all-50 end score
banned <- grep("^(overlap_|curated_|uniprot_|lr_|label)|^mmseqs_(?!ns_)|^(ll|pep|prop)_(?!ns_)|_score_end$|_score_centre$|^end_score_max$|^pair_score$|_es_arms$|_dist$|total_offset|_es_anchors$|^[nc]_score_",
               feat, value = TRUE, perl = TRUE)
if (length(banned)) stop("leakage guard: these must not be predictors: ", paste(banned, collapse = ", "))
cat(sprintf("%d features: %d held-out end, %d motif, %d MultiPep, %d MMseqs (other precursors only)\n",
            length(feat), length(grep("_ho_", feat)), length(grep("_motif_", feat)),
            length(grep("^mp_", feat)), length(grep("ns_", feat))))

## ---- folds by gene, positives spread evenly -------------------------------------
set.seed(SEED)
gpos <- cd[, .(pos = any(label %in% 1L)), by = gene]
gpos[, r := runif(.N)]
setorder(gpos, -pos, r)
gpos[, fold := ((seq_len(.N) - 1L) %% NFOLD) + 1L, by = pos]
cd[gpos, on = "gene", xgb_fold := i.fold]
cat("positives per fold:", paste(cd[label %in% 1L, .N, keyby = xgb_fold]$N, collapse = ", "), "\n")

X <- as.matrix(cd[, feat, with = FALSE]); storage.mode(X) <- "double"
y <- cd$label
params <- list(objective = "binary:logistic", eval_metric = "aucpr", eta = 0.05, max_depth = 5,
               min_child_weight = 5, subsample = 0.8, colsample_bytree = 0.8, lambda = 2, nthread = 8)
auc <- function(s, l) { r <- rank(s); n1 <- sum(l == 1); n0 <- sum(l == 0)
  if (!n1 || !n0) return(NA_real_); (sum(r[l == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
ap  <- function(s, l) { o <- order(-s); l <- l[o]; if (!sum(l)) return(NA_real_)
  sum(cumsum(l)[l == 1] / which(l == 1)) / sum(l) }
fit <- function(rows, fs) {
  yy <- y[rows]
  d  <- xgb.DMatrix(X[rows, fs, drop = FALSE], label = yy)
  ## (the R package takes its randomness from set.seed, not a `seed` param)
  xgb.train(c(params, list(scale_pos_weight = sum(yy == 0) / max(1, sum(yy == 1)))),
            d, nrounds = NROUNDS, verbose = 0)
}
## out-of-fold score over the feature set fs into column sc; the per-fold table back
cv_score <- function(fs, sc) {
  cd[, (sc) := NA_real_]
  res <- list()
  for (k in seq_len(NFOLD)) {
    tr <- which(cd$xgb_fold != k & !is.na(y)); te <- which(cd$xgb_fold == k)
    m <- fit(tr, fs)
    set(cd, i = te, j = sc, value = predict(m, xgb.DMatrix(X[te, fs, drop = FALSE])))
    ev <- te[!is.na(y[te])]
    res[[k]] <- data.table(fold = k, n = length(ev), pos = sum(y[ev] == 1),
                           auc = auc(cd[[sc]][ev], y[ev]), ap = ap(cd[[sc]][ev], y[ev]),
                           base = mean(y[ev] == 1))
    cat(sprintf("  fold %d: trained on %s rows (%d pos), scored %s; held-out ROC AUC %.3f, AP %.3f (base rate %.4f)\n",
                k, format(length(tr), big.mark = ","), sum(y[tr] == 1), format(length(te), big.mark = ","),
                res[[k]]$auc, res[[k]]$ap, res[[k]]$base))
  }
  res <- rbindlist(res)
  cat(sprintf("%s PER-FOLD held-out (gene-grouped): ROC AUC %.3f +/- %.3f, AP %.3f +/- %.3f (base rate %.4f)\n",
              sc, mean(res$auc), sd(res$auc), mean(res$ap), sd(res$ap), mean(res$base)))
  res
}
feat_nomp <- grep("^mp_", feat, value = TRUE, invert = TRUE)
## xgb_score is drawn exactly as before (CV, then the full model, one RNG
## stream); the no-MultiPep score replays that stream from the same point, so
## the only difference between the two is the features
rng0 <- .Random.seed
cat(sprintf("xgb_score: all %d features\n", length(feat)))
res <- cv_score(feat, "xgb_score")
ev <- which(!is.na(y))
full <- fit(ev, feat)
invisible(xgb.save(full, f_md))
.Random.seed <- rng0
cat(sprintf("xgb_nomp_score: %d features, no MultiPep (%d mp_* dropped)\n",
            length(feat_nomp), length(feat) - length(feat_nomp)))
res_nomp <- cv_score(feat_nomp, "xgb_nomp_score")
full_nomp <- fit(ev, feat_nomp)
invisible(xgb.save(full_nomp, f_md_nomp))
## single features, the same folds, for scale
cat("single features, per-fold mean ROC AUC / AP:\n")
for (sf in intersect(c("mp_signalling_max", "mp_peptidehormone", "n_ho_max", "c_ho_max", "mmseqs_ns_sim",
                       "mmseqs_ns_bits_per_res"), feat)) {
  a <- sapply(seq_len(NFOLD), function(k) { ev <- which(cd$xgb_fold == k & !is.na(y))
    v <- replace(cd[[sf]][ev], is.na(cd[[sf]][ev]), 0)
    c(auc(v, y[ev]), ap(v, y[ev])) })
  cat(sprintf("   %-24s %.3f / %.3f\n", sf, mean(a[1, ], na.rm = TRUE), mean(a[2, ], na.rm = TRUE)))
}
cat(sprintf("rank agreement, xgb_score vs xgb_nomp_score: Spearman %.3f over all candidates\n",
            cor(cd$xgb_score, cd$xgb_nomp_score, method = "spearman")))
## which features carry the full-data models
imp <- xgb.importance(feature_names = feat, model = full)
cat("strongest features (gain, full-data model):\n"); print(head(imp[, .(Feature, Gain = round(Gain, 4))], 15), row.names = FALSE)
imp <- xgb.importance(feature_names = feat_nomp, model = full_nomp)
cat("strongest features without MultiPep (gain, full-data model):\n"); print(head(imp[, .(Feature, Gain = round(Gain, 4))], 15), row.names = FALSE)

## ---- write back, and the tidy scored table -----------------------------------------
cd[, `:=`(xgb_score = round(xgb_score, 5), xgb_nomp_score = round(xgb_nomp_score, 5))]
keep_full <- setdiff(names(cd), c("row", grep("^ns_type_|^[nc]_motif_(?!end$|centre$)", names(cd), value = TRUE, perl = TRUE),
                                   "mmseqs_ns_log10_evalue", "mp_best_class", "mp_best_score"))
fwrite(cd[, keep_full, with = FALSE], f)
arm_cols <- function(tm) c(sprintf("%s_score_%s", tm, ARMS), sprintf("%s_ho_%s", tm, ARMS))
tidy <- c("peptide_id", "gene", "accession", "start", "end", "length", "sequence", "group_id", "group_rank",
          "xgb_score", "xgb_nomp_score", "xgb_fold", "label", "label_match", "label_type",
          ## curated knowns 8-60 aa on the same precursor
          "curated_relation", "curated_match", "curated_match_id", "curated_match_len",
          "curated_shared_res", "curated_frac_of_known", "curated_frac_of_cand",
          "curated_start_diff", "curated_end_diff",
          ## UniProt peptides 8-60 aa on the same precursor
          "uniprot_relation", "uniprot_match", "uniprot_match_id", "uniprot_match_len",
          "uniprot_shared_res", "uniprot_frac_of_known", "uniprot_frac_of_cand",
          "uniprot_start_diff", "uniprot_end_diff",
          ## end
          "n_score_end", "n_score_centre", "n_dist", "n_motif_end", "n_es_arms", "n_es_arm_list", "n_centre_arm",
          arm_cols("n"), "n_ho_max", "n_ho_arms", "n_ho_peak",
          "c_score_end", "c_score_centre", "c_dist", "c_motif_end", "c_es_arms", "c_es_arm_list", "c_centre_arm",
          arm_cols("c"), "c_ho_max", "c_ho_arms", "c_ho_peak",
          "pair_score", "end_score_max", "total_offset",
          ## MultiPep
          "mp_best_class", "mp_best_score", "mp_signalling_max", mp,
          ## MMseqs: closest known peptide of any precursor, then of another precursor
          "mmseqs_best_hit", "mmseqs_best_hit_name", "mmseqs_hit_type", "mmseqs_pident", "mmseqs_tcov",
          "mmseqs_evalue", "mmseqs_bits_per_res", "mmseqs_sim",
          "mmseqs_ns_hit", "mmseqs_ns_hit_name", "mmseqs_ns_hit_type", "mmseqs_ns_pident", "mmseqs_ns_tcov",
          "mmseqs_ns_evalue", "mmseqs_ns_bits_per_res", "mmseqs_ns_sim")
sc <- cd[order(-xgb_score), intersect(tidy, names(cd)), with = FALSE]
fwrite(sc, f_sc)
cat(sprintf("wrote xgb_score and xgb_nomp_score into %s\nwrote %s (%d columns, best first)\nwrote %s and %s\n",
            basename(f), basename(f_sc), ncol(sc), basename(f_md), basename(f_md_nomp)))
cat("\ntop 10 candidates with no known peptide at their range (label 0):\n")
print(sc[label %in% 0L][1:10, .(peptide_id, length, xgb_score, xgb_nomp_score, mp_best_class, mp_best_score,
                               mmseqs_ns_hit_name, mmseqs_ns_sim, n_ho_max, c_ho_max)], row.names = FALSE)
