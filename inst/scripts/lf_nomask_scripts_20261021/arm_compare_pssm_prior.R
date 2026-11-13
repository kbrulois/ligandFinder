## A position-specific scoring matrix (PSSM) of peptide-end windows as a POST-HOC
## PRIOR on an arm's full-proteome end scores (Kevin, 2026-11-13, on a fixed motif
## enrichment matrix: "1 sounds interesting. put a pin in that for now" -- option 1
## being a prior added after the fact, not an input channel):
##
##   combined = logit(CNN score) + lambda * window score
##
## No retraining, no TensorFlow: the arm's tracks (lf_pepend_cv5_tracks_<T>_flatten_
## <ARM>.rds), re-ranked. The question: the msa arm barely uses the cleavage motif
## (its motif-matched negatives neutralise it by design); does a fixed matrix still
## add the context residues around a real end?
##
## SITES, leak-free. UniProt peptide / propeptide ends (pepend_locus/uniprot_peptides
## .csv) plus the ligand list's, minus every cross-validated known (lf_pepend_cv5_
## knowns_<T>.csv), every never-seen end (cv5_unseen_ends*.csv) and every site in a
## gene or precursor of one of those; IDENT also drops a site whose window shares
## that share of its 36 residues with an evaluated end's window, so a paralog under
## another name cannot carry a known end in. Two site sets:
##   all  peptide + propeptide ends, as the idea was pinned
##   pep  peptide ends only. A propeptide's C end KEEPS its cleavage site (the
##        furin site ends the propeptide: R at 0, K at -1, R at -3), where a mature
##        peptide's C end -- every known -- stops before it (G at +1, K/R at
##        +1..+3); propeptides are most of the C sites
## MATRIX. Offsets -7..+28 (N) / -27..+8 (C) around the end, the windows' own; a
## residue outside the mature precursor [n_prot, c_prot] is padding and scores 0,
## as lf_pepend_labels() pads. Log2 odds the way arm_compare_motif_enrichment.R
## computes them, log2((observed + 2) / (expected + 2)), expected = the residues
## present at the offset x
##   sec  the mature secretome's composition: one matrix
##   neg  the composition AT THAT OFFSET of the arm's train-split NEGATIVE windows
##        (lf_pepend_window_groups_<T>_cv5_f<k><sfx>.npz; motif-matched by design),
##        one matrix per fold, never from validation or held-out windows. The
##        dibasic / G columns then cancel; what is left is context
## Window score = sum over the offsets of the log-odds of the residue present.
##
## LAMBDA. A fold's validation ends ARE its held-out ends: fold k's val split holds
## exactly the knowns with fold == k (checked below). So fold k's lambda is fit on
## the OTHER four folds' validation ends, each scored by the fold that held it out,
## against those folds' validation negatives -- leave one fold out -- and the
## never-seen ends (all 50 members) take the lambda fit on all five. Never on the
## ends being scored. Objective: mean log10 rank of the ends among the negatives,
## over a grid; ties go to the smaller lambda. Sensitivity (variants "..., pool
## lambda"): the same fit against each fold's held-out proteome pool.
##
## EVALUATION, as arm_compare_scan_tables.R ranks ends (CNN head):
##   rank_std  among every non-end anchor. UniProt ends count as negatives here,
##             and they are the PSSM's own sites: this ranking penalises the prior
##             for scoring real ends high (scan-pr-hold-out-other-ends)
##   rank_ho   among held-out negatives only: anchors no train window of any
##             compared arm reaches (within SPAN residues), UniProt peptide /
##             propeptide ends dropped too
## Cross-validated ends by the fold that held them out; never-seen ends (the
## intersection of cv5_unseen_ends<sfx>.csv over UNSEEN_SETS) by all 50. Reported:
## n50 / n80 (non-end anchors read for 50 / 80% recall), median and geometric-mean
## rank, recall in the top 100 / 1,000, and the motif make-up of the top TOPK
## non-end anchors (fp_motif; held_out = each fold's top TOPK/5 of the anchors it
## never trained near, all50 = the top TOPK as deployed).
##
## FOUND (2026-11-13, ARM = msa; rank_ho medians; every number in pssm_prior_summary.csv):
##  - Sites are UniProt's alone: every ligand-list gene is a known or never-seen gene. N 490
##    ends (109 peptide, 381 propeptide), C 475 (94 / 381); the closest shares 0.39 of an
##    evaluated end's window.
##  - Peptide + propeptide sites (the design as pinned): lambda 0 on every C fold, either
##    background -- the propeptide C ends put their furin site at 0 / -1 / -3. On N a small
##    median gain (527 -> 493) for fewer ends in the top 100 (28 -> 24 of 91).
##  - Negatives background (context only, motif cancelled): lambda <= 0.02 on C, <= 0.10 on N,
##    no never-seen median moves by more than 5%. The context of these sites adds nothing msa lacks.
##  - Peptide sites vs the secretome is the version that helps, and it helps as a MOTIF prior:
##    geometric-mean rank x0.76-0.93 on every set; N cross-validated 527 -> 356, C never-seen
##    peptides 20,856 -> 11,457. The gain is largest at dibasic / G + dibasic ends (x0.57-0.66);
##    N monobasic and terminus ends lose (x1.02, x1.12). Costs: C cross-validated median 341 ->
##    416, never-seen C peptides in the top 100 6 -> 4, and basic-site anchors in the top 1,000
##    held-out non-ends 17 -> 34% (N), 25 -> 34% (C): the hard negatives msa's motif-matched
##    training learned to discount, the trade cleavage-motif-prior-c-term found for dibasic.
##  - motifind (learned G / basic / P channels, original windows) beats curated_AAswap's median on
##    every set and msa's on never-seen large products (C 39,268 vs 83,880; N 11,510 vs 37,078); msa +
##    the prior matches it on never-seen C peptides (11,457 vs 12,559), not on large products.
##
##   LF_SRC=<checkout> LC_ALL=en_US.UTF-8 RETICULATE_PYTHON=~/.virtualenvs/lf-plm/bin/python \
##     Rscript arm_compare_pssm_prior.R
## env: ARM  msa (the arm the prior goes on)     REFS  curated_AAswap,motifind (compared as they are)
##      UNSEEN_SETS  _msa,_padg     UNSEEN_DIR  ~/AF2_analysis/lf_nomask_scripts_20261021
##      OUT  ~/AF2_analysis/lf_viz_20261113      SEQS  OUT/msa_arms_data/prec_seqs.rds
##      IDENT 0.5   SPAN 35   TOPK 1000
## Writes OUT/pssm_prior_summary.csv (one row per terminus, set, ranking and variant),
## OUT/pssm_prior_fp_motif.csv, OUT/pssm_prior_{matrix,lambda,recall,paired,fp_motif,
## table}.png, and the working tables in OUT/pssm_prior_data/ (sites, matrices, lambda
## curves, end ranks -- whose score is the ranking score, logit-scale for ARM and its priors,
## a probability for REFS -- top non-end anchors). ~5 min; ~4 GB at a time, one terminus.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(reticulate)})
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
source(lf_src("R/pepend_windows.R"))
A2 <- path.expand("~/AF2_analysis")
ARM <- Sys.getenv("ARM", "msa")
REFS <- strsplit(Sys.getenv("REFS", "curated_AAswap,motifind"), ",")[[1]]
US <- strsplit(Sys.getenv("UNSEEN_SETS", "_msa,_padg"), ",")[[1]]
UD <- path.expand(Sys.getenv("UNSEEN_DIR", file.path(A2, "lf_nomask_scripts_20261021")))
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113")))
SEQS <- path.expand(Sys.getenv("SEQS", file.path(OUT, "msa_arms_data", "prec_seqs.rds")))
IDENT <- as.numeric(Sys.getenv("IDENT", "0.5")); SPAN <- as.integer(Sys.getenv("SPAN", "35"))
TOPK <- as.integer(Sys.getenv("TOPK", "1000"))
DATA <- file.path(OUT, "pssm_prior_data")
dir.create(DATA, showWarnings = FALSE, recursive = TRUE)
## each arm's training windows (its window-groups suffix = its INPUT_SFX), as arm_compare_scan_tables.R
GSFX <- c(curated = "", flatten = "", curated_AAswap = "", curated_wAAid = "", motifind = "", esmonly36 = "",
          curated_AAswap_padg = "_padg", curated_AAswap_msacons_padg = "_padg", msa = "_msa", msa_gap = "_msagap")
if (length(miss <- setdiff(c(ARM, REFS), names(GSFX)))) stop("no window suffix known for ", paste(miss, collapse = ", "), call. = FALSE)
np <- import("numpy", convert = TRUE)
LEV <- c("K", "R", "H", "D", "E", "N", "Q", "S", "T", "G", "P", "C", "A", "V", "L", "I", "M", "F", "W", "Y")
OFF <- list(N = -7:28, C = -27:8)
SETS <- c(all = "UniProt peptide + propeptide ends", pep = "UniProt peptide ends only")
BGS <- c(sec = "vs the mature secretome", neg = "vs the arm's train negatives")
LAM <- list(val = seq(0, 1, by = 0.01), pool = seq(0, 1, by = 0.025))   # every optimum on 2026-11-13 was <= 0.22
rank_in <- function(sp, neg) findInterval(-sp, -sort(neg, decreasing = TRUE)) + 1L
lgt <- function(p) qlogis(pmin(pmax(p, 1e-7), 1 - 1e-7))
key_of <- function(acc, anc) paste(acc, anc)
## a variant's display name: msa+pssm_pep_sec_pool -> "msa + PSSM, peptide sites, secretome bg (pool lambda)"
vlab <- function(v) {
  m <- regmatches(v, regexec("^(.*)\\+pssm_(all|pep)_(sec|neg)(_pool)?$", v))[[1]]
  if (!length(m)) return(v)
  sprintf("%s + PSSM, %s, %s%s", m[2], c(all = "peptide + propeptide sites", pep = "peptide sites")[[m[3]]],
          c(sec = "secretome bg", neg = "negatives bg")[[m[4]]], if (nzchar(m[5])) " (pool lambda)" else "")
}
gname <- function(sfx) if (nzchar(sfx)) sfx else "orig"      # "" (the original windows) cannot name a list element
cat_of <- function(motif, near) factor(fifelse(rep_len(near, length(motif)), "next to an annotated end (+/-3)", fifelse(motif == "terminus", "precursor end (padded)", motif)),
                                       c("next to an annotated end (+/-3)", "precursor end (padded)", "dibasic", "G + dibasic", "monobasic", "other"))

## ---- sequences as letter codes (1..20 = LEV, 0 = anything else) ---------------------
if (!file.exists(SEQS)) {
  dir.create(dirname(SEQS), showWarnings = FALSE, recursive = TRUE)
  fc <- readRDS(file.path(A2, "lf_pepend_residue_cache_v2.rds"))
  saveRDS(setNames(fc$prec$seq, fc$prec$accession), SEQS); rm(fc); invisible(gc())
}
seqs <- readRDS(SEQS)
prec <- fread(file.path(A2, "pepend_locus", "prec_order.csv"))[accession %in% names(seqs)]
stopifnot(all(prec$c_prot == nchar(seqs[prec$accession])))
cd <- lapply(seqs[prec$accession], function(s) match(strsplit(s, "", fixed = TRUE)[[1]], LEV, nomatch = 0L))
CS <- unlist(cd, use.names = FALSE); ST <- c(0L, cumsum(lengths(cd)))[seq_along(cd)]
## the background: the mature secretome, every precursor from its first mature residue on
bg <- tabulate(unlist(lapply(seq_along(cd), function(r)
  if (prec$n_prot[r] <= prec$c_prot[r]) cd[[r]][prec$n_prot[r]:prec$c_prot[r]])), 20L)
bg <- setNames(bg / sum(bg), LEV)
Q_SEC <- matrix(bg, 36L, 20L, byrow = TRUE)

## the 36 letter codes of each window; padding (outside [n_prot, c_prot]) = 0
win_codes <- function(acc, anchor, tm) {
  i <- match(acc, prec$accession); stopifnot(!anyNA(i))
  M <- matrix(0L, length(i), 36L)
  for (j in 1:36) {
    p <- anchor + OFF[[tm]][j]; ok <- which(p >= prec$n_prot[i] & p <= prec$c_prot[i])
    M[ok, j] <- CS[ST[i[ok]] + p[ok]]
  }
  M
}
## the window score of every anchor under each matrix of LOs, without building the windows
win_scores <- function(acc, anchor, tm, LOs) {
  i <- match(acc, prec$accession); n0 <- prec$n_prot[i]; c0 <- prec$c_prot[i]; st <- ST[i]
  out <- matrix(0, length(i), length(LOs), dimnames = list(NULL, names(LOs)))
  for (j in 1:36) {
    p <- anchor + OFF[[tm]][j]; ok <- which(p >= n0 & p <= c0)
    code <- integer(length(i)); code[ok] <- CS[st[ok] + p[ok]]
    for (m in seq_along(LOs)) out[, m] <- out[, m] + c(0, LOs[[m]][j, ])[code + 1L]
  }
  out
}
counts <- function(W) vapply(1:20, function(a) colSums(W == a), numeric(36))
## log2((obs + 2) / (expected + 2)); expected = residues present at the offset x q (36 x 20)
pssm <- function(W, q, tm) {
  obs <- counts(W)
  L <- log2((obs + 2) / (rowSums(obs) * q + 2))
  dimnames(L) <- list(OFF[[tm]], LEV); L
}
## one fold's train and val windows: accession, anchor, known
groups <- function(tm, k, sfx) {
  g <- np$load(sprintf("%s/lf_pepend_window_groups_%s_cv5_f%d%s.npz", A2, tm, k, sfx), allow_pickle = TRUE)
  rbindlist(lapply(c("train", "val"), function(sp) data.table(
    split = sp, accession = as.character(g[[sprintf("%s__%s__accession", tm, sp)]]),
    anchor = as.integer(g[[sprintf("%s__%s__anchor", tm, sp)]]), known = as.integer(g[[sprintf("%s__%s__known", tm, sp)]]))))
}

## ---- the ends: knowns, never-seen, and every annotated end ------------------------------
up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
ll <- as.data.table(readRDS(lf_src("inst/extdata/ligand_list.rds")))[!is.na(start) & !is.na(end)]
kn_all <- rbindlist(lapply(c("N", "C"), function(tm)
  fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)), select = c("accession", "anchor", "gene.x", "fold"))[, term := tm]))
setnames(kn_all, "gene.x", "gene")
uf <- list.files(UD, "^cv5_unseen_ends.*\\.csv$", full.names = TRUE)
un_all <- unique(rbindlist(lapply(uf, function(f) fread(f, select = c("gene", "accession", "term", "anchor")))))
## one row per end: the files list an end once per isoform name (CXCL12's N end four times), and the
## join below would multiply those (16 rows); arm_compare_scan_tables.R's never-seen n still counts them
un <- Reduce(function(a, b) a[b[, .(term, accession, anchor)], on = .(term, accession, anchor), nomatch = 0L],
             lapply(US, function(s) unique(fread(file.path(UD, paste0("cv5_unseen_ends", s, ".csv"))), by = c("term", "accession", "anchor"))))
un[, set := fifelse(grepl("^never seen", set), "never seen: peptides", "never seen: large products")]
message("never-seen ends = intersection of ", paste0("cv5_unseen_ends", US, ".csv", collapse = ", "), ": ",
        paste(un[, .N, by = .(term, set)][, sprintf("%s %s %d", term, set, N)], collapse = "; "),
        "\nexcluded from the sites: every end of ", paste(basename(uf), collapse = ", "))

## ---- 1. the sites ----------------------------------------------------------------------
s0 <- rbind(up[, .(term = "N", accession, anchor = start, source = paste("UniProt", type))],
            up[, .(term = "C", accession, anchor = end, source = paste("UniProt", type))],
            ll[, .(term = "N", accession, anchor = start, source = "ligand list")],
            ll[, .(term = "C", accession, anchor = end, source = "ligand list")])
ex_acc <- unique(c(kn_all$accession, un_all$accession))
ex_gene <- unique(na.omit(c(kn_all$gene, un_all$gene, prec$gene[match(ex_acc, prec$accession)])))
ev_key <- c(kn_all[, paste(term, accession, anchor)], un_all[, paste(term, accession, anchor)])
STEPS <- list()
step <- function(d, what) STEPS[[length(STEPS) + 1L]] <<- d[, .(step = what, n = .N), by = .(term, source)]
step(s0, "1 annotated end")
s0[, `:=`(gene = prec$gene[match(accession, prec$accession)], n_prot = prec$n_prot[match(accession, prec$accession)],
          c_prot = prec$c_prot[match(accession, prec$accession)])]
s0 <- s0[!is.na(gene)];                                    step(s0, "2 in the 5,616 precursors")
s0 <- s0[anchor >= n_prot & anchor <= c_prot];             step(s0, "3 inside the mature precursor")
s0 <- s0[!paste(term, accession, anchor) %in% ev_key];     step(s0, "4 not a known or never-seen end")
s0 <- s0[!(accession %in% ex_acc | gene %in% ex_gene)];    step(s0, "5 not in their genes / precursors")
site <- unique(s0[, .(term, accession, gene, anchor, source)])
site <- site[, .(source = paste(sort(unique(source)), collapse = " + ")), by = .(term, accession, gene, anchor)]
for (tm in c("N", "C")) {
  w <- site$term == tm
  Ws <- win_codes(site$accession[w], site$anchor[w], tm)
  ev <- rbind(kn_all[term == tm, .(accession, anchor)], un_all[term == tm & accession %in% prec$accession, .(accession, anchor)])
  We <- win_codes(ev$accession, ev$anchor, tm)
  best <- numeric(nrow(Ws))
  for (e in seq_len(nrow(We))) best <- pmax(best, rowSums(Ws == rep(We[e, ], each = nrow(Ws)) & Ws > 0L) / 36)
  site[w, `:=`(ident = best, window = apply(Ws, 1, function(r) paste(c("-", LEV)[r + 1L], collapse = "")))]
}
site <- site[ident < IDENT]
step(site[, .(term, source)], sprintf("6 one row per end, window identity to every evaluated end < %.2f", IDENT))
## one site per distinct window (overlapping annotations, duplicated genes)
site <- site[, .(accession = accession[1], gene = gene[1], anchor = anchor[1], ident = max(ident), copies = .N,
                 source = paste(sort(unique(unlist(strsplit(source, " + ", fixed = TRUE)))), collapse = " + ")),
             by = .(term, window)]
step(site[, .(term, source)], "7 distinct windows")
site[, `:=`(in_all = TRUE, in_pep = grepl("UniProt peptide|ligand list", source))]
STEPS <- rbindlist(STEPS)
fwrite(site, file.path(DATA, "sites.csv"))
fwrite(STEPS, file.path(DATA, "site_steps.csv"))
options(width = 200)
print(dcast(STEPS, term + source ~ step, value.var = "n", fill = 0L))
message("sites per terminus: ", paste(site[, .(all = .N, pep = sum(in_pep)), by = term][
  , sprintf("%s all %d (peptide-only %d)", term, all, pep)], collapse = "; "),
  "; closest window to an evaluated end shares ", sprintf("%.2f", max(site$ident)), " of its residues")

## ---- 2. the matrices: per terminus, site set, background (and fold, for neg) -----------
G <- list()        # G[[tm]][[sfx]][[k]]: the fold's train / val windows
LO <- list()
for (tm in c("N", "C")) {
  for (sfx in unique(GSFX[c(ARM, REFS)])) G[[tm]][[gname(sfx)]] <- lapply(1:5, function(k) groups(tm, k, sfx))
  ## a fold's validation ends must be its held-out ends, or the leave-one-fold-out logic is wrong
  for (k in 1:5) {
    v <- G[[tm]][[gname(GSFX[[ARM]])]][[k]][split == "val" & known == 1L]
    if (!setequal(key_of(v$accession, v$anchor), kn_all[term == tm & fold == k, key_of(accession, anchor)]))
      stop(tm, " fold ", k, ": the val split's ends are not the fold's held-out knowns", call. = FALSE)
  }
  for (st in names(SETS)) {
    d <- site[term == tm & get(paste0("in_", st))]
    W <- win_codes(d$accession, d$anchor, tm)
    LO[[tm]][[paste0(st, "_sec")]] <- pssm(W, Q_SEC, tm)
    for (k in 1:5) {
      ng <- G[[tm]][[gname(GSFX[[ARM]])]][[k]][split == "train" & known == 0L & accession %in% prec$accession]
      o <- counts(win_codes(ng$accession, ng$anchor, tm))
      LO[[tm]][[sprintf("%s_neg_f%d", st, k)]] <- pssm(W, (o + 2 * Q_SEC) / (rowSums(o) + 2), tm)
    }
  }
}
MX <- rbindlist(lapply(names(LO), function(tm) rbindlist(lapply(names(LO[[tm]]), function(m) {
  d <- as.data.table(as.table(LO[[tm]][[m]])); setnames(d, c("offset", "letter", "log2_odds"))
  d[, `:=`(term = tm, matrix = m, offset = as.integer(as.character(offset)))]
}))))
fwrite(MX, file.path(DATA, "matrices.csv"))

## ---- 3. per terminus: scores, lambda, ranks --------------------------------------------
tracks_path <- function(tm, a) file.path(A2, if (a %in% c("curated", "flatten")) sprintf("lf_pepend_cv5_tracks_%s_flatten.rds", tm)
                                         else sprintf("lf_pepend_cv5_tracks_%s_flatten_%s.rds", tm, a))
read_tracks <- function(tm, a) as.data.table(readRDS(tracks_path(tm, a)))[
  , c("accession", "anchor", sprintf("cnn__f%d", 1:5), "cnn__all50"), with = FALSE]
PRIORS <- CJ(set = names(SETS), bg = names(BGS), sorted = FALSE)[, id := paste0(ARM, "+pssm_", set, "_", bg)]
ER <- list(); FP <- list(); LC <- list(); LP <- list(); PS <- list(); PM <- list()
np_ <- setNames(prec$n_prot, prec$accession)
for (tm in c("C", "N")) {
  t0 <- Sys.time()
  x <- read_tracks(tm, ARM)
  b <- x[, .(accession, anchor)]
  b[, key := key_of(accession, anchor)]
  kn <- kn_all[term == tm]; u <- un[term == tm]
  ends_std <- unique(c(key_of(ll$accession, if (tm == "N") ll$start else ll$end), key_of(kn$accession, kn$anchor),
                       un_all[term == tm, key_of(accession, anchor)]))
  ends_ho <- unique(c(ends_std, key_of(up$accession, if (tm == "N") up$start else up$end)))
  b[, `:=`(end_std = key %in% ends_std, end_ho = key %in% ends_ho)]
  b[accession %in% names(np_), motif := lf_pepend_motif_vec(seqs[[.BY$accession]], tm, anchor, np_[[.BY$accession]]), by = accession]
  near <- unique(rbindlist(lapply(c(-3:-1, 1:3), function(d)
    data.table(accession = c(ll$accession, kn$accession, up$accession, un_all[term == tm]$accession),
               anchor = c(if (tm == "N") ll$start else ll$end, kn$anchor, if (tm == "N") up$start else up$end,
                          un_all[term == tm]$anchor) + d))))
  b[, near_end := key %in% near[, key_of(accession, anchor)] & !end_ho]
  ## which folds trained near each anchor, in any compared arm (arm_compare_scan_tables.R's pools)
  seen <- matrix(FALSE, nrow(b), 5)
  for (sfx in names(G[[tm]])) for (k in 1:5) {
    iv <- unique(G[[tm]][[sfx]][[k]][split == "train", .(accession, lo = anchor - SPAN, hi = anchor + SPAN)])
    seen[unique(b[iv, on = .(accession, anchor >= lo, anchor <= hi), which = TRUE, nomatch = 0L]), k] <- TRUE
  }
  pool_k <- lapply(1:5, function(k) !b$end_ho & !seen[, k])
  pool_all <- !b$end_ho & rowSums(seen) == 0
  PS[[tm]] <- data.table(term = tm, scanned = nrow(b), non_end = sum(!b$end_std), pool_fold_mean = round(mean(vapply(pool_k, sum, 0))),
                         pool_all = sum(pool_all), sites_in_scan = sum(b$key %in% site[term == tm, key_of(accession, anchor)]))
  ## the reference bars of the fp figure: the make-up of the pools themselves
  PM[[tm]] <- rbind(rbindlist(lapply(1:5, function(k) b[pool_k[[k]], .N, by = .(cat = cat_of(motif, near_end))][, share := N / sum(N)]))[
                      , .(share = mean(share)), by = cat][, scoring := "held_out"],
                    b[!b$end_ho, .N, by = .(cat = cat_of(motif, near_end))][, .(cat, share = N / sum(N), scoring = "all50")])[, term := tm]
  ik <- match(key_of(kn$accession, kn$anchor), b$key); iu <- match(key_of(u$accession, u$anchor), b$key)
  if (anyNA(c(ik, iu))) message(tm, ": ", sum(is.na(ik)), " known and ", sum(is.na(iu)), " never-seen ends not in ", ARM, "'s scan")

  ## the window score of every anchor under every matrix, and the arm's logits
  WS <- win_scores(b$accession, b$anchor, tm, LO[[tm]])
  Lk <- lapply(1:5, function(k) lgt(x[[sprintf("cnn__f%d", k)]])); La <- lgt(x$cnn__all50)
  Sk <- function(pr, k) WS[, if (pr$bg == "sec") paste0(pr$set, "_sec") else sprintf("%s_neg_f%d", pr$set, k)]
  Sa <- function(pr) if (pr$bg == "sec") WS[, paste0(pr$set, "_sec")] else rowMeans(WS[, sprintf("%s_neg_f%d", pr$set, 1:5)])

  ## ---- lambda: each fold's validation ends against its negatives (or its held-out pool) --
  vrow <- lapply(1:5, function(k) {
    v <- G[[tm]][[gname(GSFX[[ARM]])]][[k]][split == "val"]
    r <- match(key_of(v$accession, v$anchor), b$key)
    ## an end the scan never scored (SAA1's N end at 18, one before its mature start) has no rank in
    ## any arm; it leaves the fit as it leaves the ranking
    list(pos = r[v$known == 1L & !is.na(r)], neg = r[v$known == 0L & !is.na(r)],
         lost = sprintf("%d+%d", sum(v$known == 1L & is.na(r)), sum(v$known == 0L & is.na(r))))
  })
  message(tm, ": validation ends + negatives not in the scan (left out of the fit), folds 1-5: ",
          paste(vapply(vrow, `[[`, "", "lost"), collapse = ", "))
  for (pi in seq_len(nrow(PRIORS))) {
    pr <- PRIORS[pi]
    for (k in 1:5) {
      s <- Sk(pr, k); L <- Lk[[k]]; vp <- vrow[[k]]$pos; vn <- vrow[[k]]$neg; vn <- vn[!is.na(L[vn])]
      LC[[length(LC) + 1L]] <- data.table(term = tm, prior = pr$id, fit = "val", fold = k, lambda = LAM$val, n = length(vp),
        sumlog = vapply(LAM$val, function(l) sum(log10(rank_in(L[vp] + l * s[vp], L[vn] + l * s[vn]))), 0))
      pl <- which(pool_k[[k]] & !is.na(L)); Lp <- L[pl]; sp <- s[pl]
      LC[[length(LC) + 1L]] <- data.table(term = tm, prior = pr$id, fit = "pool", fold = k, lambda = LAM$pool, n = length(vp),
        sumlog = vapply(LAM$pool, function(l) {
          cn <- Lp + l * sp
          sum(log10(1 + vapply(L[vp] + l * s[vp], function(v) sum(cn >= v), 0)))
        }, 0))
    }
  }
  lc <- rbindlist(LC)[term == tm]
  ## leave one fold out: fold k's lambda from the other four; never-seen (fold 0) from all five
  pick <- function(d, excl) { a <- d[fold != excl, .(obj = sum(sumlog) / sum(n)), by = lambda]; a$lambda[which.min(a$obj)] }
  lp <- rbindlist(lapply(split(lc, by = c("prior", "fit")), function(d)
    data.table(prior = d$prior[1], fit = d$fit[1], fold = 0:5, lambda = vapply(0:5, function(k) pick(d, k), 0))))
  LP[[tm]] <- lp[, term := tm]
  edge <- lp[, .(top = any(lambda == max(LAM[[fit]]))), by = .(prior, fit)][top == TRUE, paste(prior, fit)]
  if (length(edge)) message(tm, ": lambda at the top of the grid for ", paste(edge, collapse = ", "), " -- widen LAM")

  ## ---- ranks, as arm_compare_scan_tables.R --------------------------------------------
  e0 <- rbind(kn[, .(term = tm, set = "cross-validated", accession, anchor, gene, fold, row = ik)],
              u[, .(term = tm, set, accession, anchor, gene, fold = NA_integer_, row = iu)])
  e0[, motif := b$motif[row]]
  run_variant <- function(id, sfun, safun) {
    e <- copy(e0)[, `:=`(variant = id, score = NA_real_, rank_std = NA_integer_, rank_ho = NA_integer_)]
    for (k in 1:5) {
      s <- sfun(k); w <- which(e$set == "cross-validated" & e$fold == k & !is.na(e$row)); sc <- s[e$row[w]]
      e[w, `:=`(score = sc, rank_std = rank_in(sc, s[!b$end_std & !is.na(s)]), rank_ho = rank_in(sc, s[pool_k[[k]] & !is.na(s)]))]
      p <- which(pool_k[[k]] & !is.na(s)); o <- p[order(-s[p])[seq_len(min(TOPK %/% 5L, length(p)))]]
      FP[[length(FP) + 1L]] <<- b[o, .(term = tm, variant = id, scoring = "held_out", fold = k, rank = seq_along(o),
                                        accession, anchor, motif, near_end)]
    }
    s <- safun(); w <- which(e$set != "cross-validated" & !is.na(e$row)); sc <- s[e$row[w]]
    e[w, `:=`(score = sc, rank_std = rank_in(sc, s[!b$end_std & !is.na(s)]), rank_ho = rank_in(sc, s[pool_all & !is.na(s)]))]
    p <- which(!b$end_ho & !is.na(s)); o <- p[order(-s[p])[seq_len(min(TOPK, length(p)))]]
    FP[[length(FP) + 1L]] <<- b[o, .(term = tm, variant = id, scoring = "all50", fold = NA_integer_, rank = seq_along(o),
                                      accession, anchor, motif, near_end)]
    ER[[length(ER) + 1L]] <<- e[, row := NULL]
  }
  run_variant(ARM, function(k) Lk[[k]], function() La)
  for (pi in seq_len(nrow(PRIORS))) for (ft in names(LAM)) {
    pr <- PRIORS[pi]; lam <- lp[prior == pr$id & fit == ft]
    lam_k <- lam$lambda[match(1:5, lam$fold)]; lam_a <- lam[fold == 0]$lambda
    run_variant(if (ft == "val") pr$id else paste0(pr$id, "_pool"),
                function(k) Lk[[k]] + lam_k[k] * Sk(pr, k), function() La + lam_a * Sa(pr))
  }
  rm(x, WS, Lk, La); invisible(gc())
  for (a in REFS) {
    y <- read_tracks(tm, a)[b[, .(accession, anchor)], on = .(accession, anchor)]
    stopifnot(nrow(y) == nrow(b))
    if (any(is.na(y$cnn__all50))) message(tm, " ", a, ": ", sum(is.na(y$cnn__all50)), " of ", ARM, "'s anchors not in its scan")
    run_variant(a, function(k) y[[sprintf("cnn__f%d", k)]], function() y$cnn__all50)
    rm(y); invisible(gc())
  }
  message(sprintf("%s: %.0f s", tm, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rm(b, seen, pool_k, pool_all); invisible(gc())
}

## ---- 4. tables ------------------------------------------------------------------------
ER <- rbindlist(ER); FP <- rbindlist(FP); LC <- rbindlist(LC); LP <- rbindlist(LP)
fwrite(ER, file.path(DATA, "end_ranks.csv")); fwrite(FP, file.path(DATA, "fp_top.csv"))
fwrite(LC, file.path(DATA, "lambda_curves.csv")); fwrite(LP, file.path(DATA, "lambda.csv"))
fwrite(rbindlist(PS), file.path(DATA, "pool_summary.csv"))
VAR <- c(REFS, ARM, PRIORS$id, paste0(PRIORS$id, "_pool"))
lam_txt <- LP[, .(lambda = sprintf("%s | %s", paste(sprintf("%.2f", lambda[fold %in% 1:5][order(fold[fold %in% 1:5])]), collapse = "/"),
                                   sprintf("%.2f", lambda[fold == 0]))),
              by = .(term, variant = fifelse(fit == "val", prior, paste0(prior, "_pool")))]
EL <- melt(ER, id.vars = c("term", "set", "variant", "accession", "anchor", "gene", "fold", "motif"),
           measure.vars = c("rank_std", "rank_ho"), variable.name = "ranking", value.name = "rank")[!is.na(rank)]
SUM <- EL[, .(n = .N, median = as.double(median(rank)), n50 = sort(rank)[ceiling(0.5 * .N)], n80 = sort(rank)[ceiling(0.8 * .N)],
              gmean = round(10^mean(log10(rank)), 1), top100 = round(mean(rank <= 100), 3), top1000 = round(mean(rank <= 1000), 3)),
          by = .(term, set, ranking, variant)]
SUM <- lam_txt[SUM, on = .(term, variant)]
SUM[, label := vapply(variant, vlab, "")]
SUM[, `:=`(variant = factor(variant, VAR), set = factor(set, c("cross-validated", "never seen: peptides", "never seen: large products")),
           term = factor(term, c("N", "C")), ranking = factor(ranking, c("rank_ho", "rank_std")))]
setorder(SUM, term, set, ranking, variant)
setcolorder(SUM, c("term", "set", "ranking", "variant", "label", "lambda"))
fwrite(SUM, file.path(OUT, "pssm_prior_summary.csv"))
FPS <-FP[, .(n = .N), by = .(term, variant, scoring, cat = cat_of(motif, near_end))][, share := round(n / sum(n), 4), by = .(term, variant, scoring)]
fwrite(dcast(FPS, term + scoring + variant ~ cat, value.var = "share", fill = 0), file.path(OUT, "pssm_prior_fp_motif.csv"))
print(SUM[ranking == "rank_ho", .(term, set, variant, lambda, n, median, n80, gmean, top100)], row.names = FALSE)

## ---- 5. figures -----------------------------------------------------------------------
comma <- scales::label_comma()
SHOW <- c(REFS, ARM, PRIORS$id)                    # the pool-lambda variants stay in the tables
## the arms in arm_compare_style.R's colours (motifind grey here: its gold is too close to curated_AAswap's
## orange at line width); the priors blue (secretome bg) / green (negatives bg),
## dark = peptide + propeptide sites, light and dashed = peptide sites only
VCOL <- c(curated = "#A6A6A6", curated_AAswap = "#E08E0B", curated_wAAid = "#5D6D7E", motifind = "#7F8C8D",
          esmonly36 = "#16A085", msa = "#C0392B", msa_gap = "#8E44AD")[c(REFS, ARM)]
VCOL[is.na(VCOL)] <- "#000000"
VCOL <- c(setNames(VCOL, c(REFS, ARM)), setNames(c("#1F618D", "#117A65", "#5DADE2", "#45B39D"), PRIORS$id))[SHOW]
VLT <- setNames(ifelse(grepl("_pep_", SHOW), "22", "solid"), SHOW)
VL <- setNames(vapply(SHOW, vlab, ""), SHOW)
th <- theme_bw(base_size = 10.5) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
tlab <- c(N = "N end (offset 0 = first peptide residue)", C = "C end (offset 0 = last peptide residue)")
nsite <- site[, .(all = .N, pep = sum(in_pep)), by = term]

## 5a. the matrices: secretome background, and the fold-mean of the negatives background
mm <- MX[, .(log2_odds = mean(log2_odds)), by = .(term, offset, letter, set = sub("_.*", "", matrix), bg = sub("^[a-z]+_([a-z]+).*", "\\1", matrix))]
mm[, row := factor(sprintf("%s\n%s", SETS[set], BGS[bg]),
                   c(outer(SETS, BGS, function(s, g) sprintf("%s\n%s", s, g))))]
mm[, `:=`(letter = factor(letter, rev(LEV)), term = factor(tlab[term], tlab[c("N", "C")]))]
elim <- max(abs(mm$log2_odds))
gm <- ggplot(mm, aes(offset, letter, fill = log2_odds)) +
  geom_tile() +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey30", linewidth = 0.25, linetype = 3) +
  scale_fill_gradient2(low = "#5e3c99", mid = "#f7f7f7", high = "#e66101", midpoint = 0, limits = c(-elim, elim), name = "log2 odds") +
  facet_grid(row ~ term, scales = "free_x") +
  labs(x = "offset from the peptide end (residues)", y = NULL,
       title = "The PSSMs: independent peptide-end sites against two backgrounds",
       subtitle = paste0("Sites: UniProt peptide / propeptide ends outside every known and never-seen gene (",
                         paste(nsite[, sprintf("%s: %d, %d peptide-only", term, all, pep)], collapse = "; "),
                         "); the ligand list adds none.\nlog2((obs + 2) / (exp + 2)); negatives background = the ", ARM,
                         " arm's train-split negative windows at each offset (mean of the five folds' matrices). Padding scores 0.")) +
  th + theme(legend.position = "right", panel.grid = element_blank())
ggsave(file.path(OUT, "pssm_prior_matrix.png"), gm, width = 14, height = 12, dpi = 160, bg = "white")

## 5b. lambda: the objective along the grid, per fold, and the lambdas picked
lcd <- LC[, .(obj = sumlog / n), by = .(term, prior, fit, fold, lambda)][, d := obj - obj[lambda == 0], by = .(term, prior, fit, fold)]
lcd[, prior_lab := factor(vapply(prior, vlab, ""), VL[PRIORS$id])]
PANELS <- c(outer(c("on validation negatives", "on the held-out pool"), c("N", "C"), function(f, t) sprintf("%s terminus, lambda fit %s", t, f)))
lcd[, panel := factor(sprintf("%s terminus, lambda fit %s", term, ifelse(fit == "val", "on validation negatives", "on the held-out pool")), PANELS)]
lpd <- LP[fold %in% 1:5][lcd[, .(term, prior, fit, fold, lambda, d)], on = .(term, prior, fit, fold, lambda), nomatch = 0L]
lpd[, `:=`(prior_lab = factor(vapply(prior, vlab, ""), VL[PRIORS$id]),
           panel = factor(sprintf("%s terminus, lambda fit %s", term, ifelse(fit == "val", "on validation negatives", "on the held-out pool")), PANELS))]
gl <- ggplot(lcd, aes(lambda, d, group = fold, colour = factor(fold))) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(data = lpd, shape = 21, fill = "white", size = 2, stroke = 0.9) +
  scale_colour_brewer(palette = "Dark2", name = "fold (its own validation ends' curve; o = the lambda it was scored with, fit on the other four)") +
  facet_grid(prior_lab ~ panel, scales = "free_y", labeller = labeller(prior_lab = label_wrap_gen(28))) +
  labs(x = "lambda", y = "change in mean log10 rank of the fold's ends (below 0 = better)",
       title = "Fitting lambda: does adding the PSSM lift each fold's held-out ends above its negatives?",
       subtitle = paste0("combined = logit(CNN) + lambda x window score. Each line: one fold's validation (= held-out) ends against that fold's negatives. ",
                         "The lambda a fold is SCORED with is fit on the other four folds' curves pooled.")) +
  th
ggsave(file.path(OUT, "pssm_prior_lambda.png"), gl, width = 15, height = 11, dpi = 160, bg = "white")

## 5c. recall against anchors scanned
NG <- unique(round(10^seq(0, 6.5, length.out = 160)))
rc <- EL[variant %in% SHOW, .(N = NG, recall = vapply(NG, function(x) mean(rank <= x), 0)), by = .(term, set, ranking, variant)]
rc[, `:=`(variant = factor(variant, SHOW), set = factor(set, levels(SUM$set)),
          panel = factor(sprintf("%s terminus\n%s", term, c(rank_std = "every non-end anchor (rank_std)", rank_ho = "held-out negatives (rank_ho)")[as.character(ranking)]),
                         c(outer(c("N", "C"), c("held-out negatives (rank_ho)", "every non-end anchor (rank_std)"),
                                 function(t, d) sprintf("%s terminus\n%s", t, d))[c(1, 3, 2, 4)])))]
nn <- EL[ranking == "rank_ho" & variant == ARM, .(n = .N), by = .(term, set)]
gr <- ggplot(rc, aes(N, recall, colour = variant, linetype = variant)) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50") +
  geom_step(linewidth = 0.6) +
  scale_x_log10(labels = comma, breaks = 10^(0:6)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = VCOL, labels = VL, name = NULL) +
  scale_linetype_manual(values = VLT, labels = VL, name = NULL) +
  facet_grid(set ~ panel) +
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
  labs(x = "non-end anchors ranked above the end (log)", y = "recall of the ends",
       title = sprintf("Does the PSSM prior lift %s's ends? Recall against anchors scanned (CNN head)", ARM),
       subtitle = paste0("Cross-validated: each end by the fold that held it out, the prior weighted by a lambda fit on the other folds.\nNever seen (",
                         paste(nn[set != "cross-validated", sprintf("%s %s %d", term, sub("never seen: ", "", set), n)], collapse = ", "),
                         "): all 50 members, lambda fit on every fold. Held-out negatives: anchors no training window of any compared arm reaches;\n",
                         "UniProt ends (the PSSM's own sites) are negatives under rank_std only. Lines on top of one another = lambda 0, no change.")) +
  th
ggsave(file.path(OUT, "pssm_prior_recall.png"), gr, width = 15, height = 10.5, dpi = 160, bg = "white")

## 5d. per end: rank with the prior against without, by the end's motif
pa <- dcast(EL[ranking == "rank_ho" & variant %in% c(ARM, PRIORS$id)], term + set + accession + anchor + motif ~ variant, value.var = "rank")
pa <- melt(pa, id.vars = c("term", "set", "accession", "anchor", "motif", ARM), measure.vars = PRIORS$id, variable.name = "variant", value.name = "rank_prior")
pa[, delta := log10(rank_prior) - log10(get(ARM))]
pa[, `:=`(variant = factor(variant, PRIORS$id), set = factor(set, levels(SUM$set)), term = factor(term, c("N", "C")),
          motif = factor(fifelse(is.na(motif), "not scanned", motif), c("G + dibasic", "dibasic", "monobasic", "terminus", "other", "not scanned")))]
pm <- pa[!is.na(delta), .(med = median(delta), n = .N, up = sum(delta < 0), down = sum(delta > 0)), by = .(term, set, motif, variant)]
gp <- ggplot(pa[!is.na(delta)], aes(motif, delta, colour = variant)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 1), size = 0.9, alpha = 0.55) +
  geom_point(data = pm, aes(y = med), position = position_dodge(width = 0.8), shape = 95, size = 9, alpha = 1) +
  geom_text(data = pm[, .(n = n[1]), by = .(term, set, motif)], aes(x = motif, y = Inf, label = paste0("n = ", n)),
            inherit.aes = FALSE, vjust = 1.4, size = 2.8, colour = "grey30") +
  scale_colour_manual(values = VCOL[PRIORS$id], labels = VL[PRIORS$id], name = NULL) +
  facet_grid(term ~ set, scales = "free_x", space = "free_x") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "what flanks the end (lf_pepend_motif_vec)", y = sprintf("log10 rank_ho with the prior - without (below 0 = %s + prior ranks it higher)", ARM),
       title = sprintf("Which ends the prior moves, by cleavage motif (held-out negatives, %s's CNN head)", ARM),
       subtitle = "One point per end and prior; bar = the median. Cross-validated ends by the fold that held them out, never-seen ends by all 50.") +
  th + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(OUT, "pssm_prior_paired.png"), gp, width = 15, height = 8.5, dpi = 160, bg = "white")
fwrite(pm, file.path(DATA, "paired_by_motif.csv"))

## 5e. what the top non-end anchors are
mcol <- c("next to an annotated end (+/-3)" = "#CC79A7", "precursor end (padded)" = "#0072B2", dibasic = "#D55E00",
          "G + dibasic" = "#009E73", monobasic = "#E69F00", other = "grey80")
## references: the pool's own make-up, and the cross-validated known ends'
REFBAR <- c("known ends (for reference)", "all non-end anchors (pool)")
ref_kn <- unique(ER[set == "cross-validated" & !is.na(motif), .(term, accession, anchor, motif)])[
  , .(n = .N), by = .(term, cat = cat_of(motif, FALSE))][, .(cat, share = n / sum(n)), by = term]
fd <- rbind(FPS[variant %in% SHOW, .(term, scoring, cat, share, bar = VL[as.character(variant)])],
            rbindlist(PM)[, .(term, scoring, cat, share, bar = REFBAR[2])],
            ref_kn[, .(scoring = c("held_out", "all50")), by = .(term, cat, share)][, .(term, scoring, cat, share, bar = REFBAR[1])])
fd[, bar := factor(bar, rev(c(REFBAR, VL)))]
fd[, `:=`(term = factor(term, c("N", "C"), c("N terminus", "C terminus")),
          scoring = factor(scoring, c("held_out", "all50"), c(sprintf("held out: each fold's top %d of the anchors it never trained near", TOPK %/% 5L),
                                                             sprintf("all 50 members: the top %s non-end anchors, as deployed", comma(TOPK)))))]
fd <- fd[order(term, scoring, bar, -as.integer(cat))][, xmid := cumsum(share) - share / 2, by = .(term, scoring, bar)]
gf <- ggplot(fd, aes(share, bar, fill = cat)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  geom_text(data = fd[cat != "other" & share >= 0.05], aes(x = xmid, label = sprintf("%.0f%%", 100 * share)), size = 2.6) +
  scale_fill_manual(values = mcol, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = scales::label_percent(), breaks = c(0, 0.25, 0.5, 0.75), expand = c(0, 0)) +
  facet_grid(term ~ scoring) +
  labs(x = "share of the anchors", y = NULL,
       title = sprintf("The false-positive cost: what the top %s non-end anchors are, with and without the prior (CNN head)", comma(TOPK)),
       subtitle = "Annotated ends excluded (ligand list, knowns, never-seen sets, UniProt peptides / propeptides); 'next to' = within 3 residues of one.") +
  th
ggsave(file.path(OUT, "pssm_prior_fp_motif.png"), gf, width = 15, height = 8.5, dpi = 160, bg = "white")

## 5f. the numbers as a table: median rank_ho, top-100 recall, lambda
tb <- SUM[ranking == "rank_ho" & variant %in% SHOW]
tb[, cell := sprintf("%s  (%.0f%%)", comma(round(median)), 100 * top100)]
tw <- dcast(tb, variant ~ term + set, value.var = "cell")
tw <- tw[match(SHOW, variant)]
lam_w <- dcast(lam_txt[variant %in% SHOW], variant ~ factor(term, c("N", "C")), value.var = "lambda")
tw <- lam_w[tw, on = "variant"]
tw[, variant := VL[as.character(variant)]]
setnames(tw, names(tw), sub("_", "\n", sub("never seen: ", "never seen\n", names(tw))))
setnames(tw, c("N", "C"), c("lambda N\nfolds | never-seen", "lambda C\nfolds | never-seen"), skip_absent = TRUE)
for (j in names(tw)) set(tw, which(is.na(tw[[j]])), j, "")
tg <- gridExtra::tableGrob(tw, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 9,
  core = list(fg_params = list(hjust = 1, x = 0.95)), colhead = list(fg_params = list(fontface = "bold", hjust = 1, x = 0.95))))
ttl <- grid::textGrob(sprintf("PSSM prior on %s: median rank among held-out negatives (rank_ho), and the share of ends in the top 100", ARM),
                      x = 0.01, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 12))
sub <- grid::textGrob(paste0("Lower median = better. Cross-validated ends by the fold that held them out; never-seen ends by all 50 members. ",
                             "lambda: the five folds' (each fit on the other four) | the never-seen ends' (fit on all five).\n",
                             "Full table, rank_std too: pssm_prior_summary.csv; the pool-lambda sensitivity rows are there as well."),
                      x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 8.5, col = "grey25"))
g <- gridExtra::arrangeGrob(ttl, sub, tg, ncol = 1, heights = grid::unit.c(grid::unit(1.6, "lines"), grid::unit(2.6, "lines"), sum(tg$heights)))
png(file.path(OUT, "pssm_prior_table.png"), width = max(grid::convertWidth(sum(tg$widths), "in", valueOnly = TRUE), 10) + 0.4,
    height = 1.5 + grid::convertHeight(sum(tg$heights), "in", valueOnly = TRUE), units = "in", res = 180, bg = "white")
grid::grid.draw(g); invisible(dev.off())
message("wrote ", OUT, "/pssm_prior_*")
