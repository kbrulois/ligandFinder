#!/usr/bin/env Rscript
## ---- peptide-END windows: a separate training set -----------------------------
## Both termini of every docked known peptide (chemokines excluded for now),
## each window anchored on the peptide's own terminal residue -- first residue
## at position 8 (N), last at 28 (C) -- with DB and gap folded away: before
## the peptide is NT_cleavage_context, after it CT_cleavage_context. Geometry,
## labels and candidate anchors live in R/pepend_windows.R.
##
## Differences from the production set (9.2_add_contact_data.R):
##   * both ends of each peptide, not only the one that inserts in the receptor
##   * no insertion-depth filter (9.2 keeps lig1_end_ind < 16) and no dibasic
##     requirement (9.2 keeps win_type == "db")
##   * a C terminus shared by several peptides (CCK-8/-33/-58) is one window,
##     labelled from the SHORTEST of them; all names are kept in `pep_ids`
##   * negatives: candidate anchors (dibasic, amidation, precursor termini;
##     see lf_pepend_candidates) from genes that are not known ligands, plus a
##     share of windows anchored at random residues, so a positive whose
##     terminus is not at any candidate motif cannot be told apart by that alone
##
## Inputs are 9.2's 26 channels with 9.2's global scaling ranges, read from its
## cache, so a residue gets exactly the values the production windows give it
## (checked below against the cached production windows).
##
##   /usr/local/bin/Rscript inst/scripts/10_7_pepend_windows.R
##   ... --n-ctrl 400 --random-frac 0.25 --seed 42
##
## Output: ~/AF2_analysis/lf_pepend_nn_input.rds
##   nn_input     list(N, C) of list(train, val, all), the 9.2 contract
##   all_params3  channel order
##   class_names  LF_PEPEND_CLASSES -- pass to lf_dcnn_config(class_names = )
##   knowns       one row per known window: peptide, anchor, split, motif,
##                insertion, AA and labels as 36-long vectors
## The first run reads secretome and secretome_aa whole (~16 GB of RAM, as 9.2
## does) and caches per-precursor features to lf_pepend_residue_cache.rds;
## later runs take seconds. --rebuild-features to rebuild that cache.
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
out_path    <- path.expand(.opt("--out", "~/AF2_analysis/lf_pepend_nn_input.rds"))
n_ctrl      <- as.integer(.opt("--n-ctrl", "400"))       # negatives per split per terminus
random_frac <- as.numeric(.opt("--random-frac", "0.25")) # of those, anchored at random residues
## NEAR-END negatives: windows anchored a few residues off a REAL peptide end,
## labelled 0. Every other negative comes from a precursor with no known peptide
## at all, so the model has never been shown "almost, but not this residue" --
## and it shows: in the step-1 scan the local argmax lands on the true end only
## 30% of the time, with d=+1 scoring HIGHER on average than d=0.
## Taken OUT of n_ctrl rather than added to it, so the negative count per split
## is unchanged and only its composition differs.
near_n      <- as.integer(.opt("--near-n", "0"))          # 0 = off, the old set
near_lo     <- as.integer(.opt("--near-lo", "2"))         # |offset| range, inclusive
near_hi     <- as.integer(.opt("--near-hi", "8"))
val_frac    <- as.numeric(.opt("--val-frac", "0.34"))
## Swap the two sides after splitting. With --val-frac 0.5 this is the second
## fold of a 2-fold cross-validation: every known that was TRAIN in the first
## build is VAL here and vice versa, so training both gives each known one
## in-sample and one held-out prediction.
swap_splits <- "--swap-splits" %in% .args
sim_h       <- as.numeric(.opt("--sim-h", "0.30"))
seed        <- as.integer(.opt("--seed", "42"))
docked_ref  <- path.expand(.opt("--docked", "~/AF2_analysis/knowns.rds"))
secretome_p <- path.expand(.opt("--secretome", "~/AF2_analysis/secretome_latest.rds"))
aa_p        <- path.expand(.opt("--secretome-aa", "~/peptide_alg/build_residue_db/processed/secretome_aa.rds"))
nn_dat_p    <- path.expand(.opt("--nn-dat-cache", "~/AF2_analysis/nn_dat_cache.rds"))
feat_cache  <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
set.seed(seed)

id_map <- readRDS(file.path(ROOT, "data", "id_mapping.rds"))
e2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])
s2e <- setNames(id_map[["Entry Name"]], id_map[["Gene Names (primary)"]])
is_chemokine <- function(gene) grepl("^(CCL|CXCL|XCL|CX3CL)", gene)

## ---- 1. the docked peptides -------------------------------------------------------
## 9.2's funnel: relevant site, rank-1 model, then per peptide the receptor model
## whose inserting residue sits closest to a terminus (ties by iptm). That model
## supplies the pocket labels.
con <- readRDS(docked_ref) %>%
  filter(location == "relevant") %>%
  filter(rank == 1, .by = afpd_dir_name) %>%
  mutate(gene      = unname(e2s[p2_name]),
         accession = p2_id,
         pep_id    = paste0(gene, "_", p2_range),
         pep_start = as.integer(stringr::str_extract(p2_range, "^\\d+")),
         pep_end   = as.integer(stringr::str_extract(p2_range, "\\d+$")),
         ins_ind   = as.integer(stringr::str_extract(lig1_end, "\\d+")),
         ins_term  = stringr::str_remove(lig1_end, "\\d+")) %>%
  group_by(pep_id) %>% arrange(ins_ind, desc(iptm), .by_group = TRUE) %>% slice(1) %>% ungroup()
message(sprintf("docked peptides: %d (%d chemokines dropped)",
                sum(!is_chemokine(con$gene)), sum(is_chemokine(con$gene))))
con <- con %>% filter(!is_chemokine(gene))

## pocket residues: 9.2's clean_contacts (area > 1, dist < 6), then in_pocket
con$pocket <- lapply(con$contacts, function(ct) {
  if (is.null(ct) || !nrow(ct)) return(integer(0))
  ct <- ct %>% filter(area > 1 & dist < 6) %>%
    group_by(lig1_residue) %>% summarise(in_pocket = first(in_pocket), .groups = "drop")
  as.integer(stringr::str_extract(ct$lig1_residue[ct$in_pocket %in% TRUE], "\\d+$"))
})
con <- con %>% select(gene, accession, pep_id, pep_start, pep_end, ins_ind, ins_term, iptm, pocket)

## names from the ligand list where a peptide matches one exactly
ll <- readRDS(file.path(ROOT, "inst", "extdata", "ligand_list.rds"))
con$pep_name <- ll$final_name[match(paste(con$accession, con$pep_start, con$pep_end),
                                    paste(ll$accession, ll$start, ll$end))]

## ---- 2. precursors and their residue features (cached) --------------------------
## Sequence + signal peptide (9.2's sp_ind rule), and one (nchar, 26) feature
## matrix per precursor in 9.2's channels and scaling. Building it reads
## secretome and secretome_aa whole (~5 GB, a few minutes); later runs load
## the cache. --rebuild-features to redo it.
if (file.exists(feat_cache) && !"--rebuild-features" %in% .args) {
  message("loading ", feat_cache)
  .fc <- readRDS(feat_cache)
  prec <- .fc$prec; feats <- .fc$feats; all_params3 <- .fc$all_params3
  rm(.fc)
} else {
  message("reading ", secretome_p, " ...")
  secretome <- readRDS(secretome_p)
  prec <- tibble(accession = secretome$accession, gene = secretome$gene,
                 seq = secretome$sequence_uni,
                 sp_ind = map_int(secretome$features, function(x) {
                   n <- x %>% filter(type == "signal peptide") %>% pull(end)
                   if (!length(n) || all(is.na(n))) 1L else as.integer(max(n, na.rm = TRUE))
                 })) %>%
    filter(!is.na(seq), nchar(seq) > 0) %>%
    distinct(accession, .keep_all = TRUE) %>%
    ## 9.2: sp_ind == 1 means "no signal peptide"
    mutate(n_prot = ifelse(sp_ind == 1L, 1L, sp_ind + 1L), c_prot = nchar(seq))
  rm(secretome); invisible(gc())

  .nd <- readRDS(nn_dat_p)
  all_params3 <- .nd$all_params3; chan_range <- .nd$chan_range
  old_known   <- .nd$known_dat %>% select(peps, gene, meta_data, data)   # for the check below
  rm(.nd); invisible(gc())

  ## Must match 9.2_add_contact_data.R's aa_props_raw exactly (then min-max scaled)
  aa_props <- tibble::tribble(
    ~AA, ~AA_hydro, ~AA_charge, ~AA_mw, ~AA_pI,
    "A",  1.8,  0.0,  71.08,  6.00, "R", -4.5,  1.0, 156.19, 10.76,
    "N", -3.5,  0.0, 114.10,  5.41, "D", -3.5, -1.0, 115.09,  2.77,
    "C",  2.5,  0.0, 103.14,  5.07, "Q", -3.5,  0.0, 128.13,  5.65,
    "E", -3.5, -1.0, 129.12,  3.22, "G", -0.4,  0.0,  57.05,  5.97,
    "H", -3.2,  0.1, 137.14,  7.59, "I",  4.5,  0.0, 113.16,  6.02,
    "L",  3.8,  0.0, 113.16,  5.98, "K", -3.9,  1.0, 128.17,  9.74,
    "M",  1.9,  0.0, 131.19,  5.74, "F",  2.8,  0.0, 147.18,  5.48,
    "P", -1.6,  0.0,  97.12,  6.30, "S", -0.8,  0.0,  87.08,  5.68,
    "T", -0.7,  0.0, 101.10,  5.60, "W", -0.9,  0.0, 186.21,  5.89,
    "Y", -1.3,  0.0, 163.18,  5.66, "V",  4.2,  0.0,  99.13,  5.96,
    "U",  2.5, -1.0, 150.04,  5.47) %>%
    mutate(across(-AA, ~ (.x - min(.x)) / (max(.x) - min(.x))))
  aa_prop_cols <- c("AA_hydro", "AA_charge", "AA_mw", "AA_pI")
  scale_cols   <- setdiff(all_params3, c(aa_prop_cols, "padding"))
  stopifnot(identical(sort(colnames(chan_range)), sort(scale_cols)))

  message("reading ", aa_p, " (~3 GB) ...")
  sa <- readRDS(aa_p)
  sa <- as.data.frame(sa)[, c("accession", "index", "AA", scale_cols)]
  for (cn in scale_cols)
    sa[[cn]] <- pmin(pmax((sa[[cn]] - chan_range[1, cn]) / (chan_range[2, cn] - chan_range[1, cn]), 0), 1)
  sa_rows <- split(seq_len(nrow(sa)), sa$accession)
  invisible(gc())

  ## a residue scaffold from the sequence, with secretome_aa joined on
  ## (index, AA) as 9.2 joins it: a residue whose letter disagrees, or that the
  ## sparse table lacks, falls through to 0
  residue_matrix <- function(acc, seq) {
    aa <- strsplit(seq, "")[[1]]; n <- length(aa)
    m <- matrix(0, n, length(all_params3), dimnames = list(NULL, all_params3))
    r <- sa_rows[[acc]]
    if (length(r)) {
      s <- sa[r, , drop = FALSE]
      s <- s[which(s$index >= 1 & s$index <= n), , drop = FALSE]
      s <- s[which(s$AA == aa[s$index]), , drop = FALSE]
      s <- s[!duplicated(s$index), , drop = FALSE]
      v <- as.matrix(s[, scale_cols]); v[is.na(v)] <- 0
      m[s$index, scale_cols] <- v
    }
    p <- as.matrix(aa_props[match(aa, aa_props$AA), aa_prop_cols]); p[is.na(p)] <- 0
    m[, aa_prop_cols] <- p
    m
  }
  message("building residue features for ", nrow(prec), " precursors ...")
  feats <- setNames(map2(prec$accession, prec$seq, residue_matrix), prec$accession)
  rm(sa, sa_rows); invisible(gc())

  ## the check: rebuild cached production windows' inputs. Their coordinates
  ## are in the name; interior (unpadded) ones must match what 9.2 stored.
  chk <- old_known %>%
    mutate(wN = as.integer(stringr::str_match(peps, "_w(\\d+)-(\\d+)$")[, 2])) %>%
    filter(map_lgl(meta_data, ~ !anyNA(.x$AA))) %>%
    left_join(prec %>% select(gene, accession), by = "gene", relationship = "many-to-many") %>%
    filter(!is.na(accession)) %>% distinct(peps, .keep_all = TRUE)
  dev <- vapply(seq_len(nrow(chk)), function(i) {
    new <- feats[[chk$accession[i]]][chk$wN[i] + 0:35, , drop = FALSE]
    max(abs(new - as.matrix(chk$data[[i]][, all_params3])))
  }, numeric(1))
  message(sprintf("feature check vs %d cached production windows: max |diff| %.2g", length(dev), max(dev)))
  if (max(dev) > 1e-6) stop("features do not reproduce 9.2's windows; see the check above")
  rm(old_known); invisible(gc())

  saveRDS(list(prec = prec, feats = feats, all_params3 = all_params3), feat_cache)
  message("cached precursor features to ", feat_cache)
}
prec_seq <- setNames(prec$seq, prec$accession)

miss <- setdiff(con$accession, prec$accession)
if (length(miss)) message("dropping ", length(miss), " peptide(s) whose precursor is not in secretome: ",
                          paste(con$gene[con$accession %in% miss], collapse = ", "))
con <- con %>% inner_join(prec %>% select(accession, seq, n_prot, c_prot), by = "accession") %>%
  ## 9.2: a peptide that starts inside the annotated signal peptide moves it
  mutate(n_prot = pmin(n_prot, pep_start))

## ---- 3. one window per peptide end ---------------------------------------------------
ends <- bind_rows(con %>% mutate(term = "N", anchor = pep_start),
                  con %>% mutate(term = "C", anchor = pep_end)) %>%
  mutate(len = pep_end - pep_start + 1L) %>%
  group_by(accession, term, anchor) %>%
  arrange(len, .by_group = TRUE) %>%
  mutate(pep_ids = paste(pep_id, collapse = ";"), n_sharing = n()) %>%
  slice(1) %>% ungroup() %>%
  mutate(w_start   = lf_pepend_start(term, anchor),
         peps      = sprintf("%s_w%d-%d", gene, w_start, w_start + LF_PEPEND$seq_len - 1L),
         target    = term,
         win_type  = "known",
         motif     = lf_pepend_motif(seq, term, anchor, n_prot),
         insertion = case_when(ins_term == term & ins_ind < 3  ~ "end insertion",
                               ins_term == term & ins_ind < 16 ~ "loop insertion",
                               TRUE                            ~ "non-inserting end"))
ends$known_idx <- pmap(list(ends$w_start, ends$pep_start, ends$pep_end, ends$n_prot,
                            ends$c_prot, ends$pocket), lf_pepend_labels)
message(sprintf("known windows: %d N, %d C (from %d peptides; %d ends shared by >1 peptide)",
                sum(ends$term == "N"), sum(ends$term == "C"), nrow(con), sum(ends$n_sharing > 1)))
print(table(term = ends$term, motif = ends$motif))

## ---- 4. candidate anchors, and the negative pool --------------------------------------
gpcr_ligs <- data.table::fread(file.path(ROOT, "inst/extdata/GPCRdb_known_pairings_human_plus2more_unique.csv")) %>%
  as_tibble() %>% mutate(gene = stringr::str_remove(lig, "^h") %>% stringr::str_remove("x\\d+x\\d+"))
uni_pep <- data.table::fread(file.path(ROOT, "inst/extdata/Candidate_peps_KB_from_known_pep_unknownGPCR_run#1.txt")) %>%
  as_tibble() %>% filter(!ALso_run %in% c("no", "not_now", "??", "did we run?"))
banned <- c(gpcr_ligs$gene, uni_pep$Gene)
banned <- unique(na.omit(c(banned, unname(e2s[banned]), unname(s2e[banned]), con$gene)))
message(sprintf("banned (known-ligand) genes: %d, kept out of the negatives", length(banned)))

prec <- prec %>% filter(!is_chemokine(gene))
cand <- prec %>%
  mutate(c = map2(seq, n_prot, lf_pepend_candidates)) %>%
  select(accession, gene, c) %>% unnest(c)
## random anchors: the same number per terminus as there are db candidates,
## drawn uniformly over the mature residues of the non-ligand precursors
pool_prec <- prec %>% filter(!gene %in% banned, c_prot - n_prot >= 10)
rand <- bind_rows(lapply(c("N", "C"), function(t) {
  k <- 4L * n_ctrl
  i <- sample(nrow(pool_prec), k, replace = TRUE, prob = pool_prec$c_prot - pool_prec$n_prot + 1)
  tibble(accession = pool_prec$accession[i], gene = pool_prec$gene[i], term = t,
         anchor = pool_prec$n_prot[i] +
           as.integer(floor(runif(k) * (pool_prec$c_prot[i] - pool_prec$n_prot[i] + 1))),
         win_type = "random")
})) %>% distinct(accession, term, anchor, .keep_all = TRUE)
message("candidate windows:"); print(as.data.frame(count(cand, term, win_type)))

## near-end anchors: every known end displaced by +/- near_lo..near_hi, kept
## inside the mature range. `parent` carries the end it came from so the split
## can follow it -- a window 2 residues from a VAL positive must not be a TRAIN
## negative, which is the whole reason these are informative.
near <- if (near_n > 0) {
  ends %>%
    select(accession, gene, term, parent = anchor, n_prot, c_prot,
           pep_start, pep_end, pocket) %>%
    tidyr::crossing(d = c(-rev(near_lo:near_hi), near_lo:near_hi)) %>%
    mutate(anchor = parent + d, win_type = "near_end") %>%
    filter(anchor >= n_prot, anchor <= c_prot) %>%
    ## an offset that lands on ANOTHER known end of the same terminus is a
    ## positive, not a near miss -- polyproteins (POMC, PENK) make this common
    anti_join(ends %>% select(accession, term, anchor), by = c("accession", "term", "anchor"))
} else ends[0, ] %>% mutate(parent = integer(0), d = integer(0))
if (near_n > 0) {
  near$w_start <- lf_pepend_start(near$term, near$anchor)
  ## real per-index labels, from the real peptide coordinates. These windows DO
  ## cover peptide residues, unlike every other negative, so labelling them
  ## "none" would teach the per-index head to deny a peptide it can see.
  near$known_idx <- pmap(list(near$w_start, near$pep_start, near$pep_end,
                              near$n_prot, near$c_prot, near$pocket), lf_pepend_labels)
  message(sprintf("near-end anchors: %d (|offset| %d-%d, %d dropped as other known ends)",
                  nrow(near), near_lo, near_hi,
                  sum(duplicated(rbind(near[, c("accession","term","anchor")],
                                       ends[, c("accession","term","anchor")]))) ))
  print(as.data.frame(count(near, term)))
}

## how many known ends a candidate anchor would have found
ends$at_candidate <- paste(ends$accession, ends$term, ends$anchor) %in%
  paste(cand$accession, cand$term, cand$anchor)
message(sprintf("known ends at a candidate anchor: N %d/%d, C %d/%d",
                sum(ends$at_candidate[ends$term == "N"]), sum(ends$term == "N"),
                sum(ends$at_candidate[ends$term == "C"]), sum(ends$term == "C")))

## ---- 5. slice every window -------------------------------------------------------------
make_windows <- function(w) {
  ## w: accession, gene, term, anchor, win_type (+ n_prot, c_prot, known_idx for knowns)
  if (!"n_prot" %in% names(w))
    w <- w %>% left_join(prec %>% select(accession, n_prot, c_prot), by = "accession")
  w <- w %>% mutate(w_start = lf_pepend_start(term, anchor),
                    peps = sprintf("%s_w%d-%d", gene, w_start, w_start + LF_PEPEND$seq_len - 1L),
                    target = term)
  by_acc <- split(seq_len(nrow(w)), w$accession)
  data <- vector("list", nrow(w)); meta <- vector("list", nrow(w))
  for (acc in names(by_acc)) {
    f  <- feats[[acc]]
    aa <- strsplit(prec_seq[[acc]], "")[[1]]
    for (i in by_acc[[acc]]) {
      data[[i]] <- lf_pepend_slice(f, w$w_start[i], w$n_prot[i], w$c_prot[i])
      r  <- w$w_start[i] + 0:35
      ok <- r >= w$n_prot[i] & r <= w$c_prot[i]
      meta[[i]] <- tibble(AA = ifelse(ok, aa[pmax(pmin(r, length(aa)), 1)], NA_character_),
                          index = ifelse(ok, r, NA_integer_))
    }
  }
  w$data <- data; w$meta_data <- meta
  if (!"known_idx" %in% names(w)) w$known_idx <- rep(list(rep("none", LF_PEPEND$seq_len)), nrow(w))
  w$meta_data <- map2(w$meta_data, w$known_idx, ~ mutate(.x, known_idx = .y))
  w
}

known_w <- make_windows(ends %>% select(accession, gene, term, anchor, win_type, n_prot, c_prot, known_idx,
                                        pep_id, pep_ids, pep_name, pep_start, pep_end, len,
                                        motif, insertion, at_candidate)) %>%
  mutate(known = 1)
message("slicing ", nrow(cand), " candidate windows ...")
cand_w  <- make_windows(cand) %>% mutate(known = 0)
rand_w  <- make_windows(rand) %>% mutate(known = 0)
near_w  <- if (near_n > 0) make_windows(
  near %>% select(accession, gene, term, anchor, win_type, n_prot, c_prot,
                  known_idx, parent, d)) %>% mutate(known = 0) else NULL

## ---- 6. split and assemble ---------------------------------------------------------------
## 9.2's split_knowns: cluster windows by aligned-sequence identity so no
## cluster straddles train/val, then send every ~1/val_frac-th cluster (in
## dendrogram order) to val
split_knowns <- function(k_sub) {
  n  <- nrow(k_sub)
  aa <- do.call(rbind, lapply(k_sub$meta_data, function(m) { a <- m$AA; a[is.na(a)] <- "X"; a }))
  hc <- hclust(as.dist(ape::dist.aa(aa, scaled = TRUE)), method = "average")
  cl <- cutree(hc, h = sim_h)
  ord <- unique(cl[hc$order])
  val_cl <- ord[seq(2L, length(ord), by = max(2L, round(1 / val_frac)))]
  out <- list(train = which(!cl %in% val_cl), val = which(cl %in% val_cl))
  ## the complement, taken AFTER clustering so both folds respect the same
  ## similarity boundaries -- no cluster straddles the split in either one
  if (swap_splits) out <- list(train = out$val, val = out$train)
  out
}

draw_negatives <- function(pool, n, exclude = character(0)) {
  pool <- pool %>% filter(!peps %in% exclude)
  slice_sample(pool, n = min(n, nrow(pool)))
}

nn_input <- lapply(c(N = "N", C = "C"), function(t) {
  k   <- known_w %>% filter(term == t)
  sp  <- split_knowns(k)
  neg_c <- cand_w %>% filter(term == t, !gene %in% banned, !peps %in% k$peps)
  neg_r <- rand_w %>% filter(term == t, !peps %in% k$peps)
  n_r <- round(n_ctrl * random_frac)
  n_n <- if (is.null(near_w)) 0L else min(near_n, n_ctrl - n_r)
  ## A near-end negative INHERITS ITS PARENT'S SPLIT. Drawing it freely would put
  ## a window 2 residues from a val positive into train -- 34 of its 36 residues
  ## are the same residues, so the val end would effectively be in training.
  near_tr <- near_va <- NULL
  if (n_n > 0) {
    nw <- near_w %>% filter(term == t, !peps %in% k$peps)
    par_tr <- k$anchor[sp$train]; par_va <- k$anchor[sp$val]
    key <- function(acc, a) paste(acc, a)
    near_tr <- nw %>% filter(key(accession, parent) %in% key(k$accession[sp$train], par_tr))
    near_va <- nw %>% filter(key(accession, parent) %in% key(k$accession[sp$val],   par_va))
  }
  tr_neg <- bind_rows(draw_negatives(neg_c, n_ctrl - n_r - n_n),
                      draw_negatives(neg_r, n_r),
                      if (n_n > 0) draw_negatives(near_tr, n_n))
  va_neg <- bind_rows(draw_negatives(neg_c, n_ctrl - n_r - n_n, tr_neg$peps),
                      draw_negatives(neg_r, n_r, tr_neg$peps),
                      if (n_n > 0) draw_negatives(near_va, n_n, tr_neg$peps))
  list(train = bind_rows(k[sp$train, ], tr_neg),
       val   = bind_rows(k[sp$val, ],   va_neg),
       ## every candidate window of the terminus, knowns first so a candidate
       ## that coincides with a known end is the known copy
       all   = bind_rows(k, cand_w %>% filter(term == t)) %>% distinct(peps, .keep_all = TRUE))
})

knowns <- bind_rows(lapply(names(nn_input), function(t) {
  s <- nn_input[[t]]
  bind_rows(s$train %>% filter(known == 1) %>% mutate(split = "train"),
            s$val   %>% filter(known == 1) %>% mutate(split = "val"))
})) %>%
  mutate(AA = map(meta_data, "AA")) %>%
  select(peps, gene, accession, term, anchor, w_start, pep_id, pep_ids, pep_name, pep_start, pep_end,
         len, motif, insertion, at_candidate, split, AA, known_idx)

for (t in names(nn_input)) for (s in c("train", "val", "all")) {
  x <- nn_input[[t]][[s]]
  stopifnot(all(map_int(x$data, nrow) == 36L), all(lengths(x$known_idx) == 36L),
            all(unlist(x$known_idx) %in% LF_PEPEND_CLASSES))
  message(sprintf("%s %-5s %6d windows, %3d known", t, s, nrow(x), sum(x$known)))
}

saveRDS(list(nn_input = nn_input, all_params3 = all_params3, class_names = LF_PEPEND_CLASSES,
             knowns = knowns,
             params = list(n_ctrl = n_ctrl, random_frac = random_frac, val_frac = val_frac,
                           swap_splits = swap_splits,
                           sim_h = sim_h, seed = seed, banned = banned, built = Sys.time())),
        out_path)
message("wrote ", out_path)
