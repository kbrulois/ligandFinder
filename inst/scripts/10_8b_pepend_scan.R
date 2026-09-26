#!/usr/bin/env Rscript
## ---- score EVERY possible peptide-end window in the secretome ------------------
## 10_8_pepend_train.R trains on the peptide-END set and ranks its held-out known
## ends among that set's `all` split -- but `all` holds only the MOTIF-anchored
## candidates (dibasic, amidation, precursor termini, plus a random sample),
## 41,870 windows for C. A known end at no such motif (25 of 86 on C) can never
## be ranked there, and neither can a novel end anywhere else.
##
## This scores every mature residue of every precursor instead: 2,980,852 C
## anchors over 5,616 precursors. That cannot go through the R array contract
## (~22 GB as a list of 36x26 tibbles), so the windows are built and scored in
## chunks in python -- lf_dcnn/scan.py -- and only the scalar window score comes
## back.
##
##   Rscript inst/scripts/10_8b_pepend_scan.R --term C
##   ... --chunk-windows 200000 --refresh-residues
##
## Needs 10_8_pepend_train.R to have run first: the scan rebuilds EVERY ensemble
## member from <run>_isolated/out/members/member_NNN_<term>.weights.h5, not just
## member 0 (which is all lf_dcnn_run() restores).
##
## Ranks are on the RAW ensemble mean. The Platt calibration the training run
## fits is monotone, so it cannot change any ranking; it only rescales.
##
## Output: ~/AF2_analysis/lf_pepend_scan_<term>.rds
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
.flag <- function(flag) flag %in% .args

term     <- toupper(.opt("--term", "C"))
cache_p  <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
in_path  <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
run_p    <- path.expand(.opt("--run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
res_npz  <- path.expand(.opt("--residues-npz", "~/AF2_analysis/lf_pepend_residues.npz"))
out_p    <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
chunk_w  <- as.integer(.opt("--chunk-windows", "200000"))
batch_sz <- as.integer(.opt("--batch-size", "4096"))
top_n    <- as.integer(.opt("--top", "50"))

stopifnot(term %in% c("N", "C"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))

## the INSTALLED ligandFinder shadows this checkout in lf_dcnn_path(); pin it
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)

message("reading ", cache_p, " ...")
fc <- readRDS(cache_p)
prec <- fc$prec; feats <- fc$feats
stopifnot(identical(colnames(feats[[1]]), fc$all_params3),
          identical(names(feats), prec$accession),
          all(vapply(feats, nrow, integer(1)) == nchar(prec$seq)))

## ---- 1. residues -> one float32 array + offsets -------------------------------
if (!file.exists(res_npz) || .flag("--refresh-residues")) {
  L   <- vapply(feats, nrow, integer(1))
  off <- c(0L, cumsum(L)[-length(L)])                  # 0-based row of residue 1
  tot <- sum(L)
  message(sprintf("packing %d precursors, %s residues x %d channels ...",
                  length(feats), format(tot, big.mark = ","), length(fc$all_params3)))
  ## channel-major (C, total): its column-major layout IS numpy's row-major
  ## (total, C), so no transpose of the big array is ever needed
  M <- matrix(0, length(fc$all_params3), tot)
  for (i in seq_along(feats)) M[, (off[i] + 1L):(off[i] + L[i])] <- t(feats[[i]])
  ## transpose on the numpy side: t(M) in R would be another 636 MB copy
  np$savez(res_npz,
           feat   = np$ascontiguousarray(np$asarray(M, dtype = "float32")$T),
           offset = np$asarray(off, dtype = "int64"),
           n_prot = np$asarray(prec$n_prot, dtype = "int64"),
           c_prot = np$asarray(prec$c_prot, dtype = "int64"))
  rm(M); invisible(gc())
  message("residues -> ", res_npz)
} else {
  message("reusing ", res_npz, " (--refresh-residues to rebuild)")
}
rm(feats); invisible(gc())

## ---- 2. scan -----------------------------------------------------------------
out_dir <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out")
if (!dir.exists(file.path(out_dir, "members")))
  stop("no trained members under ", out_dir, " -- run 10_8_pepend_train.R first",
       call. = FALSE)
scan_npz <- sub("\\.rds$", ".npz", out_p)

## How many members the run trained, so a short ensemble is an error rather than
## a quietly smaller one. `python -m lf_dcnn train` used to write weights for
## member 0 only (keep_models=(member == 0)); a run from before that was fixed
## has 1 weights file for n_seeds members and must be retrained with --refresh.
sum_p   <- sub("\\.rds$", "_summary.rds", run_p)
n_seeds <- if (file.exists(sum_p)) as.integer(readRDS(sum_p)$n_seeds) else NULL
n_seeds <- as.integer(.opt("--n-seeds", n_seeds %||% NA_integer_))
n_have  <- length(Sys.glob(file.path(out_dir, "members",
                                     sprintf("member_*_%s.weights.h5", term))))
if (!is.na(n_seeds) && n_have != n_seeds)
  stop(sprintf(paste("the run trained %d members but only %d has/have weights.",
                     "Retrain with --refresh (member weights are needed to",
                     "rebuild the ensemble for the scan)."), n_seeds, n_have),
       call. = FALSE)
message(sprintf("\nscanning with %d member(s) ...", n_have))
mod$scan$scan_to_npz(residues_npz = res_npz, out_dir = out_dir, term = term,
                     out_npz = scan_npz, n_seeds = as.integer(n_have),
                     chunk_windows = as.integer(chunk_w),
                     batch_size = as.integer(batch_sz), verbose = TRUE)

z <- np$load(scan_npz)
got <- list(score = as.numeric(reticulate::py_to_r(z[["score"]])),
            sd    = as.numeric(reticulate::py_to_r(z[["sd"]])),
            prot  = as.integer(reticulate::py_to_r(z[["prot_idx"]])),
            anchor = as.integer(reticulate::py_to_r(z[["anchor"]])))

scan_t <- tibble(accession = prec$accession[got$prot + 1L],
                 gene = prec$gene[got$prot + 1L],
                 anchor = got$anchor, score = got$score, sd = got$sd) %>%
  mutate(rank_all = as.integer(rank(-score, ties.method = "first")))

message(sprintf("\nscored %s windows over %d precursors",
                format(nrow(scan_t), big.mark = ","), length(unique(scan_t$accession))))

## ---- 3. where do the known ends land, now that nothing is excluded? -----------
x <- readRDS(in_path)
kn <- x$knowns %>% filter(term == !!term) %>%
  select(accession, gene, anchor, pep_name, split, motif, at_candidate)
kn <- kn %>% left_join(scan_t %>% select(accession, anchor, score, sd, rank_all),
                       by = c("accession", "anchor"))

miss <- sum(is.na(kn$rank_all))
if (miss) message("NOTE ", miss, " known end(s) not in the scan (precursor absent)")

message(sprintf("\nknown %s ends among ALL %s windows (n = %d):", term,
                format(nrow(scan_t), big.mark = ","), sum(!is.na(kn$rank_all))))
for (grp in list(list("all knowns", rep(TRUE, nrow(kn))),
                 list("held out (val)", kn$split == "val"),
                 list("at a candidate motif", kn$at_candidate),
                 list("at NO candidate motif", !kn$at_candidate))) {
  s <- kn$rank_all[grp[[2]] & !is.na(kn$rank_all)]
  if (!length(s)) next
  message(sprintf("  %-22s n %3d   median rank %9s   top 1k %3d   top 10k %3d   top 100k %3d",
                  grp[[1]], length(s), format(median(s), big.mark = ","),
                  sum(s <= 1000), sum(s <= 10000), sum(s <= 100000)))
}

message("\n  by motif:")
kn %>% filter(!is.na(rank_all)) %>% group_by(motif) %>%
  summarise(n = n(), median_rank = median(rank_all), top10k = sum(rank_all <= 10000),
            .groups = "drop") %>% arrange(median_rank) %>%
  as.data.frame() %>% print(row.names = FALSE)

message(sprintf("\ntop %d windows overall (known ends marked):", top_n))
scan_t %>% arrange(rank_all) %>% head(top_n) %>%
  left_join(kn %>% select(accession, anchor, pep_name), by = c("accession", "anchor")) %>%
  mutate(known = ifelse(is.na(pep_name), "", pep_name)) %>%
  select(rank_all, gene, accession, anchor, score, sd, known) %>%
  as.data.frame() %>% print(row.names = FALSE, digits = 3)

saveRDS(list(term = term, scan = scan_t, knowns = kn, n_members = NA_integer_),
        out_p)
message("\nscan -> ", out_p, "  (arrays in ", basename(scan_npz), ")")
