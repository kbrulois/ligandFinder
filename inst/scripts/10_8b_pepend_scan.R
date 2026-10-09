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
## READ THE PEAK RANKS, NOT THE WINDOW RANKS. Adjacent anchors differ by one
## residue out of 36, so one real signal becomes a long run of near-identical
## windows and a few loci swallow the top of a window-ranked list (top 1,000
## windows: 59 genes; top 1,000 peaks: 447). The script collapses each run to its
## local maximum and ranks those too; that is the number to quote.
##
## Output: ~/AF2_analysis/lf_pepend_scan_<term>.rds -- `scan` (every window),
## `peaks` (local maxima, ranked), `knowns_peak` (each known end's nearest peak)
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
## per-window seed summaries, in the same currency as the candidate tables
thresh   <- as.numeric(.opt("--thresh", "0.2"))
top_k    <- as.integer(.opt("--top-k", "5"))

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
  message(sprintf("packing %d precursors, %s residues x %d channels ...",
                  length(feats), format(sum(nchar(prec$seq)), big.mark = ","),
                  length(fc$all_params3)))
  lf_dcnn_pack_residues(prec, feats, fc$all_params3, res_npz,
                        py_path = file.path(ROOT, "inst", "python"))
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
if (file.exists(scan_npz) && !.flag("--refresh-scan")) {
  message(sprintf("\nreusing %s (--refresh-scan to rescore)", basename(scan_npz)))
} else {
  message(sprintf("\nscanning with %d member(s) ...", n_have))
  mod$scan$scan_to_npz(residues_npz = res_npz, out_dir = out_dir, term = term,
                       out_npz = scan_npz, n_seeds = as.integer(n_have),
                       chunk_windows = as.integer(chunk_w),
                       batch_size = as.integer(batch_sz),
                       top_k = as.integer(top_k),
                       thresholds = reticulate::tuple(as.numeric(thresh)),
                       verbose = TRUE)
}

z <- np$load(scan_npz)
keys <- as.character(reticulate::py_to_r(z$files))
got <- list(score = as.numeric(reticulate::py_to_r(z[["score"]])),
            sd    = as.numeric(reticulate::py_to_r(z[["sd"]])),
            prot  = as.integer(reticulate::py_to_r(z[["prot_idx"]])),
            anchor = as.integer(reticulate::py_to_r(z[["anchor"]])))
## the per-window order statistic and threshold counts, whatever they were
## named, plus the insertion head's per-class means/votes when the members carry
## that head (ins_p_*, ins_votes_*, ins_sd_inserting -- absent otherwise)
extra <- setdiff(grep("^score_top|^n_seeds_gt_|^ins_", keys, value = TRUE), names(got))
for (k in extra) got[[k]] <- as.numeric(reticulate::py_to_r(z[[k]]))
if (length(extra)) message("  carrying: ", paste(extra, collapse = ", "))

scan_t <- bind_cols(
  tibble(accession = prec$accession[got$prot + 1L],
         gene = prec$gene[got$prot + 1L],
         anchor = got$anchor, score = got$score, sd = got$sd),
  as_tibble(got[extra])) %>%
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

## ---- 4. local peaks: the honest ranking --------------------------------------
## Neighbouring anchors differ by one residue out of 36, so a single real signal
## shows up as a long run of near-identical windows -- rank the windows and a
## handful of loci swallow the whole top of the list. Collapse each run to its
## local maximum (nothing scoring higher within +/- `peak_w` residues of it) and
## rank those instead.
suppressMessages(library(data.table))
peak_w <- as.integer(.opt("--peak-window", "5"))
dt <- as.data.table(scan_t); setorder(dt, accession, anchor)
shifts <- c(lapply(seq_len(peak_w), function(k) list(k, "lag")),
            lapply(seq_len(peak_w), function(k) list(k, "lead")))
dt[, rmax := do.call(pmax, c(list(score),
     lapply(shifts, function(z) shift(score, z[[1]], fill = -Inf, type = z[[2]])))),
   by = accession]
pk <- dt[score >= rmax][order(-score)][, rank_peak := .I][, rmax := NULL]

message(sprintf("\nlocal peaks (+/-%d): %s of %s windows (%.1f%%)",
                peak_w, format(nrow(pk), big.mark = ","),
                format(nrow(dt), big.mark = ","), 100 * nrow(pk) / nrow(dt)))
message(sprintf("  concentration: top 1,000 WINDOWS span %d genes, top 1,000 PEAKS %d",
                n_distinct(head(arrange(scan_t, rank_all), 1000)$gene),
                uniqueN(head(pk, 1000)$gene)))

## A known end need not be a peak itself; the nearest peak within +/-peak_w
## stands in for it. Knowns with NO peak nearby are counted, not hidden -- the
## model simply has no local opinion there.
near <- pk[, .(accession, panchor = anchor, prank = rank_peak)]
kp <- merge(as.data.table(kn)[, .(accession, anchor, pep_name, split, at_candidate)],
            near, by = "accession", allow.cartesian = TRUE)[abs(anchor - panchor) <= peak_w]
kp <- kp[order(prank)][, .SD[1], by = .(pep_name)]
message(sprintf("  known %s ends with a peak within +/-%d: %d of %d (%d have none)",
                term, peak_w, nrow(kp), nrow(kn), nrow(kn) - nrow(kp)))
for (g in list(list("all knowns", rep(TRUE, nrow(kp))),
               list("held out (val)", kp$split == "val"),
               list("at a candidate motif", kp$at_candidate),
               list("at NO candidate motif", !kp$at_candidate))) {
  v <- kp$prank[g[[2]]]
  if (!length(v)) next
  message(sprintf("  %-22s n %3d   median peak-rank %7s   top 100 %3d   top 500 %3d   top 1k %3d",
                  g[[1]], length(v), format(as.integer(median(v)), big.mark = ","),
                  sum(v <= 100), sum(v <= 500), sum(v <= 1000)))
}

message(sprintf("\ntop %d local peaks (a known end within +/-%d is named, with its offset):",
                top_n, peak_w))
merge(head(pk, top_n),
      kp[, .(rank_peak = prank, pep_name, known_anchor = anchor)],
      by = "rank_peak", all.x = TRUE)[order(rank_peak)] %>%
  transmute(rank_peak, gene, accession, anchor, score = round(score, 5),
            sd = round(sd, 5),
            known = ifelse(is.na(pep_name), "",
                           sprintf("%s (%+d)", pep_name, known_anchor - anchor))) %>%
  as.data.frame() %>% print(row.names = FALSE)

saveRDS(list(term = term, n_members = n_have, peak_w = peak_w,
             scan = scan_t, peaks = as_tibble(pk), knowns = kn, knowns_peak = as_tibble(kp)),
        out_p)
message("\nscan -> ", out_p, "  (arrays in ", basename(scan_npz), ")")
