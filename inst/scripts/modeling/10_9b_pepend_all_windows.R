#!/usr/bin/env Rscript
## ---- every peptide-end candidate window, scored -------------------------------
## One row per window of the peptide-end set's `all` split (41,870 for C), which
## unlike the production set is NOT all dibasic-anchored:
##
##   db       33,823  anchor at d-1 of a dibasic pair -- the production set's kind
##   db_amid   2,404  anchor at d-2, where d-1 is the amidation glycine
##   pep_end   5,557  the precursor's own last residue
##   known        86  a docked known peptide's actual C terminus
##
## Columns per model, as in 10_9: mean over all 20 seeds, mean over the 5 seeds
## that score THAT window highest, how many seeds clear --thresh, across-seed sd.
##
##   pepend_*  the peptide-end model, read from its run's member files
##   db_*      the production dibasic-anchored model -- joined from 10_9's output
##             on (accession, dibasic position), NA where the production set has
##             no such window. That NA pattern is the point of this file.
##
## Note a db_amid anchor is d-2 of the SAME pair a db anchor puts at d-1, so both
## join to the one production window for that dibasic: the production model scores
## the pair, and the two pepend rows score the two cleavage hypotheses at it.
## A pep_end anchor has no dibasic at all and is always NA.
##
## The db model cannot simply be RUN over these windows: its 20-seed arm was
## trained before `python -m lf_dcnn train --member K` persisted every member's
## weights, so only member 0's survive. Joining 10_9's verified per-seed
## summaries is exact where a counterpart exists; scoring member 0 alone would be
## a different, weaker number wearing the same name.
##
##   Rscript inst/scripts/10_9b_pepend_all_windows.R --term C
##
## Output: ~/AF2_analysis/lf_pepend_all_windows_C.csv
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(readr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term    <- toupper(.opt("--term", "C"))
in_path <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
run_p   <- path.expand(.opt("--run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
db_csv  <- path.expand(.opt("--db-csv", sprintf("~/AF2_analysis/lf_pepend_vs_db_windows_%s.csv", term)))
scan_p  <- path.expand(.opt("--scan", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
out_csv <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_all_windows_%s.csv", term)))
thresh  <- as.numeric(.opt("--thresh", "0.2"))
top_k   <- as.integer(.opt("--top-k", "5"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
reticulate::use_virtualenv("r-tensorflow", required = TRUE)
npc <- reticulate::import("numpy", convert = TRUE)

## ---- 1. the windows -----------------------------------------------------------
x   <- readRDS(in_path)
all <- x$nn_input[[term]]$all
n   <- nrow(all)
message(sprintf("%s peptide-end %s windows", format(n, big.mark = ","), term))
print(as.data.frame(count(all, win_type)), row.names = FALSE)

## ---- 2. the peptide-end model's per-seed scores -------------------------------
out_dir <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out")
mem <- sort(Sys.glob(file.path(out_dir, "members", sprintf("member_*_%s.npz", term))))
if (!length(mem)) stop("no member files under ", out_dir, call. = FALSE)
S <- vapply(mem, function(f) as.numeric(npc$load(f)[[paste0(term, "__global__x")]]),
            numeric(n))                                       # (n, n_seeds)
message(sprintf("pepend model: %d seeds x %s windows", ncol(S), format(nrow(S), big.mark = ",")))

## The member files carry no row labels, so prove the order matches `all` before
## trusting the join: the cached run's ensemble mean was built from these members
## over exactly this split, so the two must agree elementwise.
run <- readRDS(run_p)
pr  <- as.numeric(run$out$pred_raw)
stopifnot(length(pr) == n)
d <- max(abs(rowMeans(S) - pr))
message(sprintf("row order check vs the cached run's pred_raw: max|diff| %.3g", d))
if (d >= 1e-6) stop("member rows do not line up with the `all` split", call. = FALSE)

## ---- 3. summarise -------------------------------------------------------------
summarise_model <- function(S, tag) {
  k <- min(top_k, ncol(S))
  out <- tibble(
    mean_all = rowMeans(S),
    mean_top = apply(S, 1, function(v) mean(sort(v, decreasing = TRUE)[seq_len(k)])),
    n_gt     = apply(S, 1, function(v) sum(v > thresh)),
    sd       = apply(S, 1, sd))
  names(out) <- sprintf("%s_%s", tag,
                        c("mean_allseeds", sprintf("mean_top%d", k),
                          sprintf("n_seeds_gt_%s", format(thresh)), "sd_seeds"))
  out
}

## Which dibasic, if any, this anchor belongs to. NOT simply anchor + 1: a plain
## dibasic anchor is d-1, but an amidated one is d-2 because the glycine leaves,
## so both can point at the SAME pair. `win_type` carries the anchoring rule for
## candidates; it is "known" for the knowns, whose own sequence context (`motif`,
## populated for those rows only) says which rule applies. An anchor at a
## monobasic, at the precursor terminus or at nothing has no dibasic to key on.
res <- bind_cols(
  all %>% transmute(peps, gene, accession, term, anchor, win_type,
                    motif, at_candidate, known, pep_name,
                    pep_ids = vapply(pep_ids, function(v) paste(v, collapse = ";"), ""),
                    pep_start, pep_end, len, insertion, w_start, n_prot, c_prot,
                    db_pos = case_when(
                      win_type == "db"            ~ anchor + 1L,
                      win_type == "db_amid"       ~ anchor + 2L,
                      motif    == "dibasic"       ~ anchor + 1L,
                      motif    == "G + dibasic"   ~ anchor + 2L,
                      TRUE                        ~ NA_integer_)),
  summarise_model(S, "pepend"))
## every row must still be one window
stopifnot(nrow(res) == n, !anyDuplicated(res$peps))

## ---- 4. the db model, where a production counterpart exists -------------------
## 10_9 verified each production window's coordinates by bit-exact rebuild; join
## only those rows (aligned_ok), keyed on the dibasic position.
if (file.exists(db_csv)) {
  db <- read_csv(db_csv, show_col_types = FALSE) %>%
    filter(aligned_ok) %>%
    transmute(accession, db_pos, db_peps = peps,
              db_mean_allseeds, db_mean_top5,
              !!sprintf("db_n_seeds_gt_%s", format(thresh)) :=
                .data[[sprintf("db_n_seeds_gt_%s", format(thresh))]],
              db_sd_seeds)
  before <- nrow(res)
  res <- res %>% left_join(db, by = c("accession", "db_pos"))
  stopifnot(nrow(res) == before)
  message(sprintf("\ndb model joined onto %s of %s windows (%.1f%%)",
                  format(sum(!is.na(res$db_mean_allseeds)), big.mark = ","),
                  format(n, big.mark = ","), 100 * mean(!is.na(res$db_mean_allseeds))))
  message("  by win_type:")
  res %>% group_by(win_type) %>%
    summarise(n = n(), with_db = sum(!is.na(db_mean_allseeds)),
              pct = round(100 * with_db / n, 1), .groups = "drop") %>%
    as.data.frame() %>% print(row.names = FALSE)
} else {
  message("\n", db_csv, " not found -- db columns omitted (run 10_9 first)")
}

## ---- 5. cross-check against the exhaustive scan -------------------------------
## Same model, same 20 members, same windows built two different ways: 10_7 sliced
## them in R into the contract, lf_dcnn.scan rebuilds them in python. The scores
## must agree, which is an end-to-end check on the scan.
if (file.exists(scan_p)) {
  sc <- readRDS(scan_p)
  res <- res %>%
    left_join(sc$scan %>% select(accession, anchor, scan_score = score,
                                 scan_rank_all = rank_all),
              by = c("accession", "anchor")) %>%
    left_join(sc$peaks %>% select(accession, anchor, scan_rank_peak = rank_peak),
              by = c("accession", "anchor"))
  ok <- !is.na(res$scan_score)
  message(sprintf("\nscan cross-check on %s windows: max|pepend_mean - scan| %.3g",
                  format(sum(ok), big.mark = ","),
                  max(abs(res$pepend_mean_allseeds[ok] - res$scan_score[ok]))))
  message(sprintf("  windows that are also a local peak (+/-%d): %s",
                  sc$peak_w, format(sum(!is.na(res$scan_rank_peak)), big.mark = ",")))
}

write_csv(res, out_csv)
message(sprintf("\n%s rows x %d columns -> %s", format(nrow(res), big.mark = ","),
                ncol(res), out_csv))

## ---- 6. what the file says ----------------------------------------------------
message(sprintf("\nmedian pepend_mean_allseeds by win_type (threshold %.2f, top-%d):",
                thresh, top_k))
res %>% group_by(win_type) %>%
  summarise(n = n(),
            median_mean = round(median(pepend_mean_allseeds), 4),
            median_top5 = round(median(pepend_mean_top5), 4),
            all20_gt = sum(.data[[sprintf("pepend_n_seeds_gt_%s", format(thresh))]] == ncol(S)),
            none_gt  = sum(.data[[sprintf("pepend_n_seeds_gt_%s", format(thresh))]] == 0),
            .groups = "drop") %>%
  arrange(desc(median_mean)) %>% as.data.frame() %>% print(row.names = FALSE)

message("\ntop 15 windows overall:")
res %>% arrange(desc(pepend_mean_allseeds)) %>% head(15) %>%
  transmute(gene, anchor, win_type, motif,
            pepend = round(pepend_mean_allseeds, 4),
            top5 = round(pepend_mean_top5, 4),
            n_gt = .data[[sprintf("pepend_n_seeds_gt_%s", format(thresh))]],
            db = ifelse(is.na(db_mean_allseeds), NA, round(db_mean_allseeds, 4)),
            known = ifelse(is.na(pep_name), "", pep_name)) %>%
  as.data.frame() %>% print(row.names = FALSE)
