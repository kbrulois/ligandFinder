#!/usr/bin/env Rscript
## ---- every known peptide end, predicted twice: in-sample and held out ----------
## Two complementary folds were trained (10_7 --val-frac 0.5, and the same with
## --swap-splits), so every known end is TRAIN in exactly one fold and VAL in the
## other. This puts both predictions on one row-pair, with a `split` column
## saying which is which -- so the same end can be read as the model that fitted
## it sees it, and as a model that never did.
##
## Twice the number of predictions, as asked: 86 C ends and 92 N ends, each
## appearing once per fold.
##
##   Rscript inst/scripts/10_9v_compile_cv_predictions.R
##
## Output: ~/AF2_analysis/lf_cv_known_predictions.csv  (+ .rds)
## ------------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(tidyr); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
out_p <- path.expand(.opt("--out", "~/AF2_analysis/lf_cv_known_predictions.csv"))
A2    <- path.expand("~/AF2_analysis")

np <- reticulate::import("numpy", convert = FALSE)

## ---- 1. the knowns, and which side of each fold they fell on -------------------
## Read from the fold INPUTS rather than from a scan summary: the input is what
## defined the split, so train/val here cannot disagree with what was trained.
split_of <- function(fold, term) {
  x <- readRDS(file.path(A2, sprintf("lf_pepend_nn_input_fold%s.rds", fold)))$nn_input[[term]]
  bind_rows(
    x$train %>% filter(known == 1) %>% transmute(accession, anchor = as.integer(anchor),
                                                 pep_name, split = "train"),
    x$val   %>% filter(known == 1) %>% transmute(accession, anchor = as.integer(anchor),
                                                 pep_name, split = "val")) %>%
    mutate(fold = fold, term = term)
}
knowns <- bind_rows(lapply(c("A", "B"), function(f)
  bind_rows(lapply(c("N", "C"), function(t) split_of(f, t)))))
message(sprintf("knowns x folds: %d rows (%d distinct ends)",
                nrow(knowns), n_distinct(paste(knowns$accession, knowns$term, knowns$anchor))))

## ---- 2. the window score, and both ranks ---------------------------------------
## rank_candidates ranks among the motif-anchored candidate set the training
## summary uses; rank_scan ranks among EVERY mature anchor of every precursor.
## They answer different questions and routinely differ by orders of magnitude,
## so both travel.
win_of <- function(fold, term) {
  s <- readRDS(file.path(A2, sprintf("lf_pepend_scan_%s_fold%s.rds", term, fold)))
  sc <- s$scan %>% transmute(accession, anchor = as.integer(anchor),
                             cnn = score, cnn_sd = sd, rank_scan = rank_all)
  k <- s$knowns %>% transmute(accession, anchor = as.integer(anchor),
                              motif, at_candidate,
                              rank_candidates = if ("rank_cand" %in% names(.)) rank_cand else NA_integer_)
  left_join(sc, k, by = c("accession", "anchor")) %>%
    semi_join(k, by = c("accession", "anchor")) %>%
    mutate(fold = fold, term = term)
}
win <- bind_rows(lapply(c("A", "B"), function(f)
  bind_rows(lapply(c("N", "C"), function(t) win_of(f, t)))))

## ---- 3. the xgb window head, on the `all` split --------------------------------
## `all` holds every candidate window AND every known, in the order the fold
## input stores it, so the head's score is matched back by position.
xgb_of <- function(fold, term) {
  f <- file.path(A2, sprintf("lf_xgbhead_%s_fold%s.npz", term, fold))
  if (!file.exists(f)) return(NULL)
  z <- np$load(f, allow_pickle = TRUE)
  v <- as.numeric(reticulate::py_to_r(z[["xgb_all"]]))
  all_rows <- readRDS(file.path(A2, sprintf("lf_pepend_nn_input_fold%s.rds", fold)))$nn_input[[term]]$all
  if (length(v) != nrow(all_rows)) {
    message(sprintf("  xgb %s fold%s: %d scores vs %d `all` rows -- SKIPPED",
                    term, fold, length(v), nrow(all_rows)))
    return(NULL)
  }
  tibble(accession = all_rows$accession, anchor = as.integer(all_rows$anchor),
         xgb = v, fold = fold, term = term)
}
xgb <- bind_rows(lapply(c("A", "B"), function(f)
  bind_rows(lapply(c("N", "C"), function(t) xgb_of(f, t)))))

## ---- 4. the window-free residue models at that residue -------------------------
## THEIR SPLIT IS NOT THE WINDOW SPLIT. lf_resid divides by precursor with its
## own shuffle; the window set divides by sequence-similarity cluster over known
## windows. Fold A of one has no relationship to fold A of the other -- 36 of 66
## genes land on opposite sides -- so these columns carry their OWN split label
## and are read from whichever fold actually held that precursor out. Labelling
## them with the window split would call an in-sample value held out.
rp_split_of <- function(fold) {
  z <- np$load(file.path(A2, sprintf("lf_resid_preds_fold%s.npz", fold)), allow_pickle = TRUE)
  have <- as.character(reticulate::py_to_r(z$files))
  if (!all(c("in_val", "prot_idx", "accession") %in% have)) return(NULL)
  iv   <- as.integer(reticulate::py_to_r(z[["in_val"]]))
  it   <- as.integer(reticulate::py_to_r(z[["in_train"]]))
  prot <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
  acc  <- as.character(reticulate::py_to_r(z[["accession"]]))
  tibble(accession = acc[prot + 1L], in_val = iv, in_train = it) %>%
    group_by(accession) %>%
    summarise(resid_split = if (sum(in_val) > 0) "val" else
                            if (sum(in_train) > 0) "train" else NA_character_,
              .groups = "drop") %>%
    mutate(fold = fold)
}
rp_split <- bind_rows(lapply(c("A", "B"), rp_split_of))

rp_of <- function(fold) {
  f <- file.path(A2, sprintf("lf_resid_preds_fold%s.npz", fold))
  if (!file.exists(f)) return(NULL)
  z <- np$load(f, allow_pickle = TRUE)
  have <- as.character(reticulate::py_to_r(z$files))
  cols <- grep("^(pep_pocket|pep_other|ct_context|nt_context)__(mlp|xgb)$", have, value = TRUE)
  acc <- as.character(reticulate::py_to_r(z[["accession"]]))
  prot <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
  res  <- as.integer(reticulate::py_to_r(z[["resno"]]))
  d <- tibble(accession = acc[prot + 1L], anchor = res)
  for (cc in cols) d[[paste0("rp_", cc)]] <- as.numeric(reticulate::py_to_r(z[[cc]]))
  d$fold <- fold
  d
}
rp <- bind_rows(lapply(c("A", "B"), rp_of))

## ---- 5. one row per known end per fold ------------------------------------------
tab <- knowns %>%
  left_join(win, by = c("accession", "anchor", "fold", "term")) %>%
  left_join(xgb, by = c("accession", "anchor", "fold", "term")) %>%
  left_join(rp,  by = c("accession", "anchor", "fold")) %>%
  left_join(rp_split, by = c("accession", "fold")) %>%
  ## `split` is the WINDOW models; `resid_split` is the residue models. They
  ## disagree for most genes and must never be collapsed into one column.
  relocate(resid_split, .after = split)

## gene symbols, from the cache the models were built from
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
tab <- tab %>%
  mutate(gene = fc$prec$gene[match(accession, fc$prec$accession)], .before = 1) %>%
  arrange(term, gene, anchor, fold)

stopifnot(!any(is.na(tab$split)))
write.csv(tab, out_p, row.names = FALSE)
saveRDS(tab, sub("\\.csv$", ".rds", out_p))
message(sprintf("\n%d rows -> %s", nrow(tab), out_p))
message(sprintf("  per fold/split: %s",
                paste(capture.output(print(table(tab$fold, tab$split))), collapse = " | ")))
message(sprintf("  with a window score: %d | with xgb: %d | with residue models: %d",
                sum(!is.na(tab$cnn)), sum(!is.na(tab$xgb)),
                sum(!is.na(tab[["rp_pep_pocket__mlp"]]))))
message("  window split vs residue split (they are independent):")
print(table(window = tab$split, residue = tab$resid_split))
