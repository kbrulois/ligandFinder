#!/usr/bin/env Rscript
## ---- train a window model on the peptide-END set ------------------------------
## The production model (10_1dcnn_new6.R) trains on the dibasic-anchored set from
## 9.2_add_contact_data.R. This trains the same architecture on the peptide-END
## set built by 10_7_pepend_windows.R, which anchors on the peptide's own
## terminal residue instead and drops the DB and gap classes.
##
## ONE TERMINUS PER RUN. The anchor sits at window position 28 for C and 8 for N,
## so the position masks differ between them while `mask_matrix_cat` does not --
## Config.pepend(term) derives the masks from the anchor and pins `term_order`,
## and lf_dcnn_align_terms() refuses to widen it. Training both from one Config
## would score one terminus through the other's mask.
##
##   Rscript inst/scripts/10_8_pepend_train.R --term C
##   ... --n-seeds 20 --refresh
##
## Reports, for the held-out (validation) known peptide ends, where they rank
## among the candidate windows the set scores -- the same currency as
## 10_5_benchmark_window_model.R, so the numbers are comparable to the production
## model's. NOTE those candidates are only the motif-anchored ones, so the known
## ends sitting at no candidate motif (25 of 86 on C) are excluded from the rank
## metric; 10_8b scores EVERY position instead.
##
## Output: ~/AF2_analysis/lf_pepend_run_<term>.rds  (cache + weights beside it)
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
.flag <- function(flag) flag %in% .args

term    <- toupper(.opt("--term", "C"))
in_path <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
n_seeds <- as.integer(.opt("--n-seeds", "20"))
verbose <- as.integer(.opt("--verbose", "0"))
out_p   <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))

stopifnot(term %in% c("N", "C"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))

message("reading ", in_path, " ...")
x <- readRDS(in_path)
stopifnot(identical(as.character(x$class_names),
                    c("CT_cleavage_context", "NT_cleavage_context",
                      "pep_other", "pep_pocket", "padding", "none")))

## ONE terminus: the Config's masks are specific to it (see the header)
nn_input <- x$nn_input[term]
message(sprintf("terminus %s: train %d, val %d, all %d windows; %d channels",
                term, nrow(nn_input[[term]]$train), nrow(nn_input[[term]]$val),
                nrow(nn_input[[term]]$all), length(x$all_params3)))

## ---- train -------------------------------------------------------------------
## isolated = TRUE: the second model trained in one Keras/TF process dies
## intermittently, and an ensemble trains n_seeds of them.
## `path`: lf_dcnn_path() prefers system.file("python", package = "ligandFinder"),
## i.e. the INSTALLED copy, which silently shadows this checkout -- pin it to the
## tree this script came from, or Config.pepend() is simply not there.
run <- lf_dcnn_run(nn_input, channel_names = x$all_params3,
                   pepend = term, n_seeds = n_seeds, isolated = TRUE,
                   verbose = verbose, cache = out_p, refresh = .flag("--refresh"),
                   keep_arrays = FALSE, path = file.path(ROOT, "inst", "python"))

cfg <- run$config
message(sprintf("\ntrained: %s, %d params, %d seeds", cfg$trunk,
                as.integer(run$params[[term]]), n_seeds))
message(sprintf("masks: nt %s  mid %s  ct %s",
                paste(as.integer(cfg$nt_span), collapse = "-"),
                paste(as.integer(cfg$mid_span), collapse = "-"),
                paste(as.integer(cfg$ct_span), collapse = "-")))

## ---- validation metrics ------------------------------------------------------
suppressMessages(library(yardstick))
val_t <- tibble(score = as.numeric(run$val_scores),
                truth = factor(as.integer(run$val_labels), levels = c(0, 1)))
message(sprintf("\nvalidation (n = %d, %d positive):  PR AUC %.3f   ROC AUC %.3f",
                nrow(val_t), sum(val_t$truth == "1"),
                pr_auc_vec(val_t$truth, val_t$score, event_level = "second"),
                roc_auc_vec(val_t$truth, val_t$score, event_level = "second")))

## ---- where the held-out known ends rank among the candidates -----------------
## Same convention as 10_5: within the terminus, best first, a known's rank
## counted against the UNKNOWN candidate windows only.
all_t <- nn_input[[term]]$all %>%
  mutate(pred = as.numeric(run$pred), pred_raw = as.numeric(run$pred_raw),
         pred_sd = as.numeric(run$pred_sd))
unk <- sort(all_t$pred[all_t$known == 0], decreasing = TRUE)
rank_among_unknown <- function(p) vapply(p, function(v) sum(unk > v) + 1L, integer(1))

## `motif` / `at_candidate` are already columns of the `all` set; only the
## train/val split has to come from the knowns table
kn <- all_t %>% filter(known == 1) %>%
  mutate(rank_cand = rank_among_unknown(pred)) %>%
  left_join(x$knowns %>% select(peps, split), by = "peps")

held <- kn %>% filter(split == "val")
message(sprintf("\nheld-out known %s ends (n = %d of %d candidate windows):",
                term, nrow(held), nrow(all_t)))
message(sprintf("  median candidate rank %s;  in top 100: %d;  top 500: %d;  top 1000: %d",
                format(median(held$rank_cand)), sum(held$rank_cand <= 100),
                sum(held$rank_cand <= 500), sum(held$rank_cand <= 1000)))
message("  by motif (what sits just outside the terminus):")
held %>% group_by(motif) %>%
  summarise(n = n(), median_rank = median(rank_cand), top500 = sum(rank_cand <= 500),
            .groups = "drop") %>%
  arrange(median_rank) %>% as.data.frame() %>% print(row.names = FALSE)

saveRDS(list(term = term, n_seeds = n_seeds, knowns_ranked = kn,
             val = val_t, config = cfg$to_dict()),
        sub("\\.rds$", "_summary.rds", out_p))
message("\nsummary -> ", sub("\\.rds$", "_summary.rds", out_p))
message("weights  -> ", paste0(tools::file_path_sans_ext(out_p), "_isolated/out/members/"))
