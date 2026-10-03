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
ins_head <- "--ins-head" %in% .args
## Config default is 1.0, against loss_weight_global = 0.05 -- i.e. the insertion
## head, trained on 57 labelled windows, carries 20x the weight of the head that
## is actually shipped. --ins-weight is how that gets tested.
ins_w    <- as.numeric(.opt("--ins-weight", "NA"))
## --positives restricts the POSITIVES of train and val to one insertion class,
## by `lf_pepend_insertion_class()` -- the same definition the 3-way head was
## predicting, so "train only on inserting ends" is the head's own question moved
## out of a loss term and into the data. Negatives are untouched, and `all` (the
## scoring set) is untouched: the model still has to rank every candidate.
##
## NOT the `insertion` column 10_7 writes. That comes from the docked model's
## ins_term/ins_ind -- which terminus it enters and how deep -- while this asks
## where the pocket labels fall inside the 36-residue window. They disagree.
pos_cls  <- .opt("--positives", NA_character_)
## --plm attaches per-residue protein-language-model channels to every window
## BEFORE training, so the embedding enters at the trunk and the convolutions see
## it, rather than being bolted onto a head downstream. 32 PCA channels on top of
## the 26 hand-built ones.
plm_p    <- .opt("--plm", NA_character_)
seq_p    <- path.expand(.opt("--sequences", "~/AF2_analysis/lf_plm/sequences.parquet"))
## A different head is a different model, so it gets its own cache rather than
## relying on the fingerprint alone to tell them apart -- and so does a different
## WEIGHT on that head. This matters more than it looks: lf_dcnn_run_isolated
## reuses any member_*.npz already sitting in <cache>_isolated/out/members/ on
## its presence alone, so two arms sharing an --out would silently train one and
## report it as both.
out_p   <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_run_%s%s%s%s.rds",
                              term, if (ins_head) "_ins" else "",
                              if (ins_head && !is.na(ins_w))
                                paste0("_w", sub("\\.", "", format(ins_w, scientific = FALSE)))
                              else "",
                              if (!is.na(pos_cls)) paste0("_pos", pos_cls) else "",
                              ## a different channel set is a different model, and
                              ## shares the member-reuse trap, so it gets its own cache
                              if (!is.na(plm_p))
                                paste0("_", sub("_pca32$", "",
                                                tools::file_path_sans_ext(basename(plm_p))))
                              else "")))

stopifnot(term %in% c("N", "C"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "pepend_windows.R"))   # lf_pepend_insertion_class
source(file.path(ROOT, "R", "plm_features.R"))    # lf_plm_read / lf_plm_attach

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

## ---- PLM channels into the trunk, if asked -------------------------------------
if (!is.na(plm_p)) {
  plm_p <- path.expand(plm_p)
  message("attaching PLM channels from ", basename(plm_p), " ...")
  plm  <- lf_plm_read(plm_p)
  ## `sequences` only has to resolve gene -> accession, and the embedding ran on
  ## 5,182 proteins while this window set spans 5,616 -- so use the residue
  ## cache's own gene/accession table, which covers all of them.
  fcp  <- readRDS(path.expand(.opt("--feature-cache",
                    "~/AF2_analysis/lf_pepend_residue_cache.rds")))$prec
  seqs <- fcp[, c("gene", "accession")]

  ## lf_plm_attach is strict on purpose: a protein with no embedding is an error,
  ## not a silent zero block. Supply explicit zero blocks for the uncovered ones
  ## instead, so the strictness still holds for everything else and the count of
  ## proteins training on zeros is stated rather than hidden.
  need <- unique(unlist(lapply(nn_input, \(tm) unlist(lapply(tm, \(s) s$accession)))))
  miss <- setdiff(as.character(need), names(plm$by_acc))
  if (length(miss)) {
    for (acc in miss) {
      sq <- fcp$seq[match(acc, fcp$accession)]
      if (is.na(sq)) next
      aa <- strsplit(sq, "")[[1]]
      m  <- matrix(0, nrow = length(aa), ncol = length(plm$channels),
                   dimnames = list(as.character(seq_along(aa)), plm$channels))
      attr(m, "AA") <- setNames(aa, as.character(seq_along(aa)))
      plm$by_acc[[acc]] <- m
    }
    message(sprintf("  PLM covers %d of %d precursors in this set; %d get ZERO channels",
                    length(need) - length(miss), length(need), length(miss)))
  }
  att  <- lf_plm_attach(nn_input, plm, seqs, x$all_params3)
  nn_input       <- att$nn_input
  x$all_params3  <- att$channels
  message(sprintf("channels: %d -> %d", length(att$channels) - length(plm$channels),
                  length(att$channels)))
}

## ---- restrict the positives, if asked -----------------------------------------
if (!is.na(pos_cls)) {
  .keep <- function(d, what) {
    if (!nrow(d)) return(d)
    pos <- which(d$known == 1L)
    if (!length(pos)) return(d)
    cls <- lf_pepend_insertion_class(d$known_idx[pos], d$term[pos], d$len[pos])
    drop <- pos[!cls %in% pos_cls]
    message(sprintf("  %-5s positives %d -> %d (dropped %s)", what, length(pos),
                    length(pos) - length(drop),
                    paste(sprintf("%d %s", table(cls[!cls %in% pos_cls]),
                                  names(table(cls[!cls %in% pos_cls]))), collapse = ", ")))
    if (length(drop)) d[-drop, ] else d
  }
  message("restricting positives to insertion class '", pos_cls, "':")
  nn_input[[term]]$train <- .keep(nn_input[[term]]$train, "train")
  nn_input[[term]]$val   <- .keep(nn_input[[term]]$val,   "val")
  ## `all` is deliberately left whole -- it is the scoring set, and a model that
  ## only learns inserting ends still has to be ranked against every candidate.
}

## ---- train -------------------------------------------------------------------
## isolated = TRUE: the second model trained in one Keras/TF process dies
## intermittently, and an ensemble trains n_seeds of them.
## `path`: lf_dcnn_path() prefers system.file("python", package = "ligandFinder"),
## i.e. the INSTALLED copy, which silently shadows this checkout -- pin it to the
## tree this script came from, or Config.pepend() is simply not there.
.extra <- if (ins_head && !is.na(ins_w)) list(loss_weight_ins = ins_w) else list()
run <- do.call(lf_dcnn_run, c(list(nn_input, channel_names = x$all_params3,
                   pepend = term, ins_head = ins_head),
                   .extra, list(
                   n_seeds = n_seeds, isolated = TRUE,
                   verbose = verbose, cache = out_p, refresh = .flag("--refresh"),
                   keep_arrays = FALSE, path = file.path(ROOT, "inst", "python"))))

cfg <- run$config
message(sprintf("loss weights: global %.3g, per_index %.3g%s",
                cfg$loss_weight_global, cfg$loss_weight_per_index,
                if (ins_head) sprintf(", ins %.3g", cfg$loss_weight_ins) else ""))
message(sprintf("\ntrained: %s, %d params, %d seeds%s", cfg$trunk,
                as.integer(run$params[[term]]), n_seeds,
                if (ins_head) sprintf(" | insertion head: %s",
                                      paste(as.character(cfg$ins_class_names), collapse = "/")) else ""))
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
