## ---- 1D-CNN window scoring --------------------------------------------------
## The modeling -- architecture, losses, oversampling, training, the pooled Platt
## calibration and the embeddings -- lives in the Python package
## `inst/python/lf_dcnn`. `R/dcnn_bridge.R` is the whole boundary: the arrays
## pass in memory through reticulate, so nothing touches disk.
##
## Everything below the bridge call stays in R: the yardstick QC, nearest-known
## retrieval, the amidation motif and all the plotting.
##
## After touching either side, run the checks:
##   cd inst/python && PYTHONPATH=. python -m lf_dcnn selftest
##   Rscript inst/python/tests/roundtrip.R
library(yardstick)

if (!exists("lf_dcnn_run")) source("R/dcnn_bridge.R")

class_cols <- setNames(c("#FED439FF", "#370335FF", "#8A9197FF", "#D2AF81FF",
                         "#D5E4A2FF", "#197EC0FF", "grey85", "#075149FF"),
                       c("CT_cleavage_context", "DB", "gap", "NT_cleavage_context", "pep_other", "pep_pocket", "padding", "none")
)

## ---- train both terminus models and score every window ----------------------
## Config defaults are the production architecture: the U-Net trunk
## (unet@16-32-64, ~21.4k params) with NO position-ramp input, chosen 2026-09-18
## on inst/scripts/10_5_benchmark_window_model.R. Pass `r_exact = TRUE` to
## reproduce the pre-port R model exactly -- flat trunk, ramp on, and the two
## places the old code diverged from its own design: the DB position mask was
## overwritten with the CT mask (so DB was allowed only at the C terminus), and
## the per-epoch resample callback rebound a variable `fit` no longer read (so
## the negative sample and the augmentation noise were drawn once).
## Training runs one python process per (member, terminus) model (`isolated`,
## the default): the second model trained in one Keras/TF process dies
## intermittently, and an ensemble trains ten. See inst/python/README.md.
## n_seeds = 5: five members per terminus. pred_raw / per_index come back as the
## ensemble MEAN, with pred_sd / per_index_sd giving the spread across members --
## a single softmax cannot tell "confidently 0.5" from "the members disagree",
## and those mean opposite things when reading a per-residue profile.
## cache = : the 10 members are trained ONCE and reloaded (weights included) on
## every later run. Delete the file, or pass refresh = TRUE, to retrain.
dcnn <- lf_dcnn_run(nn_input, all_params3, seed = 42L, verbose = 1L, n_seeds = 5L,
                    cache = "~/AF2_analysis/lf_dcnn_run_cache.rds")

## Session contract for the downstream scripts (10_2_per_ind_profiles.R,
## 10_2_model_importance.R, 10_3d_embed_umap.R): they expect these names in the
## session. Read them off the Config so the R and Python sides cannot drift --
## the 1-based/0-based translation happens only here and in dcnn_bridge.R.
dcnn_cfg    <- dcnn$config
models      <- dcnn$models      # keras Models; keras3's R generics dispatch on them
nn_in_all   <- dcnn$arrays      # list[term][split] of the contract arrays
seq_len     <- as.integer(dcnn_cfg$seq_len)
n_channels  <- as.integer(dcnn_cfg$n_channels)
classes     <- setNames(seq_along(dcnn_cfg$class_names) - 1L,
                        as.character(dcnn_cfg$class_names))
real_cols   <- as.integer(dcnn_cfg$real_cols) + 1L            # python is 0-based
pi_names    <- as.character(dcnn_cfg$pi_names)                 # [6 real, none, padding]

raw_pred_comb <- dcnn$pred_raw  # uncalibrated global score (ensemble mean)
val_pred_comb <- dcnn$pred      # pooled-Platt calibrated, comparable across models

message(sprintf("lf_dcnn: %d windows scored; pooled calibrator fit on %d val windows (%d positive)",
                length(val_pred_comb), dcnn$calibrator$n_obs, dcnn$calibrator$n_pos))


metrics <- metric_set(
  roc_auc,
  pr_auc,
  accuracy,
  mcc,
  f_meas,
  precision,
  recall,
)

df <- tibble(
  truth = factor(do.call(c, lapply(nn_input, function(x) x$all$known))),
  .pred_class = ifelse(val_pred_comb > 0.6, 1, 0) %>% factor(., levels = c(0,1)),
  .pred_1 = val_pred_comb
)


qc_mets <- metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second")
qc_mets

nn_input_comb <- bind_rows(lapply(nn_input, function(x) x$all))

nn_input_comb$pred     <- val_pred_comb   # calibrated (pooled Platt, common scale)
nn_input_comb$pred_raw <- raw_pred_comb   # uncalibrated global score, for comparison

# Per-index softmax, one tidy tibble per window: 8 cols [6 real, none, padding].
# This head IS supervised (see the port README) -- padding is a trained class and
# `none` is masked out of the loss.
nn_input_comb$per_index <- dcnn$per_index_tbl

## Ensemble spread, carried alongside rather than folded into pred/per_index so
## nothing downstream re-ranks. per_index_sd matches per_index column for column.
nn_input_comb$pred_sd       <- dcnn$pred_sd
nn_input_comb$per_index_sd  <- dcnn$per_index_sd_tbl
message(sprintf("ensemble: %d members/terminus; median score sd %.4f, mean per-residue sd %.4f",
                dcnn$n_seeds, median(dcnn$pred_sd),
                mean(vapply(dcnn$per_index_sd_tbl, \(d) mean(as.matrix(d[-1])), numeric(1)))))

# --- nearest known-peptide retrieval (reuse the trained "embed" layer) ------
# The 16-d penultimate representation, L2-normalised on the Python side: for each
# window find the most similar KNOWN peptide by cosine similarity. No new loss --
# the embedding is whatever the global classifier already learned to sit on.
emb_all <- dcnn$emb                                          # rows aligned with nn_input_comb

# reference bank = known peptides only
ref_i     <- which(nn_input_comb$known == 1)
ref_emb   <- emb_all[ref_i, , drop = FALSE]
ref_names <- nn_input_comb$peps[ref_i]

sim_mat <- emb_all %*% t(ref_emb)                            # (n_windows, n_ref) cosine sims
sim_mat[outer(nn_input_comb$peps, ref_names, `==`)] <- -Inf # never match a window to itself
best <- max.col(sim_mat, ties.method = "first")

nn_input_comb$nn_closest_peptide <- ref_names[best]
nn_input_comb$nn_closest_sim     <- sim_mat[cbind(seq_along(best), best)]

nn_input_comb <- nn_input_comb %>%
  # per-model (terminus) category, and its combination with win_type
  mutate(model    = stringr::str_remove(as.character(target), "^loop_"),  # "N" or "C"
         category = paste0(model, "_", win_type)) %>%                     # e.g. "N_db", "C_chym", "C_pep_end"
  arrange(desc(pred)) %>%
  #filter(win_type == "db") %>%
  mutate(rank = row_number(), .by = known) %>%                            # existing global rank (within known/unknown)
  # rank the score WITHIN each combined (model x win_type) category, kept
  # separate for known vs candidate. Drop `known` from .by to rank the two together.
  mutate(rank_cat = row_number(), .by = c(category, known))

## ---- amidation-motif windows -----------------------------------------------
## ...X-G | K/R-K/R at the dibasic anchor: the G is the amide donor. The rule
## (and why the anchor is a fixed local position, and why both termini count)
## lives in R/amidation_motif.R so the benchmark CSVs flag the same windows.
if (!exists("lf_amidation_motif")) source("R/amidation_motif.R")
## (braced: at top level R parses `if (...) x` and a following `else` as two
## statements, so the else must not start its own line)
.ws_start <- if (exists("win_size")) {
  c(N = win_size$N$start, C = win_size$C$start)
} else {
  c(N = -5L, C = -30L)                                             # 9.2 defaults
}
nn_input_comb <- lf_amidation_motif(nn_input_comb, win_start = .ws_start)
rm(.ws_start)

uniprot_peps <- data.table::fread("~/Desktop/Peptides/uniprot_peptides.csv") %>% as_tibble()

nn_input_comb <- nn_input_comb %>%
  mutate(gene = stringr::str_extract(peps, "^[^_]+")) %>%
  {left_join(., secretome %>% distinct(gene, .keep_all = T) %>% select(gene, location), by = "gene")} %>%
  mutate(uni_pep = if_else(gene %in% uniprot_peps$gene, 1, 0))

## peptide length from the pep_id "<start>x<end>" range. Vectorized (one str_match
## over the whole column) instead of a per-row map_int -- the latter ran the regex
## ~40K times and returned NA for every non-known window. NA where no range (controls).
nn_input_comb <- nn_input_comb %>%
  mutate(length = {
    rng <- stringr::str_match(pep_id, "(\\d+)x(\\d+)")
    as.integer(rng[, 3]) - as.integer(rng[, 2])
  })

## --- annotate windows that ANCHOR a uniprot-peptide terminus -----------------
## The rule lives in R/uniprot_terminus_hits.R so the benchmark scripts
## (10_5_benchmark_window_model.R) score "top hits" identically: a window hits
## only if the peptide's matching terminus lands at the window's anchor residue
## (wN+5 for N windows, wC-5 for C), +/- anchor_tol -- not merely inside the span.
if (!exists("lf_uniprot_terminus_hits")) source("R/uniprot_terminus_hits.R")
anchor_tol <- 2L
nn_input_comb <- lf_uniprot_terminus_hits(nn_input_comb, uniprot_peps, anchor_tol = anchor_tol)

message(sprintf("uniprot-terminus hits: %d windows (N: %d, C: %d)",
                sum(nn_input_comb$pep_terminus_hit),
                sum(nn_input_comb$pep_terminus_hit & nn_input_comb$model == "N"),
                sum(nn_input_comb$pep_terminus_hit & nn_input_comb$model == "C")))

## --- violin: score of hits vs non-hits, segregated by window terminus (N/C) ---
## swap `y = pred` for `y = rank_cat` (or `rank`) to plot rank instead of score.
pep_hit_violin <- nn_input_comb %>%
  mutate(hit = factor(if_else(pep_terminus_hit, "anchors uniprot pep terminus", "no"),
                      levels = c("no", "anchors uniprot pep terminus"))) %>%
  ggplot(aes(x = hit, y = pred, fill = hit)) +
  geom_violin(scale = "width", alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(data = ~ dplyr::filter(.x, pep_terminus_hit),
              width = 0.15, size = 0.9, alpha = 0.7) +
  facet_grid(rows = vars(win_type), cols = vars(model)) +
  labs(x = NULL, y = "score (pred)",
       title = "Windows anchoring a uniprot-peptide terminus vs not, by model x win_type") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("~/AF2_analysis/uniprot_terminus_score_violin.svg", pep_hit_violin, width = 9, height = 5)

nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(peps %in% c("ANO8_w12-47", "ANO8_w36-71", "ANO8_w14-49")) %>%
  View()

nn_input_comb %>% filter(known == 0) %>%
  filter(gene == "ANO8") %>%
  View()


nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(grepl("GDF", gene)) %>%
  View()

nn_input_comb %>% filter(known == 0) %>%
  filter(grepl("ASIP", gene)) %>%
  View()


View(nn_input_comb %>%
       filter(known == 0 & location %in% c("2t", "3t", "2l", "3l", "4l")) %>%
       filter(win_type == "chym") %>%
       filter(target %in% c("C", "C_loop")))

View(nn_input_comb)

nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(grepl("BRINP3", gene) & target == "C") %>%
  View()


## ---- per-gene protein plots -------------------------------------------------
## Moved to its own script so the html files can be regenerated without
## retraining. Needs nn_input_comb (above) plus secretome, secretome_aa and
## peps_tp in the session.
##
##   source("inst/scripts/10_4_plot_proteins.R")
##
## and, once that has written the bundle, a fresh session can re-render with no
## keras and no secretome_aa:  Rscript inst/scripts/plot_only.R


p <- yardstick::roc_curve(data = df, truth = "truth", ".pred_1",
                          event_level = "second") %>%
  autoplot() +
  ggtitle("Mixed training set; Mixed test",
          subtitle = paste(paste0("roc auc: ", round(qc_mets %>% filter(.metric == "roc_auc") %>% pull(.estimate), 2)),
                           paste0("pr auc: ", round(qc_mets %>% filter(.metric == "pr_auc") %>% pull(.estimate), 2)), collapse = "\n"))


ggsave("~/AF2_analysis/all_peps_roc_AUC_mixed-train-mixed_test.svg", p)

model_stats <- metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second")

write.csv(model_stats, "~/AF2_analysis/model_stats_new.csv")




predict_combined <- function(x, category) {
  ifelse(
    category == 1,
    predict(model_A, x),
    predict(model_B, x)
  )
}














