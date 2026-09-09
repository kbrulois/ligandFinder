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
## Config defaults implement the intended architecture. Pass `r_exact = TRUE` to
## reproduce the pre-port R model exactly, including the two places the old code
## diverged from its own design: the DB position mask was overwritten with the CT
## mask (so DB was allowed only at the C terminus), and the per-epoch resample
## callback rebound a variable `fit` no longer read (so the negative sample and
## the augmentation noise were drawn once). See inst/python/README.md.
dcnn <- lf_dcnn_run(nn_input, all_params3, seed = 42L, verbose = 1L)

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

raw_pred_comb <- dcnn$pred_raw  # uncalibrated global score
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
## An amidated peptide is cut at a dibasic site with a glycine immediately 5' of
## it: ...X-G | K/R-K/R. The G is the amide donor, so a db window carrying one is
## a candidate amidation site.
##
## The anchor sits at a FIXED position in every db window, which is what makes
## this a lookup rather than a search: 9.2 sets window_origin = db_ind and then
## wN = db_ind + win_size[[t]]$start, so db_ind always lands at local position
## 1 - start (31 for C windows, 6 for N). Clamped windows are padded back out to
## seq_len at the front, so the offset holds there too. Two asymmetries matter:
##   * db_ind is the SECOND basic residue for C-target windows (the C branch of
##     9.2 adds a full lookahead offset) but the FIRST for N-target ones.
##   * BOTH termini are eligible. The motif is a property of the dibasic SITE,
##     not of the peptide you approach it from: in a polyprotein precursor one
##     dibasic pair is simultaneously the C-terminal cut of the peptide before it
##     and the start of the peptide after it, so a G sitting 5' of that pair is a
##     real amide donor regardless of which direction the window was anchored
##     from. An N-anchored window therefore reads its G at local position 5, a
##     C-anchored one at 29.
amid_targets <- c("N", "loop_N", "C", "loop_C")

## (braced: at top level R parses `if (...) x` and a following `else` as two
## statements, so the else must not start its own line)
.ws_start <- if (exists("win_size")) {
  c(N = win_size$N$start, C = win_size$C$start)
} else {
  c(N = -5L, C = -30L)                                             # 9.2 defaults
}

## local positions of the glycine and the two basic residues, per target
.motif_pos <- function(tg) {
  if (tg %in% c("C", "loop_C")) {
    a <- 1L - .ws_start[["C"]]                  # db_ind = 2nd basic
    c(g = a - 2L, b1 = a - 1L, b2 = a)
  } else {
    a <- 1L - .ws_start[["N"]]                  # db_ind = 1st basic
    c(g = a - 1L, b1 = a, b2 = a + 1L)
  }
}

.aa_at <- function(md, i) {
  aa <- as.character(md[["AA"]])
  if (i < 1L || i > length(aa)) NA_character_ else aa[[i]]
}

## Sanity check FIRST: if the anchor offset were wrong, every window would
## quietly come back FALSE and look like "no amidation motifs found". Confirm the
## two anchor positions really are basic residues before trusting the G test.
.db_i <- which(nn_input_comb$win_type == "db")
.chk  <- vapply(.db_i, function(i) {
  p <- .motif_pos(as.character(nn_input_comb$target[[i]]))
  md <- nn_input_comb$meta_data[[i]]
  isTRUE(.aa_at(md, p[["b1"]]) %in% c("K", "R") &&
         .aa_at(md, p[["b2"]]) %in% c("K", "R"))
}, logical(1))
message(sprintf("amidation: dibasic anchor confirmed at the expected offset in %d/%d db windows (%.1f%%)",
                sum(.chk), length(.chk), 100 * mean(.chk)))
if (mean(.chk) < 0.9)
  warning("amidation: the dibasic anchor is often NOT at the expected local position -- ",
          "check win_size against 9.2 before using the `amidation` column", immediate. = TRUE)

nn_input_comb$amidation <- vapply(seq_len(nrow(nn_input_comb)), function(i) {
  tg <- as.character(nn_input_comb$target[[i]])
  if (!identical(as.character(nn_input_comb$win_type[[i]]), "db")) return(FALSE)
  if (!tg %in% amid_targets) return(FALSE)
  p  <- .motif_pos(tg)
  md <- nn_input_comb$meta_data[[i]]
  isTRUE(identical(.aa_at(md, p[["g"]]), "G") &&
         .aa_at(md, p[["b1"]]) %in% c("K", "R") &&
         .aa_at(md, p[["b2"]]) %in% c("K", "R"))
}, logical(1))

message(sprintf("amidation: %d windows carry the G|dibasic motif (%d of them known peptides)",
                sum(nn_input_comb$amidation),
                sum(nn_input_comb$amidation & nn_input_comb$known == 1)))
rm(.ws_start, .motif_pos, .aa_at, .db_i, .chk)

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
## Windows are built as anchor + [-5,30] (N) / [-30,5] (C), so the putative peptide
## boundary sits at a fixed anchor residue: wN+5 for N windows, wC-5 for C windows.
## A window "hits" a uniprot peptide only if that peptide's matching terminus --
## start (N-terminus) for an N window, end (C-terminus) for a C window -- lands at
## the window's anchor residue (+/- anchor_tol), i.e. at the correct position, not
## merely somewhere inside the window span.
anchor_tol <- 2L

## per-model terminus (N/C); derive here so this block is self-contained even if
## the ranking mutate above hasn't been run on this nn_input_comb
if (!"model" %in% names(nn_input_comb))
  nn_input_comb$model <- stringr::str_remove(as.character(nn_input_comb$target), "^loop_")

wm <- stringr::str_match(nn_input_comb$peps, "_w(\\d+)-(\\d+)$")
nn_input_comb$wN         <- as.integer(wm[, 2])
nn_input_comb$wC         <- as.integer(wm[, 3])
nn_input_comb$anchor_res <- ifelse(nn_input_comb$model == "N",
                                   nn_input_comb$wN + 5L,     # expected peptide N-terminus
                                   nn_input_comb$wC - 5L)     # expected peptide C-terminus

starts_by_gene <- split(as.integer(uniprot_peps$start), uniprot_peps$gene)  # peptide N-ends
ends_by_gene   <- split(as.integer(uniprot_peps$end),   uniprot_peps$gene)  # peptide C-ends

anchor_hits <- function(gene, model, anchor, tol) {
  ter <- if (model == "N") starts_by_gene[[gene]] else ends_by_gene[[gene]]
  if (is.null(ter) || is.na(anchor)) return(FALSE)
  any(abs(ter - anchor) <= tol)
}

nn_input_comb$pep_terminus_hit <- FALSE
idx <- which(nn_input_comb$gene %in% names(starts_by_gene))   # only genes with uniprot peptides
if (length(idx) > 0) {
  nn_input_comb$pep_terminus_hit[idx] <- mapply(
    anchor_hits,
    nn_input_comb$gene[idx], nn_input_comb$model[idx], nn_input_comb$anchor_res[idx],
    MoreArgs = list(tol = anchor_tol))
}

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














