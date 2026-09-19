#!/usr/bin/env Rscript
## ---- window-model benchmark ----------------------------------------------------
## Any set of lf_dcnn variants -- Config overrides and/or extra input channels --
## trained under ONE recipe (same windows, seeds, epochs), compared side by side.
##
## A `--preset` names the arms:
##   position_ramp   the U-Net with vs without the in-graph position channel
##                   (what set the 2026-09-18 defaults)
##   t5              C terminus only: the U-Net on the 26 hand-built channels vs
##                   the same plus ProtT5 embedding channels (PCA-reduced by
##                   `python -m lf_plm reduce`, attached by R/plm_features.R)
##
## The page compares, per arm:
##   * validation ROC / PR curves, ensemble curve bold and every member thin,
##     with AUC quoted as mean +/- sd ACROSS SEEDS (the AUC of the averaged
##     score is a different, usually higher, number; the table has both), and
##     again over the top-3 members only so one collapsed seed cannot own the sd
##   * the "top hits" violin from 10_1dcnn_new6.R -- score of windows anchoring a
##     UniProt-peptide terminus vs everything else -- faceted by arm
##   * how much the candidate rankings agree: top-k overlap, rank scatter,
##     and the top windows of each arm with their rank under the others
##
## Run from the repo root, in a fresh session:
##
##   /usr/local/bin/Rscript inst/scripts/10_5_benchmark_window_model.R --preset position_ramp
##   /usr/local/bin/Rscript inst/scripts/10_5_benchmark_window_model.R --preset t5
##   ... --refresh --n-seeds 5 --terms C
##
## or source() it after 9.2 has put nn_input + all_params3 in the session.
##
## Training goes through the STANDALONE path -- lf_dcnn_export() once per input
## variant, then `python -m lf_dcnn train --isolated` per arm (one python
## process per MODEL, then a combine step), then lf_dcnn_import(). Not
## in-process like 10_1dcnn_new6.R used to: the second model trained in one
## Keras/TF process dies intermittently in a retraced tf.function. Models land
## under <cache-dir>/<arm>/members/, so an interrupted arm resumes and
## re-running only the plotting is seconds; --refresh retrains everything.
##
## Output: --out (default per preset), plus <stem>_metrics.csv,
## <stem>_overlap.csv, <stem>_windows.csv and svg/png copies of every panel.
## ------------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(yardstick)
  library(ggiraph); library(htmltools)
})

## ---- options ---------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
.flag <- function(flag) flag %in% .args

preset      <- .opt("--preset", "position_ramp")
n_seeds     <- as.integer(.opt("--n-seeds", "5"))
base_seed   <- as.integer(.opt("--seed", "42"))
nn_cache    <- path.expand(.opt("--nn-input", "~/AF2_analysis/lf_dcnn_compare_nn_input.rds"))
uniprot_csv <- path.expand(.opt("--uniprot-peps", "~/Desktop/Peptides/uniprot_peptides.csv"))
plm_parquet <- path.expand(.opt("--plm", "~/AF2_analysis/lf_plm/prot_t5_pca32.parquet"))
seq_parquet <- path.expand(.opt("--sequences", "~/AF2_analysis/lf_plm/sequences.parquet"))
refresh     <- .flag("--refresh")
verbose     <- as.integer(.opt("--verbose", "1"))   # 1 = one line per model in <cache-dir>/<arm>.log
epochs      <- .opt("--epochs", NA)     # smoke-test knob; NA = Config default (2000)
anchor_tol  <- 2L                       # as in 10_1dcnn_new6.R
top_k       <- c(10L, 25L, 50L, 100L, 250L, 500L)
n_top_table <- 25L                      # windows listed per approach x terminus
top_n_members <- 3L                     # "top3": the best members only, outliers dropped

## ---- presets -------------------------------------------------------------------
## An arm: `cfg` are Config overrides passed to the CLI as --set FIELD=JSON
## (spelled out even where they equal the defaults, so an arm means the same
## thing whatever the defaults become); `input` names the channel set it
## trains on (see `input_builders` below); `dir` is its cache directory under
## the preset's cache_dir (fixed, so a renamed arm keeps its trained models);
## `desc` is for the page; `short` suffixes the arm's columns in the wide
## windows CSV (pred_<short>, rank_cand_<short>, ...). The first arm is the
## reference the rank-agreement panels compare the others against.
## The U-Net is the compare_r_python.R harness's "unet@16-32-64": trunk "unet",
## filters 16-32-64 (two pooling levels, 36 -> 18 -> 9), dropout 0.2 per level,
## ~21.4k parameters; the production model since 2026-09-18.
unet <- list(trunk = "unet", unet_filters = c(16L, 32L, 64L), unet_dropout = c(0.2, 0.2, 0.2),
             position_ramp = FALSE)
presets <- list(
  position_ramp = list(
    title     = "U-Net window model \u2014 does it need the position input?",
    cache_dir = "~/AF2_analysis/lf_dcnn_bench",
    out       = "~/AF2_analysis/ligandFinder_v6_benchmark.html",
    terms     = NULL,                                   # all
    arms = list(
      "unet, position ramp" = list(
        cfg = c(unet[names(unet) != "position_ramp"], list(position_ramp = TRUE)),
        input = "base", dir = "unet_16-32-64", short = "ramp",
        desc = "the pooling U-Net trunk (16 &rarr; 32 &rarr; 64) with the in-graph [0,1] position channel appended to the input"),
      "unet, no position ramp" = list(
        cfg = unet, input = "base", dir = "unet_16-32-64_no_ramp", short = "noramp",
        desc = "the same U-Net trunk on the raw channels only")),
    cols = c("unet, position ramp" = "#7570B3", "unet, no position ramp" = "#E7298A")),
  t5 = list(
    title     = "U-Net window model \u2014 do ProtT5 embedding channels help? (C terminus)",
    cache_dir = "~/AF2_analysis/lf_dcnn_bench_t5",
    out       = "~/AF2_analysis/ligandFinder_v7_benchmark_t5.html",
    terms     = "C",
    arms = list(
      "unet, 26 channels" = list(
        cfg = unet, input = "base", dir = "unet_base", short = "unetonly",
        desc = "the production U-Net (no position ramp) on the 26 hand-built per-residue channels"),
      "unet, 26 + ProtT5" = list(
        cfg = unet, input = "t5", dir = "unet_t5", short = "t5",
        desc = "the same U-Net with ProtT5-XL-U50 per-residue embeddings appended, PCA-reduced to 32 channels fit over the whole secretome sample and scaled to [0,1] (inst/python/lf_plm)")),
    cols = c("unet, 26 channels" = "#E7298A", "unet, 26 + ProtT5" = "#1B9E77"))
)
if (!preset %in% names(presets))
  stop("--preset must be one of: ", paste(names(presets), collapse = ", "))
P <- presets[[preset]]
approaches    <- P$arms
approach_cols <- P$cols
cache_dir     <- path.expand(.opt("--cache-dir", P$cache_dir))
html_path     <- path.expand(.opt("--out", P$out))
terms_keep    <- { t <- .opt("--terms", ""); if (nzchar(t)) strsplit(t, ",")[[1]] else P$terms }
stopifnot(all(names(approaches) %in% names(approach_cols)),
          !anyDuplicated(vapply(approaches, `[[`, "", "dir")))

## ---- repo root + bridge -----------------------------------------------------
## The installed ligandFinder package ALSO ships inst/python, and lf_dcnn_path()
## prefers it. Pin the python package to this checkout so the arm being
## benchmarked is the code in the working tree, not whatever was last installed.
.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
if (!exists("lf_dcnn_run"))              source(file.path(ROOT, "R", "dcnn_bridge.R"))
if (!exists("lf_uniprot_terminus_hits")) source(file.path(ROOT, "R", "uniprot_terminus_hits.R"))
if (!exists("lf_plm_attach"))            source(file.path(ROOT, "R", "plm_features.R"))
invisible(lf_dcnn_python(path = file.path(ROOT, "inst", "python")))

## ---- the window data ----------------------------------------------------------
if (!exists("nn_input") || !exists("all_params3")) {
  if (file.exists(nn_cache)) {
    message("loading cached nn_input from ", nn_cache)
    .cc <- readRDS(nn_cache); nn_input <- .cc$nn_input; all_params3 <- .cc$all_params3
  } else {
    message("building nn_input via 9.2_add_contact_data.R (needs ~/AF2_analysis/nn_dat_cache.rds) ...")
    source(file.path(ROOT, "inst/scripts/putative_peptide/generate_residue_db_untested/9.2_add_contact_data.R"))
    saveRDS(list(nn_input = nn_input, all_params3 = all_params3), nn_cache)
  }
}
if (!is.null(terms_keep)) {
  missing <- setdiff(terms_keep, names(nn_input))
  if (length(missing)) stop("--terms: not in nn_input: ", paste(missing, collapse = ", "))
  nn_input <- nn_input[terms_keep]
  message("terms restricted to: ", paste(names(nn_input), collapse = ", "))
}
for (t in names(nn_input)) for (s in c("train", "val", "all"))
  message(sprintf("  %s/%-5s n=%5d  positives=%3d", t, s,
                  nrow(nn_input[[t]][[s]]), sum(nn_input[[t]][[s]]$known)))

## ---- input variants: which channels each arm trains on ------------------------
## Each builder returns list(nn_input, channels) with the SAME windows in the
## same order (only columns differ), so every arm scores identical rows.
input_builders <- list(
  base = function() list(nn_input = nn_input, channels = all_params3),
  t5   = function() {
    for (f in c(plm_parquet, seq_parquet)) if (!file.exists(f))
      stop("input 't5' needs ", f, " -- run `python -m lf_plm embed` then `reduce` (see inst/python/lf_plm)")
    plm  <- lf_plm_read(plm_parquet)
    seqs <- lf_read_parquet(seq_parquet)        # not arrow:: -- see lf_read_parquet
    lf_plm_attach(nn_input, plm, seqs, all_params3)
  }
)

uniprot_peps <- data.table::fread(uniprot_csv) %>% as_tibble()

## ---- train / reload the arms -------------------------------------------------
## One arrays.npz per input variant (arms sharing channels share it), one
## `--isolated` CLI run per arm, outputs read back with lf_dcnn_import().
## `params` and `n_seeds` come from history.json, written by the CLI.
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
in_dirs <- list()
for (inp in unique(vapply(approaches, `[[`, "", "input"))) {
  in_dir <- file.path(cache_dir, if (inp == "base") "in" else paste0("in_", inp))
  in_dirs[[inp]] <- in_dir
  if (refresh || !file.exists(file.path(in_dir, "arrays.npz"))) {
    built <- input_builders[[inp]]()
    stopifnot(identical(lapply(built$nn_input, function(x) x$all$peps),
                        lapply(nn_input, function(x) x$all$peps)))
    meta <- bind_rows(lapply(built$nn_input, function(x) x$all)) %>%
      select(any_of(c("peps", "gene", "win_type", "target", "known")))
    lf_dcnn_export(built$nn_input, built$channels, in_dir, meta = meta, seed = base_seed)
    message(sprintf("exported arrays for input '%s' (%d channels) to %s",
                    inp, length(built$channels), in_dir))
    rm(built)
  }
}

py <- file.path(dirname(reticulate::py_config()$python), "python")
py_path <- lf_dcnn_path(file.path(ROOT, "inst", "python"))

runs <- lapply(names(approaches), function(nm) {
  out_dir <- file.path(cache_dir, approaches[[nm]]$dir)
  done    <- file.path(out_dir, "outputs.npz")
  if (refresh) unlink(out_dir, recursive = TRUE)     # members too, or they would be reused
  in_dir <- in_dirs[[approaches[[nm]]$input]]
  if (refresh || !file.exists(done)) {
    message(sprintf("\n== %s  (n_seeds = %d, seed = %d): training in a subprocess ==",
                    nm, n_seeds, base_seed))
    t0 <- Sys.time()
    log_f <- file.path(cache_dir, paste0(approaches[[nm]]$dir, ".log"))
    message("   progress: tail -f ", log_f,
            "\n   curves:   tensorboard --logdir ", file.path(cache_dir, "tb"))
    ov <- approaches[[nm]]$cfg
    sets <- vapply(names(ov), function(f)
      sprintf("%s=%s", f, jsonlite::toJSON(ov[[f]], auto_unbox = TRUE)), character(1))
    argv <- c("-u", "-m", "lf_dcnn", "train", "--isolated",
              "--input-dir", in_dir, "--output-dir", out_dir,
              "--n-seeds", n_seeds, "--seed", base_seed, "--verbose", verbose,
              as.vector(rbind("--set", sets)),
              "--set", sprintf("tensorboard_dir=\"%s\"", file.path(cache_dir, "tb", approaches[[nm]]$dir)))
    if (!is.na(epochs)) argv <- c(argv, "--epochs", epochs)
    rc <- system2(py, argv, env = paste0("PYTHONPATH=", py_path), stdout = log_f, stderr = log_f)
    if (rc != 0 || !file.exists(done))
      stop(sprintf("%s: `python -m lf_dcnn train` failed (rc=%d); see %s", nm, rc, log_f))
    message(sprintf("   %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else {
    message(sprintf("== %s: reusing %s (no training) ==", nm, done))
  }
  r <- lf_dcnn_import(out_dir)
  r$params <- unlist(r$params)[names(nn_input)]
  r$n_channels <- as.integer(jsonlite::read_json(file.path(in_dir, "config.json"))$n_channels)
  ## sd across members, as run() computes it (sample sd); derive it when an
  ## outputs.npz predates the field
  if (is.null(r$pred_sd)) r$pred_sd <- apply(r$pred_raw_members, 2, sd)
  stopifnot(identical(as.character(unlist(r$term_order)), names(nn_input)),
            as.integer(r$n_seeds) == n_seeds)
  message(sprintf("   %s: %s params/terminus, %d input channels", nm,
                  paste(r$params, collapse = "/"), r$n_channels))
  r
})
names(runs) <- names(approaches)

## ---- validation scores: one long table ---------------------------------------
## Row order of val_scores / val_labels is names(nn_input) concatenated, each
## terminus in its own nn_input[[t]]$val order -- the same order lf_dcnn_arrays()
## builds the arrays in.
val_term <- rep(names(nn_input), times = vapply(nn_input, function(x) nrow(x$val), integer(1)))
val_peps <- unlist(lapply(nn_input, function(x) x$val$peps), use.names = FALSE)

val_long <- bind_rows(lapply(names(runs), function(nm) {
  r <- runs[[nm]]
  members <- r$val_scores_members                          # (n_seeds, m)
  stopifnot(ncol(members) == length(val_term))
  bind_rows(
    tibble(approach = nm, member = "ensemble", seed = NA_integer_,
           term = val_term, peps = val_peps, truth = as.integer(r$val_labels),
           score = as.numeric(r$val_scores)),
    bind_rows(lapply(seq_len(nrow(members)), function(k)
      tibble(approach = nm, member = paste0("seed ", base_seed + k - 1L),
             seed = base_seed + k - 1L,
             term = val_term, peps = val_peps, truth = as.integer(r$val_labels),
             score = as.numeric(members[k, ]))))
  )
})) %>%
  mutate(approach = factor(approach, levels = names(approaches)),
         truth_f  = factor(truth, levels = c(0L, 1L)))

## "pooled" = both termini together. The calibrator is ONE logistic fit on the
## pooled raw scores, i.e. monotone, so pooled AUC on raw == on calibrated.
val_long <- bind_rows(val_long, val_long %>% mutate(term = "pooled")) %>%
  mutate(term = factor(term, levels = c(names(nn_input), "pooled")))

## ---- AUCs: ensemble, and mean +/- sd across members ---------------------------
auc_by <- val_long %>%
  group_by(approach, term, member, seed) %>%
  summarise(roc_auc = roc_auc_vec(truth_f, score, event_level = "second"),
            pr_auc  = pr_auc_vec(truth_f, score, event_level = "second"),
            n = n(), n_pos = sum(truth), .groups = "drop")

## top-N: the N best members by that metric, so one bad seed does not dominate
## the sd. Biased upward by construction -- read it as "what the method does
## when it works", next to the all-member columns.
top_mean <- function(x, n = top_n_members) mean(sort(x, decreasing = TRUE)[seq_len(min(n, length(x)))])
top_sd   <- function(x, n = top_n_members) sd(sort(x, decreasing = TRUE)[seq_len(min(n, length(x)))])

auc_summary <- auc_by %>%
  group_by(approach, term) %>%
  summarise(
    n_val = first(n), n_pos = first(n_pos),
    roc_auc_ensemble = roc_auc[member == "ensemble"],
    roc_auc_mean = mean(roc_auc[member != "ensemble"]),
    roc_auc_sd   = sd(roc_auc[member != "ensemble"]),
    roc_auc_top3_mean = top_mean(roc_auc[member != "ensemble"]),
    roc_auc_top3_sd   = top_sd(roc_auc[member != "ensemble"]),
    pr_auc_ensemble = pr_auc[member == "ensemble"],
    pr_auc_mean = mean(pr_auc[member != "ensemble"]),
    pr_auc_sd   = sd(pr_auc[member != "ensemble"]),
    pr_auc_top3_mean = top_mean(pr_auc[member != "ensemble"]),
    pr_auc_top3_sd   = top_sd(pr_auc[member != "ensemble"]),
    .groups = "drop") %>%
  left_join(tibble(approach = factor(names(runs), levels = names(approaches)),
                   params = vapply(runs, function(r) paste(r$params, collapse = " / "), character(1))),
            by = "approach") %>%
  arrange(term, approach)

message(sprintf("\nvalidation AUC: all members mean+/-sd | top-%d members | [ensemble of the averaged score]",
                top_n_members))
print(auc_summary %>%
        transmute(term, approach, params, n_val, n_pos,
                  roc_auc = sprintf("%.3f+/-%.3f | %.3f+/-%.3f | [%.3f]", roc_auc_mean, roc_auc_sd,
                                    roc_auc_top3_mean, roc_auc_top3_sd, roc_auc_ensemble),
                  pr_auc  = sprintf("%.3f+/-%.3f | %.3f+/-%.3f | [%.3f]", pr_auc_mean,  pr_auc_sd,
                                    pr_auc_top3_mean, pr_auc_top3_sd, pr_auc_ensemble)) %>%
        as.data.frame(), row.names = FALSE)

## ---- ROC / PR curves -------------------------------------------------------------
curve_by <- function(fun) {
  val_long %>%
    group_by(approach, term, member) %>%
    group_modify(~ fun(.x, truth_f, score, event_level = "second")) %>%
    ungroup() %>%
    mutate(weight = if_else(member == "ensemble", "ensemble", "member"))
}
roc_pts <- curve_by(roc_curve)
pr_pts  <- curve_by(pr_curve) %>% filter(is.finite(.threshold) | recall > 0)

## legend labels carry the numbers so the curve panel reads on its own
auc_lab <- auc_summary %>%
  mutate(roc_lab = sprintf("%s  ROC %.2f ± %.2f", approach, roc_auc_mean, roc_auc_sd),
         pr_lab  = sprintf("%s  PR %.2f ± %.2f",  approach, pr_auc_mean,  pr_auc_sd))

curve_theme <- list(
  facet_wrap(vars(term), nrow = 1),
  scale_colour_manual(values = approach_cols, name = NULL),
  scale_linewidth_manual(values = c(ensemble = 1.1, member = 0.35), guide = "none"),
  scale_alpha_manual(values = c(ensemble = 1, member = 0.35), guide = "none"),
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)),
  theme_bw(base_size = 11),
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
)

p_roc <- ggplot(roc_pts, aes(x = 1 - specificity, y = sensitivity,
                             colour = approach, group = interaction(approach, member),
                             linewidth = weight, alpha = weight)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
  geom_path_interactive(aes(tooltip = paste0(approach, " — ", member),
                            data_id = paste(approach, member))) +
  labs(x = "1 - specificity (FPR)", y = "sensitivity (TPR)",
       title = "ROC, validation windows",
       subtitle = "bold = ensemble of the averaged score; thin = each seed") +
  curve_theme

p_pr <- ggplot(pr_pts, aes(x = recall, y = precision,
                           colour = approach, group = interaction(approach, member),
                           linewidth = weight, alpha = weight)) +
  geom_path_interactive(aes(tooltip = paste0(approach, " — ", member),
                            data_id = paste(approach, member))) +
  labs(x = "recall", y = "precision",
       title = "Precision-recall, validation windows",
       subtitle = sprintf("positives: %s",
                          paste(sprintf("%s=%d", levels(val_long$term),
                                        auc_summary$n_pos[match(levels(val_long$term), auc_summary$term)]),
                                collapse = "  "))) +
  curve_theme

## ---- every window, both approaches ---------------------------------------------
## The `all` split of both termini, in the row order lf_dcnn_run() scores it.
win_base <- bind_rows(lapply(nn_input, function(x) x$all)) %>%
  select(any_of(c("peps", "win_type", "gene", "target", "pep_id", "known"))) %>%
  mutate(model = sub("^loop_", "", as.character(target)))
if (!"gene" %in% names(win_base)) win_base$gene <- sub("_.*$", "", win_base$peps)

win_base <- lf_uniprot_terminus_hits(win_base, uniprot_peps, anchor_tol = anchor_tol)
message(sprintf("uniprot-terminus hits: %d windows (N: %d, C: %d)",
                sum(win_base$pep_terminus_hit),
                sum(win_base$pep_terminus_hit & win_base$model == "N"),
                sum(win_base$pep_terminus_hit & win_base$model == "C")))

## per window, the top-N members' mean and sd: a trimmed aggregate that drops
## the members that scored it lowest. For a window the members split on (say
## 0.88 / 0.10 / 0.51 / 0.66 / 0.98) the all-member mean sits at 0.63 with sd
## 0.35; the top-3 mean says what the members that "saw it" agree on.
top_members <- function(M, n = top_n_members) {
  n <- min(n, nrow(M))
  S <- matrix(apply(M, 2, function(v) sort(v, decreasing = TRUE)[seq_len(n)]), nrow = n)
  list(mean = colMeans(S), sd = if (n > 1) apply(S, 2, sd) else rep(NA_real_, ncol(S)))
}

windows <- bind_rows(lapply(names(runs), function(nm) {
  r <- runs[[nm]]
  stopifnot(length(r$pred) == nrow(win_base), ncol(r$pred_raw_members) == nrow(win_base))
  tm <- top_members(r$pred_raw_members)
  win_base %>%
    mutate(approach = nm,
           pred     = as.numeric(r$pred),        # pooled-Platt calibrated
           pred_raw = as.numeric(r$pred_raw),    # ensemble mean, uncalibrated
           pred_sd  = as.numeric(r$pred_sd),     # spread across members
           pred_raw_top3 = tm$mean,              # mean of the top-N member scores
           pred_sd_top3  = tm$sd)
})) %>%
  mutate(approach = factor(approach, levels = names(approaches))) %>%
  ## candidate rank: within terminus, AMONG unknown windows only (knowns are
  ## ranked separately and then dropped), best first -- the same convention as
  ## `rank` in 10_1dcnn_new6.R and the ensemble_global_predictions CSVs
  group_by(approach, model, known) %>%
  mutate(rank_cand      = as.integer(rank(-pred, ties.method = "first")),
         rank_cand_top3 = as.integer(rank(-pred_raw_top3, ties.method = "first"))) %>%
  ungroup() %>%
  mutate(rank_cand      = if_else(known == 0, rank_cand, NA_integer_),
         rank_cand_top3 = if_else(known == 0, rank_cand_top3, NA_integer_))

## ---- violin: uniprot-terminus hits vs the rest, faceted by approach ----------
## Same plot as 10_1dcnn_new6.R's pep_hit_violin, with the approach as the outer
## column facet so the two arms sit side by side per terminus.
viol_dat <- windows %>%
  mutate(hit = factor(if_else(pep_terminus_hit, "anchors uniprot\npep terminus", "no"),
                      levels = c("no", "anchors uniprot\npep terminus")))

## y = pred_raw, the ensemble-mean sigmoid the viewer plots as "prediction"
## (ligandFinder_v5 and the ensemble_global_predictions CSVs use it too). The
## pooled-Platt `pred` is so steep for these models that it pins almost every
## window to 0 and the violins collapse to lines.
p_violin <- ggplot(viol_dat, aes(x = hit, y = pred_raw, fill = hit)) +
  geom_violin(scale = "width", alpha = 0.5, quantiles = c(0.25, 0.5, 0.75),
              quantile.linetype = 1) +
  geom_jitter_interactive(
    data = ~ dplyr::filter(.x, pep_terminus_hit),
    aes(tooltip = sprintf("%s\n%s | pred %.3f (raw %.3f ± %.3f) | top%d raw %.3f ± %.3f\ncand. rank %s (top%d: %s)",
                          peps, approach, pred, pred_raw, pred_sd, top_n_members,
                          pred_raw_top3, pred_sd_top3,
                          ifelse(is.na(rank_cand), "known", rank_cand), top_n_members,
                          ifelse(is.na(rank_cand_top3), "known", rank_cand_top3)),
        data_id = peps),
    width = 0.15, size = 0.9, alpha = 0.7) +
  ## one row per arm, terminus x win_type across (nested strips when ggh4x is
  ## installed; plain facet_grid repeats the terminus label per win_type)
  (if (requireNamespace("ggh4x", quietly = TRUE))
     ggh4x::facet_nested(rows = vars(approach), cols = vars(model, win_type))
   else facet_grid(rows = vars(approach), cols = vars(model, win_type))) +
  scale_fill_manual(values = c("grey70", "#7570B3"), guide = "none") +
  labs(x = NULL, y = "prediction (pred_raw, ensemble mean)",
       title = "Windows anchoring a UniProt-peptide terminus vs not, by approach (rows) x terminus x win_type",
       subtitle = "hover a point for the window; the same window is highlighted in both approaches") +
  theme_bw(base_size = 11)

## ---- ranking agreement between the approaches ------------------------------------
## Wide table: one row per candidate window, one column block per arm.
arm_names <- names(approaches)
ref_arm   <- arm_names[[1]]
cand <- windows %>% filter(known == 0) %>%
  select(approach, model, peps, gene, win_type, pep_terminus_hit,
         pred, pred_sd, pred_raw, pred_raw_top3, pred_sd_top3, rank_cand, rank_cand_top3) %>%
  pivot_wider(id_cols = c(model, peps, gene, win_type, pep_terminus_hit),
              names_from = approach,
              values_from = c(pred, pred_sd, pred_raw, pred_raw_top3, pred_sd_top3,
                              rank_cand, rank_cand_top3),
              names_sep = "|")
col <- function(what, arm) paste0(what, "|", arm)

## every pair of arms: windows in the top k of BOTH, plus Spearman over all candidates
pairs <- combn(arm_names, 2, simplify = FALSE)
overlap <- bind_rows(lapply(pairs, function(pr) {
  bind_rows(lapply(split(cand, cand$model), function(d) {
    r1 <- d[[col("rank_cand", pr[[1]])]]; r2 <- d[[col("rank_cand", pr[[2]])]]
    tibble(arm_a = pr[[1]], arm_b = pr[[2]], model = d$model[[1]], k = top_k,
           overlap = vapply(top_k, function(k) sum(r1 <= k & r2 <= k), integer(1)),
           spearman_all = suppressWarnings(cor(d[[col("pred_raw", pr[[1]])]],
                                               d[[col("pred_raw", pr[[2]])]], method = "spearman")))
  }))
})) %>% mutate(frac = overlap / k)

message("\ntop-k overlap of candidate rankings, per pair of arms")
print(overlap %>% pivot_wider(id_cols = c(arm_a, arm_b, model, spearman_all), names_from = k,
                              values_from = overlap, names_prefix = "top") %>%
        as.data.frame(), row.names = FALSE, digits = 3)

## rank-rank scatter: every other arm against the reference arm, for everything
## that makes the top `rank_lim` under either
rank_lim <- max(top_k)
scat <- bind_rows(lapply(setdiff(arm_names, ref_arm), function(other) {
  r1 <- cand[[col("rank_cand", ref_arm)]]; r2 <- cand[[col("rank_cand", other)]]
  keep <- r1 <= rank_lim | r2 <= rank_lim
  tibble(comparison = other, model = cand$model[keep], peps = cand$peps[keep],
         pep_terminus_hit = cand$pep_terminus_hit[keep],
         rank_ref = r1[keep], rank_other = r2[keep],
         pred_ref = cand[[col("pred", ref_arm)]][keep], pred_other = cand[[col("pred", other)]][keep],
         where = case_when(r1[keep] <= rank_lim & r2[keep] <= rank_lim ~ "both",
                           r1[keep] <= rank_lim ~ paste0("only ", ref_arm),
                           TRUE ~ "only the other arm"))
})) %>% mutate(comparison = factor(comparison, levels = setdiff(arm_names, ref_arm)))

p_rank <- ggplot(scat, aes(x = rank_ref, y = rank_other)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
  geom_point_interactive(
    aes(colour = where, shape = pep_terminus_hit,
        tooltip = sprintf("%s\n%s: rank %d (pred %.3f)\n%s: rank %d (pred %.3f)%s",
                          peps, ref_arm, rank_ref, pred_ref, comparison, rank_other, pred_other,
                          ifelse(pep_terminus_hit, "\nanchors a UniProt peptide terminus", "")),
        data_id = peps),
    size = 1.3, alpha = 0.75) +
  scale_x_log10() + scale_y_log10() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("no", "yes"),
                     name = "anchors a UniProt peptide terminus") +
  scale_colour_manual(values = c("both" = "grey30",
                                 setNames(approach_cols[[ref_arm]], paste0("only ", ref_arm)),
                                 "only the other arm" = "#E6AB02"),
                      name = sprintf("in top %d under", rank_lim)) +
  facet_grid(rows = vars(comparison), cols = vars(model)) +
  coord_equal() +
  labs(x = paste0("candidate rank, ", ref_arm), y = "candidate rank, other arm (row)",
       title = sprintf("Candidate rank agreement vs \"%s\" (unknown windows in the top %d of either)",
                       ref_arm, rank_lim)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical")

## top-N table per approach x terminus, with every other arm's rank alongside
top_tbl <- bind_rows(lapply(arm_names, function(nm) {
  others <- setdiff(arm_names, nm)
  d <- cand %>% filter(.data[[col("rank_cand", nm)]] <= n_top_table) %>%
    arrange(model, .data[[col("rank_cand", nm)]])
  out <- d %>% transmute(approach = nm, model, rank = .data[[col("rank_cand", nm)]], peps, win_type,
                         pred = round(.data[[col("pred", nm)]], 3),
                         raw = sprintf("%.3f ± %.3f", .data[[col("pred_raw", nm)]], .data[[col("pred_sd", nm)]]),
                         top3 = sprintf("%.3f ± %.3f", .data[[col("pred_raw_top3", nm)]], .data[[col("pred_sd_top3", nm)]]),
                         rank_top3 = .data[[col("rank_cand_top3", nm)]],
                         uniprot_hit = pep_terminus_hit)
  for (o in others) out[[paste0("rank | ", o)]] <- d[[col("rank_cand", o)]]
  out
}))

## ---- html ----------------------------------------------------------------------
gir <- function(p, w, h) girafe(
  ggobj = p, width_svg = w, height_svg = h,
  options = list(opts_hover(css = "stroke:black;stroke-width:1.5px;"),
                 opts_hover_inv(css = "opacity:0.25;"),
                 opts_tooltip(css = "background:#fff;border:1px solid #999;padding:4px 6px;font:12px system-ui;white-space:pre;",
                              use_fill = FALSE),
                 opts_sizing(rescale = TRUE, width = 1),
                 opts_toolbar(saveaspng = TRUE)))

fmt_pm <- function(m, s) sprintf("%.3f ± %.3f", m, s)
summary_html <- auc_summary %>%
  transmute(terminus = as.character(term), approach = as.character(approach), params,
            `val n (pos)` = sprintf("%d (%d)", n_val, n_pos),
            `ROC AUC, all seeds` = fmt_pm(roc_auc_mean, roc_auc_sd),
            `ROC AUC, top3 seeds` = fmt_pm(roc_auc_top3_mean, roc_auc_top3_sd),
            `ROC AUC, ensemble` = sprintf("%.3f", roc_auc_ensemble),
            `PR AUC, all seeds`  = fmt_pm(pr_auc_mean, pr_auc_sd),
            `PR AUC, top3 seeds` = fmt_pm(pr_auc_top3_mean, pr_auc_top3_sd),
            `PR AUC, ensemble`  = sprintf("%.3f", pr_auc_ensemble))

overlap_html <- overlap %>%
  transmute(pair = paste(arm_a, "vs", arm_b), terminus = model, k,
            overlap = sprintf("%d / %d (%.0f%%)", overlap, k, 100 * frac),
            `Spearman (all candidates)` = sprintf("%.3f", spearman_all)) %>%
  pivot_wider(id_cols = c(pair, terminus, `Spearman (all candidates)`), names_from = k,
              values_from = overlap, names_prefix = "top ")

html_table <- function(d, caption = NULL, digits = 3) {
  d <- as.data.frame(d)
  th <- tags$tr(lapply(names(d), tags$th))
  rows <- lapply(seq_len(nrow(d)), function(i)
    tags$tr(lapply(d[i, , drop = TRUE], function(v) tags$td(as.character(v)))))
  tags$table(class = "bench", if (!is.null(caption)) tags$caption(caption), tags$thead(th), tags$tbody(rows))
}

css <- tags$style(HTML('
  body { font: 14px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif; color:#222;
         max-width: 1200px; margin: 0 auto; padding: 0 16px 40px; }
  h1 { font-size: 22px; margin: 18px 0 4px; }  h2 { font-size: 17px; margin: 28px 0 6px; }
  p.note { color:#555; margin: 4px 0 10px; }
  table.bench { border-collapse: collapse; margin: 8px 0 14px; font-size: 13px; }
  table.bench caption { text-align:left; font-weight:600; padding: 4px 0; }
  table.bench th, table.bench td { border: 1px solid #ddd; padding: 3px 8px; text-align: left; white-space: nowrap; }
  table.bench th { background: #f3f3f3; }
  details { margin: 8px 0; }  summary { cursor: pointer; font-weight: 600; }
  .row { display:flex; gap: 20px; flex-wrap: wrap; }  .row > div { flex: 1 1 480px; }
'))

intro <- tagList(
  tags$h1(P$title),
  tags$p(class = "note", HTML(sprintf(
    "Same windows, same training recipe, %d arms; termini: %s; %d members per
     terminus each (seeds %d–%d); scores are the ensemble mean, pooled-Platt calibrated.
     Trained %s.",
    length(arm_names), paste(names(nn_input), collapse = ", "), n_seeds, base_seed,
    base_seed + n_seeds - 1L, format(Sys.time(), "%Y-%m-%d %H:%M")))),
  tags$ul(class = "note", lapply(arm_names, function(nm)
    tags$li(HTML(sprintf("<b>%s</b> &mdash; %s (%d input channels, %s params)", nm,
                         approaches[[nm]]$desc, runs[[nm]]$n_channels,
                         paste(runs[[nm]]$params, collapse = " / ")))))),
  tags$h2("Validation AUC"),
  tags$p(class = "note", HTML(sprintf(
    "<b>All seeds</b> is the mean ± sd of the AUC of each member on its own — the
     honest spread of the method. <b>Top3 seeds</b> is the same over the %d best members only,
     so a single collapsed seed does not dominate the sd (biased upward by construction: read it
     as what the method does when it works). <b>Ensemble</b> is the single AUC of the averaged
     score, the number the pipeline actually ranks with. With 6–7 validation positives per
     terminus, PR AUC swings by tenths across seeds: read the sd before the difference.",
    top_n_members))),
  html_table(summary_html)
)

page <- tagList(
  css, intro,
  tags$h2("ROC and precision-recall"),
  tags$div(class = "row",
           tags$div(gir(p_roc, 11, 4.6)),
           tags$div(gir(p_pr,  11, 4.6))),
  tags$h2("Top hits: windows anchoring a UniProt-peptide terminus"),
  tags$p(class = "note", HTML(sprintf(
    "The violin from 10_1dcnn_new6.R with the approach as the outer column facet. Each
     point is a window whose anchor residue sits within ±%d residues of a UniProt-annotated
     peptide terminus of the matching kind (start for N windows, end for C); %d such windows.
     Hover a point: the same window lights up under both approaches.",
    anchor_tol, sum(win_base$pep_terminus_hit)))),
  gir(p_violin, 10, 1.2 + 2.4 * length(arm_names)),
  tags$h2("Do the two arms rank the same candidates?"),
  html_table(overlap_html, caption = "Unknown windows in the top k of BOTH arms of a pair, per terminus"),
  gir(p_rank, 10, 3.2 + 2.6 * (length(arm_names) - 1)),
  tags$details(
    tags$summary(sprintf("Top %d candidate windows per approach and terminus", n_top_table)),
    tags$p(class = "note", HTML(sprintf(
      "<b>raw</b> = ensemble mean ± sd across all %d members; <b>top3</b> = mean ± sd of the
       window’s %d highest member scores, <b>rank_top3</b> its candidate rank by that; the
       <b>rank | …</b> columns give the same window’s rank under each other arm.",
      n_seeds, top_n_members))),
    lapply(split(top_tbl, top_tbl$approach)[arm_names], function(d)
      html_table(d %>% select(-approach), caption = d$approach[[1]]))
  )
)

dir.create(dirname(html_path), showWarnings = FALSE, recursive = TRUE)
save_html(page, file = html_path, libdir = "dependency_files")
## same inlining the per-gene pages use, so the file stands alone
.inl <- if (exists("nn_inline_deps", mode = "function")) nn_inline_deps else
        tryCatch(get("nn_inline_deps", envir = asNamespace("ligandFinder")), error = function(e) NULL)
if (is.function(.inl)) .inl(html_path, file.path(dirname(html_path), "dependency_files"))

stem <- sub("\\.html$", "", html_path)
write.csv(auc_summary, paste0(stem, "_metrics.csv"), row.names = FALSE)
write.csv(overlap,     paste0(stem, "_overlap.csv"), row.names = FALSE)
win_cols <- c("pred", "pred_raw", "pred_sd", "rank_cand", "pred_raw_top3", "pred_sd_top3", "rank_cand_top3")
write.csv(windows %>% select(approach, model, win_type, peps, gene, known, pep_terminus_hit, all_of(win_cols)),
          paste0(stem, "_windows.csv"), row.names = FALSE)
## the same, one row per window with an arm suffix on every score column
short_of <- vapply(approaches, function(a) a$short %||% gsub("[^A-Za-z0-9]+", "_", a$dir), character(1))
windows_wide <- windows %>%
  mutate(arm = unname(short_of[as.character(approach)])) %>%
  select(model, win_type, peps, gene, known, pep_terminus_hit, arm, all_of(win_cols)) %>%
  pivot_wider(id_cols = c(peps, gene, model, win_type, known, pep_terminus_hit),
              names_from = arm, values_from = all_of(win_cols), names_glue = "{.value}_{arm}")
write.csv(windows_wide, paste0(stem, "_windows_wide.csv"), row.names = FALSE)
## static copies of every panel, svg for figures and png for a quick look
for (ext in c("svg", "png")) {
  ggsave(sprintf("%s_violin.%s", stem, ext), p_violin, width = 10, height = 1.2 + 2.4 * length(arm_names), dpi = 150)
  ggsave(sprintf("%s_roc.%s",    stem, ext), p_roc,    width = 11, height = 4.6, dpi = 150)
  ggsave(sprintf("%s_pr.%s",     stem, ext), p_pr,     width = 11, height = 4.6, dpi = 150)
  ggsave(sprintf("%s_rank.%s",   stem, ext), p_rank,   width = 10, height = 3.2 + 2.6 * (length(arm_names) - 1), dpi = 150)
}
message("\nwrote ", html_path, "\n      ", stem, "_{metrics,overlap,windows,windows_wide}.csv",
        "\n      ", stem, "_{violin,roc,pr,rank}.{svg,png}")
