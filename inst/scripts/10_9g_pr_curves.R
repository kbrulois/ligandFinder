#!/usr/bin/env Rscript
## ---- precision-recall curves for both window models, over the seeds -----------
## Both models scored over the SAME 29,599 db-anchored candidate windows (the
## production `all` split), 20 seeds each, as one PR curve per model with a band
## across seeds.
##
## The curve is built exactly, not interpolated. With 22 positives there are only
## 22 achievable recall levels, k/22, and precision at that level is k divided by
## the rank of the k-th positive -- so every seed lands on the same recall grid by
## construction and the band across seeds is a like-for-like spread of precision,
## never an artefact of resampling one curve onto another's thresholds.
##
## The band is mean +/- 1 sd of precision ACROSS SEEDS at each recall level. It is
## not a confidence interval and it is not the spread of the AUC; the inset AUCs
## are computed per seed by yardstick (exact step-wise) and summarised separately,
## so the two describe different things and are not derived from each other.
##
## Y IS LOG10. Precision runs from ~0.03 to ~0.5 over a random baseline of
## 0.00074 (22 of 29,599), so a linear axis puts the entire working range in the
## bottom fifth of the panel. The baseline is drawn.
##
## THE COMPARISON IS TILTED TOWARDS THE DB MODEL: these are its window geometry,
## most of the 22 knowns were in its training set, and the peptide-end model is
## scored one residue off its own convention (Spearman 0.996 vs aligned, so the
## cost is small). See 10_9f.
##
##   Rscript inst/scripts/10_9g_pr_curves.R
##
## Output: ~/AF2_analysis/lf_pr_curves_<term>.{csv,svg,png}
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(ggplot2); library(yardstick) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term   <- toupper(.opt("--term", "C"))
arm    <- path.expand(.opt("--arm", "~/AF2_analysis/lf_dcnn_bench_none20/unet_none_w1"))
in_dir <- path.expand(.opt("--in-dir", file.path(dirname(arm), "in")))
run_p  <- path.expand(.opt("--pepend-run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
cache  <- path.expand(.opt("--score-cache", sprintf("~/AF2_analysis/lf_pr_curve_scores_%s.rds", term)))
out_st <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pr_curves_%s", term)))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))

## ---- 1. scores: 20 seeds per model over the same windows ----------------------
if (file.exists(cache) && !("--refresh" %in% .args)) {
  cc <- readRDS(cache); y <- cc$y; db <- cc$db; pe <- cc$pe
  message("reusing ", basename(cache), " (--refresh to rescore)")
} else {
  mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
  np  <- reticulate::import("numpy", convert = FALSE)
  npc <- reticulate::import("numpy", convert = TRUE)
  z <- npc$load(file.path(in_dir, "arrays.npz"))
  X <- z[[paste0(term, "__all__x")]]
  y <- as.numeric(z[[paste0(term, "__all__y_global")]])
  dbm <- sort(Sys.glob(file.path(arm, "members", sprintf("member_*_%s.npz", term))))
  db <- vapply(dbm, function(f) as.numeric(npc$load(f)[[paste0(term, "__global__x")]]),
               numeric(length(y)))
  out_dir <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out")
  cfg <- mod$Config$from_json(file.path(dirname(out_dir), "in", "config.json"))
  models <- mod$scan$load_member_models(out_dir, cfg, term)
  xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))
  pe <- vapply(models, function(m)
    as.numeric(reticulate::py_to_r(mod$pipeline$predict_all(m, xp)[["global"]])),
    numeric(length(y)))
  saveRDS(list(y = y, db = db, pe = pe), cache)
  message("scored -> ", basename(cache))
}
np_pos <- sum(y == 1); base_prec <- mean(y == 1)
message(sprintf("%s windows, %d positive (baseline precision %.5f); %d + %d seeds",
                format(length(y), big.mark = ","), np_pos, base_prec, ncol(db), ncol(pe)))

## ---- 2. the exact curve, per seed ---------------------------------------------
## precision at recall k/P is k / (rank of the k-th positive)
curve_one <- function(sc) {
  r <- which(y[order(-sc, seq_along(sc))] == 1)   # stable: ties broken by index
  tibble(k = seq_along(r), recall = k / np_pos, precision = k / r)
}
LEV <- c("db-anchored", "peptide-end")
cur <- bind_rows(
  bind_rows(lapply(seq_len(ncol(db)), function(i) curve_one(db[, i]) %>% mutate(seed = i - 1L))) %>%
    mutate(model = LEV[1]),
  bind_rows(lapply(seq_len(ncol(pe)), function(i) curve_one(pe[, i]) %>% mutate(seed = i - 1L))) %>%
    mutate(model = LEV[2])) %>%
  mutate(model = factor(model, levels = LEV))

band <- cur %>% group_by(model, k, recall) %>%
  summarise(mean = mean(precision), sd = sd(precision), .groups = "drop") %>%
  ## the lower edge can fall below zero at recall 1/22, where the rank of the
  ## first positive swings hardest; clamp to the random baseline for the log axis
  mutate(lo = pmax(mean - sd, base_prec), hi = mean + sd)

## ---- 3. AUC per seed, summarised independently of the band --------------------
truth <- factor(as.integer(y), levels = c(0, 1))
auc <- bind_rows(
  tibble(model = LEV[1], seed = seq_len(ncol(db)) - 1L,
         pr_auc = apply(db, 2, function(s) pr_auc_vec(truth, s, event_level = "second"))),
  tibble(model = LEV[2], seed = seq_len(ncol(pe)) - 1L,
         pr_auc = apply(pe, 2, function(s) pr_auc_vec(truth, s, event_level = "second")))) %>%
  mutate(model = factor(model, levels = LEV))
sa <- auc %>% group_by(model) %>%
  summarise(mean = mean(pr_auc), sd = sd(pr_auc), .groups = "drop")
message("\nPR AUC per seed (yardstick, exact):")
print(as.data.frame(sa %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
readr::write_csv(band %>% left_join(sa %>% rename(auc_mean = mean, auc_sd = sd), by = "model"),
                 paste0(out_st, ".csv"))

## ---- 4. plot -------------------------------------------------------------------
FILL <- setNames(c(LF_VIZ$s1, LF_VIZ$s2), LEV)
## Inset legend: the swatch carries identity, the text stays in ink tokens.
## placed in the empty low-precision region, well clear of the curves, which run
## along the top of the panel for the whole recall range
ins_x <- 0.02
ins_y <- c(0.0070, 0.0042)
lab <- sprintf("%.3f +/- %.3f", sa$mean, sa$sd)

p <- ggplot(band, aes(recall, mean, colour = model, fill = model)) +
  geom_hline(yintercept = base_prec, colour = LF_VIZ$ink2, linewidth = 0.35, linetype = "22") +
  annotate("text", x = 0.99, y = base_prec, hjust = 1, vjust = -0.7, size = 2.9,
           colour = LF_VIZ$ink2, label = sprintf("random baseline %.5f", base_prec)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5, shape = 21, colour = LF_VIZ$surface, stroke = 0.4) +
  ## inset: swatch + name + mean AUC, text in ink not in the series colour
  annotate("point", x = ins_x, y = ins_y, size = 2.6, shape = 21,
           fill = FILL[LEV], colour = LF_VIZ$surface, stroke = 0.6) +
  annotate("text", x = ins_x + 0.03, y = ins_y, hjust = 0, size = 3.1, colour = LF_VIZ$ink,
           label = sprintf("%s   PR AUC %s", format(LEV, width = 11), lab)) +
  annotate("text", x = ins_x, y = 0.0026, hjust = 0, size = 2.7, colour = LF_VIZ$ink2,
           label = "mean +/- sd over 20 seeds") +
  scale_colour_manual(values = FILL, guide = "none") +
  scale_fill_manual(values = FILL, guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.02)) +
  scale_y_log10(breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1),
                labels = c("0.001", "0.003", "0.01", "0.03", "0.1", "0.3", "1")) +
  labs(title = sprintf("Precision-recall on the db-anchored windows (%s terminus)", term),
       subtitle = sprintf(paste("%s candidate windows, %d known ends, 20 seeds per model.",
                                "Band is +/- 1 sd of precision across seeds.\nEvaluated at the",
                                "%d achievable recall levels; y is log10."),
                          format(length(y), big.mark = ","), np_pos, np_pos),
       x = "recall", y = "precision (log)",
       caption = paste("Tilted towards db-anchored: its window geometry, and most of the 22 knowns were in its training set.",
                       "\nThe band is precision spread, not AUC spread; the inset AUCs are computed per seed and summarised separately.")) +
  lf_viz_theme()

ggsave(paste0(out_st, ".svg"), p, width = 7.6, height = 5.2, device = svglite::svglite)
ggsave(paste0(out_st, ".png"), p, width = 7.6, height = 5.2, dpi = 200, bg = LF_VIZ$surface)
message("\nplot  -> ", out_st, ".{svg,png}")
message("table -> ", out_st, ".csv")
