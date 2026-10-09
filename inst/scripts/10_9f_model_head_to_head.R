#!/usr/bin/env Rscript
## ---- db-anchored model vs peptide-end model, on ONE common split --------------
## The two models were trained on different sets with different validation
## splits, so their headline PR AUCs (0.868 pepend / its own 429 windows;
## whatever the db arm reports on its own 407) are NOT comparable -- different
## data, different positives, different prevalence. This scores both models,
## 20 seeds each, on a SINGLE labelled set so the numbers mean the same thing.
##
## The common set is the DB MODEL'S OWN VALIDATION SPLIT (407 windows), because
## that is the only direction that can be run at full ensemble strength:
##   * db model     -- its 20 members' scores on that split are already stored in
##                     their member files, nothing is recomputed
##   * pepend model -- its 20 member weight files exist, so its whole ensemble
##                     can be run over those windows here
## The reverse (both models on the peptide-end split) is NOT possible: the db arm
## was trained before `train --member K` saved every member's weights, so only
## member 0 survives and there is no db ensemble to run on new windows.
##
## SO THE COMPARISON IS TILTED TOWARDS THE DB MODEL, and deliberately in that
## direction. This split is the db model's early-stopping monitor set -- its
## weights are the ones that maximised PR AUC here, over up to 2,000 epochs --
## while the peptide-end model has never seen these windows in any capacity.
## Read it as: what does the peptide-end model score on the other model's home
## ground, against an opponent that selected its epoch on exactly this data.
##
## The peptide-end model also sees these windows one residue off its own
## convention (a db-anchored window pins the dibasic at 30-31, so the implied
## peptide end sits at 29, where the peptide-end model was trained to expect 28).
## 10_9 measured that: as-is and +1-aligned rank these windows at Spearman 0.996,
## so scoring as-is costs essentially nothing. It is scored as-is here.
##
##   Rscript inst/scripts/10_9f_model_head_to_head.R
##
## Output: ~/AF2_analysis/lf_model_head_to_head_<term>.{csv,svg,png}
## ------------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(ggplot2); library(ggbeeswarm); library(yardstick)
})

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term   <- toupper(.opt("--term", "C"))
arm    <- path.expand(.opt("--arm", "~/AF2_analysis/lf_dcnn_bench_none20/unet_none_w1"))
in_dir <- path.expand(.opt("--in-dir", file.path(dirname(arm), "in")))
run_p  <- path.expand(.opt("--pepend-run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
out_st <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_model_head_to_head_%s", term)))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)
npc <- reticulate::import("numpy", convert = TRUE)

## ---- 1. the common split ------------------------------------------------------
z <- npc$load(file.path(in_dir, "arrays.npz"))
X <- z[[paste0(term, "__val__x")]]
y <- as.numeric(z[[paste0(term, "__val__y_global")]])
truth <- factor(as.integer(y), levels = c(0, 1))
message(sprintf("common split (the db arm's val): %d windows, %d positive (%.1f%%)",
                length(y), sum(y == 1), 100 * mean(y == 1)))

pr <- function(s) pr_auc_vec(truth, s, event_level = "second")
ro <- function(s) roc_auc_vec(truth, s, event_level = "second")

## ---- 2. db model: its stored per-member scores on this split -------------------
dbm <- sort(Sys.glob(file.path(arm, "members", sprintf("member_*_%s.npz", term))))
db <- vapply(dbm, function(f) as.numeric(npc$load(f)[[paste0(term, "__val__x")]]),
             numeric(length(y)))
stopifnot(nrow(db) == length(y))
message(sprintf("db model    : %d seeds, scores read from %s", ncol(db), basename(arm)))

## ---- 3. peptide-end model: run its whole ensemble over the same windows --------
out_dir <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out")
cfg <- mod$Config$from_json(file.path(dirname(out_dir), "in", "config.json"))
models <- mod$scan$load_member_models(out_dir, cfg, term)
stopifnot(identical(as.integer(cfg$n_channels), dim(X)[[3]]),
          identical(as.integer(cfg$seq_len), dim(X)[[2]]))
xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))
pe <- vapply(models, function(m)
  as.numeric(reticulate::py_to_r(mod$pipeline$predict_all(m, xp)[["global"]])),
  numeric(length(y)))
message(sprintf("pepend model: %d seeds, run here over the same %d windows",
                ncol(pe), length(y)))

## ---- 3b. the same two models on the FULL candidate set ------------------------
## The val split above turns out to be saturated (7 positives in 407; both models
## sit at PR AUC ~1), so it cannot separate them. The `all` split is the harder
## set: 29,599 db-anchored windows with 22 knowns, i.e. the retrieval task the
## models actually do. CAVEAT: most of those 22 were in the db model's TRAINING
## set, so this direction is tilted towards db-anchored too.
Xa <- z[[paste0(term, "__all__x")]]
ya <- as.numeric(z[[paste0(term, "__all__y_global")]])
truth_a <- factor(as.integer(ya), levels = c(0, 1))
pra <- function(s) pr_auc_vec(truth_a, s, event_level = "second")
message(sprintf("\nfull candidate set: %s windows, %d known", format(length(ya), big.mark=","), sum(ya==1)))
db_a <- vapply(dbm, function(f) as.numeric(npc$load(f)[[paste0(term, "__global__x")]]),
               numeric(length(ya)))
xa <- np$ascontiguousarray(np$asarray(Xa, dtype = "float32"))
pe_a <- vapply(models, function(m)
  as.numeric(reticulate::py_to_r(mod$pipeline$predict_all(m, xa)[["global"]])),
  numeric(length(ya)))
## median rank of the knowns among the unknowns, per seed -- the currency the
## rest of this work reports, and far more discriminating than a saturated PR AUC
med_rank <- function(sc) {
  unk <- sort(sc[ya == 0], decreasing = TRUE)
  median(vapply(sc[ya == 1], function(v) sum(unk > v) + 1L, integer(1)))
}
d2 <- bind_rows(
  tibble(model = "db-anchored", seed = seq_len(ncol(db_a)) - 1L,
         pr_auc = apply(db_a, 2, pra), med_rank = apply(db_a, 2, med_rank)),
  tibble(model = "peptide-end", seed = seq_len(ncol(pe_a)) - 1L,
         pr_auc = apply(pe_a, 2, pra), med_rank = apply(pe_a, 2, med_rank)))
d2$model <- factor(d2$model, levels = c("db-anchored", "peptide-end"))
s2 <- d2 %>% group_by(model) %>%
  summarise(n = n(), pr_mean = mean(pr_auc), pr_sd = sd(pr_auc),
            rank_mean = mean(med_rank), rank_sd = sd(med_rank),
            rank_min = min(med_rank), rank_max = max(med_rank), .groups = "drop")
message("full candidate set, per seed:")
print(as.data.frame(s2 %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
message(sprintf("Wilcoxon on median-rank: p = %.3g",
                suppressWarnings(wilcox.test(med_rank ~ model, data = d2))$p.value))
readr::write_csv(d2, sub("$", "_allsplit.csv", out_st))

## ---- 4. per-seed PR AUC, and the ensembles ------------------------------------
d <- bind_rows(
  tibble(model = "db-anchored", seed = seq_len(ncol(db)) - 1L,
         pr_auc = apply(db, 2, pr), roc_auc = apply(db, 2, ro)),
  tibble(model = "peptide-end", seed = seq_len(ncol(pe)) - 1L,
         pr_auc = apply(pe, 2, pr), roc_auc = apply(pe, 2, ro)))
d$model <- factor(d$model, levels = c("db-anchored", "peptide-end"))
ens <- tibble(model = factor(c("db-anchored", "peptide-end"), levels = levels(d$model)),
              pr_auc = c(pr(rowMeans(db)), pr(rowMeans(pe))),
              roc_auc = c(ro(rowMeans(db)), ro(rowMeans(pe))))

s <- d %>% group_by(model) %>%
  summarise(n = n(), mean = mean(pr_auc), sd = sd(pr_auc),
            min = min(pr_auc), max = max(pr_auc),
            roc_mean = mean(roc_auc), .groups = "drop") %>%
  left_join(ens %>% select(model, ensemble = pr_auc), by = "model")
message("\nPR AUC on the common split:")
print(as.data.frame(s %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

## a paired test is not available (the seeds are not paired across models -- they
## are different trainings on different data), so this is the unpaired contrast
w <- suppressWarnings(wilcox.test(pr_auc ~ model, data = d))
message(sprintf("\ndifference in means: %+.4f (peptide-end minus db-anchored); Wilcoxon p = %.3g",
                diff(s$mean), w$p.value))
readr::write_csv(bind_rows(d %>% mutate(kind = "seed"),
                           ens %>% mutate(seed = NA_integer_, kind = "ensemble")),
                 paste0(out_st, ".csv"))

## ---- 5. plot -------------------------------------------------------------------
## Two x-axis groups: identity is position plus the axis label, so colour is
## redundant reinforcement and a legend box would only repeat the tick labels.
suppressMessages(library(patchwork))
FILL <- c("db-anchored" = LF_VIZ$s1, "peptide-end" = LF_VIZ$s2)

swarm <- function(dat, stat, yv, mlab, ttl, sub, ylab, fmt = "%.3f") {
  ggplot(dat, aes(model, .data[[yv]])) +
    geom_quasirandom(aes(fill = model), width = 0.17, shape = 21,
                     colour = LF_VIZ$surface, size = 2.7, stroke = 0.6) +
    geom_linerange(data = stat, aes(x = as.numeric(model) + 0.3,
                                    ymin = .data[[mlab]] - .data[[paste0(mlab, "_sd")]],
                                    ymax = .data[[mlab]] + .data[[paste0(mlab, "_sd")]]),
                   inherit.aes = FALSE, colour = LF_VIZ$ink, linewidth = 0.7) +
    geom_point(data = stat, aes(x = as.numeric(model) + 0.3, y = .data[[mlab]]),
               inherit.aes = FALSE, colour = LF_VIZ$ink, size = 2) +
    geom_text(data = stat, aes(x = as.numeric(model) + 0.38, y = .data[[mlab]],
                label = sprintf(paste0(fmt, "\n+/- ", sub("%\\.[0-9]f", "%s", fmt), " sd"),
                                .data[[mlab]], signif(.data[[paste0(mlab, "_sd")]], 3))),
              inherit.aes = FALSE, colour = LF_VIZ$ink, size = 3, hjust = 0, lineheight = 1.15) +
    scale_fill_manual(values = FILL, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.3, 0.72))) +
    labs(title = ttl, subtitle = sub, x = NULL, y = ylab) +
    lf_viz_theme()
}

s2p <- s2 %>% rename(rank = rank_mean, rank_sd = rank_sd)
pA <- swarm(d2, s2p, "med_rank", "rank",
            "Median rank of the 22 known ends",
            sprintf("Among %s db-anchored candidate windows, per seed. Lower is better.
Wilcoxon p = %.2f -- the two models are not separated.",
                    format(length(ya), big.mark = ","),
                    suppressWarnings(wilcox.test(med_rank ~ model, data = d2))$p.value),
            "median rank of a known (lower better)", fmt = "%.1f")

sv <- s %>% rename(prm = mean, prm_sd = sd)
pB <- swarm(d, sv, "pr_auc", "prm",
            "PR AUC on the db model's val split",
            sprintf("%d windows, %d positive. Both models sit at the ceiling:
this split cannot tell them apart.", length(y), sum(y == 1)),
            "PR AUC")

p <- pA + pB + plot_layout(widths = c(1, 1)) +
  plot_annotation(caption = paste(
    "Both directions favour db-anchored: these are db-anchored windows, most of the 22 knowns were in the db model's TRAINING set,",
    "\nand its val split is its own early-stopping monitor. The peptide-end model has seen none of it, and is scored one residue off its",
    "\nown convention (Spearman 0.996 vs aligned). The reverse comparison is impossible -- the db arm saved only member 0's weights."),
    theme = theme(plot.background = element_rect(fill = LF_VIZ$surface, colour = NA),
                  plot.caption = element_text(colour = LF_VIZ$ink2, size = 8, hjust = 0)))

ggsave(paste0(out_st, ".svg"), p, width = 10, height = 4.9, device = svglite::svglite)
ggsave(paste0(out_st, ".png"), p, width = 10, height = 4.9, dpi = 200, bg = LF_VIZ$surface)
message("\nplot  -> ", out_st, ".{svg,png}")
message("table -> ", out_st, ".csv")
