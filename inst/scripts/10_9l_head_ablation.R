#!/usr/bin/env Rscript
## ---- attention head vs plain average pooling, on the peptide-end C set --------
## Both arms read the SAME masked class softmax and hand the SAME width to the
## embedding; only the pooling between them differs:
##
##   attn  softmax -> MultiHeadAttention -> LayerNorm -> global average pool
##   gap   softmax ->                                    global average pool
##
## So the whole parameter difference is the attention block (383 on this
## vocabulary, 21,206 vs 20,823), and nothing else about the model moves. This is
## the ablation of `attn`, not a different head reading a different tensor --
## `bottleneck` is that, and it is not what is compared here.
##
## There is a PRIOR, recorded in model.py: on the dibasic-anchored set a plain
## GAP of the softmax "washed out all positional signal and tanked global AUC".
## That was a different training set and a different window geometry, so it is a
## prior and not a conclusion about this one.
##
##   Rscript inst/scripts/10_9l_head_ablation.R --n-seeds 10
##
## Output: ~/AF2_analysis/lf_head_ablation_<term>.{csv,svg,png}
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(ggplot2); library(ggbeeswarm); library(yardstick) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term    <- toupper(.opt("--term", "C"))
in_path <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
n_seeds <- as.integer(.opt("--n-seeds", "10"))
verbose <- as.integer(.opt("--verbose", "0"))
out_st  <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_head_ablation_%s", term)))
refresh <- "--refresh" %in% .args

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))

x <- readRDS(in_path)
nn_input <- x$nn_input[term]                       # ONE terminus: see Config.pepend
message(sprintf("%s terminus: train %d, val %d, all %d windows", term,
                nrow(nn_input[[term]]$train), nrow(nn_input[[term]]$val),
                nrow(nn_input[[term]]$all)))

HEADS <- c(attn = "attention", gap = "average pooling")
runs <- list()
for (h in names(HEADS)) {
  cp <- sprintf("~/AF2_analysis/lf_pepend_head_%s_%s.rds", h, term)
  message(sprintf("\n=== %s head (%s), %d seeds ===", h, HEADS[[h]], n_seeds))
  ## each arm spells out the field that defines it AND gets its own cache path:
  ## the run cache keys on the config, but a shared path would still be a
  ## needless way for two arms to collide
  runs[[h]] <- lf_dcnn_run(nn_input, channel_names = x$all_params3,
                           pepend = term, global_head = h,
                           n_seeds = n_seeds, isolated = TRUE, verbose = verbose,
                           cache = cp, refresh = refresh, keep_arrays = FALSE,
                           path = file.path(ROOT, "inst", "python"))
  message(sprintf("  %s: %d params", h, as.integer(runs[[h]]$params[[term]])))
}

## ---- per-seed metrics from the member files ------------------------------------
reticulate::use_virtualenv("r-tensorflow", required = TRUE)
npc <- reticulate::import("numpy", convert = TRUE)
all_known <- nn_input[[term]]$all$known
val_split <- x$knowns %>% filter(term == !!term) %>% select(peps, split)
all_peps  <- nn_input[[term]]$all$peps
held      <- all_peps %in% (val_split %>% filter(split == "val") %>% pull(peps))

per_seed <- bind_rows(lapply(names(HEADS), function(h) {
  d <- sprintf("~/AF2_analysis/lf_pepend_head_%s_%s_isolated/out/members", h, term)
  f <- sort(Sys.glob(path.expand(file.path(d, sprintf("member_*_%s.npz", term)))))
  stopifnot(length(f) == n_seeds)
  bind_rows(lapply(seq_along(f), function(i) {
    z  <- npc$load(f[[i]])
    vs <- as.numeric(z[[paste0(term, "__val__x")]])
    vl <- factor(as.integer(z[[paste0(term, "__val_labels__x")]]), levels = c(0, 1))
    ga <- as.numeric(z[[paste0(term, "__global__x")]])
    ## rank of each held-out known among the UNKNOWN candidate windows
    unk <- sort(ga[all_known == 0], decreasing = TRUE)
    r <- vapply(ga[held], function(v) sum(unk > v) + 1L, integer(1))
    tibble(head = h, seed = i - 1L,
           pr_auc  = pr_auc_vec(vl, vs, event_level = "second"),
           roc_auc = roc_auc_vec(vl, vs, event_level = "second"),
           med_rank = median(r), top500 = sum(r <= 500))
  }))
})) %>% mutate(head = factor(head, levels = names(HEADS)))

s <- per_seed %>% group_by(head) %>%
  summarise(n = n(), across(c(pr_auc, roc_auc, med_rank, top500),
                            list(mean = mean, sd = sd)), .groups = "drop")
message("\n== per-seed summary ==")
print(as.data.frame(s %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)

w <- function(v) suppressWarnings(wilcox.test(per_seed[[v]] ~ per_seed$head))$p.value
message(sprintf("\nWilcoxon (attn vs gap): PR AUC p = %.3g | median rank p = %.3g | top500 p = %.3g",
                w("pr_auc"), w("med_rank"), w("top500")))
message(sprintf("params: attn %d, gap %d (difference %d = the attention block)",
                as.integer(runs$attn$params[[term]]), as.integer(runs$gap$params[[term]]),
                as.integer(runs$attn$params[[term]] - runs$gap$params[[term]])))
readr::write_csv(per_seed, paste0(out_st, ".csv"))

## ---- plot ------------------------------------------------------------------------
FILL <- setNames(c(LF_VIZ$s1, LF_VIZ$s2), names(HEADS))
panel <- function(v, ylab, ttl, lower_better = FALSE) {
  st <- per_seed %>% group_by(head) %>%
    summarise(m = mean(.data[[v]]), sd = sd(.data[[v]]), .groups = "drop")
  ggplot(per_seed, aes(head, .data[[v]])) +
    geom_quasirandom(aes(fill = head), width = 0.16, shape = 21,
                     colour = LF_VIZ$surface, size = 2.7, stroke = 0.6) +
    geom_linerange(data = st, aes(x = as.numeric(head) + 0.3, ymin = m - sd, ymax = m + sd),
                   inherit.aes = FALSE, colour = LF_VIZ$ink, linewidth = 0.7) +
    geom_point(data = st, aes(x = as.numeric(head) + 0.3, y = m), inherit.aes = FALSE,
               colour = LF_VIZ$ink, size = 2) +
    geom_text(data = st, aes(x = as.numeric(head) + 0.38, y = m,
                             label = sprintf("%.3g\n+/- %.2g", m, sd)),
              inherit.aes = FALSE, colour = LF_VIZ$ink, size = 3, hjust = 0, lineheight = 1.15) +
    scale_fill_manual(values = FILL, guide = "none") +
    scale_x_discrete(labels = HEADS, expand = expansion(mult = c(0.3, 0.72))) +
    labs(title = ttl, subtitle = if (lower_better) "lower is better" else "higher is better",
         x = NULL, y = ylab) +
    lf_viz_theme()
}

suppressMessages(library(patchwork))
fig <- (panel("pr_auc", "validation PR AUC", "Validation PR AUC") |
        panel("med_rank", "median rank of a held-out known",
              "Held-out known ends", lower_better = TRUE)) +
  plot_annotation(
    title = sprintf("Global head: attention vs average pooling (%s terminus, %d seeds each)",
                    term, n_seeds),
    subtitle = paste("Both pool the SAME masked class softmax and hand the same width to the",
                     "embedding; only the pooling differs, so the whole\nparameter difference",
                     sprintf("(%d) is the attention block.",
                             as.integer(runs$attn$params[[term]] - runs$gap$params[[term]])),
                     "Validation PR AUC is also the early-stopping monitor, so it is a\nselection",
                     "score for both arms, not a held-out estimate."),
    theme = theme(plot.background = element_rect(fill = LF_VIZ$surface, colour = NA),
                  plot.title = element_text(colour = LF_VIZ$ink, face = "bold", size = 13),
                  plot.subtitle = element_text(colour = LF_VIZ$ink2, size = 9.5, lineheight = 1.15)))

ggsave(paste0(out_st, ".svg"), fig, width = 9.6, height = 5, device = svglite::svglite)
ggsave(paste0(out_st, ".png"), fig, width = 9.6, height = 5, dpi = 200, bg = LF_VIZ$surface)
message("\nplot  -> ", out_st, ".{svg,png}")
message("table -> ", out_st, ".csv")
