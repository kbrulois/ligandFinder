#!/usr/bin/env Rscript
## ---- per-seed validation PR AUC for one peptide-end run -----------------------
## Recomputes PR AUC separately for each ensemble member from the validation
## scores its own model produced, and draws the spread as a beeswarm with the
## mean and +/- 1 sd.
##
## READ THE SECOND PANEL BEFORE TRUSTING THE FIRST. The training callback is
## EarlyStopping(monitor = "val_global_pr_auc", mode = "max",
## restore_best_weights = TRUE) over up to 2,000 epochs, so each member's saved
## weights are the ones that maximised THIS metric on THIS split. Every value in
## panel A is therefore a max-over-epochs on the selection set, not a held-out
## estimate -- and this run has no third split: the contract carries only train,
## val and all, and the val split is also fit()'s validation_data and the set the
## Platt calibrator is fitted on.
##
## Panel B prices that in: the monitored val PR AUC at the selected epoch against
## the same quantity at the final epoch, per member. The gap between the two
## columns is the part of panel A that epoch selection bought.
##
## Panel A recomputes PR AUC exactly (yardstick, step-wise); the monitored value
## in panel B is keras's AUC(curve = "PR"), a Riemann sum over a fixed threshold
## grid. The two are different estimators of the same quantity, so they are kept
## in separate panels and never differenced.
##
##   Rscript inst/scripts/10_9e_pepend_seed_pr_auc.R --term C
##
## Output: ~/AF2_analysis/lf_pepend_seed_pr_auc_<term>.{csv,svg,png}
## ------------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(ggplot2); library(ggbeeswarm)
  library(yardstick); library(jsonlite); library(patchwork)
})

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term   <- toupper(.opt("--term", "C"))
run_p  <- path.expand(.opt("--run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
out_st <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_seed_pr_auc_%s", term)))

reticulate::use_virtualenv("r-tensorflow", required = TRUE)
npc <- reticulate::import("numpy", convert = TRUE)

## ---- palette -------------------------------------------------------------
## dataviz reference palette, categorical slots 1-2 on the light surface. Only
## the first three slots are documented as clearing the all-pairs floors (dot
## forms use the all-pairs list), which is exactly what this uses.
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"
GRID    <- "#e6e5e1"
S1      <- "#2a78d6"   # slot 1, blue
S2      <- "#eb6834"   # slot 2, orange

## ---- 1. per-member validation scores ------------------------------------------
mem_dir <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out", "members")
npz <- sort(Sys.glob(file.path(mem_dir, sprintf("member_*_%s.npz", term))))
jsn <- sub("\\.npz$", ".json", npz)
if (!length(npz)) stop("no member files under ", mem_dir, call. = FALSE)

lab <- NULL
rows <- lapply(seq_along(npz), function(i) {
  z <- npc$load(npz[[i]])
  s <- as.numeric(z[[paste0(term, "__val__x")]])
  y <- as.numeric(z[[paste0(term, "__val_labels__x")]])
  ## every member scored the same split; if not, nothing below is comparable
  if (is.null(lab)) lab <<- y else stopifnot(identical(lab, y))
  j <- jsonlite::fromJSON(jsn[[i]])
  h <- j$histories[[term]]$val_global_pr_auc
  tibble(member = i - 1L, seed = as.integer(j$seed), score = list(s),
         epochs = length(h), monitored_best = max(h), monitored_final = tail(h, 1))
})
d <- bind_rows(rows)
truth <- factor(as.integer(lab), levels = c(0, 1))
message(sprintf("%d members; validation split %d windows, %d positive",
                nrow(d), length(lab), sum(lab == 1)))

## ---- 2. PR AUC per member, and for the ensemble --------------------------------
pr <- function(s) pr_auc_vec(truth, s, event_level = "second")
ro <- function(s) roc_auc_vec(truth, s, event_level = "second")
d$pr_auc  <- vapply(d$score, pr, 1)
d$roc_auc <- vapply(d$score, ro, 1)
ens <- rowMeans(do.call(cbind, d$score))     # the ensemble is the mean of scores
ens_pr <- pr(ens); ens_roc <- ro(ens)

m  <- mean(d$pr_auc); s <- sd(d$pr_auc)
message(sprintf("\nPR AUC over %d seeds: mean %.4f, sd %.4f, min %.4f, max %.4f",
                nrow(d), m, s, min(d$pr_auc), max(d$pr_auc)))
message(sprintf("ensemble (mean of member scores): PR AUC %.4f, ROC AUC %.4f", ens_pr, ens_roc))
message(sprintf("ROC AUC over seeds: mean %.4f, sd %.4f", mean(d$roc_auc), sd(d$roc_auc)))
message(sprintf("\nepoch selection: monitored best %.4f +/- %.4f, final %.4f +/- %.4f (mean gap %.4f)",
                mean(d$monitored_best), sd(d$monitored_best),
                mean(d$monitored_final), sd(d$monitored_final),
                mean(d$monitored_best - d$monitored_final)))
readr::write_csv(d %>% select(-score), paste0(out_st, ".csv"))

## ---- 3. plot -------------------------------------------------------------------
base <- theme_minimal(base_size = 11) +
  theme(plot.background  = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
        axis.text  = element_text(colour = INK2),
        axis.title = element_text(colour = INK2),
        plot.title = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = INK2, size = 9.5, lineheight = 1.15),
        plot.caption  = element_text(colour = INK2, size = 8, hjust = 0))

## Panel A: one series, so no legend -- the title names it. Mean/sd and the
## ensemble are annotations in text ink, not a second hue.
## The ensemble (%.4f) and the seed mean (%.4f) differ by less than the marker
## radius, so drawing both collides and says nothing. Only the mean is drawn; the
## ensemble goes in the subtitle, where their near-equality is the actual point.
pa <- ggplot(d, aes(x = "", y = pr_auc)) +
  geom_quasirandom(width = 0.17, shape = 21, fill = S1, colour = SURFACE,
                   size = 2.7, stroke = 0.6) +
  geom_linerange(aes(ymin = m - s, ymax = m + s), x = 1.42, colour = INK, linewidth = 0.7) +
  geom_point(x = 1.42, y = m, colour = INK, size = 2) +
  annotate("text", x = 1.5, y = m, label = sprintf("mean %.3f\n+/- %.3f sd", m, s),
           colour = INK, size = 3.1, hjust = 0, lineheight = 1.15) +
  scale_x_discrete(expand = expansion(mult = c(0.3, 0.95))) +
  labs(title = sprintf("Validation PR AUC by seed (%s terminus, %d seeds)", term, nrow(d)),
       subtitle = sprintf("%d windows, %d positive | range %.3f-%.3f\nEnsemble of all %d: %.3f - no better than the average seed",
                          length(lab), sum(lab == 1), min(d$pr_auc), max(d$pr_auc),
                          nrow(d), ens_pr),
       x = NULL, y = "PR AUC") +
  base

## Panel B: two x-axis groups, so identity is position + axis label, not colour
## alone; colour is redundant reinforcement and a legend box would only repeat
## the tick labels. Lines pair the two values for the same seed.
b <- d %>%
  select(member, best = monitored_best, final = monitored_final) %>%
  tidyr::pivot_longer(c(best, final), names_to = "epoch", values_to = "v") %>%
  mutate(epoch = factor(epoch, levels = c("best", "final"),
                        labels = c("selected epoch", "final epoch")))
bs <- b %>% group_by(epoch) %>% summarise(m = mean(v), s = sd(v), .groups = "drop")

pb <- ggplot(b, aes(epoch, v)) +
  geom_line(aes(group = member), colour = GRID, linewidth = 0.35) +
  geom_linerange(data = bs, aes(x = as.numeric(epoch) + 0.28, ymin = m - s, ymax = m + s),
                 inherit.aes = FALSE, colour = INK, linewidth = 0.6) +
  geom_point(data = bs, aes(x = as.numeric(epoch) + 0.28, y = m),
             inherit.aes = FALSE, colour = INK, size = 1.8) +
  geom_point(aes(fill = epoch), shape = 21, colour = SURFACE, size = 2.7, stroke = 0.6,
             position = position_quasirandom(width = 0.13)) +
  scale_fill_manual(values = c("selected epoch" = S1, "final epoch" = S2), guide = "none") +
  labs(title = "What epoch selection bought",
       subtitle = sprintf("Monitored val PR AUC, same split, per seed  |  mean gap %.3f\nEarly stopping maximises it on this split, so panel A is a max over epochs.",
                          mean(d$monitored_best - d$monitored_final)),
       x = NULL, y = "val_global_pr_auc (monitored)") +
  base

p <- pa + pb + plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    caption = paste("The validation split is fit()'s validation_data, the early-stopping monitor and the",
                    "calibration set; there is no held-out test split.\nPanel A is yardstick's exact",
                    "step-wise PR AUC; panel B is keras's threshold-grid approximation. Different",
                    "estimators -- not differenced across panels."),
    theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA),
                  plot.caption = element_text(colour = INK2, size = 8, hjust = 0)))

ggsave(paste0(out_st, ".svg"), p, width = 9.4, height = 4.5, device = svglite::svglite)
ggsave(paste0(out_st, ".png"), p, width = 9.4, height = 4.5, dpi = 200, bg = SURFACE)
message("\nplot -> ", out_st, ".{svg,png}")
message("table -> ", out_st, ".csv")

print(as.data.frame(d %>% select(member, seed, epochs, pr_auc, roc_auc,
                                 monitored_best, monitored_final) %>%
                    mutate(across(where(is.numeric), ~round(.x, 4))) %>%
                    arrange(desc(pr_auc))), row.names = FALSE)
