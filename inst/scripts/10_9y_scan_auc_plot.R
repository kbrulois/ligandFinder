suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
tab <- read_csv(file.path(S2, "scan_auc_table.csv"), show_col_types = FALSE)
cur <- read_csv(file.path(S2, "scan_auc_curves.csv"), show_col_types = FALSE)
cols <- c("C: current" = "#D95F0E", "C: mask fixed" = "#8C2D04",
          "C: CT-padded, old mask" = "#FDBE85", "C: CT-padded, mask fixed" = "#E6550D",
          "N: current" = "#2BA84A", "N: NT-padded, mask fixed" = "#006D2C")
mk <- function(which, logx = FALSE) {
  d <- cur %>% filter(curve == which)
  prev <- unique(tab$prevalence)[1]
  lab <- tab %>% mutate(v = if (which == "ROC") roc_auc else pr_auc,
                        lo = if (which == "ROC") roc_lo else pr_lo,
                        hi = if (which == "ROC") roc_hi else pr_hi) %>%
    arrange(desc(v)) %>%
    mutate(txt = sprintf("%.4f [%.4f-%.4f]  %s", v, lo, hi, model)) %>%
    summarise(txt = paste(txt, collapse = "\n")) %>% pull(txt)
  p <- ggplot(d, aes(x, y, colour = model)) + geom_line(linewidth = 0.6)
  if (which == "ROC") {
    p <- p + geom_abline(linetype = 2, colour = "grey55", linewidth = 0.35) +
      annotate("text", x = 0.98, y = 0.02, label = lab, hjust = 1, vjust = 0,
               size = 2.4, family = "mono", lineheight = 1.05) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      labs(x = "false positive rate", y = "true positive rate",
           title = "ROC — every mature anchor a negative (2.98M)")
  } else {
    p <- p + geom_hline(yintercept = prev, linetype = 2, colour = "grey55", linewidth = 0.35) +
      annotate("text", x = 0.98, y = Inf, label = lab, hjust = 1, vjust = 1.2,
               size = 2.4, family = "mono", lineheight = 1.05) +
      scale_y_log10(labels = scales::label_number(accuracy = 0.0001)) +
      coord_cartesian(xlim = c(0, 1), expand = FALSE) +
      labs(x = "recall", y = "precision (log scale)",
           title = sprintf("Precision-recall — prevalence %.1e (dashed)", prev))
  }
  p + scale_colour_manual(values = cols, name = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          legend.key.size = unit(10, "pt"), plot.title = element_text(face = "bold", size = 9))
}
## the top of the ranking is where a screen actually operates; the PR curve
## squeezes it into the left-hand corner, so plot it directly
tn <- read_csv(file.path(S2, "scan_topn.csv"), show_col_types = FALSE)
p_top <- ggplot(tn, aes(N, recall, colour = model)) +
  geom_step(linewidth = 0.6) +
  scale_x_log10(breaks = 10^(0:6), labels = c("1", "10", "100", "1k", "10k", "100k", "1M")) +
  scale_colour_manual(values = cols, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "top N anchors of the whole scan", y = "fraction of held-out ends found",
       title = "Retrieval: recall at a top-N cutoff") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), plot.title = element_text(face = "bold", size = 9))
p <- (mk("ROC") | mk("PR") | p_top) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Peptide-end window models: ROC and PR over the FULL scan (CNN attention head)",
    subtitle = paste0("Held-out known ends only (29 C / 33 N) against every other mature anchor of every precursor; ",
                      "train positives removed, not counted as negatives.\n",
                      "Numbers are AUC [95% CI from bootstrapping the positives]. ",
                      "PR AUC here rests on a handful of extreme ranks -- the CIs overlap completely, so it does not separate these models."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey30"))) &
  theme(legend.position = "bottom")
f <- file.path(A2, "pepend_auc_fullscan.svg")
ggsave(f, p, width = 17, height = 6.4)
write_csv(tab, file.path(A2, "pepend_auc_fullscan_table.csv"))
message("wrote ", f)
