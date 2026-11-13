suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
tab <- read_csv(file.path(S2, "auc_table.csv"), show_col_types = FALSE)
cur <- read_csv(file.path(S2, "auc_curves.csv"), show_col_types = FALSE)
cur <- cur %>% separate(model, c("series", "head"), sep = " — ", remove = FALSE)
tab <- tab %>% separate(model, c("series", "head"), sep = " — ", remove = FALSE)

mk <- function(kind, which, cols, ttl) {
  d <- cur %>% filter(kind == !!kind, curve == which)
  t <- tab %>% filter(kind == !!kind)
  base <- if (which == "ROC") tibble(x = c(0, 1), y = c(0, 1)) else
    tibble(x = c(0, 1), y = rep(unique(round(t$prevalence, 4))[1], 2))
  lab <- t %>% mutate(v = if (which == "ROC") roc_auc else pr_auc,
                      lo = if (which == "ROC") roc_lo else pr_lo,
                      hi = if (which == "ROC") roc_hi else pr_hi) %>%
    arrange(desc(v)) %>%
    mutate(txt = sprintf("%.3f [%.3f-%.3f]  %s", v, lo, hi, series)) %>%
    group_by(head) %>% summarise(txt = paste(txt, collapse = "\n"), .groups = "drop")
  ggplot(d, aes(x, y, colour = series)) +
    geom_line(data = base, aes(x, y), inherit.aes = FALSE, linetype = 2,
              colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 0.6) +
    geom_text(data = lab, aes(x = if (which == "ROC") 0.98 else 0.98,
                              y = if (which == "ROC") 0.02 else 0.98, label = txt),
              inherit.aes = FALSE, hjust = 1, vjust = if (which == "ROC") 0 else 1,
              size = 2.3, lineheight = 1.05, family = "mono") +
    facet_wrap(~ head, nrow = 1) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = ttl,
         x = if (which == "ROC") "false positive rate" else "recall",
         y = if (which == "ROC") "true positive rate" else "precision") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          legend.key.size = unit(10, "pt"), plot.title = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "grey94", colour = NA),
          panel.spacing = unit(9, "pt"))
}
win_cols <- c("C: current" = "#D95F0E", "C: mask fixed" = "#8C2D04",
              "C: CT-padded, old mask" = "#FDBE85", "C: CT-padded, mask fixed" = "#E6550D",
              "N: current" = "#2BA84A", "N: NT-padded, mask fixed" = "#006D2C")
res_cols <- c("pep_pocket" = "#197EC0", "pep_other" = "#9CBE5A",
              "CT context" = "#F0B429", "NT context" = "#B07C5A")
p <- (mk("window", "ROC", win_cols, "Window models — ROC (held-out candidate windows)") |
      mk("window", "PR", win_cols, "Window models — precision-recall")) /
     (mk("no window", "ROC", res_cols, "Window-free per-residue models — ROC (held-out residues)") |
      mk("no window", "PR", res_cols, "Window-free per-residue models — precision-recall")) +
  plot_annotation(
    title = "Peptide-end models: ROC and precision-recall on each model's own validation split",
    subtitle = paste("Dashed: chance (ROC) / prevalence, the PR floor (PR).",
                     "Numbers are AUC [95% bootstrap CI], highest first.",
                     "\nWindow models score candidate windows (prevalence 6.8% C / 7.6% N);",
                     "residue models score single residues (0.45-1.9%), so their PR AUCs are NOT comparable across the two rows."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey30")))
f <- file.path(A2, "pepend_auc_curves.svg")
ggsave(f, p, width = 13, height = 10)
write_csv(tab %>% select(-series, -head), file.path(A2, "pepend_auc_table.csv"))
message("wrote ", f, " and ", file.path(A2, "pepend_auc_table.csv"))
