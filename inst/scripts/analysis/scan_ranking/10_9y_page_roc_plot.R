## ROC for exactly the series the protein pages draw. The full-scan row is
## included only when its scans exist, so this can run before they finish.
suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
term_cols <- c("N terminus" = "#2BA84A", "N terminus, NT-padded" = "#006D2C",
               "C terminus" = "#D95F0E", "C terminus, CT-padded" = "#8C2D04")
res_cols  <- c("pep_pocket" = "#197EC0", "pep_other" = "#9CBE5A",
               "CT context" = "#F0B429", "NT context" = "#B07C5A")
lab_of <- function(d, key) d %>% arrange(desc(roc_auc)) %>%
  mutate(txt = sprintf("%.3f [%.3f-%.3f]  %s", roc_auc, roc_lo, roc_hi, .data[[key]])) %>%
  summarise(txt = paste(txt, collapse = "\n")) %>% pull(txt)
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.spacing = unit(9, "pt"))
roc_panel <- function(cur, tab, key, cols, facet, ttl, sub) {
  ggplot(cur, aes(x, y, colour = .data[[key]])) +
    geom_abline(linetype = 2, colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 0.6) +
    facet_wrap(as.formula(paste("~", facet)), nrow = 1) +
    geom_text(data = tab %>% group_by(.data[[facet]]) %>%
                group_modify(~ tibble(txt = lab_of(.x, key))),
              aes(x = 0.98, y = 0.02, label = txt), inherit.aes = FALSE,
              hjust = 1, vjust = 0, size = 2.4, family = "mono", lineheight = 1.05) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(title = ttl, subtitle = sub, x = "false positive rate", y = "true positive rate") + th
}
## (1) window models, validation negatives
wv  <- read_csv(file.path(S2, "page_roc_val.csv"), show_col_types = FALSE)
wvc <- read_csv(file.path(S2, "page_roc_val_curves.csv"), show_col_types = FALSE)
p1 <- roc_panel(wvc, wv, "series", term_cols, "head",
                "Window models — validation negatives",
                ## C and N have different splits; quote both rather than one
                wv %>% distinct(series, positives, n, prevalence) %>%
                  mutate(t = substr(series, 1, 1)) %>% distinct(t, positives, n, prevalence) %>%
                  summarise(s = paste(sprintf("%s: %d held-out ends among %d candidate windows (%.1f%%)",
                                              t, positives, n, 100 * prevalence), collapse = "   |   ")) %>% pull(s))
## (2) window-free residue models (their validation set IS every residue)
rv  <- read_csv(file.path(S2, "auc_table.csv"), show_col_types = FALSE) %>%
  filter(kind == "no window") %>% separate(model, c("target", "mdl"), sep = " — ")
rvc <- read_csv(file.path(S2, "auc_curves.csv"), show_col_types = FALSE) %>%
  filter(kind == "no window", curve == "ROC") %>% separate(model, c("target", "mdl"), sep = " — ")
p2 <- roc_panel(rvc, rv, "target", res_cols, "mdl",
                "Window-free per-residue models",
                "every residue of every held-out precursor (37,015) -- already the full universe")
plots <- list(p1, p2); hs <- c(1, 1)
## (3) window models over the full scan, when the scans are there
fw <- file.path(S2, "page_roc_scan.csv")
if (file.exists(fw)) {
  sw  <- read_csv(fw, show_col_types = FALSE)
  swc <- read_csv(file.path(S2, "page_roc_scan_curves.csv"), show_col_types = FALSE)
  p3 <- roc_panel(swc, sw, "series", term_cols, "head",
                  "Window models — FULL SCAN negatives",
                  sprintf("the same held-out ends against every other mature anchor (%s anchors, prevalence ~%.0e)",
                          format(sw$n[1], big.mark = ","), mean(sw$prevalence)))
  plots <- list(p1, p3, p2); hs <- c(1, 1, 1)
}
p <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "ROC for the model tracks drawn on the protein pages",
    subtitle = paste("AUC [95% bootstrap CI], highest first. Only the four trunks the pages draw;",
                     "the comparison-only arms are not here.\nEach family is scored against its own",
                     "universe of negatives, so panels are comparable within a row, not across rows."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey30")))
f <- file.path(A2, "pepend_roc_pages.svg")
ggsave(f, p, width = 11, height = 4.1 * length(plots))
message("wrote ", f, " (", length(plots), " rows)")
