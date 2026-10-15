suppressMessages({library(tidyverse)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
d <- read_csv(file.path(S2, "page_topn.csv"), show_col_types = FALSE)
cols <- c("C terminus" = "#D95F0E", "C terminus, CT-padded" = "#8C2D04",
          "N terminus" = "#2BA84A", "N terminus, NT-padded" = "#006D2C")
npos <- tibble(head = unique(d$head))
## recall at the largest N at or below each cutoff -- the log-spaced grid does
## not land on round numbers
at <- function(rec, N, target) rec[max(which(N <= target))]
lab <- d %>% group_by(head, series) %>%
  summarise(r1k = at(recall, N, 1000), r10k = at(recall, N, 10000), .groups = "drop") %>%
  arrange(head, desc(r10k)) %>%
  mutate(t = sprintf("%5.0f%% @1k  %3.0f%% @10k   %s", 100 * r1k, 100 * r10k, series)) %>%
  group_by(head) %>% summarise(txt = paste(t, collapse = "\n"), .groups = "drop")
p <- ggplot(d, aes(N, recall, colour = series)) +
  geom_vline(xintercept = c(100, 1000, 10000), linetype = 3, colour = "grey70", linewidth = 0.3) +
  geom_step(linewidth = 0.65) +
  geom_text(data = lab, aes(x = 1.3, y = 0.97, label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 2.5, family = "mono", lineheight = 1.1) +
  facet_wrap(~ head, nrow = 1) +
  scale_x_log10(breaks = 10^(0:6), labels = c("1","10","100","1k","10k","100k","1M")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Recall at a top-N cutoff over the full scan",
       subtitle = paste("Only the four trunks the protein pages draw. Held-out known ends (29 C / 33 N) ranked against",
                        "every mature anchor of every precursor (2.98M).\nThis is the curve a screen actually operates on:",
                        "how deep you must go to recover a given share of real peptide ends."),
       x = "top N anchors of the whole scan", y = "held-out ends recovered") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.spacing = unit(10, "pt"))
f <- file.path(A2, "pepend_topn_pages.svg")
ggsave(f, p, width = 11.5, height = 5.0)
write_csv(d, file.path(A2, "pepend_topn_pages.csv"))
message("wrote ", f)
