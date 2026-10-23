suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
meta <- tribble(~tag, ~label,
  "C_near", "C terminus", "C_ctpad", "C terminus, CT-padded",
  "N_near", "N terminus", "N_ntpad", "N terminus, NT-padded")
cols <- c("C terminus" = "#D95F0E", "C terminus, CT-padded" = "#8C2D04",
          "N terminus" = "#2BA84A", "N terminus, NT-padded" = "#006D2C")
got <- meta %>% filter(file.exists(file.path(S2, paste0("grp_", tag, ".json"))))
d <- map_dfr(seq_len(nrow(got)), function(i) {
  j <- fromJSON(file.path(S2, paste0("grp_", got$tag[i], ".json")), simplifyVector = FALSE)
  base <- j$baseline
  map_dfr(names(j$groups), function(g) {
    x <- j$groups[[g]]
    tibble(label = got$label[i], group = g,
           median_rank = x$median_rank, sd = x$median_rank_sd,
           top500 = x$top500, base_rank = base$median_rank, base_top500 = base$top500,
           n_unknown = j$n_unknown, n_pos = j$n_val_pos,
           fold = x$median_rank / base$median_rank)
  })
})
stopifnot(nrow(d) > 0)
ord <- d %>% group_by(group) %>% summarise(m = median(fold), .groups = "drop") %>% arrange(m) %>% pull(group)
d <- d %>% mutate(group = factor(group, levels = ord))
rand <- d$n_unknown[1] / 2

p1 <- ggplot(d, aes(median_rank, group, colour = label)) +
  geom_vline(aes(xintercept = rand), linetype = 3, colour = "grey45", linewidth = 0.4) +
  geom_linerange(aes(xmin = pmax(median_rank - sd, 1), xmax = median_rank + sd),
                 position = position_dodge(width = 0.65), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.65), size = 1.7) +
  geom_point(aes(base_rank, group), shape = 124, size = 3, colour = "grey25",
             position = position_dodge(width = 0.65)) +
  scale_x_log10(breaks = c(100, 300, 1000, 3000, 10000, 30000),
                labels = scales::label_comma()) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Where the held-out ends rank after permuting a group",
       subtitle = "median rank among the unknown candidate windows; grey tick = that trunk's baseline, dotted = random\nhigher is worse, so a long move right means the group carried real signal",
       x = "median rank of the held-out known ends (log)", y = NULL)

p2 <- ggplot(d, aes(top500, group, colour = label)) +
  geom_point(position = position_dodge(width = 0.65), size = 1.7) +
  geom_point(aes(base_top500, group), shape = 124, size = 3, colour = "grey25",
             position = position_dodge(width = 0.65)) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Held-out ends recovered in the top 500",
       subtitle = "grey tick = baseline; of 29 (C) / 33 (N)",
       x = "ends in the top 500 candidates", y = NULL)

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), axis.text.y = element_text(size = 7.5),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.3, colour = "grey35"))
p <- ((p1 & th) | (p2 & th & theme(axis.text.y = element_blank()))) +
  plot_layout(widths = c(1.6, 1), guides = "collect") +
  plot_annotation(
    title = "Grouped permutation importance, scored on retrieval",
    subtitle = paste("A whole channel (or correlated family) is permuted across all 36 positions by shuffling window identity,",
                     "then the held-out ends are re-ranked among the candidate windows.\nMean of 3 repeats.",
                     "Correlated channels are permuted together because singly they mask one another."),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey30"))) &
  theme(legend.position = "bottom")
f <- file.path(A2, "pepend_grouped_importance.svg")
ggsave(f, p, width = 13, height = 5.6)
write_csv(d, file.path(A2, "pepend_grouped_importance.csv"))
message("wrote ", f, " (", n_distinct(d$label), " trunks)")
