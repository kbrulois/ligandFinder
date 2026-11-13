suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
meta <- tribble(~tag, ~label, ~term, ~anchor,
  "C_near",  "C terminus",             "C", 28,
  "C_ctpad", "C terminus, CT-padded",  "C", 28,
  "N_near",  "N terminus",             "N",  8,
  "N_ntpad", "N terminus, NT-padded",  "N",  8)
cols <- c("C terminus" = "#D95F0E", "C terminus, CT-padded" = "#8C2D04",
          "N terminus" = "#2BA84A", "N terminus, NT-padded" = "#006D2C")
got <- meta %>% filter(file.exists(file.path(S2, paste0("imp_", tag, ".json"))))
perm <- map_dfr(seq_len(nrow(got)), function(i) {
  d <- fromJSON(file.path(S2, paste0("imp_", got$tag[i], ".json")))
  if (is.null(d$perm_drop)) return(NULL)
  m <- d$perm_drop            # channels x positions
  expand_grid(ci = seq_len(nrow(m)), pos = seq_len(ncol(m))) %>%
    mutate(channel = d$channels[ci], drop = map2_dbl(ci, pos, ~ m[.x, .y]),
           label = got$label[i], anchor = got$anchor[i], baseline = d$baseline_roc)
})
stopifnot(nrow(perm) > 0)
## channels the question was about, plus the strongest others for context
focus <- c("relASA", "cons_rs_n", "min_afm", "mean_afm")
nice <- c(relASA = "relASA (solvent accessibility)", cons_rs_n = "conservation (cons_rs_n)",
          min_afm = "min_afm", mean_afm = "mean_afm")

p_line <- perm %>% filter(channel %in% focus) %>%
  mutate(channel = factor(nice[channel], levels = nice[focus])) %>%
  ggplot(aes(pos, drop, colour = label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(data = distinct(perm, label, anchor), aes(xintercept = anchor),
             linetype = 2, linewidth = 0.4, colour = "grey35") +
  geom_line(linewidth = 0.55) + geom_point(size = 0.7) +
  facet_wrap(~ channel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = c(1, 8, 15, 22, 28, 36)) +
  labs(title = "Input-channel importance along the window",
       subtitle = "drop in validation ROC AUC when that channel is permuted at that position (mean of 3 repeats)\ndashed line = the anchor: the peptide's own terminal residue (28 for C, 8 for N)",
       x = "position in window", y = "ROC AUC drop") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), strip.text = element_text(size = 7.5),
        strip.background = element_rect(fill = "grey94", colour = NA),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))

ord <- perm %>% group_by(channel) %>% summarise(m = max(drop), .groups = "drop") %>% arrange(m) %>% pull(channel)
p_heat <- perm %>% mutate(channel = factor(channel, levels = ord)) %>%
  ggplot(aes(pos, channel, fill = drop)) +
  geom_tile() +
  geom_vline(data = distinct(perm, label, anchor), aes(xintercept = anchor + 0.5),
             linetype = 2, linewidth = 0.35, colour = "grey20") +
  facet_wrap(~ label, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "ROC AUC drop") +
  scale_x_continuous(breaks = c(1, 8, 15, 22, 28, 36), expand = c(0, 0)) +
  labs(title = "Every input channel x window position",
       x = "position in window", y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.key.height = unit(7, "pt"), legend.key.width = unit(26, "pt"),
        axis.text.y = element_text(size = 6),
        strip.background = element_rect(fill = "grey94", colour = NA),
        plot.title = element_text(face = "bold", size = 10))

## the head's OWN features, for contrast: it never sees a channel
gain <- map_dfr(seq_len(nrow(got)), function(i) {
  d <- fromJSON(file.path(S2, paste0("imp_", got$tag[i], ".json")))
  as_tibble(d$head_gain) %>% mutate(label = got$label[i])
})
top <- gain %>% group_by(feature) %>% summarise(m = mean(gain), .groups = "drop") %>%
  slice_max(m, n = 14) %>% pull(feature)
p_gain <- gain %>% filter(feature %in% top) %>%
  mutate(feature = fct_reorder(feature, gain, mean)) %>%
  ggplot(aes(gain, feature, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = cols, name = NULL) +
  labs(title = "The head's own 48 features (pooled CNN class surface)",
       subtitle = "mean gain over the 20 per-member boosters -- no input channel appears here",
       x = "gain", y = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))

p <- (p_line | p_gain) / p_heat + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(title = "xgb window head: what it uses, and which inputs it ultimately depends on",
                  theme = theme(plot.title = element_text(face = "bold", size = 11)))
f <- file.path(A2, "pepend_xgb_importance.svg")
ggsave(f, p, width = 14, height = 11)
write_csv(perm, file.path(A2, "pepend_channel_importance.csv"))
message("wrote ", f, " (", n_distinct(perm$label), " trunks)")
