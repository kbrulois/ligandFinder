suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
meta <- tribble(~tag, ~label,
  "C_near", "C terminus", "C_ctpad", "C terminus, CT-padded",
  "N_near", "N terminus", "N_ntpad", "N terminus, NT-padded")
cols <- c("C terminus" = "#D95F0E", "C terminus, CT-padded" = "#8C2D04",
          "N terminus" = "#2BA84A", "N terminus, NT-padded" = "#006D2C")
d <- map_dfr(seq_len(nrow(meta)), function(i) {
  j <- fromJSON(file.path(S2, paste0("shap_", meta$tag[i], ".json")), simplifyVector = FALSE)
  tot <- sum(unlist(j$shapley))
  tibble(label = meta$label[i], group = names(j$shapley),
         phi = unlist(j$shapley), share = unlist(j$shapley) / tot,
         fold = exp(unlist(j$shapley)), total_fold = exp(tot),
         med_none = j$median_none, med_all = j$median_all)
})
ord <- d %>% group_by(group) %>% summarise(m = mean(phi), .groups = "drop") %>% arrange(m) %>% pull(group)
d <- d %>% mutate(group = factor(group, levels = ord))

p1 <- ggplot(d, aes(phi, group, colour = label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(position = position_dodge(width = 0.65), size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(sec.axis = sec_axis(~ exp(.), name = "as a fold improvement in rank",
                                         breaks = c(1, 1.5, 2, 3, 4))) +
  labs(title = "Shapley value by channel group",
       subtitle = "credit for the ranking improvement, averaged over every subset of the other groups\nunits are log rank; the four values sum exactly to the model's total (checked)",
       x = "Shapley value (log-rank units)", y = NULL)

p2 <- ggplot(d, aes(share, group, colour = label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(position = position_dodge(width = 0.65), size = 2) +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Share of the total", subtitle = "same values as a percentage of each trunk's total",
       x = "share of the ranking improvement", y = NULL)

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.4, colour = "grey35"))
lab <- d %>% distinct(label, med_none, med_all, total_fold) %>%
  mutate(t = sprintf("%s: %s → %s (%.0fx)", label, format(round(med_none), big.mark = ","),
                     format(round(med_all), big.mark = ","), total_fold)) %>% pull(t)
p <- ((p1 & th) | (p2 & th & theme(axis.text.y = element_blank()))) +
  plot_layout(widths = c(1.35, 1), guides = "collect") +
  plot_annotation(
    title = "Grouped Shapley: who gets credit when the channels are redundant",
    subtitle = paste0("Permutation gives a group nothing when a correlated neighbour still carries its signal; Shapley splits the shared credit instead.\n",
                      "Median rank of the held-out ends, all channels permuted → intact:  ",
                      paste(lab, collapse = "   ·   ")),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7.6, colour = "grey30"))) &
  theme(legend.position = "bottom")
f <- file.path(A2, "pepend_shapley_importance.svg")
ggsave(f, p, width = 13, height = 5.2)
write_csv(d, file.path(A2, "pepend_shapley_importance.csv"))
message("wrote ", f)
