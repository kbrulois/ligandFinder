suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
lab <- c(pep_pocket = "pocket residue", pep_any = "any peptide residue")
cols <- c("pocket residue" = "#197EC0", "any peptide residue" = "#C2663A")
j <- fromJSON(file.path(S2, "resid_shapley.json"), simplifyVector = FALSE)
d <- map_dfr(names(j), function(t) {
  x <- j[[t]]; tot <- sum(unlist(x$shapley))
  tibble(label = lab[[t]], group = names(x$shapley), phi = unlist(x$shapley),
         share = unlist(x$shapley) / tot, fold = exp(unlist(x$shapley)),
         total_fold = exp(tot), med_none = x$median_none, med_all = x$median_all)
})
d <- d %>% mutate(label = factor(label, levels = unname(lab)))
ord <- d %>% group_by(group) %>% summarise(m = mean(phi), .groups = "drop") %>% arrange(m) %>% pull(group)
d <- d %>% mutate(group = factor(group, levels = ord))

p1 <- ggplot(d, aes(phi, group, colour = label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(position = position_dodge(width = 0.65), size = 2) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(sec.axis = sec_axis(~ exp(.), name = "as a fold improvement in rank",
                                         breaks = c(1, 1.5, 2, 3, 5, 8))) +
  labs(title = "Shapley value by channel group",
       subtitle = "credit for the ranking improvement, averaged over every subset of the other groups\nunits are log rank; the seven values sum exactly to the model's total (checked)",
       x = "Shapley value (log-rank units)", y = NULL)
p2 <- ggplot(d, aes(share, group, colour = label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey55") +
  geom_point(position = position_dodge(width = 0.65), size = 2) +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Share of the total", subtitle = "same values as a percentage of each model's total",
       x = "share of the ranking improvement", y = NULL)
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.4, colour = "grey35"))
hdr <- d %>% distinct(label, med_none, med_all, total_fold) %>%
  mutate(t = sprintf("%s: %s -> %s (%.0fx)", label, format(round(med_none), big.mark = ","),
                     format(round(med_all), big.mark = ","), total_fold)) %>% pull(t)
p <- ((p1 & th) | (p2 & th & theme(axis.text.y = element_blank()))) +
  plot_layout(widths = c(1.35, 1), guides = "collect") +
  plot_annotation(
    title = "Window-free xgb: grouped Shapley over the per-residue channels",
    subtitle = paste0("Each channel carries 5 columns -- itself and +/-2 of context -- and all five move together. ",
                      "`padding` is constant outside the signal peptide, so it is dropped: 7 groups, all 128 subsets.\n",
                      "Geometric-mean rank of the held-out positives among 30,410 rows, all channels permuted -> intact:  ",
                      paste(hdr, collapse = "   .   ")),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7.6, colour = "grey30"))) &
  theme(legend.position = "bottom")
f <- file.path(A2, "pepend_resid_shapley_importance.svg")
ggsave(f, p, width = 13, height = 5.0)
ggsave(sub("[.]svg$", ".png", f), p, width = 13, height = 5.0, dpi = 190)
write_csv(d, file.path(A2, "pepend_resid_shapley_importance.csv"))
message("wrote ", f)
