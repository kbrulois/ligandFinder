suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
j <- fromJSON(file.path(S2, "seed_vs_ensemble.json"), simplifyVector = FALSE)
lab <- c(C_near = "C terminus", N_near = "N terminus")
d <- imap_dfr(j, function(x, tag) tibble(
  trunk = lab[[tag]], seed = seq_along(x$per_seed) - 1,
  median_rank = map_dbl(x$per_seed, "median"), top10k = map_dbl(x$per_seed, "top10k"),
  val_ap = unlist(x$val_ap), n_val = x$ensemble$n,
  ens_median = x$ensemble$median, ens_top10k = x$ensemble$top10k,
  best_post = x$best_post, best_val = x$best_val, best_train = x$best_train))
ens <- d %>% distinct(trunk, ens_median, ens_top10k, n_val)

## rank within each trunk: reorder() alone sorts globally and leaves the
## second facet unsorted
d <- d %>% group_by(trunk) %>% mutate(pos = rank(median_rank, ties.method = "first")) %>% ungroup()
p1 <- ggplot(d, aes(pos, median_rank)) +
  geom_hline(data = ens, aes(yintercept = ens_median), colour = "#197EC0", linewidth = 0.7) +
  geom_point(size = 1.8, colour = "grey35") +
  geom_point(data = d %>% filter(seed == best_val), colour = "#C24E5C", size = 2.6) +
  geom_text(aes(label = seed), vjust = -1.1, size = 2.2, colour = "grey45") +
  facet_wrap(~ trunk) +
  scale_x_continuous(breaks = NULL) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(title = "Every seed on its own, against the 20-seed ensemble",
       subtitle = "blue line = ensemble; red = the seed you would pick by validation PR AUC; lower is better",
       x = "the 20 seeds, each ordered by its own retrieval (label = seed number)", y = "median rank of the held-out ends (log)")

p2 <- ggplot(d, aes(val_ap, median_rank)) +
  geom_hline(data = ens, aes(yintercept = ens_median), colour = "#197EC0", linewidth = 0.6) +
  geom_point(size = 1.9, colour = "grey30") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, colour = "#C24E5C", formula = y ~ x) +
  facet_wrap(~ trunk, scales = "free") +
  scale_y_log10(labels = scales::label_comma()) +
  labs(title = "Does the selection criterion you have predict the one you want?",
       subtitle = "C: Spearman -0.76 (p < 0.001), the wrong way round.  N: +0.24 (p = 0.31), nothing.",
       x = "validation-window PR AUC (what training early-stops on)",
       y = "median rank of the held-out ends (log)")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA))
p <- (p1 & th) / (p2 & th) +
  plot_annotation(title = "One seed or twenty?",
    subtitle = "The ensemble beats 19 of 20 seeds on both trunks. Seed-to-seed spread is 17x.",
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey30")))
ggsave(file.path(A2, "pepend_seed_vs_ensemble.svg"), p, width = 10.5, height = 7.4)
write_csv(d, file.path(A2, "pepend_seed_vs_ensemble.csv"))
message("wrote ", file.path(A2, "pepend_seed_vs_ensemble.svg"))
