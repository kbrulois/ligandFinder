suppressMessages({library(tidyverse); library(jsonlite); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
f <- list.files(S2, "^lc_.*_(nms|peaks)\\.csv$", full.names = TRUE)
cur <- map_dfr(f, read_csv, show_col_types = FALSE)
js <- map_dfr(list.files(S2, "^lc_.*_(nms|peaks)\\.json$", full.names = TRUE),
              ~ as_tibble(fromJSON(.x, simplifyVector = TRUE)[c("label","score_key","method",
                  "n_loci","n_val","pct_via_locus","n_neg_all","n_neg_excl","ap_all","ap_excl")]))
stat <- c(score = "mean of 20", score_top5 = "mean of top 5", score_max = "max")
cur <- cur %>% mutate(stat = factor(stat[score_key], levels = stat),
                      label = factor(label, levels = c("C terminus","C terminus, CT-padded",
                                                       "N terminus","N terminus, NT-padded")))
js  <- js %>% mutate(stat = factor(stat[score_key], levels = stat),
                     label = factor(label, levels = levels(cur$label)),
                     lift_all = ap_all / (n_val/(n_val+n_neg_all)),
                     lift_excl = ap_excl / (n_val/(n_val+n_neg_excl)))
scols <- c("mean of 20" = "#197EC0", "mean of top 5" = "#E29C45", "max" = "#C24E5C")

## (1) PR curves, NMS only, both arms
p1 <- cur %>% filter(method == "nms") %>%
  ggplot(aes(recall, precision, colour = stat, linetype = arm)) +
  geom_line(linewidth = 0.55) +
  facet_wrap(~ label, nrow = 1) +
  scale_colour_manual(values = scols, name = NULL) +
  scale_linetype_manual(values = c("all loci counted" = "dashed",
                                   "unlabelled peptides excluded" = "solid"), name = NULL) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Precision-recall over loci (greedy non-maximum suppression)",
       subtitle = "a known end scores as its best locus within 2 residues, or as its own window score if it has none, so recall is not capped\nsolid = loci on real peptides outside this split (train ends, ligand_list, UniProt peptide/propeptide) held out rather than counted as mistakes",
       x = "recall of held-out known ends", y = "precision (log)")

## (2) recall at a budget of loci -- comparable across statistics, unlike AP
bud <- map_dfr(list.files(S2, "^lc_.*\\.json$", full.names = TRUE), function(p) {
  j <- fromJSON(p, simplifyVector = FALSE)
  map_dfr(names(j$budget), ~ tibble(label = j$label, score_key = j$score_key, method = j$method,
                                    N = as.integer(.x), found = j$budget[[.x]], n_val = j$n_val))
}) %>% mutate(stat = factor(stat[score_key], levels = stat),
              label = factor(label, levels = levels(cur$label)),
              recall = found / n_val,
              method = recode(method, nms = "greedy NMS", peaks = "find_peaks"))
p2 <- ggplot(bud, aes(N, recall, colour = stat, linetype = method)) +
  geom_line(linewidth = 0.55) + geom_point(size = 0.8) +
  facet_wrap(~ label, nrow = 1) +
  scale_x_log10(breaks = c(100, 1000, 10000, 50000), labels = c("100","1k","10k","50k")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = scols, name = NULL) +
  scale_linetype_manual(values = c("greedy NMS" = "solid", "find_peaks" = "dotted"), name = NULL) +
  labs(title = "Recall at a budget of loci",
       subtitle = "the comparable measure: AP shifts with how many loci a method emits, this does not",
       x = "loci followed up", y = "held-out ends found")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.box = "vertical",
        legend.key.size = unit(10, "pt"), legend.spacing.y = unit(1, "pt"),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.spacing = unit(8, "pt"))
p <- (p1 & th) / (p2 & th) +
  plot_annotation(title = "Peptide ends as loci: precision-recall and retrieval",
                  subtitle = "Held-out ends only. A locus counts as a hit if its centre is within 2 residues of a known end.",
                  theme = theme(plot.title = element_text(face = "bold", size = 11),
                                plot.subtitle = element_text(size = 8, colour = "grey30")))
ggsave(file.path(A2, "pepend_locus_pr.svg"), p, width = 13.5, height = 8.4)
write_csv(js, file.path(A2, "pepend_locus_summary.csv"))
write_csv(bud, file.path(A2, "pepend_locus_budget.csv"))
message("wrote ", file.path(A2, "pepend_locus_pr.svg"))
