#!/usr/bin/env Rscript
## ---- what prior knowledge of the cleavage motif is worth -------------------------
## Draws 10_9z_motif_pr.py: PR curves for six search spaces on one recall axis,
## then the comparable number underneath -- precision at a matched recall.
##
##   Rscript inst/scripts/10_9z_motif_codes.R --term C --tables <tables>
##   python  inst/scripts/10_9z_motif_pr.py   --term C --tables <tables>
##   Rscript inst/scripts/10_9z_motif_pr_plot.R
##
## Run it in a UTF-8 locale: R in a bare "C" locale drops the middle dots and
## curly quotes from the labels without a warning.
## ------------------------------------------------------------------------------
suppressMessages({library(tidyverse); library(patchwork)})
A2 <- path.expand("~/AF2_analysis"); M <- file.path(A2, "pepend_motif")
lv <- c("all anchors", "basic at +1 (any)", "dibasic", "dibasic or G-dibasic",
        "G-dibasic or protein end", "G-dibasic (amidation)")
cols <- c("all anchors" = "grey35", "basic at +1 (any)" = "#9BB7D4", "dibasic" = "#197EC0",
          "dibasic or G-dibasic" = "#5BA053", "G-dibasic or protein end" = "#E29C45",
          "G-dibasic (amidation)" = "#C24E5C")
sp <- c(heldout = "29 held-out C ends", all = "all 86 known C ends")
am <- c(counted  = "train / unlabelled ends counted as mistakes",
        excluded = "train / unlabelled ends held out")
fac <- function(d) mutate(d, universe = factor(universe, levels = lv),
                          split = factor(sp[split], levels = sp),
                          arm   = factor(am[arm], levels = am))
cur <- read_csv(file.path(M, "motif_curves.csv"), show_col_types = FALSE) %>%
  filter(precision > 0) %>% fac()
tab <- read_csv(file.path(M, "motif_table.csv"), show_col_types = FALSE) %>% fac()

## the legend carries the one number that defines each prior: how much is left to search
key <- tab %>% filter(split == sp[[1]], arm == am[[1]]) %>%
  transmute(universe, lab = sprintf("%s  — %s anchors", universe,
                                    format(n_anchors, big.mark = ",", trim = TRUE)))
labs_v <- setNames(key$lab, as.character(key$universe))

## ---- (1) PR curves, one recall axis for every prior ------------------------------
p1 <- ggplot(cur, aes(recall, precision, colour = universe)) +
  geom_hline(data = tab, aes(yintercept = prevalence, colour = universe),
             linetype = "dotted", linewidth = 0.35, show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  geom_point(data = tab, aes(x = recall_ceiling, y = prevalence), size = 1.1,
             show.legend = FALSE) +
  facet_grid(split ~ arm) +
  scale_colour_manual(values = cols, labels = labs_v, name = NULL) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.0001)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  guides(colour = guide_legend(ncol = 3, byrow = TRUE)) +
  labs(title = "Precision-recall for the C-terminal window model, by how far prior knowledge narrows the search",
       subtitle = paste0("recall is over EVERY known C end, not the ends that survive the restriction, so each prior pays for what it discards\n",
                         "right column is the fair one: an anchor that is a real end from outside the positive set is dropped from the ranking, not scored as a mistake\n",
                         "dotted line = that search space's prevalence floor (what random ranking gives); the dot where a curve ends is its recall ceiling"),
       x = "recall of known C ends", y = "precision (log)")

## ---- (2) the trade, read at a matched recall -------------------------------------
pa <- tab %>%
  select(universe, split, arm, matches("^p_at_r")) %>%
  pivot_longer(matches("^p_at_r"), names_to = "r", values_to = "precision") %>%
  mutate(r = parse_number(r),
         r = factor(sprintf("%g%% of known ends found", r),
                    levels = sprintf("%g%% of known ends found", sort(unique(r))))) %>%
  filter(!is.na(precision)) %>%
  pivot_wider(names_from = arm, values_from = precision) %>%
  rename(counted = !!am[[1]], excluded = !!am[[2]])
p2 <- ggplot(pa, aes(y = fct_rev(universe), colour = universe)) +
  geom_segment(aes(x = counted, xend = excluded, yend = fct_rev(universe)), linewidth = 0.5) +
  geom_point(aes(x = counted), size = 1.6, shape = 21, fill = "white", stroke = 0.6) +
  geom_point(aes(x = excluded), size = 2.1) +
  geom_text(aes(x = excluded, label = scales::percent(excluded, accuracy = 0.1)),
            hjust = -0.3, size = 2.3, show.legend = FALSE) +
  facet_grid(split ~ r) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0.02, 0.42))) +
  labs(title = "Precision at a matched recall",
       subtitle = paste0("the comparable number: at the same fraction of known ends recovered, how clean is the candidate list. A missing row is a prior that cannot reach that recall.\n",
                         "hollow = counting train and unlabelled ends as mistakes, filled and labelled = holding them out. The gap is how much of the apparent error was memory, not error."),
       x = "precision among the candidates inspected", y = NULL)

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), legend.text = element_text(size = 7.2),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.1, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.6),
        panel.spacing = unit(7, "pt"))
p <- (p1 & th) / (p2 & th) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(
    title = "Prior knowledge of the cleavage motif buys precision, and pays for it in recall",
    subtitle = paste0("2,980,852 C anchors over 5,616 precursors, scored by the 20-seed C_near ensemble (mean of seeds). Motif classes come from lf_pepend_motif_vec, read just outside the anchor:\n",
                      "dibasic = K/R at +1 and +2  |  G-dibasic = G at +1 then K/R at +2 and +3 (the amidation motif)  |  “basic at +1 (any)” also admits the 269,303 monobasic sites"),
    caption = paste0("AP / lift are not comparable between these search spaces -- the prevalence floor moves with every restriction. Read the matched-recall precision in the lower panel instead.\n",
                     "Ends held out in the right-hand arm: the 57 train C ends, plus ligand_list and UniProt peptide/propeptide ends -- 1,216 anchors in all. Only exact anchors are dropped; a neighbour of a real end stays a false positive."),
    theme = theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 7.6, colour = "grey30"),
            plot.caption = element_text(size = 6.8, colour = "grey45", hjust = 0)))
ggsave(file.path(A2, "pepend_motif_prior_pr.pdf"), p, width = 12.2, height = 10.4)
ggsave(file.path(A2, "pepend_motif_prior_pr.png"), p, width = 12.2, height = 10.4, dpi = 185)
message("wrote ", file.path(A2, "pepend_motif_prior_pr.png"))
