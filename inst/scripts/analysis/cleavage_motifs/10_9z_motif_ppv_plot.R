#!/usr/bin/env Rscript
## ---- positive predictive value of each cleavage motif ----------------------------
## Draws 10_9z_motif_ppv.py. Panel 1 is the motif on its own -- no model, so the
## PPV of "flag every anchor in this class". Panel 2 adds the model: the PPV of
## "flag the top k of this class".
##
##   python  inst/scripts/10_9z_motif_ppv.py --tables <tables> --term C
##   Rscript inst/scripts/10_9z_motif_ppv_plot.R
##
## Run it in a UTF-8 locale, or R drops the non-ASCII labels without a warning.
## ------------------------------------------------------------------------------
suppressMessages({library(tidyverse); library(patchwork)})
A2 <- path.expand("~/AF2_analysis"); M <- file.path(A2, "pepend_motif")
cl <- c("G + dibasic", "terminus", "dibasic", "monobasic", "other", "ALL (whole scan)")
cols <- c("G + dibasic" = "#C24E5C", "terminus" = "#E29C45", "dibasic" = "#197EC0",
          "monobasic" = "#9BB7D4", "other" = "grey72", "ALL (whole scan)" = "grey25")
gt <- c(annotated = "any annotated end (86 curated + ligand_list + UniProt)",
        knowns    = "the 86 curated ends only")
d <- read_csv(file.path(M, "motif_ppv.csv"), show_col_types = FALSE) %>%
  mutate(motif = factor(motif, levels = rev(cl)))
base <- d %>% filter(motif == "ALL (whole scan)")

## ---- (1) the motif on its own ----------------------------------------------------
p1 <- d %>%
  select(motif, n_anchors, matches("^(ppv|lo|hi)_(knowns|annotated)$")) %>%
  pivot_longer(-c(motif, n_anchors), names_to = c(".value", "truth"), names_sep = "_") %>%
  mutate(truth = factor(gt[truth], levels = gt)) %>%
  ggplot(aes(ppv, motif, colour = motif, shape = truth)) +
  geom_vline(data = tibble(ppv = c(base$ppv_annotated, base$ppv_knowns),
                           truth = factor(gt, levels = gt)),
             aes(xintercept = ppv), linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 2) +
  ## fixed accuracy rounds the curated-ends column to "0.00%"; use significant digits
  geom_text(aes(label = sprintf("%.3g%%", 100 * ppv)), hjust = -0.35, nudge_y = 0.22,
            size = 2.2, show.legend = FALSE) +
  facet_wrap(~ truth, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = cols, guide = "none") +
  scale_shape_manual(values = c(16, 1), guide = "none") +
  scale_x_log10(labels = scales::label_percent(accuracy = 0.001),
                expand = expansion(mult = c(0.08, 0.45))) +
  labs(title = "The motif on its own, with no model",
       subtitle = paste0("PPV of “flag every anchor in this class”, with its Wilson 95% interval. Dashed line = the whole scan, so distance from it is the enrichment.\n",
                         "monobasic and other sit BELOW the baseline: as rules they are worse than picking an anchor at random."),
       x = "positive predictive value (log)", y = NULL)

## ---- (2) the motif plus the model ------------------------------------------------
tk <- d %>%
  select(motif, n_heldout, matches("^ppv_top[0-9]+_(annotated|heldout)$")) %>%
  pivot_longer(starts_with("ppv_top"), names_to = "k", values_to = "ppv") %>%
  separate_wider_regex(k, c("ppv_top", k = "[0-9]+", "_", arm = ".*")) %>%
  mutate(k = as.integer(k),
         arm = factor(c(annotated = "vs any annotated end — the model trained on many of these",
                        heldout = "vs the 29 held-out ends, train / unlabelled dropped")[arm])) %>%
  filter(!is.na(ppv))
## with so few held-out ends the honest arm has a hard ceiling: min(ends, k) / k
ceil <- d %>% filter(motif != "ALL (whole scan)") %>%
  tidyr::expand_grid(k = c(10, 50, 100, 500)) %>%
  transmute(motif, k, ppv = pmin(n_heldout, k) / k,
            arm = factor(levels(tk$arm)[2], levels = levels(tk$arm)))
p2 <- ggplot(tk, aes(k, ppv, colour = motif)) +
  geom_line(data = ceil, aes(group = motif), linetype = "dotted", linewidth = 0.35) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.7) +
  facet_wrap(~ arm, nrow = 1) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_log10(breaks = c(10, 50, 100, 500)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "The motif plus the model: PPV of the top k anchors in each class",
       subtitle = paste0("dotted = the ceiling in the honest arm, min(held-out ends in the class, k) / k. With 12 held-out ends, G + dibasic cannot exceed 24% at k = 50,\n",
                         "and reaches 22% — so it is near the ceiling, not merely low."),
       x = "k, anchors flagged per class (log)", y = "PPV of the flagged set")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), legend.text = element_text(size = 7.2),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.1, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.6), panel.spacing = unit(9, "pt"))
p <- (p1 & th) / (p2 & th) +
  plot_annotation(
    title = "Positive predictive value by cleavage motif, C terminus",
    subtitle = paste0("PPV = TP / (TP + FP): of the anchors a rule flags, the fraction that really are peptide C ends. For a motif alone the rule flags the whole class, so its PPV is the class prevalence.\n",
                      "Both ground truths are lower bounds — an anchor counted as a false positive may simply never have been annotated. 2,980,852 C anchors over 5,616 precursors."),
    caption = paste0("Classes are disjoint and come from lf_pepend_motif_vec: dibasic = K/R at +1,+2 | G + dibasic = G at +1 then K/R at +2,+3 | monobasic = K/R at +1 only | terminus = the anchor is the protein's last residue.\n",
                     "terminus leads on PPV because it is one anchor per precursor and 526 of the 5,616 precursors have an annotated end there -- it is a strong rule that returns at most one candidate per protein."),
    theme = theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 7.6, colour = "grey30"),
            plot.caption = element_text(size = 6.8, colour = "grey45", hjust = 0)))
ggsave(file.path(A2, "pepend_motif_ppv.pdf"), p, width = 11.4, height = 8.2)
ggsave(file.path(A2, "pepend_motif_ppv.png"), p, width = 11.4, height = 8.2, dpi = 185)
message("wrote ", file.path(A2, "pepend_motif_ppv.png"))
