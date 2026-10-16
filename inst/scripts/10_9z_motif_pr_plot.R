#!/usr/bin/env Rscript
## ---- how cleanly the model ranks inside each cleavage-motif group ----------------
## Draws 10_9z_motif_pr.py: PR curves with recall taken WITHIN each group, the
## coverage that normalisation hides, and precision at a matched within-group
## recall.
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
lv   <- c("all anchors", "monobasic", "dibasic", "G-dibasic (amidation)")
cols <- c("all anchors" = "grey35", "monobasic" = "#9BB7D4", "dibasic" = "#197EC0",
          "G-dibasic (amidation)" = "#C24E5C")
## the coverage panel also needs the two classes no group covers
mcol <- c(cols[-1], "terminus" = "#E29C45", "other" = "grey72")
names(mcol)[names(mcol) == "G-dibasic (amidation)"] <- "G + dibasic"
sp <- c(heldout = "29 held-out C ends", all = "all 86 known C ends")
am <- c(counted  = "train / unlabelled ends counted as mistakes",
        excluded = "train / unlabelled ends held out")
fac <- function(d) mutate(d, group = factor(group, levels = lv),
                          split = factor(sp[split], levels = sp),
                          arm   = factor(am[arm], levels = am))
cur <- read_csv(file.path(M, "motif_curves.csv"), show_col_types = FALSE) %>%
  filter(precision > 0) %>% fac()
tab <- read_csv(file.path(M, "motif_table.csv"), show_col_types = FALSE) %>% fac()
cov <- read_csv(file.path(M, "motif_coverage.csv"), show_col_types = FALSE) %>%
  mutate(split = factor(sp[split], levels = sp),
         motif = factor(motif, levels = c("G + dibasic", "dibasic", "monobasic",
                                          "terminus", "other")))

key <- tab %>% filter(split == sp[[1]], arm == am[[1]]) %>%
  transmute(group, lab = sprintf("%s  — %s anchors", group,
                                 format(n_anchors, big.mark = ",", trim = TRUE)))
labs_v <- setNames(key$lab, as.character(key$group))

## ---- (1) PR within each group ----------------------------------------------------
p1 <- ggplot(cur, aes(recall, precision, colour = group)) +
  geom_hline(data = tab, aes(yintercept = prevalence, colour = group),
             linetype = "dotted", linewidth = 0.35, show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  facet_grid(split ~ arm) +
  scale_colour_manual(values = cols, labels = labs_v, name = NULL) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.0001)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Precision-recall inside each cleavage-motif group",
       subtitle = paste0("recall is within the group — denominated by the known ends that group contains — so every curve runs to 100%\n",
                         "right column is the fair one: an anchor that is a real end from outside the positive set is dropped from the ranking, not scored as a mistake\n",
                         "dotted line = that group's prevalence floor (what random ranking gives)"),
       x = "recall of the known C ends in that group", y = "precision (log)")

## ---- (2) what within-group recall hides: the groups are narrow -------------------
p2 <- ggplot(cov, aes(n_ends, fct_rev(split), fill = motif)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_ends > 0, n_ends, "")), position = position_stack(vjust = 0.5),
            size = 2.4, colour = "grey15") +
  scale_fill_manual(values = mcol, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "What each group actually holds",
       subtitle = paste0("100% recall within a group is not 100% of the known ends.\n",
                         "terminus and other sit in no group plotted above."),
       x = "known C ends", y = NULL)

## ---- (3) precision at a matched within-group recall ------------------------------
pa <- tab %>%
  select(group, split, arm, matches("^p_at_r")) %>%
  pivot_longer(matches("^p_at_r"), names_to = "r", values_to = "precision") %>%
  mutate(r = parse_number(r),
         r = factor(sprintf("%g%%", r), levels = sprintf("%g%%", sort(unique(r))))) %>%
  pivot_wider(names_from = arm, values_from = precision) %>%
  rename(counted = !!am[[1]], excluded = !!am[[2]])
p3 <- ggplot(pa, aes(y = fct_rev(group), colour = group)) +
  geom_segment(aes(x = counted, xend = excluded, yend = fct_rev(group)), linewidth = 0.5) +
  geom_point(aes(x = counted), size = 1.6, shape = 21, fill = "white", stroke = 0.6) +
  geom_point(aes(x = excluded), size = 2.1) +
  geom_text(aes(x = pmax(counted, excluded),
                label = ifelse(excluded < 0.001, "<0.1%",
                               scales::percent(excluded, accuracy = 0.1)),
                ## a label that would run off the panel goes above its point instead,
                ## clear of the segment joining the two arms
                hjust = ifelse(pmax(counted, excluded) > 0.55, 1, -0.3),
                vjust = ifelse(pmax(counted, excluded) > 0.55, -0.8, 0.5)),
            size = 2.3, show.legend = FALSE) +
  facet_grid(split ~ r) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0.03, 0.45))) +
  labs(title = "Precision at a matched within-group recall",
       subtitle = paste0("columns are the fraction of that group's own known ends recovered; at each, how clean is the candidate list\n",
                         "hollow = train and unlabelled ends counted as mistakes, filled and labelled = held out"),
       x = "precision among the candidates inspected", y = NULL)

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), legend.text = element_text(size = 7.2),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.1, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.6),
        panel.spacing = unit(7, "pt"))
p <- (p1 & th) / ((p2 & th) | (p3 & th)) +
  plot_layout(heights = c(1.15, 1), widths = c(1, 2.3)) +
  plot_annotation(
    title = "Inside the amidation motif the model is sharp; inside a plain dibasic site it is not",
    subtitle = paste0("2,980,852 C anchors over 5,616 precursors, scored by the 20-seed C_near ensemble (mean of seeds). Motif classes come from lf_pepend_motif_vec, read just outside the anchor:\n",
                      "monobasic = K/R at +1 only  |  dibasic = K/R at +1 and +2  |  G-dibasic = G at +1 then K/R at +2 and +3 (the amidation motif). The classes are disjoint; “all anchors” is the whole scan."),
    caption = paste0("AP and lift are not comparable between groups -- the prevalence floor moves with the group, from 1 in 103,000 for all anchors to 1 in 204 inside the amidation motif.\n",
                     "The held-out monobasic curve rests on 2 ends and should not be read as an estimate; the all-86 panel gives it 10. Ends held out in the right-hand arm: the 57 train C ends, plus\n",
                     "ligand_list and UniProt peptide/propeptide ends -- 1,216 anchors. Only exact anchors are dropped; a neighbour of a real end stays a false positive."),
    theme = theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 7.6, colour = "grey30"),
            plot.caption = element_text(size = 6.8, colour = "grey45", hjust = 0)))
ggsave(file.path(A2, "pepend_motif_prior_pr.pdf"), p, width = 12.2, height = 10.2)
ggsave(file.path(A2, "pepend_motif_prior_pr.png"), p, width = 12.2, height = 10.2, dpi = 185)
message("wrote ", file.path(A2, "pepend_motif_prior_pr.png"))
