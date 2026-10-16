#!/usr/bin/env Rscript
## ---- does a solvent-exposure floor add anything to the ranking? ------------------
## Draws 10_9z_motif_pr.py run with --asa-thresholds: the two groups worth
## keeping (whole scan, amidation motif) crossed with a floor on the mean relASA
## over +/-3 residues around the anchor.
##
##   Rscript inst/scripts/10_9z_relasa_window.R --term C
##   python  inst/scripts/10_9z_motif_pr.py --term C --tables <tables> \
##             --asa-thresholds 0,0.2,0.4,0.6,0.8 \
##             --only-groups "all anchors,G-dibasic (amidation)" --prefix asa
##   Rscript inst/scripts/10_9z_motif_asa_plot.R
##
## Run it in a UTF-8 locale, or R drops the non-ASCII labels without a warning.
## ------------------------------------------------------------------------------
suppressMessages({library(tidyverse); library(patchwork)})
A2 <- path.expand("~/AF2_analysis"); M <- file.path(A2, "pepend_motif")
ARM <- "excluded"   # the fair arm; the counted one is in the csv
tl  <- c("no floor", "≥ 0.2", "≥ 0.4", "≥ 0.6", "≥ 0.8")
tcol <- setNames(c("grey35", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B"), tl)
gl <- c("all anchors", "G-dibasic (amidation)")
sp <- c(heldout = "29 held-out C ends", all = "all 86 known C ends")
fac <- function(d) mutate(d,
  asa_f = factor(tl[match(asa, c(0, 0.2, 0.4, 0.6, 0.8))], levels = tl),
  motif_group = factor(motif_group, levels = gl),
  split = factor(sp[split], levels = sp))
cur <- read_csv(file.path(M, "asa_curves.csv"), show_col_types = FALSE) %>%
  filter(arm == ARM, precision > 0) %>% fac()
tab <- read_csv(file.path(M, "asa_table.csv"), show_col_types = FALSE) %>%
  filter(arm == ARM) %>% fac()

## ---- (1) PR inside each group, one curve per exposure floor ----------------------
p1 <- ggplot(cur, aes(recall, precision, colour = asa_f)) +
  geom_line(linewidth = 0.6) +
  facet_grid(split ~ motif_group) +
  scale_colour_manual(values = tcol, name = "mean relASA over ±3") +
  scale_y_log10(labels = scales::label_number(accuracy = 0.0001)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Precision-recall with a floor on solvent exposure",
       subtitle = paste0("recall is within the group, denominated by the known ends that survive the floor, so every curve runs to 100%\n",
                         "train and unlabelled real ends are held out of the ranking, not scored as mistakes"),
       x = "recall of the known C ends in that group", y = "precision (log)")

## ---- (2) what each floor costs ---------------------------------------------------
ret <- tab %>%
  select(motif_group, split, asa_f, n_anchors, n_ends) %>%
  pivot_longer(c(n_anchors, n_ends), names_to = "what", values_to = "n") %>%
  filter(!(what == "n_anchors" & split == sp[[2]])) %>%   # anchors do not depend on split
  mutate(what = case_when(what == "n_anchors" ~ "anchors",
                          split == sp[[1]] ~ "known ends, held out",
                          TRUE ~ "known ends, all 86")) %>%
  group_by(motif_group, what) %>% mutate(frac = n / n[asa_f == tl[1]]) %>% ungroup()
wcol <- c("anchors" = "grey45", "known ends, held out" = "#C24E5C",
          "known ends, all 86" = "#E29C45")
p2 <- ggplot(ret, aes(asa_f, frac, colour = what, group = what)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.7) +
  ## the two end series nearly coincide; stagger the labels or they overprint
  geom_text(aes(label = n, vjust = c("anchors" = 1.9, "known ends, held out" = -1.0,
                                     "known ends, all 86" = -2.3)[what]),
            size = 2.1, show.legend = FALSE) +
  facet_wrap(~ motif_group, nrow = 1) +
  scale_colour_manual(values = wcol, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.12)) +
  labs(title = "What the floor removes",
       subtitle = "fraction kept, relative to no floor; labels are counts",
       x = "mean relASA over ±3", y = "kept")

## ---- (3) the point: AP barely moves until the floor starts eating real ends ------
p3 <- ggplot(tab, aes(asa_f, ap, colour = motif_group, group = motif_group)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.9) +
  geom_text(aes(label = scales::number(ap, accuracy = 0.001)), vjust = -1.1,
            size = 2.1, show.legend = FALSE) +
  facet_wrap(~ split, nrow = 1) +
  scale_colour_manual(values = c("all anchors" = "grey35",
                                 "G-dibasic (amidation)" = "#C24E5C"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(title = "Average precision within the group",
       subtitle = "flat through 0.4 even though three quarters of the anchors are gone",
       x = "mean relASA over ±3", y = "AP (within group)")

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), legend.text = element_text(size = 7.2),
        plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.1, colour = "grey35"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.6),
        panel.spacing = unit(7, "pt"))
p <- (p1 & th) / ((p2 & th) | (p3 & th)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "A solvent-exposure floor removes most of the scan without changing the ranking",
    subtitle = paste0("Anchors are dropped when the mean relASA over the 7 residues centred on them falls below the floor. relASA is already one of the model's 26 input channels, and it shows:\n",
                      "the anchors a floor of 0.4 removes were all ranked below every known end already, so precision at every recall is unchanged — the prior is redundant, not wrong."),
    caption = paste0("relASA is 0 at every residue of 716 of the 5,616 precursors (missing DSSP, not burial). No known C end lies in one, so any floor above 0 discards those 351,718 anchors as free negatives.\n",
                     "AP and lift are not comparable across floors -- the prevalence floor rises as anchors are removed. Within the amidation group, lift falls 118 -> 43 from 0 to 0.4 while AP is unchanged at 0.586."),
    theme = theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 7.6, colour = "grey30"),
            plot.caption = element_text(size = 6.8, colour = "grey45", hjust = 0)))
ggsave(file.path(A2, "pepend_motif_relasa_pr.pdf"), p, width = 11.6, height = 9.4)
ggsave(file.path(A2, "pepend_motif_relasa_pr.png"), p, width = 11.6, height = 9.4, dpi = 185)
message("wrote ", file.path(A2, "pepend_motif_relasa_pr.png"))
