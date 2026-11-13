## relASA vs AlphaMissense, three panels as asa_vs_afm.png, now with each peptide
## residue filled by its insertion depth (asa_afm_depth.R).
##
## Depth is POLARITY around the extracellular surface (0): inside the receptor vs
## above it, so it takes a diverging scale -- blue = inside, red = outside, a
## neutral grey midpoint. The two arms differ in extent (-1.07 .. 0.33), so each
## gets the full ramp around 0 rather than the inside arm reading pale. The fill
## carries depth; a thin dark rim keeps a near-zero (near-white) residue visible.
## Grey points have no depth: random residues are not peptides, and a peptide
## residue with no receptor contact (area > 1, dist < 6) has no placement.
suppressMessages({library(tidyverse); library(data.table)})
A2 <- path.expand("~/AF2_analysis")
d <- fread(file.path(A2, "lf_headcmp", "asa_afm_depth.csv"))
lv <- c("100,000 random residues", "all peptide residues", "pocket residues only")
nn <- d[, .N, by = set]; nn <- setNames(nn$N, nn$set)
lab <- c(`100,000 random residues` = "100,000 random residues",
         `all peptide residues` = sprintf("all peptide residues (n = %s)", format(nn[["all peptide residues"]], big.mark = ",")),
         `pocket residues only` = sprintf("pocket residues only (n = %s)", nn[["pocket residues only"]]))
d[, set := factor(set, levels = lv)]
lo <- min(d$depth, na.rm = TRUE); hi <- max(d$depth, na.rm = TRUE)
OUT_COL <- "#e34948"; MID_COL <- "#f0efec"; IN_COL <- "#1c5cab"
nodep <- d[is.na(depth)]
dep   <- d[!is.na(depth)][order(depth)]          # the most inserted residues drawn last, on top

p <- ggplot(d, aes(mean_afm, relASA)) +
  ## no depth: HOLLOW and quieter, so it reads as background and never as a
  ## near-zero (near-white) depth
  geom_point(data = nodep, aes(size = set, stroke = set), shape = 21,
             fill = NA, colour = alpha("grey45", 0.40)) +
  geom_point(data = dep, aes(fill = depth, size = set), shape = 21, stroke = 0.3,
             colour = alpha("grey15", 0.7)) +
  geom_density_2d(colour = "#3F4A56", linewidth = 0.35, bins = 9, alpha = 0.85) +
  facet_wrap(~ set, nrow = 1, labeller = as_labeller(lab)) +
  scale_size_manual(values = c(0.45, 1.5, 2.1), guide = "none") +
  scale_discrete_manual("stroke", values = c(0.10, 0.30, 0.30), guide = "none") +
  scale_fill_gradientn(colours = c(OUT_COL, MID_COL, IN_COL),
                       values = scales::rescale(c(lo, 0, hi)), limits = c(lo, hi),
                       breaks = c(-1, -0.5, 0, 0.3),
                       labels = c("-1.0  above the receptor", "-0.5", "0  extracellular surface", "0.3  inside"),
                       name = "insertion depth\n(1 - IC_dist)") +
  guides(fill = guide_colourbar(barheight = unit(150, "pt"), barwidth = unit(9, "pt"),
                                ticks.colour = "grey20")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.50", "0.75", "1"),
                     expand = expansion(0.02)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(0.02)) +
  labs(title = "Solvent exposure against missense pathogenicity, coloured by insertion depth",
       subtitle = paste0(
         "Residues with a DSSP assignment only (the other 46% of the secretome is intracellular or transmembrane and carries zeros by design).  ",
         "Grey contours: 2D density of each panel's own points.\n",
         "Fill: how deep the residue sits in its receptor in the docked model, along the receptor's axis as a fraction of its length (0 = extracellular surface).  ",
         "Hollow grey points have no depth: random residues are not docked,\n",
         sprintf("and %s of the %s peptide residues make no receptor contact.  Pocket residues are a subset of the peptide panel.",
                 format(sum(is.na(d[set == "all peptide residues"]$depth)), big.mark = ","),
                 format(nn[["all peptide residues"]], big.mark = ","))),
       x = "mean AlphaMissense pathogenicity (scaled)",
       y = "relative solvent accessibility") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 8.6, colour = "grey35", lineheight = 1.15),
        legend.title = element_text(size = 9.5, colour = "grey15"),
        legend.text = element_text(size = 8.5, colour = "grey25"),
        panel.spacing = unit(16, "pt"))
out <- file.path(A2, "asa_vs_afm_depth.png")
ggsave(out, p, width = 16.2, height = 6.0, dpi = 200)
ggsave(sub("\\.png$", ".pdf", out), p, width = 16.2, height = 6.0)
message("wrote ", out)
