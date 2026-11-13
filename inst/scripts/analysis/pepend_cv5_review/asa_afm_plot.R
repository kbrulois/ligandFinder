suppressMessages({library(tidyverse); library(data.table)})
S2 <- commandArgs(TRUE)[1]
d <- fread(file.path(S2, "asa_afm.csv"))
d[, set := factor(set, levels = unique(set))]

## hollow points, grey fill, black rim. 100,000 of them need a thin rim and a
## small mark or the random panel fills in solid black; 2,571 need a bigger one
## or they vanish, so size and stroke are both mapped to the facet.
p <- ggplot(d, aes(mean_afm, relASA)) +
  geom_point(aes(size = set, stroke = set), shape = 21,
             fill = alpha("grey55", 0.22), colour = alpha("black", 0.55)) +
  geom_density_2d(colour = "#197EC0", linewidth = 0.4, bins = 9, alpha = 0.9) +
  facet_wrap(~ set, nrow = 1) +
  scale_size_manual(values = c(0.45, 1.6, 2.0), guide = "none") +
  scale_discrete_manual("stroke", values = c(0.10, 0.35, 0.40), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1"),
                     expand = expansion(0.02)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(0.02)) +
  labs(title = "Solvent exposure against missense pathogenicity",
       subtitle = paste("Residues with a DSSP assignment only: the other 46% of the secretome",
                        "is intracellular or transmembrane and carries zeros by design.",
                        "\nBlue contours are the 2D density of each panel's own points.",
                        "Pocket residues are a subset of the peptide panel, not a third group."),
       x = "mean AlphaMissense pathogenicity (scaled)",
       y = "relative solvent accessibility") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, colour = "grey35", lineheight = 1.15),
        panel.spacing = unit(16, "pt"))
out <- path.expand("~/Desktop/lf_presentation/asa_vs_afm.png")
ggsave(out, p, width = 14.4, height = 5.6, dpi = 200)
message("wrote ", out)
