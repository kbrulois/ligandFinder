suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
d <- read_csv(file.path(S2, "pep_afm_asa.csv"), show_col_types = FALSE) %>%
  mutate(class = factor(class, levels = c("pep_other", "pep_pocket")))
cols <- c(pep_other = "#D5E4A2", pep_pocket = "#197EC0")
labs_ <- c(pep_other = "peptide (non-inserting)", pep_pocket = "peptide (inserting)")

st <- d %>% group_by(class) %>%
  summarise(n = n(), asa = median(relASA), afm = median(min_afm), .groups = "drop")
w_asa <- wilcox.test(relASA ~ class, data = d)
w_afm <- wilcox.test(min_afm ~ class, data = d)
sub <- sprintf(paste0("%s residues from %s peptides over %s precursors  |  medians: relASA %.2f vs %.2f, ",
                      "min_afm %.3f vs %.3f (non-inserting vs inserting)"),
               format(nrow(d), big.mark = ","), n_distinct(d$pep_id),
               n_distinct(d$accession), st$asa[1], st$asa[2], st$afm[1], st$afm[2])

base <- function(p) p +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
  labs(x = "min_afm  (minimum AlphaMissense pathogenicity over all substitutions)",
       y = "relASA  (relative solvent accessibility)") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.3, colour = "grey35"),
        legend.position = "bottom", legend.key.size = unit(10, "pt"))

p1 <- base(ggplot(d, aes(min_afm, relASA)) +
  geom_point(colour = "grey35", alpha = 0.35, size = 0.9, stroke = 0) +
  geom_density_2d(colour = "grey15", linewidth = 0.25, bins = 7) +
  labs(title = "Every residue of every known peptide",
       subtitle = "no distinction drawn between inserting and non-inserting residues"))

p2 <- base(ggplot(d, aes(min_afm, relASA, colour = class)) +
  geom_point(alpha = 0.5, size = 0.9, stroke = 0) +
  geom_density_2d(linewidth = 0.3, bins = 6) +
  scale_colour_manual(values = cols, labels = labs_, name = NULL) +
  labs(title = "The same residues, split by whether they enter the pocket",
       subtitle = sprintf("Wilcoxon: relASA p = %.2g, min_afm p = %.2g", w_asa$p.value, w_afm$p.value)))

ggsave(file.path(A2, "pep_afm_asa_plain.png"),  p1, width = 5.4, height = 5.6, dpi = 190)
ggsave(file.path(A2, "pep_afm_asa_pocket.png"), p2, width = 5.4, height = 5.9, dpi = 190)
p <- (p1 | p2) + plot_annotation(
  title = "Solvent accessibility against AlphaMissense constraint, per peptide residue",
  subtitle = sub,
  caption = paste0("min_afm is the MINIMUM AlphaMissense pathogenicity across every substitution at that residue, so a high value means even the mildest substitution is damaging -- the residue is strongly constrained.\n",
                 "min_afm and relASA are two of the 26 input channels (lf_pepend_residue_cache.rds); peptide spans and pocket membership from lf_pocket_by_peptide.rds. 3 of the 145 cached peptides have no feature matrix and are absent."),
  theme = theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 7.6, colour = "grey30"),
          plot.caption = element_text(size = 6.8, colour = "grey45", hjust = 0)))
ggsave(file.path(A2, "pep_afm_asa.png"), p, width = 10.8, height = 6.2, dpi = 190)
ggsave(file.path(A2, "pep_afm_asa.pdf"), p, width = 10.8, height = 6.2)
print(st); cat(sprintf("\nrelASA  p = %.3g\nmin_afm p = %.3g\n", w_asa$p.value, w_afm$p.value))
message("wrote ", file.path(A2, "pep_afm_asa.png"))
