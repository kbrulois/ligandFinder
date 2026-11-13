## relASA against each per-residue channel, for every residue of every known
## peptide, with and without the pep_pocket split.
suppressMessages({library(tidyverse); library(patchwork)})
A2 <- path.expand("~/AF2_analysis")
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
pk <- readRDS(file.path(A2, "lf_pocket_by_peptide.rds"))
CH <- c("min_afm", "mean_afm", "cons_rs_n")

d <- pmap_dfr(pk, function(accession, gene, pep_id, pep_start, pep_end, pocket) {
  m <- fc$feats[[accession]]
  if (is.null(m)) return(NULL)
  r <- seq.int(pep_start, min(pep_end, nrow(m)))
  tibble(accession, gene, pep_id, resno = r, relASA = m[r, "relASA"],
         min_afm = m[r, "min_afm"], mean_afm = m[r, "mean_afm"],
         cons_rs_n = m[r, "cons_rs_n"],
         class = ifelse(r %in% unlist(pocket), "pep_pocket", "pep_other"))
}) %>% mutate(class = factor(class, levels = c("pep_other", "pep_pocket")),
              ## cons_rs_n, min_afm, mean_afm and relASA are ALL exactly 0 for a
              ## residue the feature pipeline never covered -- proteome-wide that
              ## is 48% of residues, and they would otherwise pile up at (0, 0).
              ## Where the channels disagree, a 0 is a real value.
              no_data = min_afm == 0 & mean_afm == 0 & cons_rs_n == 0 & relASA == 0)
n_nd <- sum(d$no_data)
d <- d %>% filter(!no_data)

## A precursor whose channel is 0 at EVERY residue has no data for it, rather
## than a uniformly zero one. Check whether any peptide sits in such a protein.
cat(sprintf("%s residues | %s peptides | %s precursors\n", format(nrow(d), big.mark=","),
            n_distinct(d$pep_id), n_distinct(d$accession)))
for (ch in c(CH, "relASA")) {
  dead <- names(which(vapply(fc$feats, function(m) all(m[, ch] == 0), logical(1))))
  cat(sprintf("  %-10s all-zero precursors %4d of %d; peptides affected here: %d; "
              , ch, length(dead), length(fc$feats), n_distinct(d$pep_id[d$accession %in% dead])))
  cat(sprintf("residues at exactly 0: %d (%.1f%%)\n", sum(d[[ch]] == 0), 100*mean(d[[ch]] == 0)))
}
write_csv(d, file.path(A2, "pep_residue_channels.csv"))
cat(sprintf("dropped %d residues with no data in ANY channel\n", n_nd))

cols  <- c(pep_other = "#D5E4A2", pep_pocket = "#197EC0")
labs_ <- c(pep_other = "peptide (non-inserting)", pep_pocket = "peptide (inserting)")
NICE  <- c(min_afm   = "min_afm  (minimum AlphaMissense pathogenicity over all substitutions)",
           mean_afm  = "mean_afm  (mean AlphaMissense pathogenicity over all substitutions)",
           cons_rs_n = "cons_rs_n  (conservation, rescaled 0-1; higher = more conserved)")
STEM  <- c(min_afm = "pep_minafm_asa", mean_afm = "pep_meanafm_asa", cons_rs_n = "pep_cons_asa")

for (ch in CH) {
  dd <- d %>% rename(xv = all_of(ch))
  st <- dd %>% group_by(class) %>%
    summarise(x = median(xv), a = median(relASA), .groups = "drop")
  wx <- wilcox.test(xv ~ class, data = dd); wa <- wilcox.test(relASA ~ class, data = dd)
  base <- function(p) p +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
    labs(x = NICE[[ch]], y = "relASA  (relative solvent accessibility)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 7.3, colour = "grey35"),
          legend.position = "bottom", legend.key.size = unit(10, "pt"))
  p1 <- base(ggplot(dd, aes(xv, relASA)) +
    geom_point(colour = "grey35", alpha = 0.35, size = 0.9, stroke = 0) +
    geom_density_2d(colour = "grey15", linewidth = 0.25, bins = 7) +
    labs(title = "Every residue of every known peptide",
         subtitle = "no distinction drawn between inserting and non-inserting residues"))
  p2 <- base(ggplot(dd, aes(xv, relASA, colour = class)) +
    geom_point(alpha = 0.5, size = 0.9, stroke = 0) +
    geom_density_2d(linewidth = 0.3, bins = 6) +
    scale_colour_manual(values = cols, labels = labs_, name = NULL) +
    labs(title = "The same residues, split by whether they enter the pocket",
         subtitle = sprintf("Wilcoxon: %s p = %.2g, relASA p = %.2g", ch, wx$p.value, wa$p.value)))
  ggsave(file.path(A2, paste0(STEM[[ch]], "_plain.png")),  p1, width = 5.4, height = 5.6, dpi = 190)
  ggsave(file.path(A2, paste0(STEM[[ch]], "_pocket.png")), p2, width = 5.4, height = 5.9, dpi = 190)
  p <- (p1 | p2) + plot_annotation(
    title = sprintf("Solvent accessibility against %s, per peptide residue", ch),
    subtitle = sprintf(paste0("%s residues from %s peptides over %s precursors  |  medians: ",
                              "%s %.3f vs %.3f, relASA %.2f vs %.2f  (non-inserting vs inserting)"),
                       format(nrow(dd), big.mark=","), n_distinct(dd$pep_id),
                       n_distinct(dd$accession), ch, st$x[1], st$x[2], st$a[1], st$a[2]),
    caption = paste0("Channels from lf_pepend_residue_cache.rds; peptide spans and pocket membership from lf_pocket_by_peptide.rds. ",
                     sprintf("3 of the 145 cached peptides have no feature matrix and are absent; a further %d residues are dropped because every channel is exactly 0, which marks a residue the feature pipeline never covered.", n_nd)),
    theme = theme_bw(base_size = 9) +
      theme(plot.title = element_text(face = "bold", size = 11.5),
            plot.subtitle = element_text(size = 7.6, colour = "grey30"),
            plot.caption = element_text(size = 6.8, colour = "grey45", hjust = 0)))
  ggsave(file.path(A2, paste0(STEM[[ch]], ".png")), p, width = 10.8, height = 6.2, dpi = 190)
  ggsave(file.path(A2, paste0(STEM[[ch]], ".pdf")), p, width = 10.8, height = 6.2)
  cat(sprintf("\n%-10s medians  other %.3f / pocket %.3f   p = %.2g   -> %s.png\n",
              ch, st$x[1], st$x[2], wx$p.value, STEM[[ch]]))
}
