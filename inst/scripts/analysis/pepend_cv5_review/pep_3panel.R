## mean_afm vs relASA: background residues, peptide residues, pocket residues
suppressMessages({library(tidyverse); library(patchwork)})
A2 <- path.expand("~/AF2_analysis"); set.seed(7)
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
pk <- readRDS(file.path(A2, "lf_pocket_by_peptide.rds"))

## every residue covered by a known peptide, so the background can exclude them
pep_span <- pk %>% transmute(accession, r = map2(pep_start, pep_end, seq.int)) %>%
  unnest(r) %>% distinct(accession, r)
in_pep <- split(pep_span$r, pep_span$accession)

## a residue with ALL FOUR channels at exactly 0 has no feature data (~48% of
## the proteome); it is not a measurement and must not join the background
keep <- function(m, r) {
  m[r, "min_afm"] != 0 | m[r, "mean_afm"] != 0 |
    m[r, "cons_rs_n"] != 0 | m[r, "relASA"] != 0
}
bg <- map_dfr(sample(names(fc$feats)), function(acc) {
  m <- fc$feats[[acc]]
  r <- setdiff(seq_len(nrow(m)), in_pep[[acc]])
  if (!length(r)) return(NULL)
  r <- r[keep(m, r)]
  if (!length(r)) return(NULL)
  r <- if (length(r) > 40) sample(r, 40) else r      # cap per protein, spread the draw
  tibble(accession = acc, resno = r, mean_afm = m[r, "mean_afm"], relASA = m[r, "relASA"])
}) %>% slice_sample(n = 10000)

pep <- read_csv(file.path(A2, "pep_residue_channels.csv"), show_col_types = FALSE)
d <- bind_rows(
  bg  %>% transmute(set = "background", mean_afm, relASA),
  pep %>% transmute(set = "peptide",    mean_afm, relASA),
  pep %>% filter(class == "pep_pocket") %>% transmute(set = "pocket", mean_afm, relASA))
LV <- c(background = "10,000 random non-peptide residues",
        peptide    = "all known-peptide residues",
        pocket     = "pocket (inserting) residues only")
d <- d %>% mutate(set = factor(LV[set], levels = LV))
cols <- setNames(c("grey45", "#A9BF63", "#197EC0"), LV)

st <- d %>% group_by(set) %>%
  summarise(n = n(), afm = median(mean_afm), asa = median(relASA), .groups = "drop")
print(st)
cat("\nWilcoxon vs background -- mean_afm:\n")
for (k in LV[2:3]) cat(sprintf("  %-34s p = %.2g\n", k,
  wilcox.test(d$mean_afm[d$set == LV[[1]]], d$mean_afm[d$set == k])$p.value))

p <- ggplot(d, aes(mean_afm, relASA, colour = set)) +
  geom_point(alpha = 0.4, size = 0.85, stroke = 0) +
  geom_density_2d(linewidth = 0.3, bins = 7, colour = "grey15") +
  facet_wrap(~ set, nrow = 1) +
  geom_text(data = st, aes(x = 0.03, y = 0.045,
            label = sprintf("n = %s\nmedian mean_afm %.3f\nmedian relASA %.2f",
                            format(n, big.mark = ","), afm, asa)),
            hjust = 0, vjust = 0, size = 2.5, colour = "grey20", lineheight = 1.15,
            inherit.aes = FALSE) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
  labs(x = "mean_afm  (mean AlphaMissense pathogenicity over all substitutions)",
       y = "relASA  (relative solvent accessibility)",
       title = "Peptide residues are more constrained than the proteome, and pocket residues more again",
       subtitle = sprintf(paste0("median mean_afm %.3f → %.3f → %.3f across the three panels. ",
                                 "Background is drawn from all %s precursors, excluding every residue ",
                                 "inside a known peptide."),
                          st$afm[1], st$afm[2], st$afm[3], format(length(fc$feats), big.mark = ",")),
       caption = paste0("A residue whose cons_rs_n, min_afm, mean_afm and relASA are ALL exactly 0 has no feature data (~48% of the proteome) and is excluded from every panel.\n",
                        "Channels from lf_pepend_residue_cache.rds; peptide spans and pocket membership from lf_pocket_by_peptide.rds. Background capped at 40 residues per precursor before sampling.")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 7.6, colour = "grey30"),
        plot.caption = element_text(size = 6.6, colour = "grey45", hjust = 0),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 8.2), panel.spacing = unit(8, "pt"))
ggsave(file.path(A2, "pep_meanafm_3panel.png"), p, width = 11.6, height = 4.9, dpi = 190)
ggsave(file.path(A2, "pep_meanafm_3panel.pdf"), p, width = 11.6, height = 4.9)
message("wrote ", file.path(A2, "pep_meanafm_3panel.png"))
