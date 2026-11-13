suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2,"cv5_pca2.rds")); m <- X$meta; ve <- X$ve
m <- m %>% mutate(term=factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
                  set=factor(set, levels=c("known (train in 4 folds, held out in 1)",
                                           "never seen (no fold trained on it)")))
lowest <- m %>% filter(!is.na(score_cv)) %>% group_by(set, term) %>% slice_min(score_cv, n=4) %>% ungroup()
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
  legend.key.size=unit(10,"pt"), plot.title=element_text(face="bold", size=9.5),
  plot.subtitle=element_text(size=7.3, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))

p1 <- ggplot(m %>% filter(!is.na(score_cv)), aes(PC1, PC2)) +
  geom_point(aes(fill=score_cv), shape=21, size=2.4, stroke=0.3, colour="grey25") +
  geom_text_repel(data=lowest, aes(label=label), size=2.1, max.overlaps=14,
                  min.segment.length=0, segment.size=0.22, colour="grey30") +
  facet_grid(term ~ set) +
  scale_fill_viridis_c(option="plasma", name="held-out score", trans="sqrt") +
  labs(title=sprintf("Every peptide-end window, coloured by its HELD-OUT score (PC1 %.1f%%, PC2 %.1f%%)", ve[1], ve[2]),
       subtitle="knowns use the fold that held them out; the never-seen set has no such fold, so its all-50 score is already honest",
       x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]))

p2 <- ggplot(m %>% filter(!is.na(score_cv)), aes(d_cen, score_cv)) +
  geom_point(aes(fill=score_cv), shape=21, size=2.2, stroke=0.25, colour="grey25") +
  geom_smooth(method="lm", se=FALSE, linewidth=0.5, colour="grey45") +
  facet_grid(term ~ set) + scale_y_log10() +
  scale_fill_viridis_c(option="plasma", guide="none", trans="sqrt") +
  labs(title="Score against how typical the window is",
       subtitle="distance to the centroid of all peptide-end windows, first 10 PCs -- essentially flat in all four panels",
       x="distance to centroid", y="held-out score (log)")

k <- m %>% filter(grepl("^known", set), !is.na(score_cv), !is.na(score_all50))
p3 <- ggplot(k, aes(score_cv, score_all50)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.35, colour="grey55") +
  geom_point(shape=21, fill="#197EC0", colour="grey25", size=2, stroke=0.25, alpha=0.8) +
  facet_wrap(~term, nrow=1) + scale_x_log10() + scale_y_log10() +
  labs(title="What training on a peptide buys it: in-sample vs held-out, same 177 known ends",
       subtitle="above the line = the all-50 score flatters it. C: 86% of ends; N: 97%. This is why all50 must never be quoted as held out for a known.",
       x="held-out score, from the fold that excluded it (log)", y="all-50 score, 40 of which trained on it (log)")

p <- ((p1 & th) / (p2 & th) / (p3 & th)) + plot_layout(heights=c(1.3,1,1)) +
  plot_annotation(title="5-fold CV: every known end is training in 4 folds and held out in 1",
    theme=theme(plot.title=element_text(face="bold", size=12)))
ggsave(file.path(A2,"pepend_cv5_pca_v2.png"), p, width=12, height=13, dpi=185)
ggsave(file.path(A2,"pepend_cv5_pca_v2.pdf"), p, width=12, height=13)
message("wrote ", file.path(A2,"pepend_cv5_pca_v2.png"))
