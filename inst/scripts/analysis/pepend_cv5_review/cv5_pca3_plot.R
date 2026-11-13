suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2,"cv5_pca2.rds")); m <- X$meta; ve <- X$ve
## fillable shapes only (21-25), so fill can carry the score
shp <- c(terminus=22, dibasic=21, monobasic=24, `G + dibasic`=25, other=23)
m <- m %>% filter(!is.na(score_cv)) %>%
  mutate(term=factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
         set=factor(set, levels=c("known (train in 4 folds, held out in 1)",
                                  "never seen (no fold trained on it)")),
         motif=factor(motif, levels=names(shp)))
lowest <- m %>% group_by(set, term) %>% slice_min(score_cv, n=4) %>% ungroup()
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
  legend.box="vertical", legend.key.size=unit(10,"pt"), legend.spacing.y=unit(1,"pt"),
  plot.title=element_text(face="bold", size=9.5), plot.subtitle=element_text(size=7.3, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))
guides_common <- guides(shape=guide_legend(order=1, override.aes=list(fill="grey70", size=2.6)),
                        fill=guide_colourbar(order=2, barwidth=unit(55,"pt"), barheight=unit(6,"pt")))

p1 <- ggplot(m, aes(PC1, PC2)) +
  geom_point(aes(fill=score_cv, shape=motif), size=2.5, stroke=0.3, colour="grey25") +
  geom_text_repel(data=lowest, aes(label=label), size=2.1, max.overlaps=14,
                  min.segment.length=0, segment.size=0.22, colour="grey30") +
  facet_grid(term ~ set) +
  scale_shape_manual(values=shp, name="end type", drop=FALSE) +
  scale_fill_viridis_c(option="plasma", name="held-out score", trans="sqrt") +
  guides_common +
  labs(title=sprintf("Peptide-end windows by end type, coloured by held-out score (PC1 %.1f%%, PC2 %.1f%%)", ve[1], ve[2]),
       subtitle="knowns use the fold that held them out; never-seen have no such fold, so their all-50 score is already honest",
       x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]))

p2 <- ggplot(m, aes(d_cen, score_cv)) +
  geom_point(aes(fill=score_cv, shape=motif), size=2.3, stroke=0.25, colour="grey25") +
  geom_smooth(method="lm", se=FALSE, linewidth=0.5, colour="grey45") +
  facet_grid(term ~ set) + scale_y_log10() +
  scale_shape_manual(values=shp, name="end type", drop=FALSE) +
  scale_fill_viridis_c(option="plasma", guide="none", trans="sqrt") +
  guides(shape=guide_legend(override.aes=list(fill="grey70", size=2.6))) +
  labs(title="Score against how typical the window is",
       subtitle="distance to the centroid of all peptide-end windows, first 10 PCs",
       x="distance to centroid", y="held-out score (log)")

p3 <- ggplot(m, aes(motif, score_cv)) +
  geom_boxplot(outlier.shape=NA, linewidth=0.3, fill="grey96", colour="grey50") +
  geom_jitter(aes(fill=score_cv, shape=motif), width=0.16, size=2, stroke=0.25, colour="grey25") +
  facet_grid(term ~ set) + scale_y_log10() +
  scale_shape_manual(values=shp, guide="none") +
  scale_fill_viridis_c(option="plasma", guide="none", trans="sqrt") +
  labs(title="Held-out score by end type", subtitle="the split that actually explains the low scorers",
       x=NULL, y="held-out score (log)") +
  theme(axis.text.x=element_text(angle=25, hjust=1, size=7))

p <- ((p1 & th) / (p2 & th) / (p3 & th)) + plot_layout(heights=c(1.3,1,1)) +
  plot_annotation(title="5-fold CV: every known end is training in 4 folds and held out in 1",
    theme=theme(plot.title=element_text(face="bold", size=12)))
ggsave(file.path(A2,"pepend_cv5_pca_v3.png"), p, width=12.5, height=14, dpi=185)
ggsave(file.path(A2,"pepend_cv5_pca_v3.pdf"), p, width=12.5, height=14)
message("wrote ", file.path(A2,"pepend_cv5_pca_v3.png"))
