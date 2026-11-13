suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2,"cv5_pca4.rds")); m <- X$meta; ve <- X$ve
shp <- c(terminus=22, dibasic=21, monobasic=24, `G + dibasic`=25, other=23)
SETS <- c("known (train in 4 folds, held out in 1)",
          "never seen (peptide, < 60 aa)",
          "large cleavage product ends")
m <- m %>% filter(!is.na(score_cv)) %>%
  mutate(term = factor(paste(term, "terminus"), levels = c("N terminus","C terminus")),
         set   = factor(set, levels = SETS),
         motif = factor(motif, levels = names(shp)))
## label the worst few per panel, plus the TGF-beta family ends that the
## length filter had been excluding
named <- m %>% filter(gene %in% c("GDF5","GDF6","BMP2","TGFB1")) %>% mutate(label=gene)
lowest <- m %>% group_by(set, term) %>% slice_min(score_cv, n = 3) %>% ungroup()
labs_df <- bind_rows(lowest, named) %>% distinct(set, term, label, PC1, PC2, d_cen, score_cv)

th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
  legend.box="vertical", legend.key.size=unit(10,"pt"), legend.spacing.y=unit(1,"pt"),
  plot.title=element_text(face="bold", size=9.5), plot.subtitle=element_text(size=7.3, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))
gc_ <- guides(shape=guide_legend(order=1, override.aes=list(fill="grey70", size=2.6)),
              fill=guide_colourbar(order=2, barwidth=unit(130,"pt"), barheight=unit(7,"pt")))
FILL <- scale_fill_viridis_c(option="plasma", name="honest score", trans="sqrt",
                             limits=c(0,1), breaks=c(0, 0.1, 0.4, 1),
                             labels=c("0", "0.1", "0.4", "1"))

p1 <- ggplot(m, aes(PC1, PC2)) +
  geom_point(aes(fill=score_cv, shape=motif), size=1.9, stroke=0.22, colour="grey25", alpha=0.95) +
  geom_text_repel(data=labs_df, aes(label=label), size=2.1, max.overlaps=20,
                  min.segment.length=0, segment.size=0.2, colour="grey30",
                  box.padding=0.45, point.padding=0.2, force=4, seed=1) +
  facet_grid(term ~ set) + scale_shape_manual(values=shp, name="end type", drop=FALSE) +
  FILL + gc_ +
  labs(title=sprintf("Peptide-end windows by end type, coloured by honest score (PC1 %.1f%%, PC2 %.1f%%)", ve[1], ve[2]),
       subtitle="one PCA over all 879 windows, so the three panels share axes; knowns use the fold that held them out, the two unseen sets their all-50 score",
       x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]))

p2 <- ggplot(m, aes(d_cen, score_cv)) +
  geom_point(aes(fill=score_cv, shape=motif), size=1.8, stroke=0.2, colour="grey25", alpha=0.95) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.5, colour="grey45") +
  facet_grid(term ~ set) + scale_y_log10() +
  scale_shape_manual(values=shp, name="end type", drop=FALSE) +
  scale_fill_viridis_c(option="plasma", guide="none", trans="sqrt", limits=c(0,1)) +
  guides(shape=guide_legend(override.aes=list(fill="grey70", size=2.6))) +
  labs(title="Score against how typical the window is",
       subtitle="distance to the centroid of all end windows, first 10 PCs",
       x="distance to centroid", y="honest score (log)")

p3 <- ggplot(m, aes(motif, score_cv)) +
  geom_boxplot(outlier.shape=NA, linewidth=0.3, fill="grey96", colour="grey50") +
  geom_jitter(aes(fill=score_cv, shape=motif), width=0.17, size=1.6, stroke=0.2,
              colour="grey25", alpha=0.9) +
  facet_grid(term ~ set) + scale_y_log10() +
  scale_shape_manual(values=shp, guide="none") +
  scale_fill_viridis_c(option="plasma", guide="none", trans="sqrt", limits=c(0,1)) +
  labs(title="Honest score by end type",
       subtitle="large products are overwhelmingly 'terminus' ends -- the mature protein's own N and C termini",
       x=NULL, y="honest score (log)") +
  theme(axis.text.x=element_text(angle=25, hjust=1, size=7))

n_lab <- m %>% count(set) %>% mutate(s=sprintf("%s: %d ends", set, n)) %>% pull(s) %>% paste(collapse="   |   ")
p <- ((p1 & th) / (p2 & th) / (p3 & th)) + plot_layout(heights=c(1.3,1,1)) +
  plot_annotation(title="5-fold CV, three end populations",
    subtitle=n_lab,
    theme=theme(plot.title=element_text(face="bold", size=12),
                plot.subtitle=element_text(size=8.5, colour="grey30")))
ggsave(file.path(A2,"pepend_cv5_pca_v4.png"), p, width=15.5, height=14.2, dpi=185)
ggsave(file.path(A2,"pepend_cv5_pca_v4.pdf"), p, width=15.5, height=14.2)
message("wrote ", file.path(A2,"pepend_cv5_pca_v4.png"))
