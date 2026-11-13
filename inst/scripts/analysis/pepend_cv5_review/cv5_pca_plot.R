suppressMessages({library(tidyverse); library(patchwork); library(ggrepel)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2,"cv5_pca.rds")); m <- X$meta; ve <- X$ve
m <- m %>% mutate(term=factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
                  set=factor(set, levels=c("known (train)","known (val)","never seen")))
ts <- m %>% filter(set=="never seen", !is.na(score)) %>%
  mutate(worst=pmin(f1,f2,f3,f4,f5), best=pmax(f1,f2,f3,f4,f5))
lowest <- ts %>% group_by(term) %>% slice_min(score, n=5) %>% ungroup()
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
  legend.key.size=unit(10,"pt"), plot.title=element_text(face="bold", size=9.5),
  plot.subtitle=element_text(size=7.4, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA))

p1 <- ggplot() +
  geom_point(data=m %>% filter(set!="never seen"), aes(PC1, PC2, shape=set),
             colour="grey72", size=1.3, alpha=0.75) +
  geom_point(data=ts, aes(PC1, PC2, fill=score), shape=21, size=3, stroke=0.4, colour="grey20") +
  geom_text_repel(data=lowest, aes(PC1, PC2, label=label), size=2.3, max.overlaps=20,
                  min.segment.length=0, segment.size=0.25, colour="grey25") +
  facet_wrap(~term, nrow=1) +
  scale_fill_viridis_c(option="plasma", name="window score\n(mean of 50)", trans="sqrt") +
  scale_shape_manual(values=c(`known (train)`=3, `known (val)`=4), name=NULL) +
  labs(title=sprintf("Peptide-end windows in feature space (PC1 %.1f%%, PC2 %.1f%%)", ve[1], ve[2]),
       subtitle="each window is its 36 x 26 channels flattened to 936 numbers; grey = the ends the models trained and validated on",
       x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]))

p2 <- ggplot(ts, aes(d_train, score)) +
  geom_point(aes(fill=score), shape=21, size=2.6, stroke=0.3, colour="grey20") +
  geom_smooth(method="lm", se=FALSE, linewidth=0.5, colour="grey45") +
  geom_text_repel(data=lowest, aes(label=label), size=2.3, max.overlaps=20,
                  min.segment.length=0, segment.size=0.25, colour="grey25") +
  facet_wrap(~term, nrow=1) + scale_y_log10() +
  scale_fill_viridis_c(option="plasma", guide="none", trans="sqrt") +
  labs(title="Is a low score explained by being unlike the training windows?",
       subtitle="distance from the centre of the training windows, in the first 10 PCs",
       x="distance to the training centroid", y="window score (log)")

p3 <- ts %>% select(term, label, f1:f5) %>%
  pivot_longer(f1:f5, names_to="fold", values_to="s") %>%
  group_by(term, label) %>% mutate(med=median(s)) %>% ungroup() %>%
  ggplot(aes(reorder(label, med), s)) +
  geom_line(aes(group=label), colour="grey70", linewidth=0.3) +
  geom_point(aes(colour=fold), size=1.3) +
  facet_wrap(~term, scales="free_x", nrow=1) + scale_y_log10() +
  scale_colour_brewer(palette="Set2", name=NULL) +
  labs(title="Every fold's score for every never-seen peptide",
       subtitle="a peptide low in ALL five is missed by the method, not by one unlucky fold",
       x=NULL, y="score (log)") +
  theme(axis.text.x=element_text(angle=60, hjust=1, size=5.4))

p <- ((p1 & th) / (p2 & th) / (p3 & th)) + plot_layout(heights=c(1.25,0.85,1)) +
  plot_annotation(title="The 37 never-seen peptides: do the ones that score badly look unusual?",
    theme=theme(plot.title=element_text(face="bold", size=12)))
ggsave(file.path(A2,"pepend_cv5_pca.png"), p, width=12.5, height=13.5, dpi=185)
ggsave(file.path(A2,"pepend_cv5_pca.pdf"), p, width=12.5, height=13.5)
write.csv(ts %>% arrange(score) %>% select(term, gene, label, anchor, score, f1:f5, PC1, PC2, d_train),
          file.path(A2,"pepend_cv5_neverseen_pca.csv"), row.names=FALSE)
message("wrote ", file.path(A2,"pepend_cv5_pca.png"))
