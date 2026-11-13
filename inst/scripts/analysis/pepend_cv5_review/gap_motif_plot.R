suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
R <- as_tibble(readRDS(file.path(S2,"gap_by_motif.rds"))) %>%
  mutate(term = factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
         set  = factor(set, levels=c("training","validation","never seen")),
         motif = factor(motif, levels=c("terminus","G + dibasic","dibasic","monobasic","other")))
cols <- c(training="#7FA66B", validation="#197EC0", `never seen`="#C24E5C")
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(),
  legend.position="bottom", legend.key.size=unit(10,"pt"),
  plot.title=element_text(face="bold", size=9.5),
  plot.subtitle=element_text(size=7.4, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))

p1 <- ggplot(R, aes(motif, 100*pct_rank, fill=set)) +
  geom_col(position=position_dodge2(width=0.8, preserve="single"), width=0.72) +
  geom_text(aes(label=sprintf("%.2f", 100*pct_rank)),
            position=position_dodge2(width=0.8, preserve="single"), vjust=-0.35, size=2.2) +
  facet_wrap(~term, nrow=1) + scale_fill_manual(values=cols, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Median percentile rank of the end, by motif class",
       subtitle="lower is better: the fraction of the ~2.98M anchors that outscore the true end",
       x=NULL, y="% of anchors ranked above")

p2 <- ggplot(R, aes(motif, 100*top10k, fill=set)) +
  geom_col(position=position_dodge2(width=0.8, preserve="single"), width=0.72) +
  geom_text(aes(label=sprintf("%.0f%%\nn=%d", 100*top10k, n)),
            position=position_dodge2(width=0.8, preserve="single"), vjust=-0.2, size=2.1, lineheight=0.9) +
  facet_wrap(~term, nrow=1) + scale_fill_manual(values=cols, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.28))) +
  labs(title="Share of ends recovered inside a 10,000-locus budget",
       subtitle="higher is better; n is the number of ends in that cell -- several are 2-5, so read those as anecdotes",
       x=NULL, y="% found")

p <- ((p1 & th) / (p2 & th)) +
  plot_annotation(title="Where the never-seen gap actually falls",
    subtitle="Current trunks, one shared negative pool. Splitting the never-seen set by motif shows the gap is NOT a composition artefact:\non C, `terminus` is the BEST never-seen class, not the worst, and the cliff sits in dibasic / monobasic / other.",
    theme=theme(plot.title=element_text(face="bold", size=12),
                plot.subtitle=element_text(size=8.2, colour="grey30")))
ggsave(file.path(A2,"pepend_gap_by_motif.png"), p, width=11.5, height=8.2, dpi=190)
ggsave(file.path(A2,"pepend_gap_by_motif.pdf"), p, width=11.5, height=8.2)
message("wrote ", file.path(A2,"pepend_gap_by_motif.png"))
