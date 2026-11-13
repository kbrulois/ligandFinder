suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
O <- as_tibble(readRDS(file.path(S2,"offset_sets.rds"))) %>%
  pivot_longer(c(top_in_prec, near_peak), names_to="measure", values_to="offset") %>%
  mutate(term = factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
         set = factor(set, levels=c("training","validation","never seen")),
         measure = recode(measure,
            top_in_prec = "offset of the protein's TOP-scoring anchor",
            near_peak   = "distance to the NEAREST local peak"),
         measure = factor(measure, levels=c("offset of the protein's TOP-scoring anchor",
                                            "distance to the NEAREST local peak")),
         o = pmax(pmin(offset, 50), -50))
cols <- c(training="#7FA66B", validation="#197EC0", `never seen`="#C24E5C")
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(),
  legend.position="bottom", legend.key.size=unit(10,"pt"),
  plot.title=element_text(face="bold", size=9.5), plot.subtitle=element_text(size=7.4, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))

p1 <- ggplot(O, aes(o, fill=set)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey35") +
  geom_histogram(binwidth=3, position="identity", alpha=0.5, colour=NA) +
  facet_grid(measure ~ term, scales="free_y") +
  scale_fill_manual(values=cols, name=NULL) +
  labs(title="Distance to the peptide end", subtitle="0 = called exactly on the end; clipped at +/-50 residues",
       x="offset (residues; negative = called too early)", y="ends")

s <- O %>% group_by(term, set, measure) %>%
  summarise(`exactly 0`=mean(offset==0, na.rm=TRUE), `within 2`=mean(abs(offset)<=2, na.rm=TRUE),
            `within 5`=mean(abs(offset)<=5, na.rm=TRUE), n=sum(!is.na(offset)), .groups="drop") %>%
  pivot_longer(c(`exactly 0`,`within 2`,`within 5`), names_to="k", values_to="v")
p2 <- ggplot(s, aes(k, 100*v, fill=set)) +
  geom_col(position=position_dodge(width=0.78), width=0.7) +
  geom_text(aes(label=sprintf("%.0f", 100*v)), position=position_dodge(width=0.78), vjust=-0.3, size=2.3) +
  facet_grid(measure ~ term) + scale_fill_manual(values=cols, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="How close, as a percentage of ends",
       subtitle="n = 57/29/37 (C) and 59/33/37 (N) for training / validation / never seen",
       x=NULL, y="% of ends")

p <- ((p1 & th) / (p2 & th)) +
  plot_annotation(title="Localising the peptide end: training vs held out vs never seen",
    subtitle="Current trunks. Top row of each panel: is the protein's single best call the end? Bottom row: is there a peak on the end at all, wherever it ranks?",
    theme=theme(plot.title=element_text(face="bold", size=12), plot.subtitle=element_text(size=8.2, colour="grey30")))
ggsave(file.path(A2,"pepend_offset_by_set.png"), p, width=11.5, height=9.6, dpi=190)
ggsave(file.path(A2,"pepend_offset_by_set.pdf"), p, width=11.5, height=9.6)
message("wrote ", file.path(A2,"pepend_offset_by_set.png"))
