suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2,"xgb_redo.rds")); R <- as_tibble(X$R); O <- as_tibble(X$O)
hl <- c(cnn="CNN attention head", xgb="xgb window head")
tl <- function(d) d %>% mutate(term=factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
                               head=factor(hl[head], levels=hl),
                               set=factor(set, levels=c("training","validation","never seen")))
R <- tl(R); O <- tl(O)
cols <- c(training="#7FA66B", validation="#197EC0", `never seen`="#C24E5C")
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
  legend.key.size=unit(10,"pt"), plot.title=element_text(face="bold", size=9.5),
  plot.subtitle=element_text(size=7.4, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))

p1 <- ggplot(R, aes(set, 100*pct_rank, fill=set)) + geom_col(width=0.65) +
  geom_text(aes(label=sprintf("%.2f%%", 100*pct_rank)), vjust=-0.35, size=2.5) +
  facet_grid(term ~ head) + scale_fill_manual(values=cols, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(title="Median percentile rank of the end (lower is better)",
       subtitle="training is IN-SAMPLE for the xgb head -- it is fitted on those rows, so read only validation and never seen",
       x=NULL, y="% of 2.98M anchors above the end")

sm <- O %>% pivot_longer(c(top_in_prec, near_peak), names_to="measure", values_to="offset") %>%
  mutate(measure = recode(measure, top_in_prec="protein's TOP call is the end",
                                   near_peak="a PEAK sits on the end"),
         measure = factor(measure, levels=c("protein's TOP call is the end","a PEAK sits on the end"))) %>%
  group_by(term, head, set, measure) %>%
  summarise(v = mean(abs(offset) <= 2, na.rm=TRUE), .groups="drop")
p2 <- ggplot(sm, aes(set, 100*v, fill=set)) + geom_col(width=0.65) +
  geom_text(aes(label=sprintf("%.0f%%", 100*v)), vjust=-0.35, size=2.5) +
  facet_grid(measure + term ~ head) + scale_fill_manual(values=cols, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.2)), limits=c(0,100)) +
  labs(title="Within 2 residues of the true end", subtitle="higher is better",
       x=NULL, y="% of ends")

p <- ((p1 & th) / (p2 & th)) + plot_layout(heights=c(1, 2)) +
  plot_annotation(title="CNN attention head vs xgb window head, by how much the model has seen",
    subtitle="Same trunks, same 2.98M anchors, both heads from lf_winxgb_fullscan_{C,N}_near.npz (its cnn column reproduces the scan to 1.8e-07).\nThe xgb head is fitted on the training rows, so its training bars are in-sample and not a performance estimate.",
    theme=theme(plot.title=element_text(face="bold", size=12), plot.subtitle=element_text(size=8.2, colour="grey30")))
ggsave(file.path(A2,"pepend_cnn_vs_xgb_by_set.png"), p, width=10.5, height=12, dpi=190)
ggsave(file.path(A2,"pepend_cnn_vs_xgb_by_set.pdf"), p, width=10.5, height=12)
write.csv(R, file.path(A2,"pepend_cnn_vs_xgb_by_set.csv"), row.names=FALSE)
message("wrote ", file.path(A2,"pepend_cnn_vs_xgb_by_set.png"))
