suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
G <- readRDS(file.path(S2,"gen_gap.rds")); R <- G$R; C <- G$C
cols <- c(training="#7FA66B", validation="#197EC0", `never seen`="#C24E5C")
tl <- function(d) d %>% mutate(term = factor(paste(term,"terminus"), levels=c("N terminus","C terminus")))
C <- tl(C); R <- tl(R)
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(),
  legend.position="bottom", legend.key.size=unit(10,"pt"),
  plot.title=element_text(face="bold", size=9.5),
  plot.subtitle=element_text(size=7.4, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA))

p1 <- ggplot(C, aes(fpr, rec, colour=set)) +
  geom_abline(slope=1, intercept=0, linewidth=0.3, colour="grey70", linetype="dashed") +
  geom_line(linewidth=0.65) + facet_wrap(~term, nrow=1) +
  scale_colour_manual(values=cols, name=NULL) +
  geom_text(data=R, aes(x=0.55, y=c(0.30,0.19,0.08)[as.integer(set)],
    label=sprintf("%s (n=%d)  AUC %.3f", set, n, roc_auc)), size=2.6, hjust=0, show.legend=FALSE) +
  labs(title="ROC by how much the model has seen of each peptide",
       subtitle="one model (the current trunks), one shared negative pool of ~2.98M anchors; only the positives differ",
       x="false positive rate", y="recall")

p2 <- ggplot(C, aes(rec, prec, colour=set)) + geom_line(linewidth=0.65) +
  facet_wrap(~term, nrow=1) + scale_colour_manual(values=cols, name=NULL) +
  scale_y_log10(labels=scales::label_number(accuracy=0.0001)) +
  labs(title="Precision-recall, same three sets",
       subtitle="prevalence is near-identical across the three (same negatives, 29-58 positives), so these are comparable",
       x="recall", y="precision (log)")

b <- R %>% select(term, set, pct_rank_median, `in top 10k`=top10k) %>%
  pivot_longer(-c(term,set)) %>%
  mutate(name = recode(name, pct_rank_median = "median percentile rank of the end",
                             `in top 10k` = "share of ends inside the top 10,000 anchors"))
p3 <- ggplot(b, aes(set, value, fill=set)) + geom_col(width=0.65) +
  geom_text(aes(label=ifelse(grepl("percentile", name), sprintf("%.2f%%", 100*value), sprintf("%.0f%%", 100*value))),
            vjust=-0.35, size=2.6) +
  facet_grid(name ~ term, scales="free_y", switch="y") +
  scale_fill_manual(values=cols, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(title="The gap that matters for a screen", subtitle="lower is better on top; higher is better below",
       x=NULL, y=NULL)

p <- ((p1 & th) / (p2 & th) / (p3 & th)) +
  plot_annotation(title="Peptide-end retrieval: training vs held out vs never seen",
    subtitle="training = ends the trunks were fitted on; validation = ends held out by the similarity-cluster split; never seen = 37 peptides whose genes contribute no training window at all and sit in the banned list, so they were never a positive and never a negative in any arm.",
    theme=theme(plot.title=element_text(face="bold", size=12),
                plot.subtitle=element_text(size=8.2, colour="grey30")))
ggsave(file.path(A2,"pepend_generalisation_gap.png"), p, width=11.5, height=12, dpi=190)
ggsave(file.path(A2,"pepend_generalisation_gap.pdf"), p, width=11.5, height=12)
write.csv(R, file.path(A2,"pepend_generalisation_gap.csv"), row.names=FALSE)
message("wrote ", file.path(A2,"pepend_generalisation_gap.png"))
