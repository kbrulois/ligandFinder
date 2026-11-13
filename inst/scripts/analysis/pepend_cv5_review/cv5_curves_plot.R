suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2,"cv5_curves.rds"))
lv <- c("all","dibasic","monobasic","G + dibasic")
fx <- function(d) d %>% mutate(
  group = factor(group, levels=lv),
  term  = factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
  scorer= factor(scorer, levels=c("mean of 50","max of folds","min of folds")))
PR <- fx(X$pr); RR <- fx(X$rr)
lab <- RR %>% distinct(term, group, n) %>% mutate(tag=sprintf("n = %d", n))
cols <- c(`mean of 50`="#197EC0", `max of folds`="#C24E5C", `min of folds`="#7FA66B")
th <- theme_bw(base_size=9) + theme(panel.grid.minor=element_blank(), legend.position="bottom",
  legend.key.size=unit(10,"pt"), plot.title=element_text(face="bold", size=9.5),
  plot.subtitle=element_text(size=7.4, colour="grey35"),
  strip.background=element_rect(fill="grey94", colour=NA), panel.spacing=unit(7,"pt"))

p1 <- ggplot(PR, aes(recall, precision, colour=scorer)) +
  geom_step(linewidth=0.6, direction="vh") +
  facet_grid(term ~ group, scales="free_y") +
  geom_text(data=lab, aes(x=0.5, y=Inf, label=tag), inherit.aes=FALSE,
            vjust=1.6, size=2.6, colour="grey30") +
  scale_colour_manual(values=cols, name=NULL) +
  scale_y_log10(labels=scales::label_number(accuracy=0.001)) +
  labs(title="Precision-recall within each motif class",
       subtitle="positives = the 37 peptides held out of every fold; negatives = anchors of the SAME class, so recall reaches 100% in every panel",
       x="recall", y="precision (log)")

p2 <- ggplot(RR, aes(N, recall, colour=scorer)) +
  geom_line(linewidth=0.6) +
  facet_grid(term ~ group) +
  geom_text(data=lab, aes(x=10, y=1.02, label=tag), inherit.aes=FALSE,
            vjust=1, hjust=0, size=2.6, colour="grey30") +
  scale_colour_manual(values=cols, name=NULL) +
  scale_x_log10(breaks=c(1,10,100,1e3,1e4,1e5,1e6),
                labels=c("1","10","100","1k","10k","100k","1M")) +
  scale_y_continuous(labels=scales::percent, limits=c(0,1.05)) +
  labs(title="Recall at a budget of anchors",
       subtitle="how many of that class you follow, against how many real ends you have found -- the comparable measure across classes",
       x="anchors of that class followed (log)", y="ends found")

p <- ((p1 & th) / (p2 & th)) +
  plot_annotation(title="Ensemble rules by motif class, on peptides no fold ever trained on",
    subtitle="Three ways to turn the 5 fold-means into one score: their average, their max at each anchor, or their min.\nG + dibasic is C-terminal only -- an N-terminal amidation motif does not exist.",
    theme=theme(plot.title=element_text(face="bold", size=12), plot.subtitle=element_text(size=8.2, colour="grey30")))
ggsave(file.path(A2,"pepend_cv5_motif_curves.png"), p, width=13, height=9.5, dpi=190)
ggsave(file.path(A2,"pepend_cv5_motif_curves.pdf"), p, width=13, height=9.5)
message("wrote ", file.path(A2,"pepend_cv5_motif_curves.png"))
