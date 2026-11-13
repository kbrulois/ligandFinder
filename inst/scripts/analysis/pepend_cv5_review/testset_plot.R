suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
F <- readRDS(file.path(S2, "testset_fair.rds")); O <- readRDS(file.path(S2, "testset_roc.rds"))$O
cols <- c(old = "#197EC0", rebuilt = "#C24E5C")
lab <- c("all negatives" = "all negatives (2.98M)",
         "complete-data negatives only" = "complete-data negatives only (1.4M)")
C <- F$C %>% mutate(pool = factor(lab[pool], levels = lab),
                    term = factor(paste(term, "terminus"), levels = c("N terminus","C terminus")))
R <- F$R %>% mutate(pool = factor(lab[pool], levels = lab),
                    term = factor(paste(term, "terminus"), levels = c("N terminus","C terminus")))
th <- theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
  legend.position = "bottom", legend.key.size = unit(10,"pt"),
  plot.title = element_text(face="bold", size=9.5),
  plot.subtitle = element_text(size=7.4, colour="grey35"),
  strip.background = element_rect(fill="grey94", colour=NA), panel.spacing = unit(8,"pt"))

p1 <- ggplot(C, aes(fpr, rec, colour = arm)) +
  geom_abline(slope=1, intercept=0, linewidth=0.3, colour="grey70", linetype="dashed") +
  geom_line(linewidth=0.6) + facet_grid(term ~ pool) +
  scale_colour_manual(values=cols, name=NULL) +
  geom_text(data = R, aes(x=0.62, y=ifelse(arm=="old", 0.18, 0.08),
            label=sprintf("%s  AUC %.3f", arm, roc_auc)), size=2.5, hjust=0, show.legend=FALSE) +
  labs(title="ROC: 37 peptides neither model has ever seen",
       subtitle="positives are the exact end anchors; every OTHER known end (3,039) is held out rather than counted as a mistake",
       x="false positive rate", y="recall of the 37 test ends")

p2 <- ggplot(C, aes(rec, prec, colour = arm)) + geom_line(linewidth=0.6) +
  facet_grid(term ~ pool) + scale_colour_manual(values=cols, name=NULL) +
  scale_y_log10(labels = scales::label_number(accuracy=0.0001)) +
  geom_text(data = R, aes(x=0.3, y=ifelse(arm=="old", 0.5, 0.15),
            label=sprintf("%s  AP %.4f (%.0fx)", arm, pr_auc, lift)), size=2.5, hjust=0, show.legend=FALSE) +
  labs(title="Precision-recall, same test set",
       subtitle="AP = mean precision at each positive; lift is against the prevalence of that pool",
       x="recall", y="precision (log)")

off <- O %>% mutate(term = factor(paste(term,"terminus"), levels=c("N terminus","C terminus")),
                    o = pmax(pmin(offset, 60), -60))
p3 <- ggplot(off, aes(o, fill = arm)) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey40") +
  geom_histogram(binwidth=4, position="identity", alpha=0.55, colour=NA) +
  facet_wrap(~ term, nrow=1) + scale_fill_manual(values=cols, name=NULL) +
  labs(title="Distance from the precursor's best-scoring anchor to the true end",
       subtitle="0 = the top-scoring anchor in that protein IS the peptide end; clipped at +/-60",
       x="offset (residues; negative = called too early)", y="peptides")

s <- off %>% group_by(arm, term) %>%
  summarise(`exactly 0` = sum(offset==0), `within 2` = sum(abs(offset)<=2),
            `within 5` = sum(abs(offset)<=5), `top-1 in gene` = sum(rank_in_gene==1),
            `top-10 in gene` = sum(rank_in_gene<=10), .groups="drop") %>%
  pivot_longer(-c(arm, term), names_to="k", values_to="n") %>%
  mutate(k = factor(k, levels=c("exactly 0","within 2","within 5","top-1 in gene","top-10 in gene")))
p4 <- ggplot(s, aes(k, n, fill=arm)) +
  geom_col(position=position_dodge(width=0.75), width=0.7) +
  geom_text(aes(label=n), position=position_dodge(width=0.75), vjust=-0.35, size=2.4) +
  facet_wrap(~ term, nrow=1) + scale_fill_manual(values=cols, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(title="How often the end is actually found", subtitle="of 37 N ends / 36 C ends",
       x=NULL, y="peptides")

p <- ((p1 & th) / (p2 & th) / ((p3 & th) | (p4 & th))) +
  plot_annotation(
    title = "Peptide-end prediction on 37 peptides held out of every training fold",
    subtitle = "Never a positive, never a negative: their genes contribute no training window and sit in the 250-gene banned list.\nLeft column is the pool as scanned; right column keeps only negatives whose 36-residue window had complete feature data, where coverage cannot decide.",
    theme = theme(plot.title = element_text(face="bold", size=12),
                  plot.subtitle = element_text(size=8.2, colour="grey30")))
ggsave(file.path(A2, "pepend_testset_roc_pr.png"), p, width=12.5, height=13.5, dpi=190)
ggsave(file.path(A2, "pepend_testset_roc_pr.pdf"), p, width=12.5, height=13.5)
write.csv(F$R, file.path(A2, "pepend_testset_summary.csv"), row.names=FALSE)
message("wrote ", file.path(A2, "pepend_testset_roc_pr.png"))
