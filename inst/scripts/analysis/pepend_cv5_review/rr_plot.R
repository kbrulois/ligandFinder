suppressMessages({library(tidyverse); library(patchwork)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
X <- readRDS(file.path(S2, "rr_curves.rds")); A <- X$A; B <- X$B
tl <- c(N = "N terminus", C = "C terminus")
A <- A %>% mutate(term = factor(tl[term], levels = tl),
                  set = factor(set, levels = c("cross-validated", "never seen: peptides",
                                               "never seen: large products")))
B <- B %>% mutate(term = factor(tl[term], levels = tl),
                  group = factor(group, levels = c("dibasic", "G + dibasic", "monobasic"),
                                 labels = c("dibasic", "G + dibasic (amidation)", "monobasic")))
cA <- c("cross-validated" = "#197EC0", "never seen: peptides" = "#C2663A",
        "never seen: large products" = "#6B7784")
cB <- c("dibasic" = "#197EC0", "G + dibasic (amidation)" = "#C2663A", "monobasic" = "#1B7837")
labA <- A %>% distinct(term, set, n) %>% mutate(t = sprintf("%s (n=%d)", set, n))
labB <- B %>% distinct(term, group, n, pool) %>%
  mutate(t = sprintf("%s (n=%d of %s)", group, n, format(pool, big.mark = ",")))
brk <- c(1, 10, 100, 1e3, 1e4, 1e5, 1e6)
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(10, "pt"), strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))
mk <- function(d, cols, key, title, sub, lab) {
  ggplot(d, aes(N, recall, colour = .data[[key]])) +
    geom_hline(yintercept = c(0.5, 0.8), linetype = 3, linewidth = 0.3, colour = "grey60") +
    geom_step(linewidth = 0.6, direction = "hv") +
    facet_wrap(~ term, nrow = 1) +
    scale_colour_manual(values = cols, name = NULL,
                        labels = \(x) lab$t[match(x, as.character(lab[[key]]))]) +
    scale_x_log10(breaks = brk, labels = c("1","10","100","1k","10k","100k","1M")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    annotation_logticks(sides = "b", short = unit(1,"pt"), mid = unit(2,"pt"),
                        long = unit(3,"pt"), linewidth = 0.25) +
    labs(title = title, subtitle = sub, x = "top N anchors retrieved (log)", y = "recall") + th
}
p1 <- mk(A, cA, "set", "Recall against how many anchors you take",
         paste("negatives = every mature anchor that is not an annotated end anywhere;",
               "a known end is ranked by the fold that HELD IT OUT, the never-seen sets by all 50"),
         labA %>% distinct(set, .keep_all = TRUE))
p2 <- mk(B, cB, "group", "The never-seen peptides, split by what flanks the end",
         paste("negatives are anchors of the SAME motif class, so recall reaches 100% in every",
               "group and the classes compare despite very different pool sizes"),
         labB %>% distinct(group, .keep_all = TRUE))
p <- (p1 / p2) + plot_annotation(
  title = "xgb window head: recall as a function of how deep you go into the ranking",
  subtitle = "dotted lines mark 50% and 80% recall",
  theme = theme(plot.title = element_text(face = "bold", size = 11),
                plot.subtitle = element_text(size = 8, colour = "grey30")))
f <- file.path(A2, "pepend_recall_rank.svg")
ggsave(f, p, width = 11, height = 8.2)
ggsave(sub("[.]svg$", ".png", f), p, width = 11, height = 8.2, dpi = 190)
## where each curve crosses 50% / 80%
cross <- function(d, key) d %>% group_by(term, .data[[key]]) %>%
  summarise(N50 = suppressWarnings(min(N[recall >= 0.5])), N80 = suppressWarnings(min(N[recall >= 0.8])), .groups="drop")
cat("\n-- anchors needed for 50% / 80% recall --\n")
print(as.data.frame(cross(A, "set")), row.names = FALSE)
print(as.data.frame(cross(B, "group")), row.names = FALSE)
message("wrote ", f)
