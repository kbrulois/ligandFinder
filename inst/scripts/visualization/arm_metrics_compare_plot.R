## Arm-vs-arm full-proteome metrics as one figure + table (Kevin, 2026-11-12:
## "compare cv model metrics with the new msa versus the old. make table or
## figure", then "add f1 score and roc auc").
##
## Reads metrics_final.R's two outputs from one HEADCMP directory -- every arm
## scored on ONE never-seen set, at the anchors the end tables carry -- and draws,
## per terminus x evaluation set x window head, four panels: normalised PR AUC,
## ROC AUC, best F1 and anchors to scan for 50% recall. The LAST arm in ARMS is
## the one under test (red), the one before it its comparator (dark grey, joined
## to it by a line), any earlier ones references (light grey).
##
## Best F1 comes from the recall-vs-rank curve, so it is the same exact-anchor
## ranking as the other panels: at a cut with N non-end anchors ranked above it
## (metrics_final.R's rank: each positive's rank among the anchors that are no
## annotated end), TP = recall x n, FP = N, FN = n - TP, so F1 = 2TP/(TP + N + n),
## maximised over the curve's rank grid. No +/- tolerance, unlike the segment and
## anchor F1 of segment_T_sweep_plot.R -- the two are not comparable -- and small
## on the never-seen sets, where ends are ~14 per million anchors.
##
##   LC_ALL=en_US.UTF-8 Rscript arm_metrics_compare_plot.R
##   env: HEADCMP  ~/AF2_analysis/lf_headcmp/unseen_padg (final_metrics.csv, final_recall_rank.csv)
##        ARMS     curated_AAswap,curated_AAswap_padg,curated_AAswap_msacons_padg
##        LABELS   one per arm, separated by ";" (labels may hold commas)
##        OUT      ~/AF2_analysis/lf_viz_20261112    NAME  cv_metrics_old_vs_new_conservation
##        TITLE    figure title
## Writes <OUT>/<NAME>.{png,pdf,csv}; the csv holds every metric for every arm.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
HC   <- path.expand(Sys.getenv("HEADCMP", "~/AF2_analysis/lf_headcmp/unseen_padg"))
ARMS <- strsplit(Sys.getenv("ARMS", "curated_AAswap,curated_AAswap_padg,curated_AAswap_msacons_padg"), ",")[[1]]
LABS <- trimws(strsplit(Sys.getenv("LABELS", paste("original windows, Aminode cons.", "window fix, OLD Aminode cons.",
                                                   "window fix, NEW MSA cons.", sep = ";")), ";")[[1]])
OUT  <- path.expand(Sys.getenv("OUT", "~/AF2_analysis/lf_viz_20261112"))
NAME <- Sys.getenv("NAME", "cv_metrics_old_vs_new_conservation")
TITLE <- Sys.getenv("TITLE", "Old (Aminode) vs new (MSA) conservation: full-proteome scan metrics")
stopifnot(length(ARMS) >= 2, length(LABS) == length(ARMS))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m <- fread(file.path(HC, "final_metrics.csv"))[trunk %in% ARMS]
r <- fread(file.path(HC, "final_recall_rank.csv"))[trunk %in% ARMS]
miss <- setdiff(ARMS, m$trunk)
if (length(miss)) stop("not in ", file.path(HC, "final_metrics.csv"), ": ", paste(miss, collapse = ", "), call. = FALSE)
## best F1 along the recall curve (header)
r <- r[m[, .(term, set, trunk, head, n)], on = .(term, set, trunk, head)]
r[, tp := recall * n][, f1 := 2 * tp / (tp + N + n)]
bf <- r[, .SD[which.max(f1)], by = .(term, set, trunk, head)][, .(term, set, trunk, head, best_f1 = f1, f1_at_N = N)]
m <- bf[m, on = .(term, set, trunk, head)]

focus <- ARMS[length(ARMS)]; comp <- ARMS[length(ARMS) - 1L]
m[, arm := factor(LABS[match(trunk, ARMS)], LABS)]
m[, head := factor(head, c("cnn", "xgb_pm"), c("CNN head", "xgb head"))]
m[, set := factor(set, c("cross-validated", "never seen: peptides", "never seen: large products"),
                  c("cross-validated\n(5-fold, held-out ends)", "never-seen genes:\npeptide ends",
                    "never-seen genes:\nlarge cleavage products"))]
m[, term := factor(term, c("N", "C"), c("N terminus", "C terminus"))]

## the table: every metric, one column per arm -- n too, since the arms' end
## sets can differ by a merge (86 vs 85 cross-validated C ends)
tab <- dcast(m[, .(term, set = gsub("\n", " ", set), head, n, trunk, norm_pr_auc, roc_auc, best_f1, n50, n80)],
             term + set + head ~ trunk, value.var = c("n", "norm_pr_auc", "roc_auc", "best_f1", "n50", "n80"))
fwrite(tab, file.path(OUT, paste0(NAME, ".csv")))

cols <- setNames(c(rep("grey70", length(ARMS) - 2L), "#7F8C8D", "#C0392B"), LABS)
panel <- function(var, title, xlab, fmt, xscale) {
  ggplot(m, aes(.data[[var]], head, colour = arm)) +
    geom_line(data = m[trunk %in% c(comp, focus)], aes(group = head), colour = "grey60", linewidth = 0.8) +
    geom_point(size = 3) +
    geom_text(data = m[trunk == focus], aes(label = fmt(.data[[var]])),
              colour = "#C0392B", vjust = -1, size = 2.9, show.legend = FALSE) +
    facet_grid(set ~ term, switch = "y") + xscale +
    scale_colour_manual(values = cols, name = NULL) +
    labs(x = xlab, y = NULL, title = title) +
    theme_minimal(base_size = 10.5) +
    theme(legend.position = "top", strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1), strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}
comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
## the ROC axis starts below the lowest value, rounded down to 0.05, so no point drops off it
roc_lo <- floor(min(m$roc_auc) * 20) / 20
p1 <- panel("norm_pr_auc", "Ranking quality: normalised PR AUC", "normalised PR AUC (higher = better)",
            function(x) sprintf("%.2f", x), scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)))
p2 <- panel("roc_auc", "ROC AUC", sprintf("ROC AUC (higher = better; axis starts at %.2f)", roc_lo),
            function(x) sprintf("%.3f", x), scale_x_continuous(limits = c(roc_lo, 1), breaks = seq(roc_lo, 1, 0.05)))
p3 <- panel("best_f1", "Best F1 over thresholds (anchor level)", "best F1 (log; higher = better)",
            function(x) sprintf("%.3f", x),
            scale_x_log10(breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
                          labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)))
p4 <- panel("n50", "Search cost", "anchors to scan for 50% recall (log; lower = better)",
            comma, scale_x_log10(labels = comma))
g <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") +
  plot_annotation(
    title = TITLE,
    subtitle = paste0("Red = ", LABS[length(LABS)], "; dark grey = ", LABS[length(LABS) - 1L], " (joined by a line)",
                      if (length(ARMS) > 2) paste0("; light grey = ", paste(LABS[seq_len(length(ARMS) - 2L)], collapse = ", "), ", for reference") else "", ".\n",
                      "Cross-validated: each known end scored by the fold that held it out. Never-seen: genes no fold trained on, scored by all 50 members. Candidates = every extracellular anchor of the proteome.\n",
                      "Best F1: exact end anchor, maximised over score thresholds (FP = non-end anchors ranked above the cut); small on never-seen sets because ends are ~14 per million anchors.")) &
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, paste0(NAME, ".png")), g, width = 16, height = 13, dpi = 200, bg = "white")
ggsave(file.path(OUT, paste0(NAME, ".pdf")), g, width = 16, height = 13)
message("wrote ", file.path(OUT, paste0(NAME, ".{png,pdf,csv}")))
