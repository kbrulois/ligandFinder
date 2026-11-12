## Segment (or anchor) classifier against threshold T, curated arm vs another
## arm, one curve per arm in each set x terminus x head panel (Kevin,
## 2026-10-24; F1 and the anchor level added 2026-11-02). Reads the
## per-threshold confusion tables of segments_roc.R (LEVEL=segment) or
## anchors_roc.R (LEVEL=anchor).
##
##   METRIC=auc_T Rscript segment_T_sweep_plot.R   # (1 + TPR - FPR) / 2, the original figure
##   METRIC=F1    Rscript segment_T_sweep_plot.R   # 2 P R / (P + R)
##   LEVEL=anchor METRIC=F1 Rscript segment_T_sweep_plot.R   # each anchor a call, not each segment
##   ARM2=esmonly36 TOL=4 (defaults)
##   ARM1=curated_AAswap ARM2=esmonly36   # two arms, neither the curated model
##
## F1 is the detection F1 of segments_roc.R's own definitions: precision over
## the positive SEGMENTS at T (centre within TOL of an evaluation end), recall
## (TPR) over the ENDS. Where T leaves no positive segment, precision is
## undefined and recall is 0: F1 = 0. Ends and segments are different units, so
## F1 is computed from P and R, not as 2 TP / (2 TP + FP + FN).
suppressMessages({library(data.table); library(ggplot2)})
METRIC <- Sys.getenv("METRIC", "auc_T"); stopifnot(METRIC %in% c("auc_T", "F1"))
LEVEL  <- Sys.getenv("LEVEL", "segment"); stopifnot(LEVEL %in% c("segment", "anchor"))
ARM2 <- Sys.getenv("ARM2", "esmonly36"); TOL <- as.integer(Sys.getenv("TOL", "4"))
LABELS <- c(esmonly36 = "ESM-C only (36 PCA)", curated_AAswap = "AAswap", curated_wAAid = "curated + AA ID")
lab2 <- if (ARM2 %in% names(LABELS)) LABELS[[ARM2]] else ARM2
## ARM1 is the first curve; "" (the default) is the curated model, whose tables carry no arm token
ARM1 <- Sys.getenv("ARM1", "")
lab1 <- if (!nzchar(ARM1)) "curated (26 hand-built)" else if (ARM1 %in% names(LABELS)) LABELS[[ARM1]] else ARM1
stopifnot(lab1 != lab2)
H <- path.expand("~/AF2_analysis/lf_headcmp")
cur <- fread(file.path(H, if (nzchar(ARM1)) sprintf("%s_confusion_%s_tol%d.csv", LEVEL, ARM1, TOL)
                          else sprintf("%s_confusion_tol%d.csv", LEVEL, TOL)))[, arm := lab1]
oth <- fread(file.path(H, sprintf("%s_confusion_%s_tol%d.csv", LEVEL, ARM2, TOL)))[, arm := lab2]
D <- rbind(cur, oth, fill = TRUE)
D[, panel := paste0(set, " — ", term, " ", head)]
D[, F1 := fifelse(is.na(precision) | precision + TPR == 0, 0, 2 * precision * TPR / (precision + TPR))]
D[, y := get(METRIC)]
opt <- D[, .SD[which.max(y)], by = .(arm, panel)]
## the lower of a panel's two optima is labelled under its ring, so labels at
## the same threshold do not overprint
opt[, vj := fifelse(y == min(y) & .N > 1, 1.9, -0.9), by = panel]
cols <- setNames(c("#3c63a0", "#c96a35"), c(lab1, lab2))
p <- ggplot(D, aes(thr, y, colour = arm)) +
  (if (METRIC == "auc_T") geom_hline(yintercept = 0.5, linewidth = .3, colour = "grey75", linetype = "22")) +
  geom_vline(xintercept = 0.5, linewidth = .3, colour = "grey85") +
  geom_line(linewidth = .7) +
  geom_point(data = opt, size = 2.2, shape = 21, fill = "white", stroke = 1.1) +
  ## on a white backing, so a steep curve running under a label does not hide it
  (if (METRIC == "F1") geom_label(data = opt, aes(label = sprintf("%.2g", y), vjust = vj), size = 2.7,
                                  fill = scales::alpha("white", 0.85), border.colour = NA,
                                  label.padding = grid::unit(0.08, "lines"), show.legend = FALSE)) +
  facet_wrap(~ panel, ncol = 4, scales = if (METRIC == "auc_T" && LEVEL == "segment") "fixed" else "free_y") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1)) +
  (if (METRIC == "auc_T" && LEVEL == "segment") coord_cartesian(ylim = c(0.45, 0.82))
   else scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)))) +
  labs(title = sprintf("%s %s against threshold, tolerance %d",
                       if (LEVEL == "segment") "Segment" else "Anchor", METRIC, TOL),
       subtitle = if (METRIC == "auc_T" && LEVEL == "segment")
         "Ringed point = each curve's own optimum. Grey line = T 0.5. auc_T = (1 + TPR - FPR) / 2"
       else if (METRIC == "auc_T")
         paste0("Ringed point = each curve's own optimum. Grey line = T 0.5. Each panel has its own y-axis. ",
                "auc_T = (1 + TPR - FPR) / 2 over anchors")
       else paste0("Ringed point = each curve's own optimum, labelled with its F1. Grey line = T 0.5. Each panel has its own y-axis.\n",
                   "F1 = 2 precision recall / (precision + recall); ",
                   if (LEVEL == "segment")
                     paste0("precision over positive segments (centre within ", TOL, " of an end), recall over ends")
                   else paste0("precision over called anchors (score >= T; correct within ", TOL,
                               " of an end), recall over ends (any called anchor within ", TOL, ")")),
       x = "threshold T", y = paste(LEVEL, METRIC)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold"))
out <- file.path(H, sprintf("%s_T_sweep%s_%s_vs_%s_tol%d.png", LEVEL,
                            if (METRIC == "auc_T") "" else paste0("_", METRIC),
                            if (nzchar(ARM1)) ARM1 else "curated", ARM2, TOL))
out <- Sys.getenv("OUT", out)
ggsave(out, p, width = 13, height = 7.5, dpi = 165, bg = "white")
cat("wrote", out, "\n")
cat(sprintf("optimum by panel (%s):\n", METRIC))
print(dcast(opt[, .(panel, arm, v = sprintf("%.3f @ T %.3f", y, thr))], panel ~ arm, value.var = "v"), row.names = FALSE)
