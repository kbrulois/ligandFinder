## The headline metrics of arm_metrics_compare_plot.R as a table (Kevin,
## 2026-11-13: "headline metrics ... check it and add a table"): normalised PR
## AUC, ROC AUC, best anchor F1 and anchors for 50% recall, per terminus, set and
## window head, one column per arm, the best arm of each row in bold.
##
## Same inputs and definitions as arm_metrics_compare_plot.R (metrics_final.R's
## final_metrics.csv and final_recall_rank.csv from ONE HEADCMP directory, so every
## arm is scored on the same never-seen set; best F1 from the recall-vs-rank curve).
##
##   LC_ALL=en_US.UTF-8 Rscript arm_compare_metrics_table.R
## env: HEADCMP ~/AF2_analysis/lf_headcmp/unseen_msa   ARMS curated_AAswap,msa,msa_gap
##      LABELS  one per arm, ";"-separated               OUT ~/AF2_analysis/lf_viz_20261113
##      NAME    cv_metrics_msa_arms_table
## Writes <OUT>/<NAME>.{png,csv}.
suppressMessages({library(data.table); library(gridExtra); library(grid)})
HC <- path.expand(Sys.getenv("HEADCMP", "~/AF2_analysis/lf_headcmp/unseen_msa"))
ARMS <- strsplit(Sys.getenv("ARMS", "curated_AAswap,msa,msa_gap"), ",")[[1]]
LABS <- trimws(strsplit(Sys.getenv("LABELS", "AAswap (site);msa;msa_gap"), ";")[[1]])
OUT <- path.expand(Sys.getenv("OUT", "~/AF2_analysis/lf_viz_20261113"))
NAME <- Sys.getenv("NAME", "cv_metrics_msa_arms_table")
stopifnot(length(LABS) == length(ARMS))

m <- fread(file.path(HC, "final_metrics.csv"))[trunk %in% ARMS]
r <- fread(file.path(HC, "final_recall_rank.csv"))[trunk %in% ARMS]
miss <- setdiff(ARMS, m$trunk)
if (length(miss)) stop("not in ", file.path(HC, "final_metrics.csv"), ": ", paste(miss, collapse = ", "), call. = FALSE)
r <- r[m[, .(term, set, trunk, head, n)], on = .(term, set, trunk, head)]
r[, f1 := 2 * recall * n / (recall * n + N + n)]
m <- r[, .(best_f1 = max(f1)), by = .(term, set, trunk, head)][m, on = .(term, set, trunk, head)]
m[, arm := factor(LABS[match(trunk, ARMS)], LABS)]
m[, set := factor(set, c("cross-validated", "never seen: peptides", "never seen: large products"),
                  c("cross-validated", "never seen: peptides", "never seen: large products"))]
m[, term := factor(term, c("N", "C"))]
m[, head := factor(head, c("cnn", "xgb_pm"), c("CNN", "xgb"))]
setorder(m, term, set, head, arm)

met <- c(norm_pr_auc = "norm. PR AUC", roc_auc = "ROC AUC", best_f1 = "best F1", n50 = "n50")
hi_better <- c(norm_pr_auc = TRUE, roc_auc = TRUE, best_f1 = TRUE, n50 = FALSE)
fmt <- list(norm_pr_auc = function(x) sprintf("%.3f", x), roc_auc = function(x) sprintf("%.4f", x),
            best_f1 = function(x) sprintf("%.3f", x), n50 = function(x) format(x, big.mark = ",", trim = TRUE))
w <- dcast(m, term + set + head ~ arm, value.var = names(met))
nmax <- m[, .(n = paste(unique(n), collapse = "/")), by = .(term, set, head)]
w <- nmax[w, on = .(term, set, head)]
fwrite(w, file.path(OUT, paste0(NAME, ".csv")))

## the displayed table, best of each row and metric in bold
disp <- w[, .(terminus = as.character(term), set = as.character(set), head = as.character(head), n)]
bold <- matrix(FALSE, nrow(w), 4 + length(met) * length(LABS))
j <- 5L
for (v in names(met)) {
  cols <- paste(v, LABS, sep = "_")
  X <- as.matrix(w[, cols, with = FALSE])
  best <- if (hi_better[[v]]) apply(X, 1, max) else apply(X, 1, min)
  for (a in seq_along(LABS)) {
    disp[[paste0(met[[v]], "\n", LABS[a])]] <- fmt[[v]](X[, a])
    bold[, j] <- X[, a] == best & rowSums(X == best) < length(LABS)
    j <- j + 1L
  }
}
ff <- ifelse(bold, "bold", "plain")
## blank repeated terminus / set labels, so the groups read down the page
disp$terminus[duplicated(w$term)] <- ""
disp$set[duplicated(w[, .(term, set)])] <- ""
th <- ttheme_minimal(base_size = 9.5,
                     core = list(fg_params = list(fontface = as.vector(ff), hjust = 1, x = 0.95),
                                 bg_params = list(fill = rep(rep(c("grey97", "white"), each = 2), length.out = nrow(disp)))),
                     colhead = list(fg_params = list(fontface = "bold", hjust = 1, x = 0.95)))
tg <- tableGrob(disp, rows = NULL, theme = th)
## a rule under the header and between the termini
tg <- gtable::gtable_add_grob(tg, segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(lwd = 1)),
                              t = 1, l = 1, r = ncol(tg))
title <- textGrob(sprintf("Full-proteome scan metrics: %s", paste(LABS, collapse = " vs ")), x = 0.01, hjust = 0,
                  gp = gpar(fontface = "bold", fontsize = 13))
sub <- textGrob(paste0("Bold = the best arm of the row (higher is better, except n50 = anchors to scan for 50% recall). ",
                       "Cross-validated: each known end scored by the fold that held it out; never seen: genes no fold trained on, all 50 members.\n",
                       "Every arm on the same ends (", basename(HC), "). Best F1 at the exact end anchor over score thresholds; ",
                       "ROC / normalised PR (1% prevalence) over every extracellular anchor of the proteome."),
                x = 0.01, hjust = 0, gp = gpar(fontsize = 9, col = "grey25"))
g <- arrangeGrob(title, sub, tg, ncol = 1, heights = unit.c(unit(1.6, "lines"), unit(2.6, "lines"), sum(tg$heights)))
wd <- max(convertWidth(sum(tg$widths), "in", valueOnly = TRUE), 11) + 0.4
png(file.path(OUT, paste0(NAME, ".png")), width = wd, height = 1.7 + convertHeight(sum(tg$heights), "in", valueOnly = TRUE),
    units = "in", res = 200, bg = "white")
grid.draw(g); invisible(dev.off())
message("wrote ", file.path(OUT, paste0(NAME, c(".png", ".csv"))), collapse = "\n")
