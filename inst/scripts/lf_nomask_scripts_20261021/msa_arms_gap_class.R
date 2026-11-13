## Does msa_gap put "gap" on the G? Mean per-residue class probabilities of the
## held-out known end windows, aligned at the end, next to their training labels
## (Kevin, 2026-11-12: "add an additional class prediction called 'gap' and
## define it as any residue between a dibasic site and the peptide (usually just
## a single G residue)").
##
## Every known end is taken from the fold that held it out (split == "val" in that
## fold's summary), so each end appears once and nothing is memorised. The class
## probabilities are the ensemble's per-residue softmax for the fold's `all`
## windows (outputs.npz `per_index`, rows in knowns_ranked order -- checked
## against `pred`), so no TensorFlow is loaded. Ends are split into C ends with a
## gap residue (30 single G, one GFL), the other C ends, and N ends (3 of 92 have
## a gap); msa, which has no gap class, is drawn on the same ends.
##
##   LC_ALL=en_US.UTF-8 Rscript msa_arms_gap_class.R
## env: ARMS  msa_gap,msa (the first must carry the gap class: its labels are drawn)
##      OUT   ~/AF2_analysis/lf_viz_20261113     NAME  gap_class_msa_arms
## Writes <OUT>/<NAME>.png and <OUT>/<NAME>.csv (mean probability per group,
## source, position and class) and <OUT>/<NAME>_at_gap.csv (per gap residue).
suppressMessages({library(data.table); library(ggplot2); library(reticulate)})
A2 <- path.expand("~/AF2_analysis")
ARMS <- strsplit(Sys.getenv("ARMS", "msa_gap,msa"), ",")[[1]]
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113")))
NAME <- Sys.getenv("NAME", "gap_class_msa_arms")
## nn_class_cols / nn_class_labels of R/plot_proteins_win_new.R (class names are
## identifiers; the figure names go through the labels)
ccol <- c(CT_cleavage_context = "#FED439", gap = "#8A9197", NT_cleavage_context = "#D2AF81",
          pep_other = "#D5E4A2", pep_pocket = "#197EC0", padding = "grey85", none = "#075149")
clab <- c(CT_cleavage_context = "CT-context", NT_cleavage_context = "NT-context", gap = "gap",
          pep_other = "peptide (non-inserting)", pep_pocket = "peptide (inserting)",
          padding = "padding", none = "none")

py_run_string("import numpy as np\ndef npz_head(f, key, n):\n    return np.load(f, allow_pickle=True)[key][:int(n)]")
P <- list(); L <- list(); have <- list()
for (tm in c("C", "N")) for (k in 1:5) for (a in ARMS) {
  base <- sprintf("%s/lf_pepend_run_%s_cv5_f%d_nomask_flatten_%s", A2, tm, k, a)
  fs <- paste0(base, "_summary.rds"); fo <- paste0(base, "_isolated/out/outputs.npz")
  if (!file.exists(fs) || !file.exists(fo)) { message("missing ", basename(base)); next }
  s <- readRDS(fs); kr <- as.data.table(s$knowns_ranked)
  cls <- unlist(s$config$class_names)
  if (is.null(cls)) cls <- unlist(jsonlite::read_json(paste0(base, "_isolated/in/config.json"))$class_names)
  pr <- py$npz_head(fo, "pred", nrow(kr))
  stopifnot(max(abs(pr - kr$pred)) < 1e-6)              # outputs rows = knowns_ranked rows
  pix <- py$npz_head(fo, "per_index", nrow(kr))          # ends x 36 x classes
  stopifnot(dim(pix)[3] == length(cls))
  for (i in which(kr$split == "val")) {
    lab <- unlist(kr$known_idx[i]); w <- length(lab)
    rel <- seq_len(w) - (kr$anchor[i] - kr$w_start[i] + 1L)   # 0 = the end residue
    id <- data.table(term = tm, fold = k, arm = a, gene = kr$gene[i], anchor = kr$anchor[i])
    P[[length(P) + 1]] <- id[, .(rel = rep(rel, length(cls)), class = rep(cls, each = w), prob = as.vector(pix[i, , ])),
                             by = .(term, fold, arm, gene, anchor)]
    if (a == ARMS[1]) L[[length(L) + 1]] <- id[, .(rel = rel, label = lab), by = .(term, fold, gene, anchor)]
  }
  have[[length(have) + 1]] <- data.table(term = tm, fold = k, arm = a)
}
P <- rbindlist(P); L <- rbindlist(L); have <- rbindlist(have)
stopifnot("gap" %in% L$label)
## an end has a gap if its gap-arm labels say so; only ends both arms scored
ends <- L[, .(has_gap = any(label == "gap"), gap_res = sum(label == "gap")), by = .(term, gene, anchor)]
both <- P[, uniqueN(arm), by = .(term, gene, anchor)][V1 == length(ARMS), .(term, gene, anchor)]
ends <- ends[both, on = .(term, gene, anchor), nomatch = 0L]
grp <- function(term, has_gap) fifelse(term == "N", "N ends",
                                       fifelse(has_gap, "C ends with a gap residue", "C ends without one"))
ends[, group := grp(term, has_gap)]
ng <- ends[, .N, by = group]
glev <- c("C ends with a gap residue", "C ends without one", "N ends")
ends[, group_n := factor(sprintf("%s (n = %d)", group, ng$N[match(group, ng$group)]),
                         sprintf("%s (n = %d)", glev, ng$N[match(glev, ng$group)]))]
P <- P[ends, on = .(term, gene, anchor), nomatch = 0L]
L <- L[ends, on = .(term, gene, anchor), nomatch = 0L]

## mean composition per position: the labels' class shares, and each arm's mean softmax
lab_m <- L[, .(prob = as.double(.N)), by = .(group_n, rel, class = label)][, prob := prob / sum(prob), by = .(group_n, rel)]
lab_m[, source := "training labels"]
arm_m <- P[, .(prob = mean(prob)), by = .(group_n, rel, class, source = paste0(arm, " (held out)"))]
M <- rbind(lab_m, arm_m)
M[, source := factor(source, c("training labels", paste0(ARMS, " (held out)")))]
M[, class := factor(class, names(ccol))]
fwrite(M[order(group_n, source, rel, class)], file.path(OUT, paste0(NAME, ".csv")))

## per gap residue: what each arm puts there, and is "gap" the argmax? Only the
## end's OWN gap (C: after the end, N: before it) -- an N window of a short
## peptide also holds that peptide's C-side G further in
gp <- L[label == "gap" & ((term == "C" & rel > 0) | (term == "N" & rel < 0)), .(term, gene, anchor, rel)]
at_gap <- P[gp, on = .(term, gene, anchor, rel)][, .(class_top = class[which.max(prob)], p_gap = sum(prob[class == "gap"]),
                                                     p_ct = sum(prob[class == "CT_cleavage_context"])),
                                                  by = .(arm, term, gene, anchor, rel)]
fwrite(at_gap, file.path(OUT, paste0(NAME, "_at_gap.csv")))
cap <- at_gap[, .(txt = sprintf("%s, %s ends' own gap residues (%d): argmax = %s; mean P(gap) %.2f", arm, term, .N,
                                paste(sprintf("%s %d", clab[names(table(class_top))], table(class_top)), collapse = ", "),
                                mean(p_gap))), by = .(term, arm)]$txt
miss <- setdiff(CJ(term = c("C", "N"), fold = 1:5, arm = ARMS)[, paste(term, fold, arm)], have[, paste(term, fold, arm)])

g <- ggplot(M, aes(rel, prob, fill = class)) +
  geom_col(width = 1, colour = NA) +
  geom_vline(data = data.table(group_n = factor(levels(M$group_n), levels(M$group_n)),
                               x = ifelse(grepl("^N", levels(M$group_n)), -0.5, 0.5)),
             aes(xintercept = x), linewidth = 0.4, linetype = 2) +
  scale_fill_manual(values = ccol, labels = clab, name = NULL, drop = TRUE) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(source ~ group_n, scales = "free_x") +
  labs(x = "residue relative to the peptide end (0 = the end residue; dashed line = the cut)",
       y = "mean class share / probability",
       title = "The gap class: what the trunks predict around the held-out known ends",
       subtitle = paste0("Top row: the 7-class training labels (share of ends). Below: each arm's mean per-residue softmax, ",
                         "every end scored by the fold that held it out.\n", paste(cap, collapse = "\n"),
                         if (length(miss)) paste0("\nDRAFT: missing ", paste(miss, collapse = ", ")) else "")) +
  theme_bw(base_size = 10.5) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        panel.grid = element_blank(), panel.spacing.x = unit(0.8, "lines"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, paste0(NAME, ".png")), g, width = 13, height = 2.3 * nlevels(M$source) + 2, dpi = 170, bg = "white")
print(at_gap[, .(n = .N, gap_top = sum(class_top == "gap"), mean_p_gap = round(mean(p_gap), 3),
                 mean_p_ct = round(mean(p_ct), 3)), by = .(arm, term)])
message("wrote ", file.path(OUT, paste0(NAME, c(".png", ".csv", "_at_gap.csv"))), collapse = "\n")
