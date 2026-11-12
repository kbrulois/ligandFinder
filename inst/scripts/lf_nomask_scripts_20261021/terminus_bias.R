## Does the window model read "padding = peptide"? The terminus-bias diagnostic
## behind the _padg training windows (Kevin, 2026-11-10).
##
## Kevin, 2026-11-10: "I believe the negative windows don't include certain
## examples such as assumed negative windows that are anchored at the N or C
## terminus of the mature precursor or other windows like mono basic anchored.
## If so, I believe this is causing the model to learn that if padding is
## present, it is a peptide window."
##
## Three measurements. 1 and 2 are what confirmed it on 2026-11-10 (and what the
## 10_7 --pad-match / --mono-match comment quotes); 3 is what showed the fix
## holding on the full proteome on 2026-11-11 (ARM_PIPELINE_README.md, "Arms on
## their own windows").
##
##   1 WINDOWS  the fold-FOLD training windows (train + val) of every window set
##     the ARMS were trained on: the share carrying padding by known / win_type,
##     and the anchor's motif and distance from the precursor terminus,
##     positives vs each kind of negative. On the old windows: pep_end negatives
##     ~7% of the negatives, positives ~20% AT a terminus and ~45% padded against
##     ~20% of the negatives. Anchors 1-7 residues from a terminus were ~5% of
##     the negatives (38 N, 37 C of 800) against 4% / 9% of the N / C positives;
##     the 2026-11-10 summary's "almost never negatives" overstated that.
##   2 BEFORE   BEFORE_ARM's full-proteome tracks, each anchor scored by the mean
##     of the five per-fold trunks: unlabelled anchors by distance from the
##     precursor terminus and by motif, the share of the top-1000 non-known
##     anchors at a terminus, the known ends at / away from one. curated_AAswap:
##     CNN > 0.5 at the terminus N 17.5%, C 16.8%, interior ~0.1-0.2%; a third
##     of the top-1000 non-known N anchors are termini; monobasic ~3x "other".
##   3 AFTER    the share scoring CNN > 0.5 by distance bin, every arm in ARMS:
##     at the terminus N 17.5% -> 5.7%, C 16.8% -> 3.3%; 1-3 residues in
##     N 9.7% -> 1.3%, C 5.2% -> 1.2%; a terminus/interior ratio of 38-99x stays,
##     because real peptide ends do concentrate there.
##
## "Unlabelled" = not a known end of the arm's own fold input (so at the arm's
## own anchor convention, arm_conv.R), on a gene outside that input's `banned`
## list (the known ligands; every known end's gene is on it). Residues excluded
## BY DESIGN -- cons_rs_n, min_afm, mean_afm and relASA all exactly 0 at the
## anchor, i.e. intracellular / TM -- are left out of 2 and 3; those zeros are
## topology, not missing data. Distance = anchor - n_prot (N), c_prot - anchor
## (C), from the residue cache; 0 is the precursor's own mature terminus. One
## training window lies before it: SAA1's known N end at 18, the cache's mature
## start being 19 (10_7 moves a known's n_prot back to its peptide's start). It
## is padded like a terminus window, so 1 counts it AT the terminus, as the
## 2026-11-10 tables did; 2 and 3 never see a negative distance.
## The per-fold mean is in-sample for the windows each fold trained on, as the
## 2026-11-10/11 numbers were; nothing here is a held-out metric.
##
##   LC_ALL=en_US.UTF-8 Rscript terminus_bias.R
##   ARMS=curated_AAswap,curated_AAswap_padg BEFORE_ARM=curated_AAswap_padg \
##     HEADCMP=/tmp/tb LC_ALL=en_US.UTF-8 Rscript terminus_bias.R
##
## Writes to HEADCMP (default ~/AF2_analysis/lf_headcmp):
##   terminus_bias_windows.csv          1: window counts by set x term x split x
##                                         group x distance bin x motif x padded
##   terminus_bias_before_<arm>.csv     2: by distance bin, CNN and xgb
##   terminus_bias.csv                  3: term x arm x distance bin
##   terminus_bias.png                  1 above 3, on the same distance bins
## ~2.5 GB of RAM and under a minute: the residue cache, two fold inputs, then
## six 150 MB track files one at a time.
suppressMessages({library(dplyr); library(data.table); library(ggplot2); library(patchwork)})
A2 <- path.expand("~/AF2_analysis")
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
## arm_var() / arm_amid_g(): each arm's training windows and C-end convention
source(file.path(Sys.getenv("LF_SCRIPTS", .lf_here), "arm_conv.R"))
source(lf_src("R/pepend_windows.R"))   # lf_pepend_motif, lf_pepend_motif_vec
source(lf_src("R/viz_theme.R"))        # LF_VIZ, lf_viz_theme

ARMS <- strsplit(Sys.getenv("ARMS", "curated_AAswap,curated_AAswap_padg,curated_AAswap_msacons_padg"), ",")[[1]]
if (length(ARMS) > 3) stop("the figure has three colour slots; pass at most three ARMS", call. = FALSE)
BEFORE_ARM <- Sys.getenv("BEFORE_ARM", "curated_AAswap")
FOLD <- as.integer(Sys.getenv("FOLD", "1"))   # any fold: all five hold the same knowns
HC <- path.expand(Sys.getenv("HEADCMP", file.path(A2, "lf_headcmp")))
dir.create(HC, showWarnings = FALSE, recursive = TRUE)
## the 2026-11-11 table's names for the three arms; any other arm keeps its own
LAB <- c(curated_AAswap = "old windows", curated_AAswap_padg = "padg",
         curated_AAswap_msacons_padg = "padg + MSA cons")
lab <- function(a) unname(ifelse(a %in% names(LAB), LAB[a], a))
BINS <- c("0 (terminus)", "1-3", "4-7", "8-15", ">15")
bin5 <- function(d) cut(d, c(-Inf, 0, 3, 7, 15, Inf), labels = BINS)

cache <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
prec  <- as.data.table(cache$prec)[, .(accession, gene, seq, n_prot, c_prot)]
excl  <- rbindlist(lapply(names(cache$feats), function(a) { m <- cache$feats[[a]]
  data.table(accession = a, anchor = seq_len(nrow(m)),
             excl = m[, "cons_rs_n"] == 0 & m[, "min_afm"] == 0 & m[, "mean_afm"] == 0 & m[, "relASA"] == 0) }))
rm(cache); invisible(gc())

## An arm's fold input (its INPUT_SFX, arms/<arm>.sh): the train + val windows,
## the known ends at the arm's own anchor convention and the banned genes. The
## inputs are ~90 MB each, so each window set is read once and only this kept.
.inputs <- list()
fold_input <- function(arm) {
  sfx <- arm_var(arm, "INPUT_SFX"); key <- paste0("f", sfx)
  if (is.null(.inputs[[key]])) {
    f <- file.path(A2, sprintf("lf_pepend_nn_input_cv5_f%d%s.rds", FOLD, sfx))
    x <- readRDS(f)
    ## an arm file pointing at another convention's windows would grade every
    ## moved C end one residue off
    if (!identical(isTRUE(x$params$amid_g), arm_amid_g(arm)))
      stop(basename(f), " has amid_g=", isTRUE(x$params$amid_g), " but arms/", arm,
           ".sh says ", arm_amid_g(arm), call. = FALSE)
    win <- rbindlist(lapply(c("N", "C"), function(t) rbindlist(lapply(c("train", "val"), function(s) {
      d <- x$nn_input[[t]][[s]]
      data.table(term = t, split = s, accession = d$accession, anchor = as.integer(d$anchor),
                 known = d$known, win_type = d$win_type,
                 pad = vapply(d$data, function(m) sum(m[, "padding"]), numeric(1)))
    }))))
    .inputs[[key]] <<- list(windows = sub("\\.rds$", "", sub("^lf_pepend_nn_input_", "", basename(f))),
                            name = if (nzchar(sfx)) paste(sub("^_", "", sfx), "windows") else "old windows",
                            amid_g = isTRUE(x$params$amid_g), win = win, banned = x$params$banned,
                            known = paste(x$knowns$accession, x$knowns$term, x$knowns$anchor))
    rm(x); invisible(gc())
  }
  .inputs[[key]]
}

## ---- 1. the training windows ---------------------------------------------------------
win_arms <- ARMS[!duplicated(vapply(ARMS, arm_var, "", var = "INPUT_SFX"))]
W1 <- rbindlist(lapply(win_arms, function(a) {
  fi <- fold_input(a)
  d <- merge(fi$win, prec[, .(accession, seq, n_prot, c_prot)], by = "accession", sort = FALSE)
  stopifnot(nrow(d) == nrow(fi$win))
  d[, `:=`(windows = fi$windows, name = fi$name,
           motif   = lf_pepend_motif(seq, term, anchor, n_prot, amid_g = fi$amid_g),
           dist    = ifelse(term == "N", anchor - n_prot, c_prot - anchor),
           group   = ifelse(known == 1, "POS", paste0("neg_", win_type)))][, seq := NULL]
}))
W1[, term := factor(term, c("N", "C"))]
for (w in unique(W1$windows)) {
  cat("\n######## training windows", w, "-- padding by known / win_type\n")
  for (t in c("N", "C")) for (s in c("train", "val")) {
    cat("==", t, s, "\n")
    print(as.data.frame(W1[windows == w & term == t & split == s,
      .(n = .N, frac_pad = round(mean(pad > 0), 3), mean_pad = round(mean(pad), 1)), keyby = .(known, win_type)]))
  }
  for (t in c("N", "C")) {
    d <- W1[windows == w & term == t]
    cat("\n=====", w, t, "motif by group (train + val)\n"); print(table(d$group, d$motif))
    cat("distance of the anchor from the precursor terminus (0 = at the terminus):\n")
    print(table(d$group, cut(d$dist, c(-Inf, 0, 3, 7, 15, Inf), labels = c("0", "1-3", "4-7", "8-15", ">15"))))
  }
}
## --pad-match / --mono-match match each SPLIT's negatives to its positives, so
## the splits are shown as well as their union (10_7's own log prints the same
## shares per split, from the window's n_prot rather than the cache's)
cat("\n######## positives vs negatives, per split and train + val\n")
## at the terminus = the "0" bin of the tables and the figure, SAA1's -1 included
hl <- function(d) d[, .(n_pos = sum(known == 1), n_neg = sum(known == 0),
                        pos_at_terminus = round(mean(dist[known == 1] <= 0), 3),
                        neg_at_terminus = round(mean(dist[known == 0] <= 0), 3),
                        neg_pep_end     = round(mean(win_type[known == 0] == "pep_end"), 3),
                        neg_1_7         = sum(known == 0 & dist >= 1 & dist <= 7),
                        pos_padded      = round(mean(pad[known == 1] > 0), 3),
                        neg_padded      = round(mean(pad[known == 0] > 0), 3),
                        pos_monobasic   = round(mean(motif[known == 1] == "monobasic"), 3),
                        neg_monobasic   = round(mean(motif[known == 0] == "monobasic"), 3)),
                    keyby = .(windows, term, split)]
print(rbind(hl(W1), hl(copy(W1)[, split := "train+val"]))[order(windows, term, match(split, c("train", "val")))],
      row.names = FALSE)
fwrite(W1[, .(n = .N), keyby = .(windows, term, split, group, dist_bin = bin5(dist), motif, padded = pad > 0)],
       file.path(HC, "terminus_bias_windows.csv"))

## ---- 2 and 3. the full-proteome tracks ---------------------------------------------------
## one arm's per-anchor scores, with everything the two tables filter on
scored <- function(t, arm) {
  f <- file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", t,
                             if (arm %in% c("", "curated", "flatten")) "" else paste0("_", arm)))
  if (!file.exists(f)) { message("no tracks for ", arm, " ", t, ": ", f); return(NULL) }
  tr <- as.data.table(readRDS(f))
  cc <- grep("^cnn__f[0-9]", names(tr), value = TRUE); xc <- grep("^xgb_pm__f[0-9]", names(tr), value = TRUE)
  tr[, cnn := rowMeans(.SD), .SDcols = cc][, xgb := rowMeans(.SD), .SDcols = xc]
  tr <- merge(tr[, .(accession, anchor, cnn, xgb)], prec, by = "accession")
  tr <- merge(tr, excl, by = c("accession", "anchor"), all.x = TRUE)
  tr[, dist := if (t == "N") anchor - n_prot else c_prot - anchor]
  tr[, known := paste(accession, t, anchor) %in% fold_input(arm)$known]
  setattr(tr, "folds", length(cc))
  tr
}

before <- function(tr, t, arm) {
  fi <- fold_input(arm)
  tr[, motif := lf_pepend_motif_vec(seq[1], t, anchor, n_prot[1], amid_g = fi$amid_g), by = accession]
  neg <- tr[!known & !gene %in% fi$banned & dist >= 0]
  neg[, dbin := cut(dist, c(-Inf, 0, 3, 7, 15, 40, Inf), labels = c("0 (terminus)", "1-3", "4-7", "8-15", "16-40", ">40"))]
  cat("\n=====", arm, t, "unlabelled anchors of non-ligand precursors, NOT excluded;",
      "score (mean of the", attr(tr, "folds"), "per-fold trunks) by distance from the precursor terminus\n")
  b <- neg[excl == FALSE, .(n = .N, cnn_med = signif(median(cnn), 3), cnn_q99 = signif(quantile(cnn, .99), 3),
                            frac_cnn_gt_0.5 = round(mean(cnn > .5), 4), xgb_med = signif(median(xgb), 3),
                            frac_xgb_gt_0.5 = round(mean(xgb > .5), 4)), keyby = dbin]
  print(b)
  cat("by motif (not excluded, dist > 15):\n")
  print(neg[excl == FALSE & dist > 15, .(n = .N, cnn_med = signif(median(cnn), 3), frac_cnn_gt_0.5 = round(mean(cnn > .5), 4),
                                         frac_xgb_gt_0.5 = round(mean(xgb > .5), 4)), keyby = motif])
  top <- tr[!excl & !known][order(-cnn)]
  cat("share of the top-1000 non-known CNN anchors at a terminus (dist 0):", round(mean(top$dist[1:1000] == 0), 3),
      " within 7:", round(mean(top$dist[1:1000] <= 7), 3), " base rate dist 0:", round(mean(top$dist == 0), 4), "\n")
  cat("known ends, at a terminus or not:\n")
  print(tr[known == TRUE, .(n = .N, cnn_med = signif(median(cnn), 3)), keyby = .(at_term = dist == 0)])
  b[, term := t]
}

B <- list(); R <- list()
for (a in union(ARMS, BEFORE_ARM)) for (t in c("N", "C")) {
  tr <- scored(t, a)
  if (is.null(tr)) next
  if (a == BEFORE_ARM) B[[t]] <- before(tr, t, a)
  if (a %in% ARMS)
    R[[paste(a, t)]] <- tr[!known & !gene %in% fold_input(a)$banned & dist >= 0 & excl == FALSE,
                           .(n = .N, n_gt_0.5 = sum(cnn > .5), frac_gt_0.5 = mean(cnn > .5)),
                           keyby = .(bin = bin5(dist))][, `:=`(term = t, arm = a)]
  tr <- NULL; invisible(gc())
}
if (length(B)) fwrite(rbindlist(B), file.path(HC, sprintf("terminus_bias_before_%s.csv", BEFORE_ARM)))
R <- rbindlist(R)
fwrite(R[, .(term, arm, label = lab(arm), bin, n, n_gt_0.5, frac_gt_0.5)], file.path(HC, "terminus_bias.csv"))

w <- dcast(R, term + arm ~ bin, value.var = "frac_gt_0.5")
tab <- w[, .(term, arm = lab(arm), term0 = `0 (terminus)`, d1_3 = `1-3`, d4_7 = `4-7`, d8_15 = `8-15`,
             interior = `>15`, ratio_term_vs_interior = round(`0 (terminus)` / `>15`, 1))]
tab <- tab[order(match(arm, lab(ARMS)), match(term, c("N", "C")))]
cat("\nShare of UNLABELLED extracellular anchors (non-ligand genes) scoring CNN > 0.5,",
    "by distance from the precursor terminus:\n")
print(tab[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], row.names = FALSE)

## ---- the figure: where the training windows sit, then where the model fires ---------------
## Colour follows the arm; a window set's negatives (upper row) wear the colour
## of the first arm in ARMS trained on that set. The positives are the
## reference, in ink.
sets <- unique(W1$name)
arm_set <- vapply(ARMS, function(a) fold_input(a)$name, "")
pal <- setNames(c(LF_VIZ$s1, LF_VIZ$s2, LF_VIZ$s3)[seq_along(ARMS)], lab(ARMS))
set_col <- setNames(pal[lab(ARMS[match(sets, arm_set)])], paste("negatives,", sets))
pos_set <- tail(sets, 1)       # the knowns of the newest set: they differ only by the C G shift
c1 <- W1[known == 0, .(n = .N), keyby = .(name, term, bin = bin5(dist))]
c1 <- rbind(c1[, .(series = paste("negatives,", name), term, bin, n)],
            W1[known == 1 & name == pos_set, .(n = .N), keyby = .(term, bin = bin5(dist))][, series := "known ends"])
c1 <- c1[CJ(series = unique(series), term = factor(c("N", "C"), c("N", "C")), bin = factor(BINS, BINS)),
         on = .(series, term, bin)][is.na(n), n := 0L][, share := n / sum(n), by = .(series, term)]
c1[, series := factor(series, c("known ends", names(set_col)))]
col1 <- c("known ends" = LF_VIZ$ink, set_col)
p1 <- ggplot(c1, aes(bin, share, colour = series, group = series)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_wrap(~ term, labeller = as_labeller(c(N = "N windows", C = "C windows"))) +
  scale_colour_manual(values = col1, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_x_discrete(expand = expansion(add = c(1, 0.4))) +     # bins above the lower row's
  labs(title = sprintf("Training windows (fold %d, train + val): where the anchor sits", FOLD),
       subtitle = "share of each group's windows by distance of the anchor from the precursor's mature terminus",
       x = NULL, y = "share of windows") +
  lf_viz_theme(10) + theme(legend.position = "top", legend.justification = "left")
R[, `:=`(term = factor(term, c("N", "C")), label = factor(lab(arm), lab(ARMS)))]
lbl <- R[bin == BINS[1]][, txt := scales::percent(frac_gt_0.5, accuracy = 0.1)]
p2 <- ggplot(R, aes(bin, frac_gt_0.5, colour = label, group = label)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  ggrepel::geom_text_repel(data = lbl, aes(label = txt), size = 3, colour = LF_VIZ$ink2,
                           nudge_x = -0.45, direction = "y", hjust = 1, segment.colour = LF_VIZ$grid,
                           min.segment.length = 0, seed = 1) +
  facet_wrap(~ term, labeller = as_labeller(c(N = "N anchors", C = "C anchors"))) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_y_log10(labels = function(x) paste0(signif(100 * x, 2), "%")) +
  scale_x_discrete(expand = expansion(add = c(1, 0.4))) +
  labs(title = "Full-proteome scan: unlabelled anchors the model calls (CNN > 0.5)",
       subtitle = paste0("non-ligand precursors, residues excluded by design left out; ",
                         "each anchor scored by the mean of the five per-fold trunks"),
       x = "distance of the anchor from the precursor's mature terminus (residues)",
       y = "share scoring CNN > 0.5 (log)",
       caption = if ("padg windows" %in% sets)
         paste0("padg windows: 10_7 --pad-match --mono-match --amid-g-anchor -- negatives topped up to ",
                "the positives' padding and monobasic shares per split,\nC ends followed by G + dibasic ",
                "anchored on the G. Known ends drawn from the ", pos_set, ".")) +
  lf_viz_theme(10) + theme(legend.position = "top", legend.justification = "left")
png_f <- file.path(HC, "terminus_bias.png")
ggsave(png_f, p1 / p2 + plot_layout(heights = c(1, 1.15)), width = 9, height = 8.5, dpi = 200, bg = LF_VIZ$surface)
cat("\nwrote", file.path(HC, c("terminus_bias_windows.csv", sprintf("terminus_bias_before_%s.csv", BEFORE_ARM),
                               "terminus_bias.csv", "terminus_bias.png")), sep = "\n  ")
