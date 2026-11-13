## The msa-arm figures drawn from arm_compare_scan_tables.R's tables (Kevin,
## 2026-11-13: are msa / msa_gap better than the site's curated_AAswap, and why):
##
##   recall_rank_<TAG>_<head>   recall of the known ends against the number of non-end
##                        anchors one would have to scan, per set and terminus;
##                        metrics_final.R's ranking (left) and the same ends against
##                        held-out negatives only (right)
##   fp_motif_<scoring>   what the top-scoring non-end anchors of the proteome are:
##                        the share at a dibasic / G + dibasic / monobasic site, a
##                        precursor end, or next to an annotated end
##   tracks_examples      end-score tracks along example precursors, as proteome
##                        ranks, with each known end's held-out score marked
##
##   LC_ALL=en_US.UTF-8 Rscript arm_compare_scan_plots.R
## env: DATA   ~/AF2_analysis/lf_viz_20261113/msa_arms_data   OUT  ~/AF2_analysis/lf_viz_20261113
##      ARMS   the arms to draw, in order (default: every arm in the tables)
##      TRACK_ARMS  curated_AAswap,msa,msa_gap      HEAD  cnn      TAG  msa_arms (in the file names)
##      LABELS arm=label pairs separated by ";"
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
## the shared arm colours (arm_compare_style.R, beside this script)
.here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.here, "arm_compare_style.R"))
A2 <- path.expand("~/AF2_analysis")
DATA <- path.expand(Sys.getenv("DATA", file.path(A2, "lf_viz_20261113/msa_arms_data")))
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113")))
HEAD <- Sys.getenv("HEAD", "cnn")
TAG <- Sys.getenv("TAG", "msa_arms")                   # in every output name
er <- fread(file.path(DATA, "end_ranks.csv"))
ARMS <- strsplit(Sys.getenv("ARMS", paste(unique(er$arm), collapse = ",")), ",")[[1]]
TRACK_ARMS <- strsplit(Sys.getenv("TRACK_ARMS", "curated_AAswap,msa,msa_gap"), ",")[[1]]
lab_env <- Sys.getenv("LABELS", paste("curated=curated (original)", "curated_AAswap=curated_AAswap (site)",
                                      "curated_AAswap_msacons_padg=msacons_padg (previous MSA arm)",
                                      "msa=msa", "msa_gap=msa_gap", sep = ";"))
kv <- strsplit(trimws(strsplit(lab_env, ";")[[1]]), "=")
LAB <- setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
lab <- function(a) unname(ifelse(a %in% names(LAB), LAB[a], a))
acol <- arm_cols(ARMS, lab(ARMS))
comma <- scales::label_comma()
hd_lab <- c(cnn = "CNN head", xgb_pm = "xgb head")[[HEAD]]
th <- theme_bw(base_size = 10.5) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ps <- fread(file.path(DATA, "pool_summary.csv"))

## ---- 2. recall against anchors scanned --------------------------------------------
NG <- unique(round(10^seq(0, 6.5, length.out = 160)))
cv <- er[head == HEAD & arm %in% ARMS]
cv <- rbind(cv[, .(def = "metrics_final.R ranking", rank = rank_std, term, set, arm)],
            cv[, .(def = "held-out negatives only", rank = rank_ho, term, set, arm)])
cv <- cv[!is.na(rank)]
nn <- cv[def == "held-out negatives only", .(n = .N), by = .(term, set, arm)][, .(n = paste(unique(n), collapse = "/")), by = .(term, set)]
rc <- cv[, .(N = NG, recall = vapply(NG, function(x) mean(rank <= x), 0)), by = .(def, term, set, arm)]
rc[, arm := factor(lab(arm), lab(ARMS))]
rc[, set := factor(set, c("cross-validated", "never seen: peptides", "never seen: large products"))]
rc[, panel := factor(sprintf("%s terminus\n%s", term, def),
                     c(outer(c("N", "C"), c("metrics_final.R ranking", "held-out negatives only"),
                             function(t, d) sprintf("%s terminus\n%s", t, d))[c(1, 3, 2, 4)]))]
nn[, set := factor(set, levels(rc$set))]
g2 <- ggplot(rc, aes(N, recall, colour = arm)) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50") +
  geom_step(linewidth = 0.7) +
  geom_text(data = nn[, .(panel = factor(levels(rc$panel)[grepl(term, substr(levels(rc$panel), 1, 1))][1:2], levels(rc$panel)),
                          lab = paste0("n = ", n)), by = .(term, set)], aes(x = 1.3, y = 0.97, label = lab),
            inherit.aes = FALSE, hjust = 0, size = 2.9, colour = "grey30") +
  scale_x_log10(labels = comma, breaks = 10^(0:6)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = acol, name = NULL) +
  facet_grid(set ~ panel) +
  labs(x = "non-end anchors ranked above the end (log) = how far down the proteome list one must read",
       y = "recall of the known ends",
       title = sprintf("Recall against anchors scanned, %s", hd_lab),
       subtitle = paste0("Cross-validated: each known end ranked by the fold that held it out. Never seen: genes no fold trained on, all 50 members; ",
                         "the same ends for every arm.\nHeld-out negatives: every anchor scored only by folds no training window of which ",
                         "reaches it, in any arm (per fold for the cross-validated ends, all folds for never-seen; ",
                         paste(sprintf("%s: %s of %s anchors", ps$term, comma(ps$pool_all), comma(ps$scanned)), collapse = ", "),
                         " reached by none); UniProt peptide/propeptide ends excluded too.")) +
  th
ggsave(file.path(OUT, sprintf("recall_rank_%s_%s.png", TAG, HEAD)), g2, width = 15, height = 10, dpi = 170, bg = "white")

## ---- 4. what the top non-end anchors are --------------------------------------------
fp <- fread(file.path(DATA, "fp_top.csv"))[head == HEAD & arm %in% ARMS]
pm <- fread(file.path(DATA, "pool_motif.csv"))
cat_of <- function(motif, near) factor(fifelse(rep_len(near, length(motif)), "next to an annotated end (+/-3)",
                                               fifelse(motif == "terminus", "precursor end (padded)", motif)),
                                       c("next to an annotated end (+/-3)", "precursor end (padded)", "dibasic",
                                         "G + dibasic", "monobasic", "other"))
mcol <- c("next to an annotated end (+/-3)" = "#CC79A7", "precursor end (padded)" = "#0072B2", dibasic = "#D55E00",
          "G + dibasic" = "#009E73", monobasic = "#E69F00", other = "grey80")
KS <- c(100L, 1000L, 10000L)
## held-out: each fold's top K/5 among the anchors it never trained near; all50: the top K
fpk <- rbindlist(lapply(KS, function(k) fp[rank <= fifelse(scoring == "held_out", k %/% 5L, k),
                                           .(K = k, n = .N), by = .(term, arm, scoring, cat = cat_of(motif, near_end))]))
fpk[, share := n / sum(n), by = .(term, arm, scoring, K)]
fpk[, bar := lab(arm)]
## references: the whole pool, and the known ends themselves
kn <- unique(er[set == "cross-validated", .(term, accession, anchor, motif)])
ref_kn <- kn[, .(n = .N), by = .(term, cat = cat_of(motif, FALSE))][, .(term, cat, share = n / sum(n)), by = term][, -1]
blev <- c("known ends (for reference)", "all non-end anchors (pool)", lab(ARMS))
for (sc in c("held_out", "all50")) {
  ref_pool <- pm[pool == sc, .(share = sum(share)), by = .(term, cat = cat_of(motif, near_end))]
  refs <- rbind(ref_pool[, bar := "all non-end anchors (pool)"], ref_kn[, bar := "known ends (for reference)"])
  d <- rbind(fpk[scoring == sc, .(term, K, bar, cat, share)],
             refs[, .(K = KS), by = .(term, bar, cat, share)][, .(term, K, bar, cat, share)])
  d[, bar := factor(bar, rev(blev))]
  d[, Kf := factor(sprintf("top %s non-end anchors", comma(K)), sprintf("top %s non-end anchors", comma(KS)))]
  d[, term := factor(term, c("N", "C"), c("N terminus", "C terminus"))]
  ## label positions by hand: ggplot stacks the first factor level furthest from zero
  d <- d[order(term, Kf, bar, -as.integer(cat))][, xmid := cumsum(share) - share / 2, by = .(term, Kf, bar)]
  g4 <- ggplot(d, aes(share, bar, fill = cat)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    geom_text(data = d[cat %in% c("dibasic", "G + dibasic", "monobasic", "precursor end (padded)") & share >= 0.06],
              aes(x = xmid, label = sprintf("%.0f%%", 100 * share)), size = 2.6) +
    scale_fill_manual(values = mcol, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = scales::label_percent(), breaks = c(0, 0.25, 0.5, 0.75), expand = c(0, 0)) +
    facet_grid(term ~ Kf) +
    labs(x = "share of the anchors", y = NULL,
         title = sprintf("What the model scores highest outside the annotated ends (%s, %s scoring)", hd_lab,
                         if (sc == "held_out") "held-out" else "all 50 members, as deployed"),
         subtitle = paste0(if (sc == "held_out")
                             "Each anchor scored only by folds no training window of which reaches it, so a trained-on motif negative cannot score low by memory.\n"
                           else "All 50 members: an anchor a member trained on as a negative is a memorised low score (the msa arms trained on most G + dibasic sites).\n",
                           "Annotated ends (ligand list, knowns, never-seen sets, UniProt peptides/propeptides) excluded; ",
                           "the two reference bars give the make-up of every such anchor and of the known ends.")) +
    th + theme(legend.position = "top")
  ggsave(file.path(OUT, sprintf("fp_motif_%s_%s_%s.png", TAG, HEAD, sc)), g4, width = 15, height = 6.5, dpi = 170, bg = "white")
  fwrite(dcast(d, term + Kf + bar ~ cat, value.var = "share", fill = 0),
         file.path(OUT, sprintf("fp_motif_%s_%s_%s.csv", TAG, HEAD, sc)))
}

## ---- 6. tracks along example precursors ------------------------------------------------
tx <- fread(file.path(DATA, "tracks_examples.csv"))
rcols <- sprintf("%s__%s__all50_rank", TRACK_ARMS, HEAD)
if (all(rcols %in% names(tx))) {
  tl <- melt(tx, id.vars = c("term", "accession", "gene", "anchor"), measure.vars = rcols, variable.name = "arm", value.name = "rank")
  tl[, arm := sub(sprintf("__%s__all50_rank$", HEAD), "", arm)]
  ## break the line wherever the scan skips residues (intracellular / TM, excluded by design)
  tl <- tl[order(term, gene, arm, anchor)]
  tl[, seg := cumsum(c(1L, diff(anchor) > 1L)), by = .(term, gene, arm)]
  kc <- er[set == "cross-validated" & head == HEAD & arm %in% TRACK_ARMS & gene %in% unique(tx$gene)]
  tl[, arm := factor(lab(arm), lab(TRACK_ARMS))]; kc[, arm := factor(lab(arm), lab(TRACK_ARMS))]
  tcol <- acol[lab(TRACK_ARMS)]; tcol[is.na(tcol)] <- "black"; names(tcol) <- lab(TRACK_ARMS)
  ymax <- max(tl$rank, na.rm = TRUE)
  one <- function(gn) {
    d <- tl[gene == gn]; k <- kc[gene == gn]
    ggplot(d, aes(anchor, rank, colour = arm)) +
      geom_vline(data = unique(k[, .(term, anchor)]), aes(xintercept = anchor), colour = "grey40", linetype = 2, linewidth = 0.3) +
      geom_line(aes(group = interaction(arm, seg)), linewidth = 0.35, alpha = 0.85) +
      geom_point(data = k, aes(y = rank_std), shape = 21, fill = "white", size = 2, stroke = 0.9) +
      scale_y_continuous(trans = scales::compose_trans("log10", "reverse"), labels = comma,
                         breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6), limits = c(ymax, 1)) +
      scale_colour_manual(values = tcol, name = NULL, drop = FALSE) +
      facet_grid(term ~ ., labeller = labeller(term = c(N = "N end", C = "C end"))) +
      labs(title = gn, x = NULL, y = NULL) +
      th + theme(legend.position = "none", plot.title = element_text(size = 10.5),
                 panel.grid.major.x = element_blank(), strip.text = element_text(size = 8))
  }
  gl <- lapply(sort(unique(tl$gene)), one)
  g6 <- wrap_plots(gl, ncol = 2) + plot_layout(guides = "collect") +
    plot_annotation(title = sprintf("End-score tracks along example precursors (%s)", hd_lab),
                    subtitle = paste0("Line: every scanned anchor's rank among the proteome's non-end anchors, all 50 members ",
                                      "(as on the site; up = better, log).\nDashed: the known ends. Open circles: each known end's ",
                                      "HELD-OUT rank, from the fold that did not train on it. Gaps = residues the scan excludes by design ",
                                      "(intracellular / TM)."),
                    theme = theme(plot.title = element_text(face = "bold", size = 13))) &
    theme(legend.position = "top") & guides(colour = guide_legend(override.aes = list(linewidth = 1.2)))
  ggsave(file.path(OUT, sprintf("tracks_examples_%s_%s.png", TAG, HEAD)), g6, width = 14, height = 2.4 * ceiling(length(gl) / 2) + 1.4,
         dpi = 170, bg = "white")
} else message("tracks_examples.csv lacks ", paste(setdiff(rcols, names(tx)), collapse = ", "), ": no track figure")
message("wrote the figures to ", OUT)
