## Each held-out known end's rank under one arm against another, coloured by the
## motif flanking the end (Kevin, 2026-11-13: "whether they are better, and why"
## for the msa / msa_gap arms of 2026-11-12).
##
## One point per known end, ranked by the fold that held it out. Below the
## diagonal = the second arm ranks the end higher (better). The ends that moved
## most either way are named. A second figure gives the median rank per motif and
## arm, which is the "why": motif-matched negatives should cost the dibasic /
## G + dibasic ends something and buy the padded / "other" ends something.
##
##   RANKS=~/AF2_analysis/lf_viz_20261113/msa_arms_data/end_ranks.csv \
##   PAIRS=curated_AAswap:msa,msa:msa_gap LC_ALL=en_US.UTF-8 Rscript arm_compare_paired_ranks.R
##
## env: RANKS   csv with term, arm, gene, anchor, rank [, head, set]: either
##              trunk_val_compare.R's trunk_val_ranks.csv (rank among the ~40k
##              candidate windows) or arm_compare_scan_tables.R's end_ranks.csv (rank
##              among every non-end anchor of the proteome)
##      RANK_COL  rank_std (default when present), rank_ho, or rank
##      HEAD    the head to plot when RANKS has a head column (cnn)
##      SET     the set to plot when RANKS has a set column (cross-validated)
##      PAIRS   x:y arm pairs, comma-separated; y is the arm under test
##      ARMS    the arms of the per-motif figure, in order (default: those in PAIRS)
##      LABELS  arm=label pairs separated by ";" (optional)
##      RANK_DESC  axis wording for the rank
##      OUT     ~/AF2_analysis/lf_viz_20261113     NAME  paired_ranks_msa_arms
## Writes <OUT>/<NAME>.png, <OUT>/<NAME>_by_motif.png and <OUT>/<NAME>.csv
## (one row per end and pair, with both ranks).
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
## the shared arm colours (arm_compare_style.R, beside this script)
.here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.here, "arm_compare_style.R"))
A2 <- path.expand("~/AF2_analysis")
RANKS <- path.expand(Sys.getenv("RANKS", file.path(A2, "lf_viz_20261113/msa_arms_data/end_ranks.csv")))
PAIRS <- strsplit(strsplit(Sys.getenv("PAIRS", "curated_AAswap:msa,msa:msa_gap"), ",")[[1]], ":")
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113")))
NAME <- Sys.getenv("NAME", "paired_ranks_msa_arms")
RANK_DESC <- Sys.getenv("RANK_DESC", "rank among every non-end anchor of the proteome")
lab_env <- Sys.getenv("LABELS", "")
LAB <- if (nzchar(lab_env)) {
  kv <- strsplit(trimws(strsplit(lab_env, ";")[[1]]), "=")
  setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
} else character(0)
lab <- function(a) unname(ifelse(a %in% names(LAB), LAB[a], a))
ARMS <- strsplit(Sys.getenv("ARMS", paste(unique(unlist(PAIRS)), collapse = ",")), ",")[[1]]

r <- fread(RANKS)
if ("head" %in% names(r)) r <- r[head == Sys.getenv("HEAD", "cnn")]
if ("set" %in% names(r)) r <- r[set == Sys.getenv("SET", "cross-validated")]
## end_ranks.csv carries two rankings (rank_std = metrics_final.R's, rank_ho =
## against held-out negatives); trunk_val_ranks.csv one (rank)
RANK_COL <- Sys.getenv("RANK_COL", if ("rank_std" %in% names(r)) "rank_std" else "rank")
r[, `:=`(arm = as.character(arm), rank = as.double(get(RANK_COL)))]
r <- r[!is.na(rank)]

## the motif flanking each end, from the shared knowns tables (terminus = the
## end is the precursor's own end, so its window is padded)
mot <- rbindlist(lapply(c("N", "C"), function(tm)
  fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)), select = c("accession", "anchor", "gene.x", "motif"))[
    , term := tm]))
setnames(mot, "gene.x", "gene")
mot[, motif := factor(fifelse(motif == "terminus", "padded (precursor end)", motif),
                      c("padded (precursor end)", "monobasic", "dibasic", "G + dibasic", "other"))]
key <- if ("accession" %in% names(r)) c("term", "accession", "anchor") else c("term", "gene", "anchor")
## the _padg arms anchored the 30 amidated C ends on their G (mature end + 1;
## removed everywhere on 2026-11-12): put those back on the mature anchor
mk <- mot[, paste(term, if ("accession" %in% key) accession else gene, anchor)]
rk <- r[, paste(term, if ("accession" %in% key) accession else gene, anchor)]
rk1 <- r[, paste(term, if ("accession" %in% key) accession else gene, anchor - 1L)]
g_anch <- !(rk %in% mk) & r$term == "C" & rk1 %in% mot[motif == "G + dibasic", mk[motif == "G + dibasic"]]
if (any(g_anch)) {
  message(sum(g_anch), " C ends of ", paste(unique(r$arm[g_anch]), collapse = ", "),
          " sit on the G: moved to the mature anchor")
  r[g_anch, anchor := anchor - 1L]
}
r <- mot[, c(key, setdiff(c("gene", "accession", "motif"), key)), with = FALSE][r, on = key]
if (r[, anyNA(motif)]) warning(r[is.na(motif), .N], " ends have no motif in the knowns tables")
mcol <- c("padded (precursor end)" = "#0072B2", monobasic = "#E69F00", dibasic = "#D55E00",
          "G + dibasic" = "#009E73", other = "#7F7F7F")

## ---- paired scatter --------------------------------------------------------
pp <- rbindlist(lapply(PAIRS, function(p) {
  x <- r[arm == p[1], unique(c(key, "gene", "accession", "motif", "rank")), with = FALSE]
  y <- r[arm == p[2], c(key, "rank"), with = FALSE]
  w <- x[y, on = key, nomatch = 0L]
  setnames(w, c("rank", "i.rank"), c("rank_x", "rank_y"))
  w[, `:=`(arm_x = p[1], arm_y = p[2], pair = sprintf("%s  ->  %s", lab(p[1]), lab(p[2])))]
}))
if (!nrow(pp)) stop("no end is ranked by both arms of any pair", call. = FALSE)
pp[, lr := log2((rank_y + 1) / (rank_x + 1))]
pp[, pair := factor(pair, unique(pair))]
pp[, term := factor(term, c("N", "C"), c("N terminus", "C terminus"))]
## the six biggest movers each way per panel
pp[, mover := frank(lr, ties.method = "first") <= 6 | frank(-lr, ties.method = "first") <= 6, by = .(pair, term)]
pp[, label := fifelse(mover & abs(lr) >= 1, sprintf("%s %d", gene, anchor), NA_character_)]
st <- pp[, .(n = .N, better = sum(rank_y < rank_x), worse = sum(rank_y > rank_x),
             med_x = median(rank_x), med_y = median(rank_y)), by = .(pair, term)]
st[, txt := sprintf("%d ends: %d better, %d worse\nmedian rank %s -> %s", n, better, worse,
                    format(med_x, big.mark = ","), format(med_y, big.mark = ","))]
lims <- range(c(pp$rank_x, pp$rank_y))
g <- ggplot(pp, aes(rank_x, rank_y)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2) +
  geom_point(aes(colour = motif), size = 2, alpha = 0.85) +
  geom_text_repel(aes(label = label), size = 2.7, min.segment.length = 0, max.overlaps = Inf,
                  seed = 1, na.rm = TRUE, colour = "grey15") +
  geom_text(data = st, aes(x = 1, y = lims[2], label = txt), hjust = 0, vjust = 1, size = 3,
            inherit.aes = FALSE) +
  scale_x_log10(labels = scales::label_comma(), limits = c(1, lims[2])) +
  scale_y_log10(labels = scales::label_comma(), limits = c(1, lims[2])) +
  scale_colour_manual(values = mcol, name = "flanking motif", drop = FALSE) +
  facet_grid(pair ~ term, switch = "y") + coord_equal() +
  labs(x = "rank under the first arm (log; lower = better)", y = "rank under the second arm (log)",
       title = "Held-out known ends: rank under one arm against another",
       subtitle = paste0("Each end ranked by the fold that held it out; ", RANK_DESC,
                         ".\nBelow the dashed diagonal = the second arm ranks the end higher. ",
                         "Named: the six biggest movers each way (at least 2-fold).")) +
  theme_bw(base_size = 10.5) +
  theme(legend.position = "top", strip.placement = "outside", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, paste0(NAME, ".png")), g, width = 11, height = 5.2 * length(PAIRS) + 1, dpi = 170, bg = "white")
fwrite(pp[order(pair, term, lr), .(pair, term, accession, gene, anchor, motif, rank_x, rank_y, log2_ratio = round(lr, 3))],
       file.path(OUT, paste0(NAME, ".csv")))

## ---- per motif: median rank per arm ----------------------------------------
bm <- r[arm %in% ARMS & !is.na(motif), .(n = .N, med = median(rank), q25 = quantile(rank, 0.25),
                                          q75 = quantile(rank, 0.75)), by = .(term, motif, arm)]
bm[, arm := factor(lab(arm), lab(ARMS))]
bm[, term := factor(term, c("N", "C"), c("N terminus", "C terminus"))]
bm[, motif_n := sprintf("%s (n=%d)", motif, n)]
g2 <- ggplot(bm, aes(med, motif, colour = arm)) +
  geom_pointrange(aes(xmin = q25, xmax = q75), position = position_dodge(width = 0.6), size = 0.35) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = arm_cols(ARMS, lab(ARMS)), name = NULL) +
  facet_wrap(~term) +
  labs(x = "held-out rank of the known ends (median, bar = quartiles; log; lower = better)", y = NULL,
       title = "Where the gain comes from: held-out rank by the motif flanking the end",
       subtitle = paste0(RANK_DESC, "; n per motif in the knowns tables (the same ends for every arm)")) +
  theme_bw(base_size = 10.5) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(OUT, paste0(NAME, "_by_motif.png")), g2, width = 10, height = 4.6, dpi = 170, bg = "white")
message("wrote ", file.path(OUT, paste0(NAME, c(".png", "_by_motif.png", ".csv"))), collapse = "\n")
