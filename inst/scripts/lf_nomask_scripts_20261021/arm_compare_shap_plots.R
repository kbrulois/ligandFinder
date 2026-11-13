## Per-position and global attribution of the window trunks, arm by arm (Kevin,
## 2026-11-13: "per window position shap values for AAid, relASA etc. ... a
## per-window position heatmap for AAid only. and global shap for each model").
##
## Reads shap_positions.py's tables -- expected gradients on each member's logit,
## the Shapley value of every (position, channel) cell against the arm's own
## candidate windows as background, every held-out known end explained by the
## fold that did not train on it -- and draws, for the held-out known ends
## (set == "positive"), fold means weighted by the ends (or windows) per fold:
##
##   shap_positions_by_family_<TAG>.png  mean |SHAP| along the window per channel
##       family (AA identity, relASA, conservation, AFM, ...), one line per arm;
##       _signed.png the same, signed
##   shap_aa_heatmap_<TAG>.png  AA identity only: when residue X sits at offset p,
##       the signed SHAP of its one-hot cell (cells with < MIN_WINDOWS windows grey)
##   shap_global_<TAG>.png  global SHAP per model: each family's and the top
##       channels' share of the arm's total |SHAP|
##
## Units are logit contributions relative to an average candidate window. Offset
## 0 = the end residue; negative = inside the peptide, positive = past the end.
##
##   LC_ALL=en_US.UTF-8 Rscript arm_compare_shap_plots.R
## env: ARMS       curated,curated_AAswap,msa,msa_gap ("curated" = unsuffixed files)
##      SHAP_DIRS  where shap_*_<arm>.csv are looked for, ":"-separated, first hit
##                 wins (~/AF2_analysis/lf_shap:~/AF2_analysis/lf_headcmp)
##      OUT  ~/AF2_analysis/lf_viz_20261113   TAG  msa_arms   MIN_WINDOWS  4   TOP  12
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
## the shared arm colours (arm_compare_style.R, beside this script)
.here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.here, "arm_compare_style.R"))
A2 <- path.expand("~/AF2_analysis")
ARMS <- strsplit(Sys.getenv("ARMS", "curated,curated_AAswap,msa,msa_gap"), ",")[[1]]
DIRS <- path.expand(strsplit(Sys.getenv("SHAP_DIRS", paste(file.path(A2, "lf_shap"), file.path(A2, "lf_headcmp"),
                                                                 sep = ":")), ":")[[1]])
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113")))
TAG <- Sys.getenv("TAG", "msa_arms")
MINW <- as.integer(Sys.getenv("MIN_WINDOWS", "4")); TOP <- as.integer(Sys.getenv("TOP", "12"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

find_tab <- function(kind, a) {
  f <- sprintf("shap_%s%s.csv", kind, if (a %in% c("curated", "flatten")) "" else paste0("_", a))
  hit <- file.path(DIRS, f); hit <- hit[file.exists(hit)]
  if (length(hit)) fread(hit[1])[, arm := a][] else NULL
}
have <- ARMS[vapply(ARMS, function(a) !is.null(find_tab("channel", a)), TRUE)]
if (length(miss <- setdiff(ARMS, have)))
  message("no attribution yet for ", paste(miss, collapse = ", "), " (run shap_positions.py --arm <arm>)")
if (!length(have)) stop("no attribution tables for any of ", paste(ARMS, collapse = ", "), call. = FALSE)
ARMS <- have
acol <- arm_cols(ARMS)
th <- theme_bw(base_size = 10.5) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
termf <- function(t) factor(t, c("N", "C"), c("N end (offset 0 = first peptide residue)", "C end (offset 0 = last peptide residue)"))

## ---- per position, per family ---------------------------------------------------
pf <- rbindlist(lapply(ARMS, find_tab, kind = "pos_family"), use.names = TRUE)[set == "positive"]
pf <- pf[, .(abs = weighted.mean(sum_of_abs, n), signed = weighted.mean(signed, n)), by = .(arm, term, family, offset)]
fam_rank <- pf[family != "padding", .(tot = sum(abs)), by = family][order(-tot)]$family
pf <- pf[family %in% fam_rank]
pf[, `:=`(family = factor(family, fam_rank), arm = factor(arm, ARMS), term = termf(term))]
for (v in c("abs", "signed")) {
  g <- ggplot(pf, aes(offset, .data[[v]], colour = arm)) +
    geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
    { if (v == "signed") geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) } +
    geom_line(linewidth = 0.7) + geom_point(size = 0.8) +
    scale_colour_manual(values = acol, name = NULL) +
    facet_grid(family ~ term, scales = "free_y", switch = "y") +
    labs(x = "offset from the peptide end (residues)",
         y = if (v == "abs") "mean |SHAP| of the family's cells (logit)" else "mean signed SHAP (logit)",
         title = sprintf("Where along the window each channel family matters (%s SHAP, held-out known ends)",
                         if (v == "abs") "absolute" else "signed"),
         subtitle = paste0("Expected-gradients SHAP on each member's logit against the arm's own candidate windows; every known end\n",
                           "explained by the fold that did not train on it. Families ordered by total |SHAP|; per-panel y scales.")) +
    th + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1))
  ggsave(file.path(OUT, sprintf("shap_positions_by_family_%s%s.png", TAG, if (v == "signed") "_signed" else "")), g,
         width = 13, height = 1.55 * length(fam_rank) + 1.6, dpi = 170, bg = "white")
}
fwrite(pf, file.path(OUT, sprintf("shap_positions_by_family_%s.csv", TAG)))

## ---- AA identity heatmap: letter x offset ----------------------------------------
pl <- rbindlist(lapply(ARMS, find_tab, kind = "pos_letter"), use.names = TRUE)
if (nrow(pl)) {
  pl <- pl[set == "positive", .(signed = weighted.mean(signed, n_windows), n_windows = sum(n_windows)),
           by = .(arm, term, letter, offset)]
  ## physico-chemical order: basic, acidic, amide / hydroxyl, G / P, C, aliphatic, aromatic
  lev <- c("K", "R", "H", "D", "E", "N", "Q", "S", "T", "G", "P", "C", "U", "A", "V", "L", "I", "M", "F", "W", "Y")
  pl <- pl[letter %in% lev]
  pl[, `:=`(letter = factor(letter, rev(lev)), arm = factor(arm, intersect(ARMS, unique(arm))), term = termf(term))]
  pl[n_windows < MINW, signed := NA_real_]
  lim <- quantile(abs(pl$signed), 0.99, na.rm = TRUE)
  g <- ggplot(pl, aes(offset, letter, fill = pmax(pmin(signed, lim), -lim))) +
    geom_tile() +
    geom_vline(xintercept = c(-0.5, 0.5), colour = "grey30", linewidth = 0.25, linetype = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         na.value = "grey88", name = "signed SHAP\n(logit)") +
    facet_grid(arm ~ term, scales = "free_x") +
    labs(x = "offset from the peptide end (residues)", y = NULL,
         title = "Amino-acid identity: what residue X at offset p contributes (held-out known ends)",
         subtitle = sprintf(paste0("Signed SHAP of the one-hot cell of the residue that is there, mean over the windows that have it ",
                                   "(fold means weighted by windows); grey = fewer than %d windows. Red = raises the end score."), MINW)) +
    th + theme(panel.grid = element_blank(), legend.position = "right")
  ggsave(file.path(OUT, sprintf("shap_aa_heatmap_%s.png", TAG)), g, width = 13, height = 2.9 * uniqueN(pl$arm) + 1.3,
         dpi = 170, bg = "white")
  fwrite(pl, file.path(OUT, sprintf("shap_aa_heatmap_%s.csv", TAG)))
}

## ---- global SHAP per model -------------------------------------------------------
ch <- rbindlist(lapply(ARMS, find_tab, kind = "channel"), use.names = TRUE)[set == "positive" & channel != "padding"]
ch <- ch[, .(abs = weighted.mean(sum_of_abs, n), signed = weighted.mean(signed, n)), by = .(arm, term, channel, family)]
ch[, share := abs / sum(abs), by = .(arm, term)]
fwrite(ch[order(arm, term, -share)], file.path(OUT, sprintf("shap_global_%s.csv", TAG)))
fs <- ch[, .(share = sum(share)), by = .(arm, term, family)]
fs[, family := factor(family, rev(c(fam_rank, setdiff(unique(fs$family), fam_rank))))]
fs[, `:=`(arm = factor(arm, rev(ARMS)), term = factor(term, c("N", "C"), c("N end", "C end")))]
fcols <- c("AA identity (one-hot)" = "#E69F00", relASA = "#0072B2", conservation = "#009E73", AFM = "#56B4E9",
           "secondary structure" = "#CC79A7", dihedrals = "#D55E00", "H-bond energies" = "#999999",
           "AA properties" = "#F0E442", "PLM (ESM-C)" = "#6A3D9A", "motif indicators" = "#B15928", other = "black")
p1 <- ggplot(fs, aes(share, arm, fill = family)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(data = fs[share >= 0.07], aes(label = sprintf("%.0f%%", 100 * share)),
            position = position_stack(vjust = 0.5), size = 2.7) +
  scale_fill_manual(values = fcols, name = NULL) +
  scale_x_continuous(labels = scales::label_percent(), expand = c(0, 0)) +
  facet_wrap(~term, ncol = 1) +
  labs(x = "share of the model's total |SHAP|", y = NULL, title = "By channel family") + th +
  guides(fill = guide_legend(ncol = 3, reverse = TRUE))
tc <- ch[order(-share), head(.SD, TOP), by = .(arm, term)]
tc[, lab := sub("^plm_", "", channel)]
tc[, `:=`(arm = factor(arm, ARMS), term = factor(term, c("N", "C"), c("N end", "C end")))]
tc[, key := paste(arm, term, lab, sep = "|")][, key := factor(key, rev(unique(key)))]
p2 <- ggplot(tc, aes(share, key, fill = family)) +
  geom_col(width = 0.7) +
  scale_y_discrete(labels = function(k) sub("^.*\\|", "", k)) +
  scale_fill_manual(values = fcols, guide = "none") +
  scale_x_continuous(labels = scales::label_percent()) +
  facet_wrap(term ~ arm, scales = "free_y", nrow = 2) +
  labs(x = "share of the model's total |SHAP|", y = NULL, title = sprintf("Top %d channels", TOP)) + th
g <- (p1 | p2) + plot_layout(widths = c(1, 2.2)) +
  plot_annotation(title = "Global SHAP per model: what each trunk's held-out end scores rest on",
                  subtitle = paste0("Sum over the window of |SHAP| per channel, mean over the held-out known ends (fold means weighted by ends); ",
                                    "padding left out.\nThe AA one-hot spreads over 21 channels, so its family total is large ",
                                    "while each letter is small."),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUT, sprintf("shap_global_%s.png", TAG)), g, width = 4.5 + 2.6 * length(ARMS), height = 8.5,
       dpi = 170, bg = "white")
message("wrote ", file.path(OUT, sprintf("shap_{positions_by_family,aa_heatmap,global}_%s.png", TAG)))
