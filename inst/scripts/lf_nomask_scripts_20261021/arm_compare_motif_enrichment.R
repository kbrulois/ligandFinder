## Amino-acid enrichment around the known peptide ends, on the same letter x
## offset grid as the AA-identity SHAP heatmap, and drawn beneath it (Kevin,
## 2026-11-13: "make a motif enrichment plot from pepend_windows_aligned_gap_C ...
## and add it below the shap values from shap_aa_heatmap_msa_vs_msa_gap").
##
## The data are what 10_7b_plot_pepend_windows.R draws one row per window of:
## every known end (lf_pepend_cv5_knowns_<T>.csv, the anchors every arm is scored
## at) aligned at its end residue, offsets -27..+8 for C and -7..+28 for N.
## Enrichment of residue X at offset p = log2((observed + 2) / (expected + 2)),
## expected = (ends with a residue there) x (X's share of the mature secretome:
## every precursor from its first mature residue on). The +2 on both sides keeps
## a letter seen 0 or 1 times in ~90 ends from reading as a strong depletion.
## The top rows are arm_compare_shap_plots.R's heatmap (its CSV), so the figure
## says side by side what is enriched at real ends and what each model uses.
##
##   LC_ALL=en_US.UTF-8 Rscript arm_compare_motif_enrichment.R
## env: TAG   the SHAP heatmap to draw above (OUT/shap_aa_heatmap_<TAG>.csv)  msa_vs_msa_gap
##      OUT   ~/AF2_analysis/lf_viz_20261113     SEQS  a named-vector rds of the
##            precursor sequences (default OUT/msa_arms_data/prec_seqs.rds, cut from
##            the v2 feature cache when absent)
## Writes <OUT>/motif_enrichment_<TAG>.png (SHAP above, enrichment below),
## <OUT>/motif_enrichment_known_ends.png (enrichment alone) and its csv.
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
A2 <- path.expand("~/AF2_analysis")
OUT <- path.expand(Sys.getenv("OUT", file.path(A2, "lf_viz_20261113")))
TAG <- Sys.getenv("TAG", "msa_vs_msa_gap")
SEQS <- path.expand(Sys.getenv("SEQS", file.path(OUT, "msa_arms_data", "prec_seqs.rds")))
if (!file.exists(SEQS)) {
  dir.create(dirname(SEQS), showWarnings = FALSE, recursive = TRUE)
  fc <- readRDS(file.path(A2, "lf_pepend_residue_cache_v2.rds"))
  saveRDS(setNames(fc$prec$seq, fc$prec$accession), SEQS); rm(fc); invisible(gc())
}
seqs <- readRDS(SEQS)
prec <- fread(file.path(A2, "pepend_locus", "prec_order.csv"))
lev <- c("K", "R", "H", "D", "E", "N", "Q", "S", "T", "G", "P", "C", "A", "V", "L", "I", "M", "F", "W", "Y")   # no U: too rare to read
OFF <- list(C = -27:8, N = -7:28)

## background: the mature secretome's composition
mat <- prec[accession %in% names(seqs)]
bg <- table(factor(unlist(strsplit(substring(seqs[mat$accession], pmax(mat$n_prot, 1L)), "")), lev))
bg <- bg / sum(bg)

en <- rbindlist(lapply(c("N", "C"), function(tm) {
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)), select = c("accession", "anchor"))
  kn <- kn[accession %in% names(seqs)]
  aa <- rbindlist(lapply(seq_len(nrow(kn)), function(i) {
    s <- seqs[[kn$accession[i]]]; p <- kn$anchor[i] + OFF[[tm]]
    ok <- p >= 1 & p <= nchar(s)
    data.table(offset = OFF[[tm]][ok], letter = substring(s, p[ok], p[ok]))
  }))
  n_off <- aa[, .(n = .N), by = offset]
  obs <- aa[letter %in% lev, .(obs = .N), by = .(offset, letter)]
  g <- CJ(offset = OFF[[tm]], letter = lev)[n_off, on = "offset"]
  g <- obs[g, on = .(offset, letter)][is.na(obs), obs := 0L]
  g[, expected := n * as.numeric(bg[letter])]
  g[, `:=`(term = tm, n_ends = nrow(kn), log2_enrich = log2((obs + 2) / (expected + 2)))]
}))
fwrite(en[order(term, offset, letter)], file.path(OUT, "motif_enrichment_known_ends.csv"))

termf <- function(t) factor(t, c("N", "C"), c("N end (offset 0 = first peptide residue)", "C end (offset 0 = last peptide residue)"))
en[, `:=`(letter = factor(letter, rev(lev)), panel = sprintf("known ends\n(N %d, C %d)", en[term == "N"]$n_ends[1],
                                                               en[term == "C"]$n_ends[1]), term = termf(term))]
elim <- max(abs(en$log2_enrich))
th <- theme_bw(base_size = 10.5) +
  theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"), panel.grid = element_blank(),
        legend.position = "right")
pe <- ggplot(en, aes(offset, letter, fill = log2_enrich)) +
  geom_tile() +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey30", linewidth = 0.25, linetype = 3) +
  scale_fill_gradient2(low = "#5e3c99", mid = "#f7f7f7", high = "#e66101", midpoint = 0, limits = c(-elim, elim),
                       name = "log2 enrichment\nvs the mature\nsecretome") +
  facet_grid(panel ~ term, scales = "free_x") +
  labs(x = "offset from the peptide end (residues)", y = NULL) + th
ggsave(file.path(OUT, "motif_enrichment_known_ends.png"),
       pe + labs(title = "Amino-acid enrichment around the known peptide ends",
                 subtitle = "log2((observed + 2) / (expected + 2)); expected from the mature secretome's composition. Orange = enriched at real ends."),
       width = 13, height = 4.2, dpi = 170, bg = "white")

## SHAP heatmap above, enrichment below
sh_p <- file.path(OUT, sprintf("shap_aa_heatmap_%s.csv", TAG))
if (file.exists(sh_p)) {
  sh <- fread(sh_p)
  arms <- unique(sh$arm)
  ## the csv carries the facet label ("N end (offset 0 = ...)"); back to N / C
  sh[, `:=`(letter = factor(letter, rev(lev)), arm = factor(arm, arms), term = termf(sub(" end.*", "", term)))]
  lim <- quantile(abs(sh$signed), 0.99, na.rm = TRUE)
  ps <- ggplot(sh, aes(offset, letter, fill = pmax(pmin(signed, lim), -lim))) +
    geom_tile() +
    geom_vline(xintercept = c(-0.5, 0.5), colour = "grey30", linewidth = 0.25, linetype = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         na.value = "grey88", name = "signed SHAP\n(logit)") +
    facet_grid(arm ~ term, scales = "free_x") +
    labs(x = NULL, y = NULL) + th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  g <- ps / pe + plot_layout(heights = c(length(arms), 1)) +
    plot_annotation(
      title = "What each model uses (SHAP, top) against what is enriched at real ends (bottom)",
      subtitle = paste0("Top: signed SHAP of the residue's one-hot cell at each offset, held-out known ends (grey = fewer than 4 windows). ",
                        "Bottom: the same ends' amino-acid enrichment\nagainst the mature secretome, log2((obs + 2) / (exp + 2)). ",
                        "Red / orange = raises the score / enriched at real ends."),
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(OUT, sprintf("motif_enrichment_%s.png", TAG)), g, width = 13, height = 2.9 * length(arms) + 4.2,
         dpi = 170, bg = "white")
  message("wrote ", file.path(OUT, sprintf("motif_enrichment_%s.png", TAG)))
} else message("no ", sh_p, ": drew the enrichment alone")
print(en[order(-log2_enrich)][1:12, .(term, offset, letter, obs, expected = round(expected, 1), log2_enrich = round(log2_enrich, 2))])

## ---- the same as sequence logos (Kevin, 2026-11-13: "in the style of letters
## with different height") -------------------------------------------------------
## A Kullback-Leibler logo: letter height = q * log2(q / background), q the
## letter's frequency among the ends at that offset with a 2-letter pseudocount
## spread by the background. Enriched letters stack above the axis, depleted ones
## below, and a column's net height is its information against the secretome; a
## rare letter, which ~90 ends say little about, stays small. Under the SHAP
## heatmap the two share the offset axis, column for column.
suppressMessages(library(ggseqlogo))
en[, tm := ifelse(grepl("^N", as.character(term)), "N", "C")]
en[, q := (obs + 2 * as.numeric(bg[as.character(letter)])) / (n + 2)]
en[, kl := q * log2(q / as.numeric(bg[as.character(letter)]))]
## zoom = TRUE (Kevin: "expand the scale so [the] dibasic site or G-dibasic
## residues are cut off and the surround residues are full height"): the y axis
## stops at the tallest column OUTSIDE the cleavage motif (N: -2, -1; C: +1..+3),
## so the motif columns are clipped -- their full height printed on them -- and
## the context letters fill the panel
logo <- function(t, xlab = TRUE, zoom = FALSE) {
  m <- dcast(en[tm == t], letter ~ offset, value.var = "kl")
  M <- as.matrix(m[, -1]); rownames(M) <- as.character(m$letter)
  off <- as.integer(colnames(M)); cut <- which(off == 0) + if (t == "C") 0.5 else -0.5
  suppressMessages(
    ggseqlogo(M, method = "custom", seq_type = "aa") +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
      geom_vline(xintercept = cut, colour = "grey30", linewidth = 0.3, linetype = 3) +
      scale_x_continuous(breaks = which(off %% 5 == 0), labels = off[off %% 5 == 0], expand = c(0, 0),
                         limits = c(0.5, length(off) + 0.5)) +
      labs(x = if (xlab) "offset from the peptide end (residues)" else NULL, y = "bits (KL)",
           title = sprintf("known %s ends (n = %d)", t, en[tm == t]$n_ends[1])) +
      theme_bw(base_size = 10.5) +
      theme(legend.position = "none", panel.grid = element_blank(), plot.title = element_text(size = 10, face = "bold")) +
      if (zoom) {
        cs <- en[tm == t, .(pos = sum(pmax(kl, 0)), neg = sum(pmin(kl, 0))), by = offset]
        motif <- if (t == "N") c(-2L, -1L) else 1:3
        cap <- 1.06 * max(cs[!offset %in% motif]$pos); lo <- 1.06 * min(cs$neg)
        cl <- cs[pos > cap]
        list(coord_cartesian(ylim = c(lo, cap), expand = FALSE),
             annotate("label", x = match(cl$offset, off), y = cap, vjust = 1.15, label = sprintf("%.1f", cl$pos),
                      size = 2.5, fontface = "bold", label.size = 0, label.padding = unit(0.08, "lines"), fill = "white"))
      })
}
gl <- (logo("N") | logo("C")) +
  plot_annotation(title = "Sequence logos of the known peptide ends against the mature secretome",
                  subtitle = paste0("Kullback-Leibler logo: height = q log2(q / background), q the frequency at that offset (2-letter pseudocount). ",
                                    "Above the axis = enriched at real ends, below = depleted. Dotted line = the cut."),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUT, "motif_logo_known_ends.png"), gl, width = 15, height = 3.6, dpi = 170, bg = "white")
glz <- (logo("N", zoom = TRUE) | logo("C", zoom = TRUE)) +
  plot_annotation(title = "Sequence logos of the known peptide ends, zoomed past the cleavage motif",
                  subtitle = paste0("As motif_logo_known_ends.png, with the y axis stopped at the tallest column outside the motif ",
                                    "(N: -2, -1; C: +1..+3): those columns are clipped, their full height in bits printed on them."),
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
ggsave(file.path(OUT, "motif_logo_known_ends_zoom.png"), glz, width = 15, height = 3.6, dpi = 170, bg = "white")

if (exists("sh") && nrow(sh)) {
  ## the SHAP heatmap one terminus per column, so patchwork lines each column up
  ## with its logo; x is the same offset axis, no expansion on either
  shap_col <- function(t) {
    d <- sh[grepl(paste0("^", t), as.character(term))]
    ggplot(d, aes(offset, letter, fill = pmax(pmin(signed, lim), -lim))) +
      geom_tile() +
      geom_vline(xintercept = if (t == "C") 0.5 else -0.5, colour = "grey30", linewidth = 0.3, linetype = 3) +
      scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                           na.value = "grey88", name = "signed SHAP\n(logit)") +
      scale_x_continuous(expand = c(0, 0)) +
      facet_grid(arm ~ term) +
      labs(x = NULL, y = NULL) + th +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  g <- (shap_col("N") | shap_col("C")) / (logo("N") | logo("C")) +
    plot_layout(heights = c(length(arms), 0.9), guides = "collect") +
    plot_annotation(
      title = "What each model uses (SHAP, top) against what is enriched at real ends (logos, bottom)",
      subtitle = paste0("Top: signed SHAP of the residue's one-hot cell, held-out known ends (grey = fewer than 4 windows; red = raises the score). ",
                        "Bottom: Kullback-Leibler logos of the same ends\nagainst the mature secretome (letters above the axis enriched, below depleted). ",
                        "Columns line up offset for offset; dotted line = the cut."),
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(OUT, sprintf("motif_logo_%s.png", TAG)), g, width = 15, height = 2.9 * length(arms) + 3.6,
         dpi = 170, bg = "white")
  gz <- (shap_col("N") | shap_col("C")) / (logo("N", zoom = TRUE) | logo("C", zoom = TRUE)) +
    plot_layout(heights = c(length(arms), 0.9), guides = "collect") +
    plot_annotation(
      title = "What each model uses (SHAP, top) against what is enriched at real ends (logos zoomed past the motif, bottom)",
      subtitle = paste0("As motif_logo_", TAG, ".png, with the logos' y axis stopped at the tallest column outside the cleavage motif ",
                        "(N: -2, -1; C: +1..+3);\nthose columns are clipped, their full height in bits printed on them. ",
                        "Columns line up offset for offset; dotted line = the cut."),
      theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(OUT, sprintf("motif_logo_%s_zoom.png", TAG)), gz, width = 15, height = 2.9 * length(arms) + 3.6,
         dpi = 170, bg = "white")
  message("wrote ", file.path(OUT, sprintf("motif_logo_%s.png", TAG)))
}
