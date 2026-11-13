#!/usr/bin/env Rscript
## ---- the MSA conservation channel along example precursors --------------------
## Kevin, 2026-11-12: "plots new conservation score for a few example proteins".
## The msacons arms replaced the Aminode channel (cons_rs_n) with plm_cons_msa,
## a per-residue score from each precursor's ColabFold MSA. Two figures:
##
##   msacons_examples     six precursors, both channels on one 0-1 axis, the
##                        known peptides and their ends (G + dibasic C ends on
##                        the G, as the _padg windows anchor them)
##   msacons_end_summary  both channels at residues near known ends against
##                        every other scored residue, AUC with a gene bootstrap
##
## EXAMPLES, BY RULE (lf_msacons_select in R/plot_msa_conservation.R):
##   3 classics -- of ADCYAP1, GCG, POMC, NPY, CALCA, PENK: the most G + dibasic
##     trained C ends (where the _padg convention moved the end onto the G),
##     then the most trained ends, then the deeper MSA;
##   1 never-seen gene -- the most ends in the never-seen peptide set (< 60 aa),
##     ties to the deeper MSA;
##   1 more never-seen peptide gene, so the excluded-by-design band has an
##     example -- the most ends among those with excluded residues, ties to the
##     shorter precursor;
##   1 shallow MSA -- the lowest neff of any gene carrying a trained or
##     never-seen end. It is also drawn WITHOUT the evidence prior (dashed),
##     recomputed by msa_conservation.py's own score_one, after prior = 10 has
##     reproduced the parquet.
## On 2026-11-12 that is POMC, ADCYAP1, GCG | INS | TGFA | MT-RNR2.
##
##   LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/msa_cons/plot_msa_conservation.R
##   ... --genes POMC,INS (skip the rule)  --k 3  --out-dir ~/AF2_analysis/lf_viz_20261112
##
## Output, in --out-dir: msacons_examples.{png,pdf}, msacons_end_summary.{png,pdf}
## and msacons_end_summary.csv. ~4 GB of RAM at peak (the residue cache read twice: once
## to build, once to spot-check what is drawn), ~20 s.
## ------------------------------------------------------------------------------

suppressMessages({ library(data.table); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
out_dir <- path.expand(.opt("--out-dir", "~/AF2_analysis/lf_viz_20261112"))
K       <- as.integer(.opt("--k", "3"))          # residues either side of an end
genes   <- .opt("--genes", "")
dpi     <- as.integer(.opt("--dpi", "170"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..", "..")) else getwd()
## pepend_windows.R: lf_pepend_motif_vec (amid_g); plot_proteins_win_new.R:
## nn_class_cols; viz_theme.R: LF_VIZ + lf_viz_theme; pepend_tracks.R: LF_TERM_COLS
for (f in c("pepend_windows.R", "plot_proteins_win_new.R", "viz_theme.R", "pepend_tracks.R",
            "plot_msa_conservation.R"))
  source(file.path(ROOT, "R", f))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cols <- lf_msacons_cols()

## ---- 1. both channels on every cached residue ------------------------------------
x <- lf_msacons_residues()
res <- x$res; prec <- x$prec; per <- x$per
message(sprintf("%s residues, %s excluded by design; plm_cons_msa = min-max of raw [%.5f, %.4f] over %s scored residues (max |dev| %.1e)",
                format(nrow(res), big.mark = ","), format(sum(res$excluded), big.mark = ","),
                x$scale$lo, x$scale$hi, format(x$scale$n_covered, big.mark = ","), x$scale$max_dev))
ends <- lf_msacons_ends(prec)
message("known ends: ", paste(sprintf("%s %d", levels(ends$set), tabulate(ends$set)), collapse = ", "))

## ---- 2. the examples ---------------------------------------------------------------
sel <- lf_msacons_select(ends, res, prec, per)
if (nzchar(genes)) {
  g <- strsplit(genes, ",")[[1]]
  acc <- prec$accession[match(g, prec$gene)]
  if (anyNA(acc)) stop("not in the residue cache: ", paste(g[is.na(acc)], collapse = ", "))
  sel <- rbind(sel[accession %in% acc],
               data.table(accession = acc[!acc %in% sel$accession], gene = g[!acc %in% sel$accession],
                          group = "asked for", why = "asked for (--genes)"), fill = TRUE)
  sel <- sel[match(acc, accession)]
}
sel[, `:=`(depth = per$depth[match(accession, per$accession)], neff = per$neff[match(accession, per$accession)],
           source = per$source[match(accession, per$accession)])]
print(sel[, .(gene, accession, group, depth, neff = round(neff, 1), source, why)])

ll <- readRDS(file.path(ROOT, "inst", "extdata", "ligand_list.rds"))
spans <- lf_msacons_spans(sel$accession, ll)

## every example rescored by the score's own code: prior = 10 must reproduce the
## parquet (else nothing is drawn); prior = 0 says what the prior does. Only the
## shallow example gets the no-prior line -- on a deep MSA it would sit on the
## solid one -- and the others' largest shift goes in the caption instead
script_dir <- file.path(ROOT, "inst", "scripts", "msa_cons")
noprior <- list(); prior_shift <- c()
for (a in sel$accession) {
  s <- prec[accession == a]$seq; src <- per[accession == a]$source
  r10 <- lf_msacons_rescore(a, s, src, 10, script_dir)
  dev <- max(abs(r10$jsd - res[accession == a][order(index)]$jsd))
  if (dev > 1e-6) stop("score_one(prior = 10) does not reproduce msa_cons_raw.parquet for ", a)
  r0 <- lf_msacons_rescore(a, s, src, 0, script_dir)
  prior_shift[a] <- max(abs(r0$jsd - r10$jsd))
  message(sprintf("%-8s score_one(prior 10) vs parquet max |diff| %.1e; mean raw %.3f with the prior, %.3f without, largest shift %.3f",
                  sel$gene[sel$accession == a], dev, mean(r10$jsd), mean(r0$jsd), prior_shift[a]))
  if (a %in% sel[group == "shallow MSA"]$accession) {
    noprior[[a]] <- data.frame(index = r0$index,
                               value = pmin(pmax((r0$jsd - x$scale$lo) / (x$scale$hi - x$scale$lo), 0), 1))
    sel[accession == a, why := sprintf("%s; mean raw score %.2f with the prior, %.2f without it (dashed)",
                                       why, mean(r10$jsd), mean(r0$jsd))]
  }
}
deep <- setdiff(sel$accession, names(noprior))

## ---- 3. figure 1: one panel per example ------------------------------------------
## wrap prose to the figure's width: a paragraph per element of ... (NULL = none)
wrap <- function(width, ...) paste(vapply(Filter(Negate(is.null), list(...)),
                                          function(t) paste(strwrap(t, width), collapse = "\n"),
                                          character(1)), collapse = "\n")
plots <- lapply(seq_len(nrow(sel)), function(i) {
  a <- sel$accession[i]
  lf_msacons_protein_plot(
    r = res[accession == a], spans = spans[accession == a],
    ends = ends[accession == a & set %in% LF_MSACONS_SETS],
    meta = list(gene = sel$gene[i], accession = a, depth = sel$depth[i], neff = sel$neff[i],
                source = sel$source[i], why = sel$why[i]),
    noprior = noprior[[a]], x_title = if (i == nrow(sel)) "residue (precursor numbering)",
    cols = cols)
})
names(plots) <- sel$accession

## spot-check: what the figure DRAWS against the parquets, read again here and
## independently of lf_msacons_residues -- the channel straight from the arm's
## parquet, the old one straight from the residue cache
chk_new <- as.data.table(arrow::read_parquet(path.expand(LF_MSACONS_PATHS$channel),
                                             col_select = c("accession", "index", "AA", "plm_cons_msa")))[accession %in% sel$accession]
fc <- readRDS(path.expand(LF_MSACONS_PATHS$cache))
set.seed(7)
spot <- rbindlist(lapply(sel$accession, function(a) {
  p <- plots[[a]]
  li <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))[1]
  ld <- as.data.table(layer_data(p, li))
  drawn_new <- ld[colour == cols$new & linetype == "solid"]
  drawn_old <- ld[colour == cols$old]
  r <- res[accession == a & excluded == FALSE]
  ## a random scored residue, the channel's maximum, and an end anchor
  idx <- unique(c(sample(r$index, 1), r$index[which.max(r$new)],
                  intersect(ends[accession == a]$anchor, r$index)[1]))
  idx <- idx[!is.na(idx)]
  rbindlist(lapply(idx, function(i) data.table(
    gene = sel$gene[sel$accession == a], index = i,
    AA = chk_new[accession == a & index == i]$AA,
    drawn_new = drawn_new[x == i]$y, parquet_new = chk_new[accession == a & index == i]$plm_cons_msa,
    drawn_old = if (nrow(drawn_old[x == i])) drawn_old[x == i]$y else NA_real_,
    cache_old = unname(fc$feats[[a]][i, "cons_rs_n"]))))
}))
rm(fc)
print(spot)
## the old line is drawn only where cons_rs_n is not 0; a 0 is drawn as a cross
ok_old <- ifelse(is.na(spot$drawn_old), spot$cache_old == 0, abs(spot$drawn_old - spot$cache_old) < 1e-9)
stopifnot(length(spot$drawn_new) == nrow(spot), max(abs(spot$drawn_new - spot$parquet_new)) < 1e-6, all(ok_old))

## the end marks as DRAWN: residue, glyph, and grey = in neither end table
marks <- rbindlist(lapply(sel$accession, function(a) {
  p <- plots[[a]]
  li <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint") && "mark" %in% names(l$data),
                     logical(1)))
  if (!length(li)) return(NULL)
  ld <- as.data.table(layer_data(p, li[1]))
  ld[, .(gene = sel$gene[sel$accession == a],
         drawn = paste(sprintf("%s%d%s", shape, as.integer(x), ifelse(colour == cols$other_end, "(grey)", "")),
                       collapse = " "))]
}))
print(marks)

## the old channel's highest value per example: its ceiling is the precursor's
## species_limit, how deep in time its Aminode alignment reaches (an invariant
## residue reads 1.00 in INS, 0.74 in GCG; see the header of
## R/plot_msa_conservation.R), so the subtitle says so with the drawn numbers,
## and with how low that ceiling sits across all precursors
top_old <- res[accession %in% sel$accession & excluded == FALSE, .(top = max(old)), by = accession]
top_old <- top_old[match(sel$accession, accession)][, gene := sel$gene]
top_txt <- paste(c(top_old[top > 0, sprintf("%s %.2f", gene, top)],
                   if (any(top_old$top == 0)) sprintf("%s none", paste(top_old[top == 0]$gene, collapse = ", "))),
                 collapse = ", ")
top_all <- res[excluded == FALSE & old > 0, .(top = max(old)), by = accession]$top
## the listed peptides drawn with grey ends, by gene in figure order, for the caption
unl <- lf_msacons_unlisted(spans, ends, prec)
unl[, gene := factor(prec$gene[match(accession, prec$accession)], levels = sel$gene)]
unl_txt <- unl[order(gene), .(n = paste(name, collapse = ", ")), by = gene][, paste(gene, n, collapse = "; ")]

hts <- vapply(plots, function(p) attr(p, "y_top"), numeric(1))
unit_in <- 1.2                                   # inches per channel unit
key_in  <- 1.0                                   # the legend row
fig_h <- sum(hts) * unit_in + 0.72 * length(plots) + key_in + 2.8
key <- patchwork::wrap_elements(full = lf_msacons_legend(cols)) +
  theme(plot.background = element_rect(fill = LF_VIZ$surface, colour = NA),
        plot.margin = margin(10, 4, 0, 4))
fig <- patchwork::wrap_plots(c(list(key), plots),
                             ncol = 1, heights = c(key_in / unit_in, hts)) +
  patchwork::plot_annotation(
    title = sprintf("The new conservation channel along %d precursors: MSA conservation (blue) against the Aminode channel it replaced (grey)",
                    nrow(sel)),
    subtitle = wrap(185,
      paste0(
        "NEW = plm_cons_msa, exactly what the msacons arms read: the raw MSA score min-max scaled over all ",
        format(x$scale$n_covered, big.mark = ","), sprintf(" scored residues (raw %.5f -> 0, %.3f -> 1). ", x$scale$lo, x$scale$hi),
        "Raw score = Capra & Singh (2007) Jensen-Shannon divergence between each column of the precursor's ColabFold MSA ",
        "(Henikoff-weighted, insertions deleted) and the BLOSUM62 background, after mixing in 10 background ",
        "pseudo-sequences against an effective count of neff x (1 - gap fraction), then multiplied by (1 - gap fraction). ",
        "Higher = more conserved; per residue, no smoothing."),
      ## Reviewer, 2026-11-12: the old channel's ceiling is species_limit (9_expand_by_residue.R),
      ## the depth in time of the Aminode alignment -- say so, rather than "set per precursor"
      paste0(
        "OLD = cons_rs_n: the Aminode evolutionary rate, reversed and rescaled to 0-1 across genes, then multiplied by ",
        "the precursor's species_limit, how far back in time its Aminode alignment reaches (1 = the deepest species; ",
        sprintf("%.0f%% of the %s precursors with an Aminode value never exceed 0.5). ",
                100 * mean(top_all < 0.5), format(length(top_all), big.mark = ",")),
        "Linear in the rate within a panel, but its ceiling is that depth, so read its shape within a panel, not its ",
        "level across panels. Highest value drawn: ", top_txt, "; its flat stretches are invariant residues at that ceiling.")),
    caption = wrap(205,
      paste0(
        if (nzchar(genes)) "Examples asked for with --genes; one the rule also picks keeps the rule's reason in its subtitle. "
        else paste0(
          "Examples by rule. CLASSICS: of ADCYAP1, GCG, POMC, NPY, CALCA and PENK, the three with the most G + dibasic ",
          "trained C ends (where the _padg windows moved the end onto the G), then the most trained ends, then neff. ",
          "NEVER-SEEN (genes no fold trained on): the most ends in the never-seen peptide set (< 60 aa), ties to the deeper ",
          "MSA; and, so the excluded-by-design band has an example, the most never-seen peptide ends among precursors with ",
          "excluded residues, ties to the shorter one. SHALLOW: the lowest neff of any gene with a known end. "),
        if (length(noprior)) paste0(
          "Dashed = the same score with the prior switched off (msa_conservation.py's score_one at prior 0, same min-max)",
          if (length(deep)) sprintf("; on the other %d the prior moves no residue's raw score by more than %.3f, so they get no dashed line",
                                    length(deep), max(prior_shift[deep])), ".")),
      paste0("Ends: trained = lf_pepend_cv5_knowns_{N,C}_padg.csv; never-seen = cv5_unseen_ends_padg.csv, C ends moved onto ",
             "the G where G + dibasic follows. Spans: ligand_list.rds (mature peptides). MSA depth = sequences incl. the ",
             "query; neff = exp(entropy of the Henikoff weights)."),
      if (nzchar(unl_txt)) paste0(
        "Grey ends (", unl_txt, ") belong to listed peptides that neither end table holds: the trained ends are 10_7's ",
        "DOCKED peptides (knowns.rds), so on a trained gene a listed peptide that was never docked has none."),
      sprintf("Built %s by inst/scripts/msa_cons/plot_msa_conservation.R.", format(Sys.Date()))),
    theme = lf_viz_theme(base_size = 9) +
      theme(plot.title = element_text(size = 11, face = "bold", colour = cols$ink),
            plot.subtitle = element_text(size = 8.5, colour = cols$ink2, lineheight = 1.15),
            plot.caption = element_text(size = 7.5, colour = cols$ink2, hjust = 0, lineheight = 1.15)))
f1 <- file.path(out_dir, "msacons_examples")
ggsave(paste0(f1, ".png"), fig, width = 12, height = fig_h, dpi = dpi, device = ragg::agg_png,
       limitsize = FALSE, bg = LF_VIZ$surface)
ggsave(paste0(f1, ".pdf"), fig, width = 12, height = fig_h, device = cairo_pdf, bg = LF_VIZ$surface)
message("wrote ", f1, ".{png,pdf}")

## ---- 4. figure 2: both channels near known ends vs elsewhere ----------------------
es <- lf_msacons_end_stats(res, ends, k = K)
print(es$stats[, .(set, channel, n_ends, n_res, n_genes, auc = round(auc, 3), lo = round(lo, 3), hi = round(hi, 3),
                   auc_same = round(auc_same, 3), auc_cov = round(auc_cov, 3),
                   median_near = round(median_near, 3), median_bg = round(median_bg, 3))])
fwrite(es$stats, file.path(out_dir, "msacons_end_summary.csv"))
st <- es$stats
a_of <- function(s, ch, col = "auc") st[set == s & channel == ch][[col]]
## a channel's contrast at known ends is not the model's gain: what the CNN
## trained with each channel scored (the two _padg arms, normalised PR AUC =
## AP at 1% prevalence over the full-proteome scan, metrics_final.R)
fm_p <- path.expand(LF_MSACONS_PATHS$arm_metrics)
model_txt <- if (file.exists(fm_p)) {
  fm <- fread(fm_p)[head == "cnn"]
  arms <- c(old = "curated_AAswap_padg", new = "curated_AAswap_msacons_padg")
  npr <- function(s, tm) vapply(arms, function(tr) {
    v <- fm[trunk == tr & set == s & term == tm]$norm_pr_auc; if (length(v) == 1) v else NA_real_
  }, numeric(1))
  sets <- c(`cross-validated` = "cross-validated", `never-seen peptides` = "never seen: peptides",
            `large products` = "never seen: large products")
  v <- lapply(sets, function(s) rbind(C = npr(s, "C"), N = npr(s, "N")))
  if (anyNA(unlist(v))) { message("arm metrics incomplete in ", fm_p, "; caption line skipped"); NULL }
  else paste0("A channel's contrast is not the model's gain. The CNN trained with each channel (arms ", arms[["old"]],
              " vs ", arms[["new"]], "; normalised PR AUC, old vs new): ",
              paste(vapply(names(v), function(s) sprintf("%s C %.3f vs %.3f, N %.3f vs %.3f", s,
                                                         v[[s]]["C", "old"], v[[s]]["C", "new"],
                                                         v[[s]]["N", "old"], v[[s]]["N", "new"]), character(1)),
                    collapse = "; "),
              " (lf_headcmp/unseen_padg/final_metrics.csv).")
}
p2 <- lf_msacons_summary_plot(es, cols) +
  labs(title = sprintf("Each channel at residues within %d of a known end, against every other scored residue", K),
       subtitle = wrap(150, paste0(
         "Same residues in both panels: scored = not excluded by design (", format(es$n_scored, big.mark = ","),
         "); \"other\" = not within ", K, " of any end shown. AUC = P(a residue near an end outranks an \"other\" ",
         "residue), ties half: 0.5 = no contrast. [ ] = 95% interval, resampling genes. Compare within a panel only: ",
         "the channels sit on different scales.")),
       caption = wrap(175,
         paste0(sprintf("The old panel keeps its %s scored residues with no Aminode value (0, as the model saw them); ",
                        format(es$n_old0, big.mark = ",")),
                sprintf("dropped from BOTH channels, AUC old/new = %.2f/%.2f (trained), %.2f/%.2f (never-seen peptide), %.2f/%.2f (large product). ",
                        a_of("trained", "old", "auc_cov"), a_of("trained", "new", "auc_cov"),
                        a_of("never-seen peptide", "old", "auc_cov"), a_of("never-seen peptide", "new", "auc_cov"),
                        a_of("never-seen large product", "old", "auc_cov"), a_of("never-seen large product", "new", "auc_cov")),
                sprintf("Against only the other residues of the SAME precursors: %.2f/%.2f, %.2f/%.2f, %.2f/%.2f.",
                        a_of("trained", "old", "auc_same"), a_of("trained", "new", "auc_same"),
                        a_of("never-seen peptide", "old", "auc_same"), a_of("never-seen peptide", "new", "auc_same"),
                        a_of("never-seen large product", "old", "auc_same"), a_of("never-seen large product", "new", "auc_same"))),
         sprintf(paste0("Over all scored residues the channels agree weakly (Spearman %.2f), and %.0f%% of the old ",
                        "channel's variance is a per-protein level (%.0f%% of the new one's; the old one's ceiling is ",
                        "how deep in time the precursor's Aminode alignment reaches): the old channel mostly says ",
                        "WHICH PROTEIN, the new one which residue."),
                 es$spearman, 100 * es$between[["old"]], 100 * es$between[["new"]]),
         model_txt,
         paste0("Ends at the windows' anchors (G + dibasic C ends on the G). Large cleavage products = the never-seen set's ",
                "longer chains: chemokines, interleukins, TNF family, WNTs, BMPs, growth factors ..."),
         sprintf("Built %s by inst/scripts/msa_cons/plot_msa_conservation.R; numbers in msacons_end_summary.csv.",
                 format(Sys.Date())))) +
  theme(plot.title = element_text(size = 10.5),
        plot.subtitle = element_text(size = 8.5),
        plot.caption = element_text(size = 7.5, lineheight = 1.15))
f2 <- file.path(out_dir, "msacons_end_summary")
ggsave(paste0(f2, ".png"), p2, width = 10, height = 5.6, dpi = dpi, device = ragg::agg_png, bg = LF_VIZ$surface)
ggsave(paste0(f2, ".pdf"), p2, width = 10, height = 5.6, device = cairo_pdf, bg = LF_VIZ$surface)
message("wrote ", f2, ".{png,pdf,csv}")
