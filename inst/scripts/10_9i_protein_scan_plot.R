#!/usr/bin/env Rscript
## ---- one protein, every step-1 scan window, in sequence context ----------------
## The repo's protein page (make_protein_plot_win) draws each candidate window as
## a segment on a stacked layer. That works when a gene has a handful of windows
## -- NPY has TWO in the production set -- but the step-1 scan gives one window
## per mature residue, and on a 97-residue precursor all 69 of them overlap, so
## assign_overlap_layers would stack 69 layers and the page would be unreadable.
##
## A scan is a continuous track, so it is drawn as one: score against the residue
## the window is ANCHORED on (its C-terminal residue), with the across-seed band,
## over the sequence and its cleavage motifs.
##
## Each x position is a whole 36-residue window, not a residue property: the
## window anchored at r covers r-27..r+8. The score belongs to the cleavage
## hypothesis "the peptide ends at r", which is why the known peptide's C
## terminus is the thing to look for, not its span.
##
##   Rscript inst/scripts/10_9i_protein_scan_plot.R --gene NPY
##
## Output: ~/AF2_analysis/lf_scan_protein_<gene>_<term>.{svg,png}
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
gene   <- .opt("--gene", "NPY")
term   <- toupper(.opt("--term", "C"))
scan_p <- path.expand(.opt("--scan", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
cache_p<- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
in_p   <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
ll_p   <- .opt("--ligand-list", "inst/extdata/ligand_list.rds")
out_st <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_scan_protein_%s_%s", gene, term)))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "viz_theme.R"))

ANCHOR <- c(N = 8L, C = 28L)[[term]]     # where the peptide end sits in a window

## ---- 1. the protein, its scan, its annotations ---------------------------------
fc <- readRDS(cache_p)
p  <- fc$prec %>% filter(gene == !!gene)
if (!nrow(p)) stop("no precursor called ", gene, call. = FALSE)
if (nrow(p) > 1) { message("note: ", nrow(p), " accessions for ", gene, "; using ", p$accession[1]); p <- p[1, ] }
L <- nchar(p$seq); aa <- strsplit(p$seq, "")[[1]]
message(sprintf("%s (%s): %d residues, mature %d-%d", gene, p$accession, L, p$n_prot, p$c_prot))

sc <- readRDS(scan_p)
tr <- sc$scan %>% filter(accession == p$accession) %>% arrange(anchor) %>%
  mutate(motif = lf_pepend_motif_vec(p$seq, term, anchor, p$n_prot),
         motif3 = case_when(motif %in% c("dibasic", "G + dibasic") ~ "dibasic",
                            motif == "monobasic" ~ "monobasic",
                            TRUE ~ "other / terminus"))
pk <- sc$peaks %>% filter(accession == p$accession) %>% arrange(anchor)
message(sprintf("scanned %d positions; %d local peaks; score %.3f-%.3f",
                nrow(tr), nrow(pk), min(tr$score), max(tr$score)))

kn <- sc$knowns %>% filter(accession == p$accession, !is.na(rank_all))
peps <- tryCatch(readRDS(file.path(ROOT, ll_p)) %>% filter(accession == p$accession) %>%
                   select(final_name, start, end) %>% distinct() %>% arrange(start),
                 error = function(e) tibble(final_name = character(), start = integer(), end = integer()))
## isoforms of one peptide share a C terminus and differ by a residue or two at
## the N side, so drawn on one row they sit on top of each other -- give each
## its own row
peps <- peps %>% mutate(row = dplyr::row_number(),
                        y = -0.085 - 0.075 * (row - 1),
                        y_lab = y - 0.048)
Y_LO <- if (nrow(peps)) min(peps$y_lab) - 0.045 else -0.06
if (nrow(kn)) { message("known ", term, " ends:"); print(as.data.frame(kn %>% select(pep_name, anchor, score, rank_all)), row.names = FALSE) }

## ---- 2. the track ---------------------------------------------------------------
MOT <- c(dibasic = LF_VIZ$s2, monobasic = LF_VIZ$s3, `other / terminus` = LF_VIZ$grid)
sig <- tibble(xmin = 0.5, xmax = p$n_prot - 0.5)          # signal peptide: not scanned

top <- tr %>% slice_max(score, n = 1)
lab_kn <- kn %>% left_join(tr %>% select(anchor, score), by = "anchor", suffix = c("", ".tr"))

p_track <- ggplot(tr, aes(anchor, score)) +
  {if (p$n_prot > 1) annotate("rect", xmin = sig$xmin, xmax = sig$xmax, ymin = -Inf, ymax = Inf,
                              fill = LF_VIZ$grid, alpha = 0.55) } +
  {if (p$n_prot > 1) annotate("text", x = (sig$xmin + sig$xmax) / 2, y = 1.04, size = 2.8,
                              colour = LF_VIZ$ink2, label = "signal peptide\n(not scanned)",
                              lineheight = 0.95) } +
  ## annotated peptides: the span, with the C terminus that the model is scoring
  {if (nrow(peps)) geom_segment(data = peps, inherit.aes = FALSE,
      aes(x = start, xend = end, y = y, yend = y),
      colour = LF_VIZ$ink2, linewidth = 1.6, lineend = "round") } +
  {if (nrow(peps)) geom_text(data = peps, inherit.aes = FALSE,
      aes(x = (start + end) / 2, y = y_lab, label = final_name),
      colour = LF_VIZ$ink2, size = 2.7) } +
  geom_ribbon(aes(ymin = pmax(score - sd, 0), ymax = pmin(score + sd, 1)),
              fill = LF_VIZ$s1, alpha = 0.2) +
  geom_line(colour = LF_VIZ$s1, linewidth = 0.8) +
  geom_point(data = pk, colour = LF_VIZ$s1, fill = LF_VIZ$surface,
             shape = 21, size = 2.2, stroke = 0.8) +
  ## the known end: a rule at the anchor, labelled
  {if (nrow(lab_kn)) geom_segment(data = lab_kn, inherit.aes = FALSE,
      aes(x = anchor, xend = anchor, y = 0, yend = score),
      colour = LF_VIZ$ink, linewidth = 0.4, linetype = "22") } +
  {if (nrow(lab_kn)) geom_point(data = lab_kn, inherit.aes = FALSE,
      aes(anchor, score), colour = LF_VIZ$ink, size = 2.1) } +
  {if (nrow(lab_kn)) geom_text(data = lab_kn, inherit.aes = FALSE,
      aes(anchor, score, label = sprintf("%s\nends here (%d)", pep_name, anchor)),
      colour = LF_VIZ$ink, size = 2.9, hjust = 1.08, vjust = 0.15, lineheight = 1.1) } +
  scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                     breaks = scales::breaks_width(10)) +
  scale_y_continuous(limits = c(Y_LO, 1.09), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = sprintf("%s (%s): every step-1 %s-terminus window", gene, p$accession, term),
       subtitle = sprintf(paste("One window per mature residue -- %d of them, each a 36-residue window",
                                "anchored on that residue.\nLine is the 20-seed mean, band +/- 1 sd.",
                                "Circles mark local maxima (+/-5)."), nrow(tr)),
       x = NULL, y = "window score") +
  lf_viz_theme() +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())

## ---- 3. sequence + motif strip ---------------------------------------------------
seqd <- tibble(pos = seq_len(L), aa = aa) %>%
  left_join(tr %>% select(pos = anchor, motif3), by = "pos") %>%
  mutate(motif3 = ifelse(is.na(motif3), "other / terminus", motif3))

p_seq <- ggplot(seqd, aes(pos, 1)) +
  geom_tile(aes(fill = motif3), height = 0.42, colour = NA) +
  geom_text(aes(label = aa), size = 1.85, colour = LF_VIZ$ink, vjust = 0.5) +
  scale_fill_manual(values = MOT, name = "what follows the anchor",
                    breaks = names(MOT)) +
  scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                     breaks = scales::breaks_width(10)) +
  labs(x = "residue (the window's anchor = the peptide's last residue)", y = NULL) +
  lf_viz_theme() +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.justification = "left",
        legend.key.size = unit(9, "pt"),
        legend.text = element_text(colour = LF_VIZ$ink2, size = 8),
        legend.title = element_text(colour = LF_VIZ$ink2, size = 8))

fig <- p_track / p_seq + plot_layout(heights = c(5, 1.5)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = LF_VIZ$surface, colour = NA)))

ggsave(paste0(out_st, ".svg"), fig, width = 11, height = 5.9, device = svglite::svglite)
ggsave(paste0(out_st, ".png"), fig, width = 11, height = 5.9, dpi = 200, bg = LF_VIZ$surface)
message("\nplot -> ", out_st, ".{svg,png}")

message("\ntop 8 positions:")
print(as.data.frame(tr %>% arrange(desc(score)) %>% head(8) %>%
        mutate(aa = aa[anchor], next3 = substring(p$seq, anchor + 1, anchor + 3)) %>%
        select(anchor, aa, next3, motif, score, sd, rank_all)), row.names = FALSE, digits = 3)
