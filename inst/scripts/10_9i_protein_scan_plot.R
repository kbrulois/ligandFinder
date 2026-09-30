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
##   ... --genes NPY,KNG1,CXCL14 --group train --no-png
##
## Several genes in one call share the loaded caches -- the residue cache and the
## scan are ~100 MB and 2.98M rows, so loading them per gene is the whole cost.
##
## Output: ~/AF2_analysis/lf_scan_protein_<gene>_<term>.svg (+ .png unless
## --no-png), and a metadata rds naming what each figure holds.
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(ggplot2); library(patchwork) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
genes  <- strsplit(.opt("--genes", .opt("--gene", "NPY")), ",")[[1]] |> trimws()
group  <- .opt("--group", NA_character_)
term   <- toupper(.opt("--term", "C"))
scan_p <- path.expand(.opt("--scan", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
cache_p<- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
ll_p   <- .opt("--ligand-list", "inst/extdata/ligand_list.rds")
out_dir<- path.expand(.opt("--out-dir", "~/AF2_analysis"))
meta_p <- path.expand(.opt("--meta", "~/AF2_analysis/lf_scan_protein_meta.rds"))
no_png <- "--no-png" %in% .args

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "viz_theme.R"))

ANCHOR <- c(N = 8L, C = 28L)[[term]]     # where the peptide end sits in a window
RES <- c(`dibasic pair (KK/KR/RK/RR)` = LF_VIZ$s2,
         `single K / R`               = LF_VIZ$s3,
         other                        = LF_VIZ$grid)

message("loading caches (once for all ", length(genes), " gene(s)) ...")
fc <- readRDS(cache_p)
sc <- readRDS(scan_p)
LL <- tryCatch(readRDS(file.path(ROOT, ll_p)), error = function(e) NULL)

## ---- one figure ------------------------------------------------------------------
make_one <- function(gene) {
  p <- fc$prec %>% filter(gene == !!gene)
  if (!nrow(p)) { message("  ", gene, ": not in the secretome -- skipped"); return(NULL) }
  if (nrow(p) > 1) message("  note: ", nrow(p), " accessions for ", gene, "; using ", p$accession[1])
  p <- p[1, ]
  L <- nchar(p$seq); aa <- strsplit(p$seq, "")[[1]]

  tr <- sc$scan %>% filter(accession == p$accession) %>% arrange(anchor) %>%
    ## a property of the WINDOW, not of the residue: what sits just outside the
    ## anchor. It belongs on the anchor axis, never as a colour on the residue.
    mutate(anchor_ctx = lf_pepend_motif_vec(p$seq, term, anchor, p$n_prot))
  if (!nrow(tr)) { message("  ", gene, ": nothing scanned -- skipped"); return(NULL) }
  pk <- sc$peaks %>% filter(accession == p$accession) %>% arrange(anchor)
  kn <- sc$knowns %>% filter(accession == p$accession, !is.na(rank_all))

  peps <- if (is.null(LL)) tibble(final_name = character(), start = integer(), end = integer())
          else LL %>% filter(accession == p$accession) %>%
               select(final_name, start, end) %>% distinct() %>% arrange(start)
  ## isoforms share a C terminus and differ by a residue or two at the N side, and
  ## KNG1 has six of them ending within a few residues -- so each gets its own row
  ## and its label goes to the RIGHT of the span's end, into the axis a C-terminal
  ## peptide leaves empty.
  peps <- peps %>% mutate(row = dplyr::row_number(), y = -0.085 - 0.085 * (row - 1))
  Y_LO <- if (nrow(peps)) min(peps$y) - 0.075 else -0.06

  ## A 97-residue precursor and a 1,232-residue one cannot share a canvas width or
  ## a letter size. Width scales with length (SVG is vector, so wide zooms fine);
  ## the letter is sized so its width fits one residue, and is dropped entirely
  ## when that would take it below legibility.
  PLOT_W  <- max(11, min(42, 3.5 + L * 0.034))
  per_res <- (PLOT_W - 2.2) * 25.4 / L                 # mm of canvas per residue
  SEQ_SZ  <- max(0.85, min(1.9, per_res / 0.62))
  SHOW_AA <- per_res / 0.62 >= 0.8
  BRK_W   <- if (L > 800) 100 else if (L > 400) 50 else if (L > 200) 25 else 10

  lab_kn <- kn %>% left_join(tr %>% select(anchor, score), by = "anchor", suffix = c("", ".tr"))

  p_track <- ggplot(tr, aes(anchor, score)) +
    {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5, ymin = -Inf, ymax = Inf,
                                fill = LF_VIZ$grid, alpha = 0.55) } +
    {if (p$n_prot > 1) annotate("text", x = p$n_prot / 2, y = 1.04, size = 2.8,
                                colour = LF_VIZ$ink2, label = "signal peptide\n(not scanned)",
                                lineheight = 0.95) } +
    {if (nrow(peps)) geom_segment(data = peps, inherit.aes = FALSE,
        aes(x = start, xend = end, y = y, yend = y),
        colour = LF_VIZ$ink2, linewidth = 1.6, lineend = "round") } +
    {if (nrow(peps)) geom_text(data = peps, inherit.aes = FALSE,
        aes(x = end + L * 0.006, y = y, label = final_name),
        colour = LF_VIZ$ink2, size = 2.6, hjust = 0) } +
    geom_ribbon(aes(ymin = pmax(score - sd, 0), ymax = pmin(score + sd, 1)),
                fill = LF_VIZ$s1, alpha = 0.2) +
    geom_line(colour = LF_VIZ$s1, linewidth = 0.8) +
    geom_point(data = pk, colour = LF_VIZ$s1, fill = LF_VIZ$surface,
               shape = 21, size = 2.2, stroke = 0.8) +
    ## the window property, marked on the anchor axis rather than on the residues
    geom_point(data = tr %>% filter(anchor_ctx %in% c("dibasic", "G + dibasic")),
               aes(anchor, -0.03), shape = 25, size = 1.6,
               colour = LF_VIZ$s2, fill = LF_VIZ$s2) +
    {if (nrow(lab_kn)) geom_segment(data = lab_kn, inherit.aes = FALSE,
        aes(x = anchor, xend = anchor, y = 0, yend = score),
        colour = LF_VIZ$ink, linewidth = 0.4, linetype = "22") } +
    {if (nrow(lab_kn)) geom_point(data = lab_kn, inherit.aes = FALSE,
        aes(anchor, score), colour = LF_VIZ$ink, size = 2.1) } +
    {if (nrow(lab_kn)) geom_text(data = lab_kn, inherit.aes = FALSE,
        aes(anchor, score, label = sprintf("%s\nends here (%d)", pep_name, anchor)),
        colour = LF_VIZ$ink, size = 2.9, hjust = 1.08, vjust = 0.15, lineheight = 1.1) } +
    scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                       breaks = scales::breaks_width(BRK_W)) +
    scale_y_continuous(limits = c(Y_LO, 1.09), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(title = sprintf("%s (%s): every step-1 %s-terminus window", gene, p$accession, term),
         subtitle = sprintf(paste("One window per mature residue -- %d of them, each a 36-residue",
                                  "window anchored on that residue.\nLine is the 20-seed mean, band",
                                  "+/- 1 sd. Circles mark local maxima (+/-5); triangles mark",
                                  "anchors followed by a dibasic site."), nrow(tr)),
         x = NULL, y = "window score") +
    lf_viz_theme() +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())

  ## the sequence strip annotates each residue by WHAT IT IS -- an earlier version
  ## coloured each anchor by what FOLLOWED it, which painted 10 of 11 non-basic
  ## residues as "dibasic"/"monobasic"
  is_basic <- aa %in% c("K", "R")
  in_pair  <- rep(FALSE, L)
  pr <- which(head(is_basic, -1) & tail(is_basic, -1))
  in_pair[c(pr, pr + 1L)] <- TRUE
  amid <- pr[pr > 1L & aa[pmax(pr - 1L, 1L)] == "G"] - 1L
  seqd <- tibble(pos = seq_len(L), aa = aa,
                 res = ifelse(in_pair, names(RES)[1],
                       ifelse(is_basic, names(RES)[2], names(RES)[3])))

  p_seq <- ggplot(seqd, aes(pos, 1)) +
    geom_tile(aes(fill = res), height = 0.42, colour = NA) +
    {if (SHOW_AA) geom_text(aes(label = aa), size = SEQ_SZ, colour = LF_VIZ$ink, vjust = 0.5) } +
    {if (length(amid)) annotate("point", x = amid, y = 0.735, shape = 17, size = 1.7,
                                colour = LF_VIZ$ink) } +
    {if (length(amid)) annotate("text", x = amid[1], y = 0.60, label = "amidation G",
                                size = 2.4, colour = LF_VIZ$ink2, hjust = 0.5) } +
    scale_fill_manual(values = RES, name = "residue", breaks = names(RES)) +
    scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                       breaks = scales::breaks_width(BRK_W)) +
    labs(x = "residue (the window's anchor = the peptide's last residue)", y = NULL) +
    lf_viz_theme() +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom", legend.justification = "left",
          legend.key.size = unit(9, "pt"),
          legend.text = element_text(colour = LF_VIZ$ink2, size = 8),
          legend.title = element_text(colour = LF_VIZ$ink2, size = 8))

  fig <- p_track / p_seq + plot_layout(heights = c(5, 1.5)) +
    plot_annotation(theme = theme(plot.background =
                                    element_rect(fill = LF_VIZ$surface, colour = NA)))
  st <- file.path(out_dir, sprintf("lf_scan_protein_%s_%s", gene, term))
  ggsave(paste0(st, ".svg"), fig, width = PLOT_W, height = 5.9, device = svglite::svglite)
  if (!no_png) ggsave(paste0(st, ".png"), fig, width = PLOT_W, height = 5.9,
                      dpi = 200, bg = LF_VIZ$surface)

  best <- tr %>% slice_max(score, n = 1, with_ties = FALSE)
  message(sprintf("  %-8s %4d res, %4d scanned, %2d peaks, best %.3f @ %d, %d known, %.0f in",
                  gene, L, nrow(tr), nrow(pk), best$score, best$anchor, nrow(kn), PLOT_W))
  tibble(gene = gene, accession = p$accession, group = group, term = term,
         n_res = L, n_scanned = nrow(tr), n_peaks = nrow(pk),
         best_score = best$score, best_anchor = best$anchor,
         width_in = PLOT_W, svg = paste0(st, ".svg"),
         n_known = nrow(kn),
         knowns = list(kn %>% select(pep_name, anchor, score, rank_all, split, motif)),
         splits = paste(sort(unique(kn$split)), collapse = "/"))
}

meta <- bind_rows(lapply(genes, make_one))
if (!nrow(meta)) stop("nothing rendered", call. = FALSE)

## accumulate across calls so train/val/other can be built separately
old <- if (file.exists(meta_p)) readRDS(meta_p) else NULL
if (!is.null(old)) meta <- bind_rows(old %>% filter(!gene %in% meta$gene), meta)
saveRDS(meta, meta_p)
message(sprintf("\n%d figure(s) this call; %d in %s", length(genes), nrow(meta), basename(meta_p)))
