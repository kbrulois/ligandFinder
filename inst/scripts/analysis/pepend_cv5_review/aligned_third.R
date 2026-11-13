#!/usr/bin/env Rscript
## Reproduces ~/AF2_analysis/pepend_windows_aligned_20260925_160838.svg at 1/3
## the rows. That figure groups rows by MOTIF (the current 10_7b groups by
## insertion class instead, a later change), and colours the row label by the
## train/val split rather than drawing a dot beside it -- both restored here.
##
## Rows are dropped at random WITHIN each term x motif cell, so every category
## survives and the group proportions are preserved.
##
## The split is shown as a coloured DOT, not a coloured row label. ggtext 0.1.2
## emits old S3 theme elements while ggplot2 4.x expects S7 ones, so
## element_markdown is silently ignored and the <span> markup prints literally
## -- which is exactly what happened in the September figure being reproduced.
suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })
.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i+1L]] }
in_path <- path.expand(.opt("--in", "~/AF2_analysis/lf_pepend_nn_input.rds"))
out_f   <- path.expand(.opt("--out", "~/AF2_analysis/pepend_windows_aligned_third.svg"))
FRAC    <- as.numeric(.opt("--frac", "0.3333"))
SEED    <- as.integer(.opt("--seed", "42"))
ROOT <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
TOL <- 2L

.x <- readRDS(in_path); k <- .x$knowns; ct_pad <- isTRUE(.x$params$ct_pad); rm(.x)
k <- bind_cols(k, lf_pepend_insertion_class(k$known_idx, k$term, k$len, tol = TOL, detail = TRUE))

class_cols <- nn_class_cols[c("NT_cleavage_context","pep_other","pep_pocket",
                              "CT_cleavage_context","padding")]
class_labs <- nn_class_labels[names(class_cols)]
motif_lv  <- c("dibasic","G + dibasic","monobasic","terminus","other")
ins_short <- c("end insertion"="end ins","loop insertion"="loop ins","non-inserting end"="non-ins")
split_cols <- c(train = "grey60", val = "#2166ac")

n_before <- nrow(k)
set.seed(SEED)
k <- k %>%
  mutate(term = factor(term, levels = c("N","C")), motif = factor(motif, levels = motif_lv)) %>%
  group_by(term, motif) %>%
  ## group_modify, not slice_sample(n = n() * FRAC): slice_sample needs a
  ## CONSTANT n, so the per-group count has to be computed inside the group.
  ## max(1) keeps every category on the figure even where the cell is tiny.
  group_modify(~ slice_sample(.x, n = max(1L, round(nrow(.x) * FRAC)))) %>%
  ungroup()
message(sprintf("kept %d of %d rows (%.0f%%)", nrow(k), n_before, 100*nrow(k)/n_before))

k <- k %>%
  mutate(term = factor(as.character(term), levels = c("N","C"),
                       labels = c("N terminus: first residue at 8","C terminus: last residue at 28")),
         name = coalesce(pep_name, sub("^[^_]+_", "", pep_id)),
         row  = sprintf("%s %s %d-%d (%d) | %s | %s", gene, name,
                        pep_start, pep_end, len, motif, ins_short[insertion])) %>%
  arrange(term, motif, len, gene) %>%
  group_by(term) %>% mutate(y = rev(row_number())) %>% ungroup()

tiles <- k %>% select(term, y, row, split, AA, known_idx) %>%
  unnest(c(AA, known_idx)) %>%
  group_by(term, y) %>% mutate(pos = row_number()) %>% ungroup() %>%
  mutate(known_idx = factor(known_idx, levels = names(class_cols)))
grp <- k %>% group_by(term, motif) %>% summarise(lo = min(y), hi = max(y), n = n(), .groups = "drop")
cut <- tibble(term = factor(levels(k$term), levels = levels(k$term)), x = c(7.5, 28.5))

mk <- function(t) {
  d <- filter(tiles, term == t); g <- filter(grp, term == t); kk <- filter(k, term == t)
  ggplot(d, aes(pos, y)) +
    geom_tile(aes(fill = known_idx), colour = "white", linewidth = 0.15, show.legend = TRUE) +
    geom_text(aes(label = AA), size = 1.7, na.rm = TRUE) +
    geom_point(data = kk, aes(x = -0.2, y = y, colour = split), size = 1.1) +
    geom_vline(data = filter(cut, term == t), aes(xintercept = x), linewidth = 0.6) +
    geom_hline(yintercept = head(sort(g$hi), -1) + 0.5, linewidth = 0.4, colour = "grey30") +
    annotate("text", x = 37, y = (g$lo + g$hi)/2,
             label = sprintf("%s\n(%d)", g$motif, g$n), hjust = 0, size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = class_cols, labels = class_labs, name = "label", drop = FALSE) +
    scale_colour_manual(values = split_cols, name = "split") +
    guides(colour = guide_legend(override.aes = list(fill = NA))) +
    scale_x_continuous(breaks = c(1,8,15,22,28,36), expand = expansion(add = c(0.3, 6))) +
    scale_y_continuous(breaks = kk$y, labels = kk$row, expand = expansion(add = 0.6)) +
    coord_cartesian(clip = "off") +
    labs(title = t, x = "position in window", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 5.2),
          plot.title = element_text(face = "bold"), legend.position = "bottom")
}
p <- patchwork::wrap_plots(lapply(levels(k$term), mk), nrow = 1, guides = "collect") +
  patchwork::plot_annotation(
    title = sprintf("Peptide-end training windows, 1 in 3: %d N and %d C of %d, from %d docked peptides (chemokines excluded)",
                    sum(k$term == levels(k$term)[1]), sum(k$term == levels(k$term)[2]), n_before,
                    length(unique(unlist(strsplit(k$pep_ids, ";"))))),
    subtitle = paste0("A random third of the rows, sampled WITHIN each motif group so every category and its proportions survive. ",
                      "The dot left of position 1 is the split: grey = train, blue = val.\n",
                      "Labels: before the peptide = NT context, after it = CT context; pocket = in the receptor pocket in the docked model. ",
                      "Groups: what flanks this end.")) &
  theme(legend.position = "bottom")
n_rows <- max(table(k$term))
ggsave(out_f, p, width = 20, height = 2 + 0.115 * n_rows, limitsize = FALSE)
message("wrote ", out_f, sprintf("  (%d rows tall)", n_rows))
print(as.data.frame(grp %>% select(term, motif, n)), row.names = FALSE)
