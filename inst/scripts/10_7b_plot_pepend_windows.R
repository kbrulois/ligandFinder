#!/usr/bin/env Rscript
## ---- every known peptide-end window, aligned, with its labels ------------------
## One row per window from 10_7_pepend_windows.R, N and C side by side: each
## residue tile is coloured by its training label and carries its letter. The
## line marks the cut -- between positions 7 and 8 (N: the peptide starts at 8)
## and 28 and 29 (C: it ends at 28).
##
## Rows are grouped by what flanks the terminus (dibasic, G + dibasic,
## monobasic, precursor terminus, other), then by peptide length. The row label
## gives the peptide (ligand-list name where one matches exactly), its range and
## how this end sits in the docked receptor model; the dot left of position 1
## is the train/val split.
##
##   /usr/local/bin/Rscript inst/scripts/10_7b_plot_pepend_windows.R
##
## Output: ~/AF2_analysis/pepend_windows_aligned_<stamp>.svg
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
in_path <- path.expand(.opt("--in", "~/AF2_analysis/lf_pepend_nn_input.rds"))
out_dir <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
stamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
## nn_class_cols / nn_class_labels: the shared class palette and figure names
source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
TOL <- as.integer(.opt("--tol", "2"))

k <- readRDS(in_path)$knowns
## where the pocket sits relative to THIS end -- the classification under review
k <- dplyr::bind_cols(k, lf_pepend_insertion_class(k$known_idx, k$term, k$len,
                                                   tol = TOL, detail = TRUE))

## the five classes this figure draws, taking colour AND label from the shared
## table rather than restating either -- a second copy of the palette here is
## how a class ends up a different colour on two pages of the same report
class_cols <- nn_class_cols[c("NT_cleavage_context", "pep_other", "pep_pocket",
                              "CT_cleavage_context", "padding")]
class_labs <- nn_class_labels[names(class_cols)]
motif_lv <- c("dibasic", "G + dibasic", "monobasic", "terminus", "other")
cls_lv   <- c("inserting", "loop", "non_inserting")
cls_lab  <- c(inserting = "INSERTING\npocket at\nthis end",
              loop = "LOOP\npocket mid-peptide,\nnear neither end",
              non_inserting = "NON-INSERTING\npocket at the far end,\nor not in the window")
ins_short <- c("end insertion" = "end ins", "loop insertion" = "loop ins",
               "non-inserting end" = "non-ins")

k <- k %>%
  mutate(term  = factor(term, levels = c("N", "C"),
                        labels = c("N terminus: first residue at 8", "C terminus: last residue at 28")),
         motif = factor(motif, levels = motif_lv),
         ins_class = factor(ins_class, levels = cls_lv),
         name  = coalesce(pep_name, sub("^[^_]+_", "", pep_id)),
         ## the two distances the call turned on, so the rule can be audited
         ## from the figure rather than taken on trust
         row   = sprintf("%s %s %d-%d (%d) | %s | end %s far %s | %s%s",
                         gene, name, pep_start, pep_end, len, motif,
                         ifelse(is.na(gap_end), "-", gap_end),
                         ifelse(is.na(gap_far), "-", gap_far),
                         ins_short[insertion],
                         ifelse(grepl(";", pep_ids), " | shared", ""))) %>%
  arrange(term, ins_class, len, gene) %>%
  group_by(term) %>% mutate(y = rev(row_number())) %>% ungroup()

tiles <- k %>%
  select(term, y, row, split, AA, known_idx) %>%
  unnest(c(AA, known_idx)) %>%
  group_by(term, y) %>% mutate(pos = row_number()) %>% ungroup() %>%
  mutate(known_idx = factor(known_idx, levels = names(class_cols)))

## motif group separators and labels
grp <- k %>% group_by(term, ins_class) %>% summarise(lo = min(y), hi = max(y), n = n(), .groups = "drop")
cut <- tibble(term = factor(levels(k$term), levels = levels(k$term)), x = c(7.5, 28.5))

split_cols <- c(train = "grey60", val = "#2166ac")

mk <- function(t) {
  d <- tiles %>% filter(term == t)
  g <- grp %>% filter(term == t)
  ggplot(d, aes(pos, y)) +
    geom_tile(aes(fill = known_idx), colour = "white", linewidth = 0.15) +
    geom_text(aes(label = AA), size = 1.7, na.rm = TRUE) +
    geom_point(data = k %>% filter(term == t), aes(x = -0.2, y = y, colour = split), size = 1.1) +
    geom_vline(data = cut %>% filter(term == t), aes(xintercept = x), linewidth = 0.6) +
    geom_hline(yintercept = head(sort(g$hi), -1) + 0.5, linewidth = 0.4, colour = "grey30") +
    annotate("text", x = 37, y = (g$lo + g$hi) / 2, label = sprintf("%s\n(%d)", cls_lab[as.character(g$ins_class)], g$n),
             hjust = 0, size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = class_cols, labels = class_labs, name = "label",
                      drop = FALSE) +
    scale_colour_manual(values = split_cols, name = "split") +
    scale_x_continuous(breaks = c(1, 8, 15, 22, 28, 36), expand = expansion(add = c(0.3, 6))) +
    scale_y_continuous(breaks = k$y[k$term == t], labels = k$row[k$term == t],
                       expand = expansion(add = 0.6)) +
    coord_cartesian(clip = "off") +
    labs(title = t, x = "position in window", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 5.2),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

p <- patchwork::wrap_plots(lapply(levels(k$term), mk), nrow = 1, guides = "collect") +
  patchwork::plot_annotation(
    title = sprintf("Peptide-end training windows: %d N and %d C, from %d docked peptides (chemokines excluded)",
                    sum(k$term == levels(k$term)[1]), sum(k$term == levels(k$term)[2]),
                    length(unique(unlist(strsplit(k$pep_ids, ";"))))),
    subtitle = paste0("Labels: before the peptide = NT context, after it = CT context (DB and gap merged in); ",
                      "pocket = in the receptor pocket in the docked model. Groups: what flanks this end.\n",
                      "Row: gene, peptide, range (length) | how this end sits in the docked model ",
                      "(end/loop insertion or non-inserting) | shared = end shared by several peptides, ",
                      "labelled from the shortest.\n",
                      sprintf(paste0("GROUPS are the classification under review, from the pocket labels in THIS window ",
                                     "(tolerance %d residues): INSERTING = a pocket residue within %d of this end; ",
                                     "LOOP = pocket present but near neither end;\nNON-INSERTING = pocket within %d of ",
                                     "the peptide's other terminus, or no pocket in the window. Row shows end/far = the ",
                                     "two distances the call turned on. `inserting` wins when the pocket reaches both ends ",
                                     "(a 4-residue peptide inserts at each)."), TOL, TOL, TOL))) &
  theme(legend.position = "bottom")

n_rows <- max(table(k$term))
f <- file.path(out_dir, sprintf("pepend_windows_aligned_%s.svg", stamp))
ggsave(f, p, width = 20, height = 2 + 0.115 * n_rows, limitsize = FALSE)
message("wrote ", f)
