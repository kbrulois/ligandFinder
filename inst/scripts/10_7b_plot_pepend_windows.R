#!/usr/bin/env Rscript
## ---- every known peptide-end window, aligned, with its labels ------------------
## One row per window from 10_7_pepend_windows.R, N and C side by side: each
## residue tile is coloured by its training label and carries its letter. The
## line marks the cut -- between positions 7 and 8 (N: the peptide starts at 8)
## and 28 and 29 (C: it ends at 28).
##
## Rows are grouped by where the docked pocket sits relative to this end
## (inserting, loop, non-inserting), then by peptide length. The row label
## gives the peptide (ligand-list name where one matches exactly), its range,
## what flanks the end and how it sits in the docked receptor model; the dot
## left of position 1 is the train/val split.
##
## Kevin, 2026-11-12: the _padg builds need two more things read off the rows.
##   * A cross-validated build (--n-folds) records only its own fold's split,
##     so its dot would say "held out in fold 1 or not". The digit there says
##     WHICH fold holds the end out, read from the sibling fold inputs (~5 s
##     each); without them it falls back to the dot.
##   * A --amid-g-anchor build ends a G + dibasic peptide ON its G. The moved G
##     is outlined (the C anchor, and inside the same peptide's N window where
##     it falls there), C rows whose end moved carry a diamond, and their range
##     reads start-<mature end>+G.
## A build with neither (no --n-folds, no --amid-g-anchor) draws as it always
## did; an old cv5 build gets the fold digits and no G marks.
##
##   /usr/local/bin/Rscript inst/scripts/10_7b_plot_pepend_windows.R
##   ... --in ~/AF2_analysis/lf_pepend_nn_input_cv5_f1_padg.rds --out-dir ~/AF2_analysis/lf_viz_20261112
##
## Output: ~/AF2_analysis/pepend_windows_aligned_[ctpad_][amidg_]<stamp>.{svg,png}
## 10_7c_plot_gshift_windows.R puts the moved C ends beside their old windows.
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
source(file.path(ROOT, "R", "viz_theme.R"))
## lf_pepend_cv_folds / lf_pepend_gshift_pos / lf_gshift_cols
source(file.path(ROOT, "R", "pepend_gshift_plots.R"))
TOL <- as.integer(.opt("--tol", "2"))

.x <- readRDS(in_path)
k <- .x$knowns
params <- .x$params
## a --ct-pad set (10_7) cuts every C window at its anchor; say so on the figure
ct_pad <- isTRUE(params$ct_pad)
## a --amid-g-anchor set ends G + dibasic peptides ON the G; mark what moved
amid_g <- isTRUE(params$amid_g)
rm(.x); invisible(gc())
## where the pocket sits relative to THIS end -- the classification under review
k <- dplyr::bind_cols(k, lf_pepend_insertion_class(k$known_idx, k$term, k$len,
                                                   tol = TOL, detail = TRUE))
## which fold holds each end out (NULL: not cross-validated, or a sibling fold
## input missing -- then the dot shows this build's own split, as before)
fold <- lf_pepend_cv_folds(in_path, k, params)
k$fold <- if (is.null(fold)) NA_integer_ else fold
## the moved G's window position, NA on every other row (and every row of an
## old build); g_amid is per PEPTIDE, so it marks both of its windows
k$g_pos <- lf_pepend_gshift_pos(k, amid_g)
k$moved <- amid_g & k$term == "C" & (if ("g_amid" %in% names(k)) k$g_amid else FALSE)
if (amid_g) message(sprintf("--amid-g-anchor build: %d C ends moved onto the G; G outlined in %d N windows",
                            sum(k$moved), sum(!is.na(k$g_pos) & k$term == "N")))
## Kevin, 2026-11-12 (review): the groups below count the pocket's distance
## from each end, and a moved peptide's C end is now its G. The docked models
## are of the mature peptides, so the G is outside every pocket and that
## distance grows by one -- enough to move 5 rows of the _padg build to another
## group than on the old plot (10_7c names them). Said on the figure only if
## true of the Gs drawn.
g_out <- amid_g && !any(vapply(which(!is.na(k$g_pos)), function(i)
  k$known_idx[[i]][[k$g_pos[[i]]]] == "pep_pocket", NA))

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

## a G-anchored peptide's range reads start-<mature end>+G: the window counts
## the G as peptide, the mature (amidated) peptide stops before it
rng <- if (amid_g) ifelse(k$g_amid, sprintf("%d-%d+G", k$pep_start, k$pep_end_mature),
                          sprintf("%d-%d", k$pep_start, k$pep_end)) else
  sprintf("%d-%d", k$pep_start, k$pep_end)
k <- k %>%
  mutate(term  = factor(term, levels = c("N", "C"),
                        labels = c("N terminus: first residue at 8",
                                   if (amid_g) "C terminus: last residue at 28 (on G + dibasic ends, the G)"
                                   else "C terminus: last residue at 28")),
         motif = factor(motif, levels = motif_lv),
         ins_class = factor(ins_class, levels = cls_lv),
         name  = coalesce(pep_name, sub("^[^_]+_", "", pep_id)),
         ## the two distances the call turned on, so the rule can be audited
         ## from the figure rather than taken on trust
         row   = sprintf("%s %s %s (%d) | %s | end %s far %s | %s%s",
                         gene, name, rng, len, motif,
                         ifelse(is.na(gap_end), "-", gap_end),
                         ifelse(is.na(gap_far), "-", gap_far),
                         ins_short[insertion],
                         ifelse(grepl(";", pep_ids), " | shared", ""))) %>%
  arrange(term, ins_class, len, gene) %>%
  group_by(term) %>% mutate(y = rev(row_number())) %>% ungroup()

tiles <- k %>%
  select(term, y, row, split, g_pos, AA, known_idx) %>%
  unnest(c(AA, known_idx)) %>%
  group_by(term, y) %>% mutate(pos = row_number()) %>% ungroup() %>%
  mutate(known_idx = factor(known_idx, levels = names(class_cols)))

## motif group separators and labels
grp <- k %>% group_by(term, ins_class) %>% summarise(lo = min(y), hi = max(y), n = n(), .groups = "drop")
cut <- tibble(term = factor(levels(k$term), levels = levels(k$term)), x = c(7.5, 28.5))

split_cols <- c(train = "grey60", val = "#2166ac")
cv <- !all(is.na(k$fold))
g_col <- lf_gshift_cols()[["G"]]
## room left of position 1 for the fold digit and the moved-end diamond
x_lft <- if (amid_g) 2.1 else if (cv) 0.9 else 0.3

mk <- function(t) {
  d <- tiles %>% filter(term == t)
  g <- grp %>% filter(term == t)
  kt <- k %>% filter(term == t)
  ggplot(d, aes(pos, y)) +
    ## show.legend = TRUE draws a key for a class the panel lacks (CT context in
    ## a --ct-pad C panel), so both panels' legends match and collect into one
    geom_tile(aes(fill = known_idx), colour = "white", linewidth = 0.15, show.legend = TRUE) +
    geom_text(aes(label = AA), size = 1.7, na.rm = TRUE) +
    (if (cv) geom_text(data = kt, aes(x = -0.2, y = y, label = fold), size = 1.7, colour = LF_VIZ$ink2)
     else geom_point(data = kt, aes(x = -0.2, y = y, colour = split), size = 1.1)) +
    ## the moved G, wherever it falls in a window, and a diamond on the C rows
    ## whose end moved onto it
    (if (amid_g) geom_tile(data = d %>% filter(pos == g_pos), fill = NA, colour = g_col,
                           linewidth = 0.55, width = 0.9, height = 0.86)) +
    (if (any(kt$moved)) geom_point(data = kt %>% filter(moved), aes(x = -1.3, y = y, shape = "moved"),
                                   colour = g_col, size = 1.5)) +
    geom_vline(data = cut %>% filter(term == t), aes(xintercept = x), linewidth = 0.6) +
    geom_hline(yintercept = head(sort(g$hi), -1) + 0.5, linewidth = 0.4, colour = "grey30") +
    annotate("text", x = 37, y = (g$lo + g$hi) / 2, label = sprintf("%s\n(%d)", cls_lab[as.character(g$ins_class)], g$n),
             hjust = 0, size = 2.6, lineheight = 0.9) +
    scale_fill_manual(values = class_cols, labels = class_labs, name = "label",
                      drop = FALSE) +
    ## ...but that also paints a tile behind the split dots and the diamond;
    ## clear it there. (Each scale only where its layer is: an unused manual
    ## scale warns that it shares no levels with the data.)
    (if (!cv) list(scale_colour_manual(values = split_cols, name = "split"),
                   guides(colour = guide_legend(override.aes = list(fill = NA))))) +
    (if (any(kt$moved)) list(scale_shape_manual(values = c(moved = 18), name = NULL,
                                         labels = c(moved = "C end moved onto the G (G outlined)")),
                      guides(shape = guide_legend(override.aes = list(fill = NA, colour = g_col,
                                                                      size = 2.5))))) +
    scale_x_continuous(breaks = c(1, 8, 15, 22, 28, 36), expand = expansion(add = c(x_lft, 6))) +
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
    title = sprintf("Peptide-end training windows: %d N and %d C, from %d docked peptides (chemokines excluded)%s",
                    sum(k$term == levels(k$term)[1]), sum(k$term == levels(k$term)[2]),
                    length(unique(unlist(strsplit(k$pep_ids, ";")))),
                    if (amid_g) " -- G + dibasic C ends anchored ON the G" else ""),
    subtitle = paste0(if (ct_pad) "CT-PADDED: every C window is cut at its anchor, so positions 29-36 are padding whatever follows the peptide.\n" else "",
                      if (amid_g) sprintf(paste0(
                        "G-ANCHORED (--amid-g-anchor): the %d C ends followed by G + KK/KR/RK/RR end ON the G (diamond), ",
                        "one past the mature amidated end, so the range reads start-<mature end>+G. The G is outlined, ",
                        "and is peptide in the same peptide's N window too (outlined where it falls inside).\n"),
                        sum(k$moved)) else "",
                      if (g_out) paste0(
                        "No docked model has the G in its pocket, so on these peptides' rows the pocket's distance from ",
                        "the C end (end on C rows, far on N rows) is one more than from the mature end, and can put a ",
                        "row in another group than under the old anchor.\n") else "",
                      "Labels: before the peptide = NT context, after it = CT context (DB and gap merged in); ",
                      "pocket = in the receptor pocket in the docked model.",
                      if (cv) sprintf(" Digit left of position 1: the cross-validation fold (of %d) that holds this end out.",
                                      params$n_folds) else "",
                      "\nRow: gene, peptide, range (length) | what flanks this end | end/far (below) | ",
                      "how this end sits in the docked model (end/loop insertion or non-inserting) | ",
                      "shared = end shared by several peptides, labelled from the shortest.\n",
                      sprintf(paste0("GROUPS are the classification under review, from the pocket labels in THIS window ",
                                     "(tolerance %d residues): INSERTING = a pocket residue within %d of this end; ",
                                     "LOOP = pocket present but near neither end;\nNON-INSERTING = pocket within %d of ",
                                     "the peptide's other terminus, or no pocket in the window. Row shows end/far = the ",
                                     "two distances the call turned on. `inserting` wins when the pocket reaches both ends ",
                                     "(a 4-residue peptide inserts at each)."), TOL, TOL, TOL))) &
  theme(legend.position = "bottom")

n_rows <- max(table(k$term))
## (+ the subtitle lines a G-anchored build adds, so its rows keep their height)
h <- 2 + 0.115 * n_rows + 0.2 * (amid_g + g_out)
f <- file.path(out_dir, sprintf("pepend_windows_aligned_%s%s%s.svg", if (ct_pad) "ctpad_" else "",
                                if (amid_g) "amidg_" else "", stamp))
ggsave(f, p, width = 20, height = h, limitsize = FALSE)
## and a PNG to look at; 200 dpi keeps the 1.7 mm residue letters readable
ggsave(sub("\\.svg$", ".png", f), p, width = 20, height = h, dpi = 200,
       bg = "white", limitsize = FALSE)
message("wrote ", f, " (+ .png)")
