## ---- G-anchored peptide-end windows: what moved, and the figures that show it --
## 10_7 --amid-g-anchor (Kevin, 2026-11-10) ends a peptide whose mature C
## terminus is followed by G + KK/KR/RK/RR ON THE G: the glycine is the amide
## donor and the residue the cleavage machinery actually leaves at the end.
## `pep_end` moves one residue for every use, so the C window is anchored on the
## G (position 28, basic pair at 29-30, where a plain dibasic end has its pair)
## and the G reads as peptide in the same peptide's N window too.
##
## Two figures read this:
##   10_7b_plot_pepend_windows.R  every known window, the moved G outlined
##   10_7c_plot_gshift_windows.R  the moved C ends before/after, with counts
##
## Kevin, 2026-11-12: "new plot of aligned known peptide windows showing new
## positioning of G dibasic". `lf_gshift_pairs()` is also the check: it stops
## unless the two builds differ by exactly the shift and the merge it implies.
##
## Needs LF_PEPEND and lf_pepend_insertion_class (pepend_windows.R),
## nn_class_cols / nn_class_labels (plot_proteins_win_new.R) and LF_VIZ
## (viz_theme.R).
## ------------------------------------------------------------------------------

#' Which cross-validation fold holds each known end out
#'
#' A `--n-folds K --fold k` build records only its own split, so a train/val
#' mark on its knowns says whether FOLD k holds an end out -- a fifth of the
#' answer. The K fold inputs share one file-name pattern (`..._f<k>_...rds`)
#' and hold the same knowns; reading the siblings' `knowns` finds, for every
#' end, the one fold whose validation set it is in.
#'
#' Costs one `readRDS()` per sibling: ~5 s and ~0.8 GB for a 90 MB fold input,
#' released before the next is read.
#'
#' @param in_path the fold input whose knowns are being drawn.
#' @param knowns its `$knowns` (`term`, `accession`, `anchor`, `split`).
#' @param params its `$params`; `n_folds` and `fold` locate the siblings.
#' @return integer vector, the held-out fold of each row of `knowns`; `NULL`
#'   (with a message) when the build is not cross-validated, a sibling is
#'   missing, or the folds do not partition the ends -- the caller then falls
#'   back to the build's own train/val split.
#' @export
lf_pepend_cv_folds <- function(in_path, knowns, params) {
  K <- params$n_folds; k0 <- params$fold
  if (!isTRUE(K > 0) || length(k0) != 1L || is.na(k0)) return(NULL)
  b   <- basename(in_path)
  tok <- sprintf("_f%d(?=[_.])", k0)
  if (lengths(regmatches(b, gregexpr(tok, b, perl = TRUE))) != 1L) {
    message("cv folds: cannot find _f", k0, " in ", b, "; showing this build's split")
    return(NULL)
  }
  sib <- file.path(dirname(in_path), vapply(seq_len(K), function(j)
    sub(tok, sprintf("_f%d", j), b, perl = TRUE), ""))
  if (!all(file.exists(sib))) {
    message("cv folds: missing ", paste(basename(sib[!file.exists(sib)]), collapse = ", "),
            "; showing this build's split")
    return(NULL)
  }
  key <- function(d) paste(d$term, d$accession, d$anchor)
  keys <- key(knowns)
  held <- vapply(seq_len(K), function(j) {
    if (j == k0) { kj <- knowns; pj <- params } else {
      x <- readRDS(sib[[j]]); kj <- x$knowns; pj <- x$params
      rm(x); invisible(gc())
    }
    ## a sibling from another build (other K, other fold) would partition
    ## nothing; neither would one whose knowns differ
    if (!identical(as.integer(pj$n_folds), as.integer(K)) ||
        !identical(as.integer(pj$fold), as.integer(j)) ||
        !setequal(key(kj), keys)) return(rep(NA, length(keys)))
    keys %in% key(kj)[kj$split == "val"]
  }, logical(length(keys)))
  held <- matrix(held, ncol = K)
  n <- rowSums(held)
  if (anyNA(n) || any(n != 1L)) {
    message(sprintf("cv folds: %d of %d ends are not held out in exactly one of %d folds; showing this build's split",
                    sum(is.na(n) | n != 1L), length(n), K))
    return(NULL)
  }
  as.integer(max.col(held, ties.method = "first"))
}

#' Where the moved G sits in each known window
#'
#' Under 10_7 `--amid-g-anchor` a peptide flagged `g_amid` ends on its G, and
#' `pep_end` moved for every use: the G is the C window's anchor (position 28)
#' and the last peptide residue of the same peptide's N window, at
#' `len + 7` -- inside the window when the shifted peptide is at most 29 long.
#'
#' @param knowns `$knowns` of a 10_7 build, `term` still `"N"`/`"C"`.
#' @param amid_g whether the build used the flag (`$params$amid_g`).
#' @return integer window position of the moved G per row; `NA` where there is
#'   none in the window, and for every row of a build without the flag. Stops if
#'   a row flagged as moved does not show a G there -- that is a build fault,
#'   not something to draw around.
#' @export
lf_pepend_gshift_pos <- function(knowns, amid_g) {
  n <- nrow(knowns)
  if (!isTRUE(amid_g) || !"g_amid" %in% names(knowns)) return(rep(NA_integer_, n))
  term <- as.character(knowns$term)
  pos <- ifelse(term == "C", LF_PEPEND$anchor[["C"]],
                LF_PEPEND$anchor[["N"]] + as.integer(knowns$len) - 1L)
  pos[!knowns$g_amid | pos > LF_PEPEND$seq_len] <- NA_integer_
  aa <- vapply(seq_len(n), function(i)
    if (is.na(pos[i])) NA_character_ else knowns$AA[[i]][[pos[i]]], "")
  bad <- which(!is.na(pos) & (is.na(aa) | aa != "G"))
  if (length(bad)) stop("g_amid rows without a G where the moved end should be: ",
                        paste(knowns$peps[bad], collapse = ", "))
  as.integer(pos)
}

#' The G and the basic pair's colours, on the tiles and in the strip
#'
#' LF_VIZ orange and aqua, not its blue: a blue outline sat beside
#' `pep_pocket`'s #197EC0 on the same figure and read as "inserting". A
#' function, not a constant, because the package collates this file before
#' viz_theme.R.
#' @keywords internal
lf_gshift_cols <- function() c(G = LF_VIZ$s2, basic = LF_VIZ$s3)

## residue at window position `p` of every window in a list of 36-long vectors
.lf_at <- function(v, p) vapply(seq_along(v), function(i) {
  q <- if (length(p) == 1L) p else p[[i]]
  if (is.na(q) || q < 1L || q > length(v[[i]])) NA_character_ else v[[i]][[q]]
}, "")
.lf_basic <- function(x) !is.na(x) & x %in% c("K", "R")

#' Pair the G-anchored build's C ends with the old build's, and check the shift
#'
#' The pairing is by accession + MATURE end: the new build keeps the mature end
#' in `pep_end_mature`, the old one anchors on it. Everything the figures state
#' is checked here, and the function stops if any of it fails:
#'
#' * a `g_amid` end has the G at 28 and a basic pair at 29-30 in the new build,
#'   the G at 29 and the pair at 30-31 in the old, and its new window is the
#'   old one shifted by one residue (`AA[1:35]` = old `AA[2:36]`); its labels
#'   are the old ones shifted, except the G, which goes CT-context -> peptide;
#' * every other C end has the same window and the same labels in both builds;
#' * an old C end with no counterpart sits exactly on a moved end's G, i.e. a
#'   moved end MERGED with a peptide that already ended on that G (GLP-1);
#' * the N windows are the same windows, and differ in their labels only at the
#'   moved G, where it falls inside the window.
#'
#' @param new,old `$knowns` of the G-anchored build and of the old-convention
#'   build of the same fold (only `split` differs between a build's folds).
#' @return a list: `pairs` (one row per moved C end; old and new `AA`,
#'   `known_idx`, and the G's position in each), `merged` (old ends that became
#'   one window with a moved end; `same_window` says whether that window is
#'   their old one, labels included), `other` (the unchanged C ends, new
#'   build), `regrouped` (N and C ends whose [lf_pepend_insertion_class()]
#'   group differs between the builds) and `counts`.
#' @export
lf_gshift_pairs <- function(new, old) {
  stopifnot(all(c("g_amid", "pep_end_mature") %in% names(new)))
  pep_cls <- c("pep_other", "pep_pocket")
  same <- function(a, b) mapply(identical, a, b)
  nc <- new[new$term == "C", ]; oc <- old[old$term == "C", ]
  kn <- paste(nc$accession, nc$pep_end_mature); ko <- paste(oc$accession, oc$anchor)
  stopifnot(!anyDuplicated(kn), !anyDuplicated(ko))
  m <- match(kn, ko)
  if (anyNA(m)) stop("new C ends with no old end at their mature end: ",
                     paste(nc$peps[is.na(m)], collapse = ", "))
  stopifnot(all(nc$anchor == nc$pep_end),
            all(nc$anchor == nc$pep_end_mature + as.integer(nc$g_amid)))

  ## the moved ends
  g <- which(nc$g_amid); o <- oc[m[g], ]; n <- nc[g, ]
  aC <- LF_PEPEND$anchor[["C"]]
  stopifnot(all(.lf_at(n$AA, aC) == "G"), all(.lf_basic(.lf_at(n$AA, aC + 1L))),
            all(.lf_basic(.lf_at(n$AA, aC + 2L))),
            all(.lf_at(o$AA, aC + 1L) == "G"), all(.lf_basic(.lf_at(o$AA, aC + 2L))),
            all(.lf_basic(.lf_at(o$AA, aC + 3L))),
            all(.lf_at(n$known_idx, aC) %in% pep_cls),
            all(.lf_at(o$known_idx, aC + 1L) == "CT_cleavage_context"),
            all(mapply(function(a, b) identical(a[1:35], b[2:36]), n$AA, o$AA)))
  relab <- mapply(function(a, b) which(a[1:35] != b[2:36]), n$known_idx, o$known_idx, SIMPLIFY = FALSE)
  stopifnot(all(vapply(relab, identical, NA, aC)))

  ## the unchanged ones
  u <- which(!nc$g_amid)
  stopifnot(all(same(nc$AA[u], oc$AA[m[u]])), all(same(nc$known_idx[u], oc$known_idx[m[u]])))

  ## old ends left over: each must be a moved end's G
  left <- oc[setdiff(seq_len(nrow(oc)), m), ]
  at <- match(paste(left$accession, left$anchor), paste(n$accession, n$anchor))
  if (anyNA(at)) stop("old C ends that neither pair nor merge: ",
                      paste(left$peps[is.na(at)], collapse = ", "))
  ## (keyed by pep_id: a merged end and the end it became are ONE window, so
  ## they share `peps`)
  merged <- left; merged$into <- n$pep_id[at]
  ## Kevin, 2026-11-12 (review): which peptide labels a merged window is a tie
  ## (both are the same length once one gains its G), so say what is CHECKED
  ## instead -- whether the window it became is its old window, residues and
  ## labels. Same anchor, so the residues always are; the labels are unless
  ## the two docked models put different residues in the pocket.
  merged$same_window <- mapply(function(a, b, x, y) identical(a, b) && identical(x, y),
                               left$AA, n$AA[at], left$known_idx, n$known_idx[at])

  ## the peptides followed by G + dibasic are the OLD windows' pep_ids; the
  ## new windows' also carry each merged peptide, which ends ON the G
  pep_g <- unique(unlist(strsplit(o$pep_ids, ";")))
  stopifnot(setequal(c(pep_g, unlist(strsplit(left$pep_ids, ";"))),
                     unlist(strsplit(n$pep_ids, ";"))))

  ## N windows: same windows, relabelled only at the moved G
  nn <- new[new$term == "N", ]; on <- old[old$term == "N", ]
  mn <- match(paste(nn$accession, nn$anchor), paste(on$accession, on$anchor))
  stopifnot(nrow(nn) == nrow(on), !anyNA(mn), all(same(nn$AA, on$AA[mn])))
  g_n <- lf_pepend_gshift_pos(nn, TRUE)
  dif <- mapply(function(a, b) which(a != b), nn$known_idx, on$known_idx[mn], SIMPLIFY = FALSE)
  stopifnot(all(mapply(function(d, p) if (is.na(p)) !length(d) else identical(d, p), dif, g_n)),
            all(.lf_at(on$known_idx[mn], g_n)[!is.na(g_n)] == "CT_cleavage_context"),
            all(.lf_at(nn$known_idx, g_n)[!is.na(g_n)] %in% pep_cls))

  ## 10_7b groups every end by where the pocket sits relative to it
  ## (lf_pepend_insertion_class), and a G-anchored end is measured from its G
  ## -- a residue no docked model of the mature peptide has in the pocket. So
  ## the shift can move an end to another group; keep which, to say so.
  grp <- function(k) lf_pepend_insertion_class(k$known_idx, k$term, k$len)
  on_m <- on[mn, ]
  rg <- function(a, b, i) tibble::tibble(
    term = a$term[i], gene = a$gene[i], name = dplyr::coalesce(a$pep_name[i], a$pep_id[i]),
    old = grp(b)[i], new = grp(a)[i])
  regrouped <- dplyr::bind_rows(rg(n, o, which(grp(n) != grp(o))),
                                rg(nn, on_m, which(grp(nn) != grp(on_m))))

  pairs <- tibble::tibble(
    id = n$pep_id, peps = n$peps, gene = n$gene, accession = n$accession,
    name = dplyr::coalesce(n$pep_name, sub("^[^_]+_", "", n$pep_id)),
    pep_start = n$pep_start, pep_end_mature = n$pep_end_mature, anchor = n$anchor,
    pep_ids = n$pep_ids, motif_old = o$motif, motif_new = n$motif,
    pair = paste0(.lf_at(n$AA, aC + 1L), .lf_at(n$AA, aC + 2L)),
    AA_old = o$AA, idx_old = o$known_idx, AA_new = n$AA, idx_new = n$known_idx,
    g_old = aC + 1L, g_new = aC)
  merged <- tibble::tibble(
    id = merged$pep_id, peps = merged$peps, into = merged$into, gene = merged$gene,
    name = dplyr::coalesce(merged$pep_name, sub("^[^_]+_", "", merged$pep_id)),
    pep_start = merged$pep_start, pep_end = merged$pep_end, anchor = merged$anchor,
    motif_old = merged$motif, AA_old = merged$AA, idx_old = merged$known_idx,
    same_window = merged$same_window, g_old = aC)
  stopifnot(all(.lf_at(merged$AA_old, aC) == "G"))
  list(pairs = pairs, merged = merged, other = nc[u, ], regrouped = regrouped,
       counts = list(old_C = nrow(oc), new_C = nrow(nc), moved = nrow(pairs),
                     unchanged = length(u), merged = nrow(merged),
                     peptides = length(pep_g),
                     pairs_by = table(pairs$pair),
                     ## what the G is labelled in its new windows (pep_other
                     ## or pep_pocket; checked above)
                     g_class = table(.lf_at(n$known_idx, aC)),
                     ## unchanged ends with a G after them: no basic PAIR
                     ## follows that G, so they keep the mature anchor
                     other_G = sum(.lf_at(nc$AA[u], aC + 1L) %in% "G"),
                     N_gamid = sum(nn$g_amid), N_relabelled = sum(!is.na(g_n))))
}

#' Before/after tiles of the moved C ends
#'
#' One row per moved end, its old window (mature end at 28, G at 29) beside its
#' new one (G at 28, basic pair at 29-30): tiles by training class, residue
#' letters, the cut line, the G and the basic pair outlined. An old end that
#' merged into a moved one gets a row of its own on the left and a note on the
#' right, under the end it became.
#'
#' @param gp the list from [lf_gshift_pairs()].
#' @param labels named `c(before = , after = )` panel titles.
#' @return a ggplot.
#' @export
lf_gshift_window_plot <- function(gp, labels = c(before = "BEFORE", after = "AFTER")) {
  aC <- LF_PEPEND$anchor[["C"]]; L <- LF_PEPEND$seq_len
  p <- gp$pairs %>% dplyr::arrange(gene, pep_end_mature)
  ## rows top to bottom; a merged old end sits directly under its partner
  ord <- unlist(lapply(p$id, function(x) c(x, gp$merged$id[gp$merged$into == x])))
  lab_p <- sprintf("%s %s %d-%d", p$gene, p$name, p$pep_start, p$pep_end_mature)
  lab_m <- sprintf("%s %s %d-%d (ends ON the G)", gp$merged$gene, gp$merged$name,
                   gp$merged$pep_start, gp$merged$pep_end)
  rows <- tibble::tibble(id = c(p$id, gp$merged$id), lab = c(lab_p, lab_m),
                         merged = c(rep(FALSE, nrow(p)), rep(TRUE, nrow(gp$merged))))
  rows <- rows[match(ord, rows$id), ]
  rows$y <- rev(seq_len(nrow(rows)))
  y_of <- stats::setNames(rows$y, rows$id)

  win <- function(id, panel, AA, idx, g) {
    pos <- rep(seq_len(L), length(id)); off <- pos - rep(g, each = L)
    tibble::tibble(id = rep(id, each = L), panel = panel, pos = pos,
                   AA = unlist(AA), class = unlist(idx),
                   emph = ifelse(off == 0L, "G", ifelse(off %in% 1:2, "basic", "none")))
  }
  tiles <- dplyr::bind_rows(
    win(p$id, "before", p$AA_old, p$idx_old, p$g_old),
    win(gp$merged$id, "before", gp$merged$AA_old, gp$merged$idx_old, gp$merged$g_old),
    win(p$id, "after", p$AA_new, p$idx_new, p$g_new)) %>%
    dplyr::mutate(y = unname(y_of[id]),
                  panel = factor(panel, levels = c("before", "after")),
                  class = factor(class, levels = names(nn_class_cols)))
  stopifnot(all(tiles$AA[tiles$emph == "G"] == "G"), all(.lf_basic(tiles$AA[tiles$emph == "basic"])))
  ## one outline per G, one around each basic PAIR
  g_box <- tiles %>% dplyr::filter(emph == "G") %>%
    dplyr::transmute(panel, y, xmin = pos - 0.5, xmax = pos + 0.5, what = "G")
  b_box <- tiles %>% dplyr::filter(emph == "basic") %>% dplyr::group_by(panel, y) %>%
    dplyr::summarise(xmin = min(pos) - 0.5, xmax = max(pos) + 0.5, .groups = "drop") %>%
    dplyr::mutate(what = "basic")
  boxes <- dplyr::bind_rows(g_box, b_box) %>%
    dplyr::mutate(what = factor(what, levels = c("G", "basic")))
  top <- max(rows$y) + 0.5
  vlines <- tibble::tibble(panel = factor(c("before", "after", "after"), levels = levels(tiles$panel)),
                           x = c(aC + 0.5, aC + 0.5, aC - 0.5),
                           kind = c("cut", "cut", "mature"),
                           note = c("cut: mature end", "cut: the G", "mature end"))
  merge_note <- tibble::tibble(
    panel = factor("after", levels = levels(tiles$panel)), y = unname(y_of[gp$merged$id]),
    label = sprintf("no window of its own: its end is the G of %s above (merged)",
                    p$name[match(gp$merged$into, p$id)]))
  ## a class the panels lack still gets a key, so this legend matches 10_7b's
  cls <- c("NT_cleavage_context", "pep_other", "pep_pocket", "CT_cleavage_context", "padding")

  ggplot2::ggplot(tiles, ggplot2::aes(pos, y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = class), colour = "white", linewidth = 0.25,
                       show.legend = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = AA, fontface = ifelse(emph == "none", "plain", "bold")),
                       size = 2.3, colour = LF_VIZ$ink, na.rm = TRUE) +
    ggplot2::geom_segment(data = vlines, ggplot2::aes(x = x, xend = x, y = 0.5, yend = top,
                                                      linetype = kind),
                          inherit.aes = FALSE, linewidth = 0.55, colour = LF_VIZ$ink,
                          show.legend = FALSE) +
    ## the outlines over the cut lines, which run along their edges
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - 0.45,
                                                  ymax = y + 0.45, colour = what),
                       inherit.aes = FALSE, fill = NA, linewidth = 0.75) +
    ggplot2::geom_text(data = vlines, ggplot2::aes(x = x, y = top + 0.35, label = note,
                                                   hjust = ifelse(kind == "mature", 1.04, -0.04)),
                       inherit.aes = FALSE, vjust = 0, size = 2.4, colour = LF_VIZ$ink2) +
    ggplot2::geom_text(data = merge_note, ggplot2::aes(x = 0.7, y = y, label = label),
                       inherit.aes = FALSE, hjust = 0, size = 2.3, fontface = "italic",
                       colour = LF_VIZ$ink2) +
    ggplot2::facet_wrap(~panel, nrow = 1, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_fill_manual(values = nn_class_cols[cls], breaks = cls,
                               labels = nn_class_label_fn, name = "label", drop = FALSE) +
    ggplot2::scale_colour_manual(values = lf_gshift_cols(),
                                 labels = c(G = "the G", basic = "basic pair (KK/KR/RK/RR)"),
                                 name = NULL) +
    ggplot2::scale_linetype_manual(values = c(cut = "solid", mature = "22")) +
    ## (show.legend on the tiles would also key the cut lines' linetype)
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(fill = NA, linewidth = 0.9),
                                                   order = 2),
                    fill = ggplot2::guide_legend(order = 1), linetype = "none") +
    ggplot2::scale_x_continuous(breaks = c(1, 8, 15, 22, 28, 36), expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::scale_y_continuous(breaks = rows$y, labels = rows$lab,
                                expand = ggplot2::expansion(add = c(0.6, 1.6))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "position in window", y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = LF_VIZ$surface, colour = NA),
                   strip.text = ggplot2::element_text(face = "bold", size = 9.5, hjust = 0,
                                                      colour = LF_VIZ$ink),
                   axis.text.y = ggplot2::element_text(size = 7, colour = LF_VIZ$ink),
                   axis.text.x = ggplot2::element_text(colour = LF_VIZ$ink2),
                   axis.title.x = ggplot2::element_text(colour = LF_VIZ$ink2),
                   panel.spacing.x = grid::unit(1.2, "lines"),
                   legend.position = "bottom")
}

#' How often each position holds the G and a basic residue, around the C anchor
#'
#' Three groups on one position axis: the moved ends at their old anchor and at
#' their new one, and every other C end (the same windows in both builds). The
#' bars are the share of ends with G, or with K/R, at each position; a padded
#' position counts as neither.
#'
#' @param gp the list from [lf_gshift_pairs()].
#' @param pos window positions to show.
#' @param groups named `c(old = , new = , other = )` facet titles.
#' @return a ggplot.
#' @export
lf_gshift_freq_plot <- function(gp, pos = 22:34,
                                groups = c(old = "moved ends, OLD anchor",
                                           new = "moved ends, NEW anchor",
                                           other = "all other C ends")) {
  aC <- LF_PEPEND$anchor[["C"]]
  ## (`ps`, not `pos`: tibble() evaluates its columns in order, so a `pos`
  ## column would shadow the argument in the columns after it)
  share <- function(AA, grp, ps = pos) {
    m <- do.call(rbind, AA)[, ps, drop = FALSE]
    tibble::tibble(grp = grp, pos = rep(ps, 2), res = rep(c("G", "K/R"), each = length(ps)),
                   share = c(colMeans(matrix(m %in% "G", nrow(m))),
                             colMeans(matrix(m %in% c("K", "R"), nrow(m)))))
  }
  d <- dplyr::bind_rows(share(gp$pairs$AA_old, "old"), share(gp$pairs$AA_new, "new"),
                        share(gp$other$AA, "other")) %>%
    dplyr::mutate(grp = factor(grp, levels = names(groups)), res = factor(res, levels = c("G", "K/R")))
  ## direct labels only where a bar is (all but) every end, placed by hand: a
  ## dodge over the labelled subset would centre a lone label on its position
  ## instead of over its own bar
  dw <- 0.84
  top <- d %>% dplyr::filter(share >= 0.95) %>%
    dplyr::mutate(x = pos + ifelse(res == "G", -dw / 4, dw / 4))
  ggplot2::ggplot(d, ggplot2::aes(pos, share, fill = res)) +
    ## the cut stops at 100%: run to the top, it crossed the "100%" over the
    ## old anchor's G bar, which sits just right of it
    ggplot2::annotate("segment", x = aC + 0.5, xend = aC + 0.5, y = 0, yend = 1,
                      linewidth = 0.55, colour = LF_VIZ$ink) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = dw), width = 0.8) +
    ggplot2::geom_text(data = top, ggplot2::aes(x = x, label = sprintf("%d%%", round(100 * share))),
                       vjust = -0.45, size = 2.3, colour = LF_VIZ$ink) +
    ggplot2::facet_wrap(~grp, nrow = 1, labeller = ggplot2::as_labeller(groups)) +
    ggplot2::scale_fill_manual(values = stats::setNames(lf_gshift_cols(), c("G", "K/R")),
                               labels = c(G = "G", `K/R` = "K or R"), name = "residue") +
    ggplot2::scale_x_continuous(breaks = pos, expand = ggplot2::expansion(add = 0.4)) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                                breaks = c(0, 0.5, 1), limits = c(0, 1.17),
                                expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::labs(x = "position in window (C anchor at 28; the line is the cut)",
                  y = "share of ends") +
    lf_viz_theme(base_size = 9) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 9, hjust = 0,
                                                      colour = LF_VIZ$ink),
                   axis.text.x = ggplot2::element_text(size = 7),
                   panel.spacing.x = grid::unit(1.4, "lines"),
                   legend.position = "right")
}

#' The before/after figure of the G-anchored C ends, with its counts
#'
#' [lf_gshift_window_plot()] over [lf_gshift_freq_plot()], titled with the
#' counts [lf_gshift_pairs()] checked.
#'
#' @param gp the list from [lf_gshift_pairs()].
#' @param old_name,new_name how the two builds are named on the figure.
#' @return a patchwork.
#' @export
lf_gshift_figure <- function(gp, old_name = "old build", new_name = "G-anchored build") {
  cn <- gp$counts
  mg <- gp$merged
  into <- gp$pairs[match(mg$into, gp$pairs$id), ]
  ## one line per merge: the full names are long (GLP-1's are 30+ characters)
  merge_txt <- if (nrow(mg)) paste0("\n", paste(sprintf(
    "MERGED: %s %s (%d-%d) moved onto G%d, where %s (%d-%d) already ended -- one window now, %s.",
    into$gene, into$name, into$pep_start, into$pep_end_mature, into$anchor, mg$name, mg$pep_start,
    mg$pep_end, ifelse(mg$same_window, "the same residues and labels as that end's old window",
                       "that end's old residues with other pocket labels")), collapse = "\n")) else ""
  ## what the G became, from the data: "peptide (non-inserting) in all 30"
  gl <- cn$g_class
  g_lab <- if (length(gl) == 1L) sprintf("%s in all %d", nn_class_label_fn(names(gl)), as.integer(gl)) else
    paste(sprintf("%s in %d", nn_class_label_fn(names(gl)), as.integer(gl)), collapse = ", ")
  ## Kevin, 2026-11-12 (review): 10_7b's groups count the pocket's distance
  ## from each end, and the C end is now the G -- outside every docked pocket
  ## -- so a few rows sit in another group than on the plot he knows. Name them.
  rg <- gp$regrouped
  rg_txt <- if (nrow(rg)) {
    st <- function(x) gsub("_", "-", x)
    by_t <- rg %>%
      dplyr::mutate(term = factor(term, levels = c("C", "N")), step = sprintf("%s -> %s", st(old), st(new))) %>%
      dplyr::group_by(term, step) %>%
      dplyr::summarise(who = paste(name, collapse = ", "), .groups = "drop") %>%
      dplyr::group_by(term) %>%
      dplyr::summarise(txt = paste(sprintf("%s (%s)", who, step), collapse = ", "), .groups = "drop")
    ## a sentence, then the names on a line of their own: wrapped as one
    ## paragraph, a break could fall inside "(inserting -> loop)"
    sprintf("\n10_7b groups ends by the pocket's distance from them, and the C end is now the G%s: %d end%s change group.\n%s.",
            if (identical(names(gl), "pep_other")) ", in no docked pocket" else "",
            nrow(rg), if (nrow(rg) == 1L) "" else "s",
            paste(sprintf("%s ends: %s", by_t$term, by_t$txt), collapse = "; "))
  } else ""
  p1 <- lf_gshift_window_plot(gp, labels = c(
    before = "BEFORE: mature end at 28, the G at 29, basic pair at 30-31",
    after  = "AFTER (--amid-g-anchor): the G at 28, basic pair at 29-30"))
  p2 <- lf_gshift_freq_plot(gp, groups = c(
    old = sprintf("the %d moved ends, OLD anchor", cn$moved),
    new = sprintf("the %d moved ends, NEW anchor", cn$moved),
    other = sprintf("the other %d C ends (unchanged)", cn$unchanged)))
  ## (free the strip's y title: aligned, it would sit out past the tile
  ## panel's long row labels, far from its own axis)
  patchwork::wrap_plots(p1, patchwork::free(p2, type = "label", side = "l"), ncol = 1,
                        heights = c(nrow(gp$pairs) + nrow(mg), 9.5)) +
    patchwork::plot_annotation(
      title = sprintf("G + dibasic C ends now end ON the G: %d of %d known C ends moved one residue%s",
                      cn$moved, cn$old_C,
                      if (cn$merged) sprintf("; %d merged, leaving %d", cn$merged, cn$new_C) else ""),
      subtitle = paste0(
        ## (the peptides are the OLD windows' pep_ids: a merged peptide ends
        ## ON the G, so it is not one of them)
        sprintf(paste0("%d C ends (%d peptides; some share an end) are followed by G + KK/KR/RK/RR. --amid-g-anchor ",
                       "anchors each ON the G, one past the mature (amidated) end.\n",
                       "Each new window is its old one shifted by one residue; its labels are the old ones shifted, ",
                       "except the G: CT-context -> %s. The other %d C windows are identical in both builds."),
                cn$moved, cn$peptides, g_lab, cn$unchanged),
        merge_txt,
        sprintf(paste0("\nN windows: the G is peptide there too. %d of these peptides' %d N windows reach it and ",
                       "are relabelled at it (outlined in 10_7b); nothing else in the N set changes."),
                cn$N_relabelled, cn$N_gamid),
        rg_txt),
      caption = paste0(
        if (cn$other_G) sprintf(paste0("Right strip: %d of the other %d C ends are followed by a G with no basic pair after it ",
                                       "(G + monobasic or other), so they keep the mature anchor -- the G bar at 29.\n"),
                                cn$other_G, cn$unchanged) else "",
        sprintf(paste0("Before: %s. After: %s. Paired by accession + mature end. ",
                       "Basic pair after the G: %s. Every statement above is computed and checked by ",
                       "lf_gshift_pairs() (R/pepend_gshift_plots.R) before the figure is drawn."),
                old_name, new_name,
                paste(sprintf("%s %d", names(cn$pairs_by), as.integer(cn$pairs_by)), collapse = ", "))),
      theme = lf_viz_theme(base_size = 9) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 13),
                       plot.subtitle = ggplot2::element_text(size = 9.5, lineheight = 1.2))) &
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = LF_VIZ$surface, colour = NA))
}
