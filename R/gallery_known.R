## Known peptides on a gallery index: one chip per end per arm, showing and
## coloured by the score at that end's anchor (Kevin, 2026-11-06: "instead of
## TP/FN color by the actual mean score for the anchor"). A TP / FN chip says
## which side of one threshold an end fell on; this says how far.
##
## Built on R/gallery_tpfn.R: an arm is an lf_tpfn_arm(), and the cards go through
## make_gallery_index.R the same way. Scores come from the arms' full-proteome
## tracks, the rds end_calls.R reads, so no scan or page is needed.

#' The arms of the two-arm gallery
#'
#' `head` is pinned to the xgb window head: it is what the end segments and the
#' candidates are cut from, so it is the score a page's candidates rest on.
#' `thr` is each arm's segment cut, here only to be quoted in the tooltips.
#' @keywords internal
lf_known_arms <- function() list(
  lf_tpfn_arm("swap", "AAswap",     "_curated_AAswap", "0.5",  short = "AAswap", head = "xgb_pm"),
  lf_tpfn_arm("esm",  "ESM-C only", "_esmonly36",      "0.65", short = "ESM&#8209;C", head = "xgb_pm")
)

#' Every known peptide with both its ends
#'
#' Two sources, disjoint by gene: the training peptides and the never-seen ones
#' (whose gene no fold trained on). Large cleavage products are left out; they
#' are chains, not peptides.
#'
#' A training end shared by several peptides is ONE row of `knowns`, named for
#' just one of them, so the peptides come from `pep_ids`. Reading `pep_id` alone
#' loses a peptide or leaves it with one end (PACAP-38 shares PACAP-27's N end).
#' @return a data.table of gene, accession, name, start, end and set.
#' @keywords internal
lf_known_peptides <- function(a2 = "~/AF2_analysis") {
  a2 <- path.expand(a2)
  k <- data.table::as.data.table(readRDS(file.path(a2, "lf_pepend_nn_input_foldA.rds"))$knowns)
  ids <- unique(unlist(strsplit(k$pep_ids, ";", fixed = TRUE)))
  m <- regmatches(ids, regexec("^(.+)_(-?[0-9]+)x(-?[0-9]+)$", ids))
  stopifnot(all(lengths(m) == 4L))
  tr <- data.table::data.table(pep_id = ids, gene = vapply(m, `[`, "", 2L),
                               start = as.integer(vapply(m, `[`, "", 3L)),
                               end = as.integer(vapply(m, `[`, "", 4L)))
  tr$accession <- k$accession[match(tr$gene, k$gene)]
  tr$name <- k$pep_name[match(tr$pep_id, k$pep_id)]
  tr$name <- ifelse(is.na(tr$name), tr$gene, tr$name)
  tr$set <- "training"

  un <- data.table::fread(file.path(a2, "lf_nomask_scripts_20261021", "cv5_unseen_ends.csv"))
  un <- un[grepl("^never", un$set), ]
  un$start <- ifelse(un$term == "N", un$anchor, un$anchor - un$len + 1L)
  un$end   <- ifelse(un$term == "N", un$anchor + un$len - 1L, un$anchor)
  ns <- unique(un[, c("gene", "accession", "name", "start", "end")])
  ns$set <- "never seen"
  stopifnot(!any(ns$gene %in% tr$gene))
  data.table::rbindlist(list(tr[, c("gene", "accession", "name", "start", "end", "set")], ns),
                        use.names = TRUE)
}

#' Each arm's scores at each end of each peptide
#'
#' `score` is the number a chip shows. For a never-seen peptide that is the
#' all-50 mean: no fold trained on its gene. For a training end the all-50 is
#' 40/50 in-sample (median 0.91 against 0.75 held out on AAswap), so `training =
#' "heldout"` takes the mean of the fold that held the end out, as end_calls.R
#' does; `"all50"` takes the page's heavy line instead. An end of a training
#' peptide that was never itself a training end has no fold and keeps the all-50.
#' @return a data.table with one row per peptide end per arm.
#' @keywords internal
lf_known_scores <- function(peps, arms = lf_known_arms(), a2 = "~/AF2_analysis",
                            training = c("heldout", "all50")) {
  training <- match.arg(training)
  a2 <- path.expand(a2)
  ends <- unique(data.table::rbindlist(lapply(c("N", "C"), function(tm) {
    x <- data.table::copy(peps); x$term <- tm
    x$anchor <- if (tm == "N") x$start else x$end
    x
  })))
  folds <- data.table::rbindlist(lapply(c("N", "C"), function(tm) {
    d <- data.table::fread(file.path(a2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)),
                           select = c("accession", "anchor", "fold"))
    d$term <- tm
    d
  }))
  ends$fold <- folds$fold[match(paste(ends$accession, ends$term, ends$anchor),
                                paste(folds$accession, folds$term, folds$anchor))]
  data.table::rbindlist(lapply(arms, function(a) data.table::rbindlist(lapply(c("N", "C"), function(tm) {
    p <- file.path(a2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, a$suffix))
    if (!file.exists(p)) stop("no ", basename(p), " for arm '", a$label, "'", call. = FALSE)
    b <- data.table::as.data.table(readRDS(p))
    e <- ends[ends$term == tm, ]
    i <- match(paste(e$accession, e$anchor), paste(b$accession, b$anchor))
    fold_val <- function(hd) vapply(seq_along(i), function(r)
      if (is.na(i[r]) || is.na(e$fold[r])) NA_real_
      else as.numeric(b[[sprintf("%s__f%d", hd, e$fold[r])]][i[r]]), 0)
    e$all50 <- b[[paste0(a$head, "__all50")]][i]
    e$held  <- fold_val(a$head)
    e$cnn_all50 <- b$cnn__all50[i]
    e$cnn_held  <- fold_val("cnn")
    e$arm <- a$key
    e$score <- if (training == "heldout") ifelse(is.na(e$held), e$all50, e$held) else e$all50
    e
  }))))
}

## turbo over 0..1: the candidate track's scale, so a colour reads the same on the
## index and on the page
.LF_KNOWN_RAMP <- substr(viridisLite::turbo(33), 1L, 7L)

#' Fill and ink for a score
#' @keywords internal
.lf_known_col <- function(v) {
  rgb <- grDevices::colorRamp(.LF_KNOWN_RAMP)(pmin(pmax(v, 0), 1)) / 255
  lum <- function(x) { x <- ifelse(x <= 0.03928, x / 12.92, ((x + 0.055) / 1.055)^2.4)
                       0.2126 * x[, 1] + 0.7152 * x[, 2] + 0.0722 * x[, 3] }
  L <- lum(rgb)
  list(fill = grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3]),
       ## whichever of the two inks reads better on this fill
       ink = ifelse((L + 0.05) / (lum(matrix(11 / 255, 1, 3)) + 0.05) >= 1.05 / (L + 0.05),
                    "#0b0b0b", "#ffffff"))
}

#' The per-gene card HTML, the CSS it needs, and each card's link
#'
#' @param peps,sc from [lf_known_peptides()] and [lf_known_scores()].
#' @param page_of named character vector, accession to the page's gene name.
#' @param base put before each link, for an index that is not beside its pages.
#' @return a list with `cards` (gene, html), `href` (gene, href) and `style`.
#' @keywords internal
lf_known_cards <- function(peps, sc, page_of, arms = lf_known_arms(), base = "") {
  AL <- stats::setNames(vapply(arms, `[[`, "", "label"), vapply(arms, `[[`, "", "key"))
  AT <- stats::setNames(vapply(arms, `[[`, "", "thr"), vapply(arms, `[[`, "", "key"))
  tip <- function(r) {
    what <- if (is.na(r$held)) "all-50 mean" else sprintf("held-out fold %d mean of 10 seeds", r$fold)
    paste0(sprintf("%s xgb head, %s %d: %.3f (%s)", AL[[r$arm]], r$term, r$anchor, r$score, what),
           if (!is.na(r$held)) sprintf("; all-50 mean %.3f, as drawn on the page", r$all50),
           sprintf("; linear head %.3f; this arm's segment cut %s",
                   if (is.na(r$held)) r$cnn_all50 else r$cnn_held, AT[[r$arm]]))
  }
  chip <- function(r) {
    if (!nrow(r) || is.na(r$score)) return('<span class="kc na" title="no window at this anchor">&middot;</span>')
    col <- .lf_known_col(r$score)
    sprintf('<span class="kc" style="background:%s;color:%s" title="%s">%.2f</span>',
            col$fill, col$ink, .lf_tpfn_esc(tip(r)), r$score)
  }
  hdr <- paste0(vapply(arms, function(a) sprintf("<span>%s</span>", a$short), ""), collapse = "")
  pg <- unname(page_of[peps$accession])
  genes <- sort(unique(stats::na.omit(pg)))
  rows <- data.table::rbindlist(lapply(genes, function(g) {
    p <- peps[which(pg == g), ]
    p <- p[order(p$start, p$end), ]
    body <- vapply(seq_len(nrow(p)), function(i) {
      x <- p[i, ]
      e <- sc[sc$accession == x$accession & sc$start == x$start & sc$end == x$end, ]
      ends <- vapply(c("N", "C"), function(tm) {
        y <- e[e$term == tm, ]
        sprintf('<div class="knr"><span>%s %d</span>%s</div>', tm, y$anchor[1],
                paste0(vapply(arms, function(a) chip(y[y$arm == a$key, ]), ""), collapse = ""))
      }, "")
      paste0(sprintf('<div class="knp">%s %d&ndash;%d <span class="kt">%s</span></div>',
                     .lf_tpfn_esc(x$name), x$start, x$end, x$set), paste(ends, collapse = ""))
    }, "")
    data.table::data.table(gene = g, accession = p$accession[1],
      href = sprintf("%s%s.html#%sx%dx%d", base, g, g, p$start[1], p$end[1]),
      html = paste0('<div class="kn"><div class="knh"><span>known peptides</span>', hdr, '</div>',
                    paste(body, collapse = ""), '</div>'))
  }))
  style <- paste0('<style>',
    ## as in lf_tpfn_cards(): flex + order lifts the block above the card's text
    '.card{display:flex;flex-direction:column}',
    '.kn{order:-1;margin:0 0 8px;padding-bottom:7px;border-bottom:1px solid var(--grid);font-size:11px}',
    sprintf('.knh,.knr{display:grid;grid-template-columns:1fr%s;align-items:center;column-gap:3px}',
            strrep(" 44px", length(arms))),
    '.knh{color:var(--ink2);font-size:10px;margin-bottom:3px}',
    '.knh span+span{text-align:center;white-space:nowrap}',
    '.knp{color:var(--ink);font-weight:600;margin-top:4px;overflow-wrap:anywhere}',
    '.kt{color:var(--ink2);font-weight:400}.knr{color:var(--ink2);margin-top:1px}',
    '.kc{display:block;text-align:center;border-radius:3px;font-size:10.5px;font-weight:650;',
    'line-height:16px;font-variant-numeric:tabular-nums}',
    '.kc.na{color:var(--ink2);font-weight:400}',
    '.kramp{display:inline-block;width:160px;height:10px;border-radius:2px;vertical-align:-1px;margin:0 6px;background:',
    sprintf('linear-gradient(to right,%s)', paste(.LF_KNOWN_RAMP, collapse = ",")), '}',
    '</style>')
  list(cards = as.data.frame(rows[, c("gene", "html")]), href = as.data.frame(rows[, c("gene", "href")]),
       accession = as.data.frame(rows[, c("gene", "accession")]), style = style)
}

#' Build the known-peptide index of a gallery
#'
#' Writes `out` into `dir`, a card for each gallery page carrying a known peptide.
#' The gallery's own `index.html` is left alone.
#' @param dir the gallery directory of `<GENE>.html` pages.
#' @param out the file to write: inside `dir`, unless an absolute path.
#' @param training how a training end is scored, see [lf_known_scores()].
#' @param base put before each link, for an index that is not beside its pages.
#' @param script path to `make_gallery_index.R`.
#' @return invisibly, the scores, for checking.
#' @keywords internal
lf_gallery_known_index <- function(dir, out = "index_known.html", arms = lf_known_arms(),
                                   training = c("heldout", "all50"),
                                   title = "LigandFinder — known peptides, AAswap and ESM-C only",
                                   a2 = "~/AF2_analysis", workdir = tempfile("lf_known_"),
                                   script, base = "") {
  training <- match.arg(training)
  dir <- path.expand(dir)
  peps <- lf_known_peptides(a2)
  sc <- lf_known_scores(peps, arms, a2, training)
  ## a known gene's name can be two genes' ("DEFB4A; DEFB4B"), a page's cannot
  first <- sub(";.*$", "", peps$gene)
  have <- file.exists(file.path(dir, paste0(first, ".html")))
  page_of <- stats::setNames(first[have], peps$accession[have])
  page_of <- page_of[!duplicated(names(page_of))]
  cc <- lf_known_cards(peps, sc, page_of, arms, base)
  miss <- unique(peps$gene[!have])
  if (length(miss)) message("no page, left off: ", paste(miss, collapse = ", "))

  shown <- sc[sc$accession %in% names(page_of), ]
  n_pep <- nrow(unique(shown[, c("accession", "start", "end")]))
  th <- paste(vapply(arms, function(a) sprintf("%s %s", a$label, a$thr), ""), collapse = ", ")
  blurb <- paste0(cc$style, sprintf(paste0(
    '<p>The %d genes with a known peptide: %d peptides, the training set and the never-seen peptides of ',
    'under 60 aa. Larger cleavage products of never-seen genes are left out. Each card lists its peptides ',
    'with both ends, and each chip is one arm&rsquo;s <b>xgb window head</b> score at the anchor on that ',
    'end, both number and colour: <span>0</span><span class="kramp"></span><span>1</span></p>\n',
    '<p><b>Which mean.</b> <i>never seen</i>: the all-50 mean, the heavy squares on the page; no fold ',
    'trained on these genes. <i>training</i>: %s. Hover a chip for both, the linear head and the fold.</p>\n',
    '<p>Read a column on its own: the arms&rsquo; scales differ, as their segment cuts do (%s). ',
    'A dot means the anchor has no window.</p>'),
    length(unique(page_of)), n_pep,
    if (training == "heldout")
      paste0('the mean of the fold that held that end out, 10 seeds. The page&rsquo;s heavy line would ',
             'overstate these: 40 of its 50 members trained on the end')
    else 'the all-50 mean, the heavy squares on the page, though 40 of those 50 members trained on the end',
    th))

  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  f <- function(x, n) { p <- file.path(workdir, n); saveRDS(x, p); p }
  blurb_p <- file.path(workdir, "blurb.html"); writeLines(blurb, blurb_p)
  ord <- data.frame(gene = cc$cards$gene, rank = seq_len(nrow(cc$cards)))
  out <- path.expand(out)
  args <- c("--dir", dir, "--out", if (startsWith(out, "/")) out else file.path(dir, out), "--title", title,
            "--card-html", f(cc$cards, "cards.rds"), "--card-href", f(cc$href, "href.rds"),
            "--card-accession", f(cc$accession, "acc.rds"), "--card-order", f(ord, "order.rds"),
            "--only-ranked", "--no-stock-desc", "--blurb", blurb_p)
  rc <- system2("/usr/local/bin/Rscript", shQuote(c(script, args)))
  if (rc != 0L) stop("make_gallery_index.R failed, rc = ", rc, call. = FALSE)
  invisible(sc)
}
