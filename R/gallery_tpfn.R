## TP / FN chips on a gallery index, one column per model arm.
##
## The per-gene pages are self-contained and `make_gallery_index.R` turns a
## directory of them into an index. This adds the layer the peptide-end
## galleries want on top: every never-seen end called TP or FN by each arm, at
## the top of that gene's card.
##
## It was a loose script in ~/AF2_analysis/lf_nomask_scripts_20261021 that had
## to be rewritten for each new arm (Kevin, 2026-10-26: "wire this index.html
## generation code into the package so i don't have to keep asking"). The arm
## set is now DATA -- a list of `lf_tpfn_arm()` -- so adding an arm is one more
## entry, not another copy of the script.

#' Describe one model arm for a TP / FN gallery index
#'
#' @param key short id used internally and in the chip tooltips.
#' @param label the arm's full name, e.g. `"ESM-C only"`.
#' @param suffix the arm's file-name suffix, `""` for the curated baseline and
#'   e.g. `"_esmonly36"` otherwise. This is what names its `end_calls_*` and
#'   `segment_confusion_*` files.
#' @param thr the segment threshold this arm is called at, as a STRING, because
#'   it appears verbatim in the file name (`"0.65"`, not `0.65`).
#' @param short an abbreviation for the column head; defaults to `label`.
#' @param head optionally pin the window head (`"cnn"` or `"xgb_pm"`); the
#'   default `NA` picks it with [lf_tpfn_pick_head()].
#' @return a list of class `lf_tpfn_arm`.
#' @export
lf_tpfn_arm <- function(key, label, suffix, thr, short = label, head = NA_character_) {
  stopifnot(is.character(key), length(key) == 1L, nzchar(key),
            is.character(label), length(label) == 1L,
            is.character(suffix), length(suffix) == 1L,
            is.character(thr), length(thr) == 1L,      # a string: it is a file name
            is.na(head) || head %in% c("cnn", "xgb_pm"))
  structure(list(key = key, label = label, suffix = suffix, thr = thr,
                 short = short, head = head), class = "lf_tpfn_arm")
}

#' The arms the peptide-end galleries currently draw
#'
#' Each at its own operating point: the curated model at the threshold its
#' gallery has always used, the others at the optimum measured for them. The
#' TOLERANCE is deliberately not here -- it belongs to the evaluation, not to a
#' model, so letting it differ per arm would put a second, unasked-for
#' difference into the chips.
#' @export
lf_tpfn_default_arms <- function() list(
  lf_tpfn_arm("cur", "curated",      "",               "0.2",  short = "curated"),
  lf_tpfn_arm("esm", "ESM-C only",   "_esmonly36",     "0.65", short = "ESM&#8209;C"),
  lf_tpfn_arm("waa", "curated+AAid", "_curated_wAAid", "0.5",  short = "+AA id")
)

#' Choose an arm's window head on segment auc_T
#'
#' Both heads of an arm sit at the SAME threshold, so they can be compared
#' directly -- but **not** on TP count, which rewards whichever head calls more
#' segments rather than the better one. At T = 0.2 the curated xgb head finds
#' exactly the same never-seen peptide ends as the linear head while calling
#' about a third more false-positive segments. `auc_T` = (1 + TPR - FPR)/2
#' charges for both, so it is what decides here.
#'
#' The confusion table holds FOUR rows per (set, threshold): two termini by two
#' heads. They are averaged. Reading one terminus picks a head on half the
#' evidence, and the two termini do disagree -- which is a mistake this function
#' exists to stop anyone repeating.
#'
#' @param arm an [lf_tpfn_arm()].
#' @param tol segment tolerance in residues.
#' @param headcmp directory holding the `segment_confusion_*` tables.
#' @param eval_set which evaluation set decides. Named `eval_set`, not `set`:
#'   inside `[.data.table` a bare symbol resolves to the COLUMN of that name
#'   first, so an argument called `set` would silently compare the column with
#'   itself and match every row.
#' @return `arm` with `head`, `auc`, `margin`, `fp` and `fpr` filled in.
#' @export
lf_tpfn_pick_head <- function(arm, tol = 4, headcmp = "~/AF2_analysis/lf_headcmp",
                              eval_set = "never seen: peptides") {
  stopifnot(inherits(arm, "lf_tpfn_arm"))
  f <- file.path(path.expand(headcmp),
                 sprintf("segment_confusion%s_tol%s.csv", arm$suffix, tol))
  if (!file.exists(f))
    stop("no ", basename(f), " for arm '", arm$label, "'.\n",
         "  Run segments_roc.R for that arm (ARM=", arm$suffix, " TOL=", tol, ").",
         call. = FALSE)
  d <- data.table::fread(f)
  ## the mask is computed OUTSIDE `[` on purpose: evaluated inside, a bare
  ## symbol that matches a column name resolves to the column.
  keep <- d$set == eval_set & abs(d$thr - as.numeric(arm$thr)) < 1e-9
  d <- d[keep, ]
  if (!nrow(d))
    stop("no rows in ", basename(f), " for set '", eval_set, "' at T = ", arm$thr,
         ".\n  Re-run segments_roc.R, or correct the arm's `thr`.", call. = FALSE)
  if (nrow(d) != 4L || data.table::uniqueN(d$head) != 2L || data.table::uniqueN(d$term) != 2L)
    stop("expected 4 rows (2 termini x 2 heads) in ", basename(f), " at T = ",
         arm$thr, ", got ", nrow(d), call. = FALSE)
  g <- d[, list(auc_T = mean(.SD$auc_T), FP = sum(.SD$FP), FPR = mean(.SD$FPR)), by = "head"]
  if (!is.na(arm$head)) {
    pick <- g[g$head == arm$head, ]
    if (!nrow(pick)) stop("arm '", arm$label, "' pins head '", arm$head,
                          "', which is not in ", basename(f), call. = FALSE)
  } else {
    pick <- g[which.max(g$auc_T), ]
  }
  other <- g[g$head != pick$head, ]
  arm$head   <- pick$head
  arm$auc    <- pick$auc_T
  arm$margin <- pick$auc_T - other$auc_T[1]
  arm$fp     <- pick$FP
  arm$fpr    <- pick$FPR
  arm
}

.lf_tpfn_esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE)
}

.LF_TPFN_HN <- c(cnn = "linear head", xgb_pm = "xgb head")
.LF_TPFN_HS <- c(cnn = "lin", xgb_pm = "xgb")

#' Read one arm's end calls
#' @keywords internal
.lf_tpfn_calls <- function(arm, tol, headcmp) {
  f <- sprintf("end_calls_T%s%s_tol%s.csv", arm$thr, arm$suffix, tol)
  p <- file.path(path.expand(headcmp), f)
  if (!file.exists(p))
    stop("no ", f, " for arm '", arm$label, "'.\n",
         "  Run end_calls.R for it: THRESH=", arm$thr, " ARM=", arm$suffix,
         " TOL=", tol, " Rscript end_calls.R", call. = FALSE)
  d <- data.table::fread(p)
  d <- d[grepl("^never", d$set), ]
  ## the all-50 profile only -- the heavy squares on the page. A per-fold row
  ## would silently multiply every end by five.
  if ("fold" %in% names(d)) d <- d[is.na(d$fold), ]
  d <- d[d$head == arm$head, ]
  d$arm <- arm$key
  d
}

#' The per-gene card HTML, and the CSS it needs
#'
#' @param arms list of arms, already through [lf_tpfn_pick_head()].
#' @param genes gene symbols to build cards for.
#' @param tol,headcmp as in [lf_tpfn_pick_head()].
#' @return a list with `cards` (a data.frame of `gene` and `html`), `style`
#'   (a `<style>` block) and `calls` (the pooled call table, for the blurb).
#' @export
lf_tpfn_cards <- function(arms, genes, tol = 4, headcmp = "~/AF2_analysis/lf_headcmp") {
  ec <- data.table::rbindlist(lapply(arms, .lf_tpfn_calls, tol = tol, headcmp = headcmp),
                              fill = TRUE)
  if (!nrow(ec)) stop("no never-seen end calls for any arm", call. = FALSE)
  ec$start <- ifelse(ec$term == "N", ec$anchor, ec$anchor - ec$len + 1L)
  ec$end   <- ifelse(ec$term == "N", ec$anchor + ec$len - 1L, ec$anchor)
  AL <- stats::setNames(vapply(arms, `[[`, "", "label"), vapply(arms, `[[`, "", "key"))

  tip <- function(r) {
    pk <- if (!is.na(r$peak_centre))
      sprintf("nearest peak %.2f at %d", r$peak_height, r$peak_centre) else
      "no peak >= T on this precursor"
    sprintf("%s %s, %s %d: %s", AL[[r$arm]], .LF_TPFN_HN[[r$head]], r$term, r$anchor,
            switch(r$reason,
              "peak on the end" = sprintf("TP -- peak %.2f centred at %d",
                                          r$peak_height, r$peak_centre),
              "below T"         = sprintf("FN, below T -- end score %.2f; %s", r$score, pk),
              "off-centre peak" = sprintf(paste("FN, inside a peak centred %d residues away",
                                                "(%.2f at %d); end score %.2f"),
                                          abs(r$peak_centre - r$anchor), r$peak_height,
                                          r$peak_centre, r$score),
              sprintf("FN, %s", r$reason)))
  }
  chip <- function(r) if (is.null(r) || !nrow(r) || is.na(r$call))
    '<span class="c na">&middot;</span>' else
    sprintf('<span class="c %s" title="%s">%s</span>', tolower(r$call),
            .lf_tpfn_esc(tip(r)), r$call)

  ## arm on one line, head on the next, each unbreakable: in a 38px column
  ## "ESM-C lin" otherwise wraps AT THE HYPHEN and reads as "ESM-" / "C lin".
  hdr <- paste0(vapply(arms, function(a)
    sprintf('<span>%s<br>%s</span>', a$short, .LF_TPFN_HS[[a$head]]), ""), collapse = "")

  cards <- data.table::rbindlist(lapply(genes, function(g) {
    e <- ec[ec$gene == g, ]
    if (!nrow(e)) return(NULL)
    pep <- unique(e[, c("set", "name", "start", "end")])
    pep <- pep[order(pep$start, -pep$end), ]
    body <- vapply(seq_len(nrow(pep)), function(i) {
      p <- pep[i, ]
      h <- sprintf('<div class="nsp">%s %d&ndash;%d%s</div>', .lf_tpfn_esc(p$name),
                   p$start, p$end,
                   if (grepl("large", p$set)) ' <span class="lp">large product</span>' else "")
      ends <- e[e$set == p$set & e$name == p$name & e$start == p$start & e$end == p$end, ]
      rows <- vapply(c("N", "C"), function(tm) {
        x <- ends[ends$term == tm, ]
        if (!nrow(x)) return("")
        sprintf('<div class="nsr"><span>%s %d</span>%s</div>', tm, x$anchor[1],
                paste0(vapply(arms, function(a) chip(x[x$arm == a$key, ][1, ]), ""),
                       collapse = ""))
      }, character(1))
      paste0(h, paste(rows, collapse = ""))
    }, character(1))
    data.table::data.table(gene = g, html = paste0(
      '<div class="ns"><div class="nsh"><span>never seen</span>', hdr, '</div>',
      paste(body, collapse = ""), '</div>'))
  }), fill = TRUE)

  n <- length(arms)
  style <- paste0('<style>',
    ## the card is display:block in make_gallery_index.R; flex + order:-1 lifts
    ## the calls to the TOP of the card without editing that shared renderer.
    '.card{display:flex;flex-direction:column}',
    '.ns{order:-1;margin:0 0 8px;padding-bottom:7px;border-bottom:1px solid var(--grid);font-size:11px}',
    sprintf('.nsh,.nsr{display:grid;grid-template-columns:1fr%s;align-items:center;column-gap:3px}',
            strrep(" 38px", n)),
    '.nsh{color:var(--ink2);font-size:10px;margin-bottom:3px}',
    '.nsh span+span{text-align:center;white-space:nowrap;line-height:1.15}',
    '.nsp{color:var(--ink);font-weight:600;margin-top:3px;overflow-wrap:anywhere}',
    '.nsr{color:var(--ink2);margin-top:1px}.lp{color:var(--ink2);font-weight:400}',
    '.c{display:inline-block;text-align:center;border-radius:3px;padding:0 4px;',
    'font-size:10.5px;font-weight:650;line-height:16px}',
    '.c.tp{background:#e3f1e5;color:#1B5E20}.c.fn{background:#f8e1e1;color:#8E1B1B}',
    '.c.na{background:transparent;color:var(--ink2);font-weight:400}',
    '</style>')
  list(cards = as.data.frame(cards), style = style, calls = ec)
}

#' The blurb paragraphs describing the chips
#'
#' Kept separate from [lf_tpfn_cards()] because it is the part that must stay
#' honest: the arms sit at DIFFERENT thresholds, so their TP counts are not a
#' ranking, and the false-positive cost is what says so.
#' @param arms,calls,genes,tol as elsewhere.
#' @export
lf_tpfn_blurb <- function(arms, calls, genes, tol = 4) {
  u <- unique(calls[, c("arm", "set", "accession", "term", "anchor", "call")])
  on_pg <- u[u$accession %in% unique(calls$accession[calls$gene %in% genes]), ]
  pep <- function(d) d[d$set == "never seen: peptides", ]
  cnt <- function(d, k) { x <- d[d$arm == k, ]
    sprintf("%d of %d", sum(x$call == "TP"), nrow(x)) }
  j <- function(f) paste(vapply(arms, f, ""), collapse = ", ")

  sprintf(paste0(
    '<p>Peptides whose gene no fold trained on are drawn above the sequence as a ',
    '<b>magenta &ldquo;never seen&rdquo; line</b> (hover for name and span). Each card lists them ',
    'at the top with every end called on the all-50 profile (the heavy squares), one column per ',
    'arm: %s &mdash; all at tolerance %s residues. <span class="c tp">TP</span> when a run of ',
    'anchors &ge; T has its maximum within tolerance of the end, <span class="c fn">FN</span> ',
    'otherwise; hover a call for the arm, the head and why.</p>\n',
    '<p><b>Which head.</b> Each arm shows the better of its two heads, chosen on auc_T &mdash; ',
    '(1 + TPR &minus; FPR)/2 &mdash; over never-seen peptides at that arm&rsquo;s own threshold, ',
    'so both heads are compared at the same operating point and the choice is charged for false ',
    'positives as well as true ones. Choosing on TP alone would mislead: a head that calls more ',
    'segments finds more ends and more false ones. Margins here are %s, so read the column heads ',
    'as a tie-break, not a verdict.</p>\n',
    '<p>On these %d pages: %s never-seen peptide ends TP. Over the whole never-seen peptide set: ',
    '%s.</p>\n',
    '<p><b>Do not read those counts as a ranking of the arms.</b> They sit at different ',
    'thresholds, so they are not the same test &mdash; a lower threshold buys ends with false ',
    'positives. False-positive segments over never-seen peptides, each arm at its own cut: %s. ',
    'The anchor-ranking metrics, which do not depend on a threshold, are the ones to compare ',
    'arms on.</p>'),
    j(function(a) sprintf('<b>%s</b> at T = %s (%s, auc_T %.3f)',
                          a$label, a$thr, .LF_TPFN_HN[[a$head]], a$auc)),
    tol,
    j(function(a) sprintf("%s %+.3f", a$label, a$margin)),
    length(genes),
    j(function(a) sprintf("%s <b>%s</b>", a$label, cnt(pep(on_pg), a$key))),
    j(function(a) sprintf("%s %s", a$label, cnt(pep(u), a$key))),
    j(function(a) sprintf("%s <b>%s</b> (FPR %.2f)", a$label,
                          format(a$fp, big.mark = ","), a$fpr)))
}

#' Build a gallery index with per-arm TP / FN chips
#'
#' Writes `index.html` into `dir`, with every never-seen end called TP or FN by
#' each arm at the top of that gene's card.
#'
#' Inputs it does NOT make for you, because they are an evaluation and not a
#' page: `end_calls_T<thr><arm>_tol<tol>.csv` and
#' `segment_confusion<arm>_tol<tol>.csv` under `headcmp`, one pair per arm. The
#' error names the command to run when one is missing.
#'
#' @param dir the gallery directory of `<GENE>.html` pages.
#' @param arms list of [lf_tpfn_arm()]; default [lf_tpfn_default_arms()].
#' @param genes gene symbols; default every `<GENE>.html` in `dir`.
#' @param tol segment tolerance in residues.
#' @param headcmp directory holding the evaluation tables.
#' @param title index title; default names the arms it drew.
#' @param blurb_template,footer optional HTML; `blurb_template` has `anchor`
#'   replaced by the generated section.
#' @param anchor the line of `blurb_template` to replace. It must match exactly
#'   one line.
#' @param notes optional per-gene notes, a data frame of `gene` and `note`.
#' @param workdir where the intermediate rds / html are written.
#' @param script path to `make_gallery_index.R`.
#' @return invisibly, the arms with their chosen heads and evidence.
#' @export
lf_gallery_tpfn_index <- function(dir,
                                  arms = lf_tpfn_default_arms(),
                                  genes = NULL,
                                  tol = 4,
                                  headcmp = "~/AF2_analysis/lf_headcmp",
                                  title = NULL,
                                  blurb_template = NULL,
                                  footer = NULL,
                                  anchor = "magenta",
                                  notes = NULL,
                                  workdir = tempdir(),
                                  script = NULL) {
  dir <- path.expand(dir)
  if (!dir.exists(dir)) stop("no such gallery directory: ", dir, call. = FALSE)
  if (is.null(genes)) {
    genes <- sub("\\.html$", "", list.files(dir, pattern = "\\.html$"))
    genes <- setdiff(genes, c("index", "index_before_tpfn"))
  }
  if (!length(genes)) stop("no <GENE>.html pages in ", dir, call. = FALSE)
  arms <- lapply(arms, lf_tpfn_pick_head, tol = tol, headcmp = headcmp)
  for (a in arms)
    message(sprintf("%-14s T=%-5s -> %-7s (auc_T %.3f, margin %+.3f; FP %s, FPR %.2f)",
                    a$label, a$thr, a$head, a$auc, a$margin,
                    format(a$fp, big.mark = ","), a$fpr))

  cc <- lf_tpfn_cards(arms, genes, tol = tol, headcmp = headcmp)
  if (is.null(title))
    title <- paste0("LigandFinder — peptide-end pages: ",
                    paste(vapply(arms, `[[`, "", "label"), collapse = ", "))

  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  cards_p <- file.path(workdir, "lf_tpfn_cards.rds")
  saveRDS(cc$cards, cards_p)
  args <- c("--dir", dir, "--title", title, "--card-html", cards_p, "--no-stock-desc")

  if (!is.null(blurb_template)) {
    tpl <- paste(readLines(path.expand(blurb_template), warn = FALSE), collapse = "\n")
    hit <- grep(anchor, strsplit(tpl, "\n")[[1]], value = TRUE)
    if (length(hit) != 1L)
      stop("the blurb anchor ", sQuote(anchor), " matches ", length(hit),
           " lines of ", basename(blurb_template), "; it must match exactly one.",
           call. = FALSE)
    blurb_p <- file.path(workdir, "lf_tpfn_blurb.html")
    writeLines(sub(hit, paste0(cc$style, lf_tpfn_blurb(arms, cc$calls, genes, tol)),
                   tpl, fixed = TRUE), blurb_p)
    args <- c(args, "--blurb", blurb_p)
  }
  if (!is.null(footer)) args <- c(args, "--footer", path.expand(footer))
  if (!is.null(notes)) {
    notes_p <- file.path(workdir, "lf_tpfn_notes.rds")
    saveRDS(notes, notes_p); args <- c(args, "--notes", notes_p)
  }
  if (is.null(script)) {
    script <- system.file("scripts", "make_gallery_index.R", package = "ligandFinder")
    if (!nzchar(script)) stop("make_gallery_index.R not found; pass `script`", call. = FALSE)
  }
  rc <- system2("/usr/local/bin/Rscript", shQuote(c(script, args)))
  if (rc != 0L) stop("make_gallery_index.R failed, rc = ", rc, call. = FALSE)
  message(sprintf("%d of %d cards carry never-seen calls", nrow(cc$cards), length(genes)))
  invisible(arms)
}

#' Per-gene notes naming each gene's training ends and the fold that held it out
#'
#' The card says whether a page is wholly out-of-sample or carries ends four of
#' five folds trained on -- those two are not interchangeable, and that belongs
#' on the card rather than in a readme.
#'
#' @param genes gene symbols.
#' @param cv5 directory holding `lf_pepend_cv5_knowns_{N,C}.csv`.
#' @return a data frame of `gene` and `note`, for `lf_gallery_tpfn_index(notes=)`.
#' @export
lf_tpfn_fold_notes <- function(genes, cv5 = "~/AF2_analysis") {
  kn <- data.table::rbindlist(lapply(c("N", "C"), function(t) {
    p <- file.path(path.expand(cv5), sprintf("lf_pepend_cv5_knowns_%s.csv", t))
    if (!file.exists(p)) stop("no ", basename(p), call. = FALSE)
    d <- data.table::fread(p)
    data.table::data.table(gene = sub("_w-?[0-9]+--?[0-9]+$", "", d$peps),
                           term = t, anchor = d$anchor, fold = d$fold)
  }))
  note <- vapply(genes, function(g) {
    k <- kn[kn$gene == g, ]
    k <- k[order(k$term, k$anchor), ]
    if (!nrow(k)) "no known end in any fold: all tracks out-of-sample" else
      paste0("training ends (4 of 5 folds saw each; held out by): ",
             paste(sprintf("%s %d (fold %d)", k$term, k$anchor, k$fold), collapse = ", "))
  }, character(1))
  data.frame(gene = genes, note = unname(note), stringsAsFactors = FALSE)
}
