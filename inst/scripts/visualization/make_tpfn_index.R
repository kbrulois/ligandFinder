#!/usr/bin/env Rscript
## ---- a gallery index with per-arm TP / FN chips -------------------------------
## The thin command-line face of lf_gallery_tpfn_index(). All the work is in
## R/gallery_tpfn.R, so a new arm is one more entry in --arms and not another
## copy of a script.
##
##   Rscript inst/scripts/make_tpfn_index.R --dir ~/AF2_analysis/pepend_pages_arms
##
## Options
##   --dir       the gallery of <GENE>.html pages        (required in practice)
##   --arms      key:label:suffix:threshold, comma-separated; default is the
##               three arms lf_tpfn_default_arms() describes. Use "-" for the
##               curated baseline's empty suffix, e.g.
##                 cur:curated:-:0.2,esm:ESM-C only:_esmonly36:0.65
##   --tol       segment tolerance, residues             4
##   --headcmp   where the evaluation tables live        ~/AF2_analysis/lf_headcmp
##   --title     index title                             names the arms it drew
##   --blurb     blurb template; its --anchor line is replaced by the generated
##               section (style + the chip paragraphs)
##   --anchor    the template line to replace            magenta
##   --footer    footer html
##   --no-notes  skip the per-gene training-end notes
##
## It does NOT compute the calls: end_calls_T<thr><arm>_tol<tol>.csv and
## segment_confusion<arm>_tol<tol>.csv must exist under --headcmp, one pair per
## arm. A missing one is an error naming the command that makes it.
## ------------------------------------------------------------------------------
suppressMessages({ library(data.table) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }

## Prefer THIS checkout's definitions over an installed copy: lf_dcnn_path()'s
## trap applies here too -- system.file() finds the installed package, which
## silently shadows the tree this script came from.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT  <- normalizePath(file.path(dirname(.self), "..", ".."))
if (file.exists(file.path(ROOT, "R", "gallery_tpfn.R"))) {
  source(file.path(ROOT, "R", "gallery_tpfn.R"))
  SCRIPT <- file.path(ROOT, "inst", "scripts", "make_gallery_index.R")
} else {
  suppressMessages(library(ligandFinder)); SCRIPT <- NULL
}

dir_p   <- .opt("--dir", "~/AF2_analysis/pepend_pages_arms")
tol     <- as.integer(.opt("--tol", "4"))
headcmp <- .opt("--headcmp", "~/AF2_analysis/lf_headcmp")
title   <- .opt("--title", NA_character_)
blurb   <- .opt("--blurb", NA_character_)
anchor  <- .opt("--anchor", "magenta")
footer  <- .opt("--footer", NA_character_)
arms_s  <- .opt("--arms", NA_character_)

arms <- if (is.na(arms_s)) lf_tpfn_default_arms() else
  lapply(strsplit(arms_s, ",", fixed = TRUE)[[1]], function(s) {
    p <- strsplit(s, ":", fixed = TRUE)[[1]]
    if (length(p) != 4L)
      stop("--arms entry must be key:label:suffix:threshold, got ", sQuote(s), call. = FALSE)
    lf_tpfn_arm(p[1], p[2], if (p[3] == "-") "" else p[3], p[4], short = p[2])
  })

genes <- sub("\\.html$", "", list.files(path.expand(dir_p), pattern = "\\.html$"))
genes <- setdiff(genes, c("index", "index_before_tpfn"))
notes <- if ("--no-notes" %in% .args) NULL else lf_tpfn_fold_notes(genes)

lf_gallery_tpfn_index(
  dir = dir_p, arms = arms, genes = genes, tol = tol, headcmp = headcmp,
  title = if (is.na(title)) NULL else title,
  blurb_template = if (is.na(blurb)) NULL else blurb,
  anchor = anchor,
  footer = if (is.na(footer)) NULL else footer,
  notes = notes,
  script = SCRIPT)
