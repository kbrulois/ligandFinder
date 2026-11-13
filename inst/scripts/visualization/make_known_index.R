#!/usr/bin/env Rscript
## ---- a gallery index of the known peptides, with a score chip per arm ---------
## The thin command-line face of lf_gallery_known_index(). See R/gallery_known.R.
##
##   Rscript inst/scripts/make_known_index.R --dir ~/AF2_analysis/pepend_pages_allsec_arms2
##
## Options
##   --dir       the gallery of <GENE>.html pages
##   --out       the file to write, inside --dir unless absolute  index_known.html
##   --base      put before each link, for an index not beside its pages
##   --training  how a training end is scored: heldout, the mean of the fold that
##               held it out, or all50, the page's heavy line          heldout
##   --title     index title
## ------------------------------------------------------------------------------
suppressMessages({ library(data.table) })
.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
## this checkout's definitions, not an installed copy's (see make_tpfn_index.R)
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT  <- normalizePath(file.path(dirname(.self), "..", ".."))
source(file.path(ROOT, "R", "gallery_tpfn.R"))
source(file.path(ROOT, "R", "gallery_known.R"))
a <- list(dir = .opt("--dir", "~/AF2_analysis/pepend_pages_allsec_arms2"),
          out = .opt("--out", "index_known.html"), training = .opt("--training", "heldout"),
          base = .opt("--base", ""),
          script = file.path(ROOT, "inst", "scripts", "make_gallery_index.R"))
ttl <- .opt("--title", NA_character_)
if (!is.na(ttl)) a$title <- ttl
do.call(lf_gallery_known_index, a)
