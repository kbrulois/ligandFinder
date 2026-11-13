#!/usr/bin/env Rscript
## ---- a gallery's candidate peptides as one clickable scatter ------------------
## The thin command-line face of lf_gallery_candidate_plot(). See
## R/gallery_candidate_plot.R; the drawing is inst/templates/candidate_plot.html.
##
##   Rscript inst/scripts/make_candidate_plot.R --dir ~/AF2_analysis/pepend_pages_allsec_arms2
##
## Options
##   --dir   the gallery of <GENE>.html pages
##   --out   the file to write, inside --dir unless absolute  index_candidates.html
##   --base  put before each link, for a page not beside the gene pages
##   --known a known-peptide index to carry inside, opened at #known
##   --key   the candidate set                                xgb_pm_T0.5_tol4_arms3
##   --thr   the end-score cut of the end score type colouring                  0.6
## ------------------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
## this checkout's definitions, not an installed copy's (see make_tpfn_index.R)
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT  <- normalizePath(file.path(dirname(.self), "..", ".."))
source(file.path(ROOT, "R", "pepend_windows.R"))      # lf_pepend_motif_vec
source(file.path(ROOT, "R", "gallery_candidate_plot.R"))
lf_gallery_candidate_plot(dir = .opt("--dir", "~/AF2_analysis/pepend_pages_allsec_arms2"),
                          template = file.path(ROOT, "inst", "templates", "candidate_plot.html"),
                          key = .opt("--key", "xgb_pm_T0.5_tol4_arms3"),
                          thr = as.numeric(.opt("--thr", "0.6")),
                          out = .opt("--out", "index_candidates.html"),
                          base = .opt("--base", ""), known = .opt("--known", NULL))
