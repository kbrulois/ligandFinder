## =============================================================================
## Bundle entries for genes plot_bundle.rds lacks (Kevin, 2026-10-28).
##
## plot_bundle.rds holds the genes the OLD dibasic-window model scored: 10_4's
## default gene list is nn_input_comb's, and 9.2 made a window only at a
## KK/KR/RK/RR after the signal peptide. A precursor without one -- CD24, DEFB1,
## PRAP1, TCTA and ~430 others -- has no entry, so final_pages.R cannot draw it,
## although every peptide-end input (residue cache, scans, candidates) covers it.
##
## This builds the_input / pep_input entries for the named genes with
## plot_all.R's and 10_4_plot_proteins.R's own code, run up to (not including)
## 10_4's render, and writes them to their own rds. final_pages.R appends that
## file when LF_EXTRA_BUNDLE points at it. plot_bundle.rds is never read or
## written. Run from the repo root:
##
##   LF_PLOT_DIR=<scratch dir> /usr/local/bin/Rscript inst/scripts/plot_bundle_extra.R CD24 DEFB1 PRAP1 TCTA
##   (writes LF_EXTRA_OUT, default ~/AF2_analysis/plot_bundle_extra.rds)
##
## It loads what plot_all.R loads (secretome.rds and secretome_aa_1.rds, ~3 GB
## on disk): minutes, and well over 10 GB of RAM -- not beside page renders.
## LF_PLOT_DIR: 10_4 creates its plot_dir even though nothing is drawn.
## =============================================================================
if (!file.exists("inst/scripts/10_4_plot_proteins.R"))
  stop("run this from the repo root (", getwd(), " is not it)", call. = FALSE)
.extra_genes <- commandArgs(trailingOnly = TRUE)
if (!length(.extra_genes)) stop("name the genes to build entries for", call. = FALSE)
.extra_out <- path.expand(Sys.getenv("LF_EXTRA_OUT", "~/AF2_analysis/plot_bundle_extra.rds"))

## a script's text up to (not including) the first line matching `pat`
.upto <- function(f, pat) {
  s <- readLines(f); i <- grep(pat, s)[1]
  if (is.na(i)) stop("no line matching ", pat, " in ", f, call. = FALSE)
  s[seq_len(i - 1L)]
}
## plot_all.R without its last step, the source() of 10_4 (which renders). Its
## own gene-argument handling makes plot_genes these genes and write_bundle FALSE.
eval(parse(text = .upto("inst/scripts/plot_all.R", '^source\\("inst/scripts/10_4_plot_proteins.R"\\)')),
     envir = globalenv())
stopifnot(identical(plot_genes, .extra_genes), isFALSE(write_bundle))
## 10_4 up to its render section: the_input / pep_input, built and paired as
## for every bundled gene
eval(parse(text = .upto("inst/scripts/10_4_plot_proteins.R", "^## ---- render")), envir = globalenv())

.got <- vapply(the_input, function(x) as.character(x$gene)[1], "")
if (length(setdiff(.extra_genes, .got)))
  warning("no entry built for: ", paste(setdiff(.extra_genes, .got), collapse = ", "), call. = FALSE)
## (final_pages.R refuses an extra entry for a gene the bundle already has)
saveRDS(list(the_input = the_input, pep_input = pep_input), .extra_out)
message(sprintf("wrote %d entr%s (%s) -> %s", length(.got), if (length(.got) == 1) "y" else "ies",
                paste(.got, collapse = ", "), .extra_out))
