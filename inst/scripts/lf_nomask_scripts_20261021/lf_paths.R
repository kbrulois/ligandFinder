## ---------------------------------------------------------------------------
## Where the ligandFinder package's own code, data, scripts and python live --
## resolved ONCE, here, instead of once per script.
##
## These scripts read R/ plotting code, inst/extdata data files, inst/scripts
## helpers and inst/python from a ligandFinder SOURCE TREE, not from the
## installed package. That is deliberate, not laziness: the install lags badly
## (built 2026-09-10 as this was written) and its inst/python carries neither
## lf_winxgb nor lf_resid at all, so system.file() would break every window
## scan and would silently render pages from six-week-old plotting code. The
## three data files the pages read through system.file() -- data/id_mapping.rds,
## extdata/species_dat.rds -- are byte-identical in the install and are left
## alone; everything that used to come from a throwaway git worktree now comes
## from here.
##
## The source tree is found as:
##   1. $LF_SRC, if set (validated, so a typo fails loudly)
##   2. the checkout this script directory sits in -- these scripts are
##      committed at inst/scripts/lf_nomask_scripts_20261021/, so a copy that
##      travels with the repository needs no configuration whatsoever
## A copy of this directory kept OUTSIDE a checkout (the working copy in
## ~/AF2_analysis, say) therefore needs LF_SRC, and says so when it does not
## have it:
##   export LF_SRC=~/R_projects/ligandFinder
##
## Use it as lf_src("R/pepend_tracks.R"), lf_src("inst/extdata/ligand_list.rds"),
## lf_src("inst/python") -- paths relative to the top of the source tree.
##
## Sourcing this file is a two-line preamble, the same in every script:
##   .lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
##   source(file.path(.lf_here, "lf_paths.R"))
## `.lf_here` is how this file knows where it was sourced from; without it the
## search starts at getwd() instead.
## ---------------------------------------------------------------------------

.lf_is_src <- function(d) {
  f <- file.path(d, "DESCRIPTION")
  isTRUE(file.exists(f)) && dir.exists(file.path(d, "R")) &&
    any(grepl("^Package:[[:space:]]*ligandFinder[[:space:]]*$",
              readLines(f, warn = FALSE)))
}

lf_src_root <- local({
  root <- NULL
  function() {
    if (!is.null(root)) return(root)
    e <- Sys.getenv("LF_SRC", "")
    if (nzchar(e)) {
      d <- normalizePath(path.expand(e), mustWork = FALSE)
      if (!.lf_is_src(d))
        stop("LF_SRC is not a ligandFinder source tree (no DESCRIPTION / R/): ", e,
             call. = FALSE)
      root <<- d
      return(root)
    }
    from <- if (exists(".lf_here", inherits = TRUE)) get(".lf_here", inherits = TRUE) else getwd()
    d <- normalizePath(from, mustWork = FALSE)
    repeat {
      if (.lf_is_src(d)) { root <<- d; return(root) }
      up <- dirname(d)
      if (identical(up, d)) break
      d <- up
    }
    stop("no ligandFinder source tree at or above ", from,
         "\n  these scripts read R/, inst/extdata, inst/scripts and inst/python from a checkout.",
         "\n  Run them from inst/scripts/lf_nomask_scripts_20261021/ in a checkout, or set LF_SRC:",
         "\n      export LF_SRC=~/R_projects/ligandFinder", call. = FALSE)
  }
})

## one file or directory inside the source tree, checked -- a missing path here
## means the tree is not the one the script expects, which is worth stopping for
lf_src <- function(...) {
  p <- file.path(lf_src_root(), ...)
  if (!file.exists(p))
    stop("not in the ligandFinder source tree ", lf_src_root(), ": ",
         file.path(...), call. = FALSE)
  p
}
