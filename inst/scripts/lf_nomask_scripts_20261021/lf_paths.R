## ---------------------------------------------------------------------------
## A forwarder. The real thing is one directory up, at inst/scripts/lf_paths.R,
## because a second set of scripts (10_9z_*) wants it too and the resolution
## logic should exist once.
##
## The 14 scripts here are untouched: they still open with
##   .lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
##   source(file.path(.lf_here, "lf_paths.R"))
## and land on this file, which hands straight over.
##
## Two ways to find the real file, matching the two ways lf_paths.R itself
## resolves a source tree -- this directory kept inside a checkout, or copied
## out of one (the working copy in ~/AF2_analysis, say) with LF_SRC set. The
## second case is why this cannot be a bare source("../lf_paths.R"): a copied
## directory has no parent inst/scripts to walk into.
## ---------------------------------------------------------------------------

local({
  here <- if (exists(".lf_here", inherits = TRUE)) get(".lf_here", inherits = TRUE) else getwd()
  cand <- file.path(here, "..", "lf_paths.R")
  if (!file.exists(cand)) {
    e <- Sys.getenv("LF_SRC", "")
    cand <- if (nzchar(e)) file.path(path.expand(e), "inst", "scripts", "lf_paths.R") else ""
  }
  if (!nzchar(cand) || !file.exists(cand))
    stop("cannot find inst/scripts/lf_paths.R from ", here,
         "\n  this file only forwards to it. Run these scripts from a checkout,",
         "\n  or point LF_SRC at one:    export LF_SRC=~/R_projects/ligandFinder",
         call. = FALSE)
  ## into the global environment, as a direct source() of the real file would:
  ## the callers go on to use lf_src() at top level
  source(normalizePath(cand), local = FALSE)
})
