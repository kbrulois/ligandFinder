## SUPERSEDED -- this is now a shim.
##
## The TP/FN gallery index moved INTO THE PACKAGE on 2026-10-26 (Kevin: "wire
## this index.html generation code into the package so i don't have to keep
## asking you to generate it that way"). The code lives in:
##
##   R/gallery_tpfn.R                  lf_gallery_tpfn_index() and friends
##   inst/scripts/make_tpfn_index.R    its command-line face
##
## An arm is now DATA -- lf_tpfn_arm() entries in lf_tpfn_default_arms() -- so
## adding one is a single line there rather than another copy of this script.
## Add it there, not here.
##
## This file stays only so the old path keeps working. It forwards its arguments
## unchanged, supplying this directory's blurb and footer when none are given.
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
S <- path.expand("~/AF2_analysis/lf_nomask_scripts_20261021")
a <- commandArgs(trailingOnly = TRUE)
if (!"--blurb" %in% a) a <- c(a, "--blurb", file.path(S, "final_blurb.html"))
if (!"--footer" %in% a) a <- c(a, "--footer", file.path(S, "final_footer.html"))
if (!"--dir" %in% a) a <- c(a, "--dir", path.expand("~/AF2_analysis/pepend_pages_arms"))
message("index_tpfn_arms.R is a shim; running inst/scripts/make_tpfn_index.R")
q(status = system2("/usr/local/bin/Rscript",
                   shQuote(c(lf_src("inst/scripts/make_tpfn_index.R"), a))))
