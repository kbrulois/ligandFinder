## arm_var(): one setting of an arm, read out of arms/<arm>.sh (Kevin, 2026-11-10).
##
## The evaluators and the gallery scripts need a few of an arm's settings -- its
## training windows (INPUT_SFX), ARM_LABEL, PLM, DROP, INPUT_FLAGS -- without
## SOURCING the arm file, which would run its training-time guards. So the
## `VAR=value` line is read as text: quotes and a trailing comment dropped, the
## first such line wins, and a variable the file does not set reads "".
##
##   source(lf_src("inst/scripts/modeling/arm_vars.R"))     # after lf_paths.R
##   sfx <- arm_var(ARM, "INPUT_SFX")
##
## `arm` may carry the leading underscore the scripts' ARM env var uses
## ("_esmonly36"); "" / "curated" / "flatten" are the curated baseline, which
## has no arm file, so every variable of it reads "".

## The arm files are inst/scripts/modeling/arms/<arm>.sh in the source tree that
## lf_paths.R found ($LF_SRC, or the checkout the calling script sits in). No
## arms/ there would read every arm as the curated baseline, so that stops here.
.lf_arms <- file.path(lf_src_root(), "inst", "scripts", "modeling", "arms")
if (!dir.exists(.lf_arms))
  stop("arm_vars.R: no arms/ directory at ", .lf_arms,
       "\n  source lf_paths.R first, from a checkout or with LF_SRC set", call. = FALSE)

arm_var <- function(arm, var) {
  a <- sub("^_", "", arm)
  f <- file.path(.lf_arms, paste0(a, ".sh"))
  if (!nzchar(a) || !file.exists(f)) return("")
  l <- grep(paste0("^", var, "="), readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return("")
  trimws(gsub("[\"']", "", sub("#.*$", "", sub(paste0("^", var, "="), "", l[[1]]))))
}
