## arm_var(): one setting of an arm, read out of arms/<arm>.sh (Kevin, 2026-11-10).
##
## The evaluators and the gallery scripts need a few of an arm's settings -- its
## training windows (INPUT_SFX), ARM_LABEL, PLM, DROP, INPUT_FLAGS -- without
## SOURCING the arm file, which would run its training-time guards. So the
## `VAR=value` line is read as text: quotes and a trailing comment dropped, the
## first such line wins, and a variable the file does not set reads "".
##
##   source(file.path(Sys.getenv("LF_SCRIPTS", .lf_here), "arm_vars.R"))
##   sfx <- arm_var(ARM, "INPUT_SFX")
##
## `arm` may carry the leading underscore the scripts' ARM env var uses
## ("_esmonly36"); "" / "curated" / "flatten" are the curated baseline, which
## has no arm file, so every variable of it reads "".

## LF_SCRIPTS = the scripts directory: this file and arms/<arm>.sh. Unset, it is
## the directory of the script that sourced this one (`.lf_here`, as lf_paths.R
## finds it). No arms/ there would read every arm as the curated baseline, so
## that stops here instead.
.lf_scripts <- path.expand(Sys.getenv("LF_SCRIPTS",
  if (exists(".lf_here", inherits = TRUE)) get(".lf_here", inherits = TRUE) else getwd()))
if (!dir.exists(file.path(.lf_scripts, "arms")))
  stop("arm_vars.R: no arms/ directory in ", .lf_scripts,
       "\n  source it from the scripts directory, or set LF_SCRIPTS to that directory", call. = FALSE)

arm_var <- function(arm, var) {
  a <- sub("^_", "", arm)
  f <- file.path(.lf_scripts, "arms", paste0(a, ".sh"))
  if (!nzchar(a) || !file.exists(f)) return("")
  l <- grep(paste0("^", var, "="), readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return("")
  trimws(gsub("[\"']", "", sub("#.*$", "", sub(paste0("^", var, "="), "", l[[1]]))))
}
