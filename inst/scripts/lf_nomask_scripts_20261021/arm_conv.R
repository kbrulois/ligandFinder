## Per-arm end conventions for the evaluation scripts (Kevin, 2026-11-10).
##
## Arms trained on windows built with 10_7 --amid-g-anchor (INPUT_FLAGS in
## arms/<arm>.sh contains amid_g=TRUE) put a C end that is followed by G + a
## dibasic pair ON THE G, one residue past the mature (amidated) end. Every
## shared end table the evaluators read -- lf_pepend_cv5_knowns_<T>.csv,
## cv5_unseen_ends*.csv, ligand_list, uniprot peptides -- is in the MATURE
## convention, so such an arm must be scored at anchor + 1 for those ends or it
## is graded one residue away from what it was trained to call.
##
## The CV folds of the _padg inputs are pinned to the old ones (10_7
## split_by_mature), so lf_pepend_cv5_knowns_<T>.csv + arm_anchor() IS the padg
## arms' known set; no separate knowns file is needed.
##
##   source(file.path(Sys.getenv("LF_SCRIPTS", .lf_here), "arm_conv.R"))
##   kn[, anchor := arm_anchor(ARM, accession, "C", anchor)]
##
## `arm` may carry the leading underscore the scripts' ARM env var uses
## ("_curated_AAswap_padg"); "" / "curated" / "flatten" are the curated baseline.

## LF_SCRIPTS = the scripts directory: this file and arms/<arm>.sh. Unset, it is
## the directory of the script that sourced this one (`.lf_here`, as lf_paths.R
## finds it). No arms/ there would read every arm as the mature convention, so
## that stops here instead of scoring a _padg arm one residue off.
.lf_scripts <- path.expand(Sys.getenv("LF_SCRIPTS",
  if (exists(".lf_here", inherits = TRUE)) get(".lf_here", inherits = TRUE) else getwd()))
if (!dir.exists(file.path(.lf_scripts, "arms")))
  stop("arm_conv.R: no arms/ directory in ", .lf_scripts,
       "\n  source it from the scripts directory, or set LF_SCRIPTS to that directory", call. = FALSE)

arm_var <- function(arm, var) {
  a <- sub("^_", "", arm)
  f <- file.path(.lf_scripts, "arms", paste0(a, ".sh"))
  if (!nzchar(a) || !file.exists(f)) return("")
  l <- grep(paste0("^", var, "="), readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return("")
  trimws(gsub("[\"']", "", sub("#.*$", "", sub(paste0("^", var, "="), "", l[[1]]))))
}

## TRUE for an arm whose training windows anchor amidated C ends on the G
arm_amid_g <- function(arm) grepl("(^| )amid_g=TRUE( |$)", arm_var(arm, "INPUT_FLAGS"))

.arm_prec_seq <- NULL
arm_prec_seq <- function() {
  if (is.null(.arm_prec_seq)) {
    p <- data.table::fread(path.expand("~/AF2_analysis/lf_plm/precursor_seqs.tsv"))
    .arm_prec_seq <<- stats::setNames(p$seq, p$accession)
  }
  .arm_prec_seq
}

## MATURE-convention anchors -> the arm's convention. Identity for every arm
## without amid_g; for an amid_g arm, a C end followed by G + KK/KR/RK/RR moves
## onto the G. Vectorised; `term` may be a scalar.
arm_anchor <- function(arm, accession, term, anchor) {
  anchor <- as.integer(anchor)
  if (!arm_amid_g(arm)) return(anchor)
  s <- unname(arm_prec_seq()[as.character(accession)])
  nxt <- substr(s, anchor + 1L, anchor + 3L)
  shift <- rep_len(term, length(anchor)) == "C" & !is.na(s) & grepl("^G(KK|KR|RK|RR)", nxt)
  anchor + as.integer(shift)
}

## After arm_anchor(), two knowns can land on ONE residue: GLP-1(7-36)amide ends
## at 127 (G at 128) and GLP-1(7-37) at 128, and the padg inputs merged them into
## one end -- held out in fold 5 (127's fold), TRAINED in fold 4. Keeping both
## rows would grade fold 4's members on an end they trained on. For an arm with
## its own inputs (INPUT_SFX), the fold of every known is read from those inputs
## (the end is `val` in exactly one of the five) and duplicates collapse to it.
## Identity for every other arm. `kn` needs accession, anchor (already mapped),
## fold; extra columns are kept from the first row of each key.
.arm_folds <- list()
arm_fix_knowns <- function(arm, kn, term) {
  sfx <- arm_var(arm, "INPUT_SFX")
  if (!nzchar(sfx)) return(kn)
  key <- paste(sfx, term)
  if (is.null(.arm_folds[[key]])) {
    ## A CACHED three-column csv, because the fold inputs are 92 MB each and a
    ## render process would otherwise read ten of them (~1.3 GB, ~45 s) before
    ## its first gene. arm_knowns_cache.R writes it; delete it to rebuild.
    cf <- path.expand(sprintf("~/AF2_analysis/lf_pepend_cv5_knowns_%s%s.csv", term, sfx))
    .arm_folds[[key]] <<- if (file.exists(cf)) {
      data.table::fread(cf)[, .(accession, anchor = as.integer(anchor), fold_in = as.integer(fold))]
    } else data.table::rbindlist(lapply(1:5, function(k) {
      x <- readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_nn_input_cv5_f%d%s.rds", k, sfx)))$knowns
      x <- x[x$term == term & x$split == "val", ]
      data.table::data.table(accession = x$accession, anchor = as.integer(x$anchor), fold_in = k)
    }))
  }
  f <- .arm_folds[[key]]
  stopifnot(!anyDuplicated(f[, .(accession, anchor)]))
  kn <- data.table::as.data.table(kn)
  kn <- kn[!duplicated(kn[, .(accession, anchor)])]
  kn <- merge(kn, f, by = c("accession", "anchor"), all.x = TRUE, sort = FALSE)
  if (anyNA(kn$fold_in))
    stop(sum(is.na(kn$fold_in)), " mapped knowns are not a val end in any ", sfx, " input")
  kn[, fold := fold_in][, fold_in := NULL][]
}
