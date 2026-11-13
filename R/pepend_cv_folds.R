## ---- the cross-validation fold of each known end ---------------------------------
## Kevin, 2026-11-12: a cross-validated build's train/val mark says only whether
## its own fold holds an end out; 10_7b_plot_pepend_windows.R prints instead
## WHICH fold does, read from the sibling fold inputs.
## ------------------------------------------------------------------------------

#' Which cross-validation fold holds each known end out
#'
#' A `--n-folds K --fold k` build records only its own split, so a train/val
#' mark on its knowns says whether FOLD k holds an end out -- a fifth of the
#' answer. The K fold inputs share one file-name pattern (`..._f<k>_...rds`)
#' and hold the same knowns; reading the siblings' `knowns` finds, for every
#' end, the one fold whose validation set it is in.
#'
#' Costs one `readRDS()` per sibling: ~5 s and ~0.8 GB for a 90 MB fold input,
#' released before the next is read.
#'
#' @param in_path the fold input whose knowns are being drawn.
#' @param knowns its `$knowns` (`term`, `accession`, `anchor`, `split`).
#' @param params its `$params`; `n_folds` and `fold` locate the siblings.
#' @return integer vector, the held-out fold of each row of `knowns`; `NULL`
#'   (with a message) when the build is not cross-validated, a sibling is
#'   missing, or the folds do not partition the ends -- the caller then falls
#'   back to the build's own train/val split.
#' @export
lf_pepend_cv_folds <- function(in_path, knowns, params) {
  K <- params$n_folds; k0 <- params$fold
  if (!isTRUE(K > 0) || length(k0) != 1L || is.na(k0)) return(NULL)
  b   <- basename(in_path)
  tok <- sprintf("_f%d(?=[_.])", k0)
  if (lengths(regmatches(b, gregexpr(tok, b, perl = TRUE))) != 1L) {
    message("cv folds: cannot find _f", k0, " in ", b, "; showing this build's split")
    return(NULL)
  }
  sib <- file.path(dirname(in_path), vapply(seq_len(K), function(j)
    sub(tok, sprintf("_f%d", j), b, perl = TRUE), ""))
  if (!all(file.exists(sib))) {
    message("cv folds: missing ", paste(basename(sib[!file.exists(sib)]), collapse = ", "),
            "; showing this build's split")
    return(NULL)
  }
  key <- function(d) paste(d$term, d$accession, d$anchor)
  keys <- key(knowns)
  held <- vapply(seq_len(K), function(j) {
    if (j == k0) { kj <- knowns; pj <- params } else {
      x <- readRDS(sib[[j]]); kj <- x$knowns; pj <- x$params
      rm(x); invisible(gc())
    }
    ## a sibling from another build (other K, other fold) would partition
    ## nothing; neither would one whose knowns differ
    if (!identical(as.integer(pj$n_folds), as.integer(K)) ||
        !identical(as.integer(pj$fold), as.integer(j)) ||
        !setequal(key(kj), keys)) return(rep(NA, length(keys)))
    keys %in% key(kj)[kj$split == "val"]
  }, logical(length(keys)))
  held <- matrix(held, ncol = K)
  n <- rowSums(held)
  if (anyNA(n) || any(n != 1L)) {
    message(sprintf("cv folds: %d of %d ends are not held out in exactly one of %d folds; showing this build's split",
                    sum(is.na(n) | n != 1L), length(n), K))
    return(NULL)
  }
  as.integer(max.col(held, ties.method = "first"))
}
