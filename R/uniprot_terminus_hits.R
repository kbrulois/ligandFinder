## Flag windows whose anchor residue sits at a UniProt-annotated peptide terminus.
##
## Pulled out of 10_1dcnn_new6.R so the benchmark scripts score "top hits" with
## exactly the rule the main pipeline uses.

#' Annotate windows that anchor a UniProt-peptide terminus
#'
#' Windows are built as anchor + `[-5, 30]` (N) / `[-30, 5]` (C), so the putative
#' peptide boundary sits at a fixed anchor residue: `wN + 5` for N windows,
#' `wC - 5` for C windows. A window "hits" a UniProt peptide only if that
#' peptide's matching terminus -- start (N-terminus) for an N window, end
#' (C-terminus) for a C window -- lands at the window's anchor residue
#' (+/- `anchor_tol`), i.e. at the correct position, not merely somewhere inside
#' the window span.
#'
#' @param windows a data frame with `peps` (`<gene>_w<N>-<C>`), `target`
#'   (`N`, `C`, `loop_N`, `loop_C`) and `gene`. A `model` column (`N`/`C`) is
#'   used if present and derived from `target` otherwise.
#' @param uniprot_peps data frame with `gene`, `start`, `end` -- one row per
#'   annotated peptide, as in `~/Desktop/Peptides/uniprot_peptides.csv`.
#' @param anchor_tol residues of slack either side of the anchor.
#' @return `windows` with `model`, `wN`, `wC`, `anchor_res` and the logical
#'   `pep_terminus_hit` added.
#' @export
lf_uniprot_terminus_hits <- function(windows, uniprot_peps, anchor_tol = 2L) {
  if (!"model" %in% names(windows))
    windows$model <- sub("^loop_", "", as.character(windows$target))

  wm <- stringr::str_match(windows$peps, "_w(\\d+)-(\\d+)$")
  windows$wN         <- as.integer(wm[, 2])
  windows$wC         <- as.integer(wm[, 3])
  windows$anchor_res <- ifelse(windows$model == "N",
                               windows$wN + 5L,     # expected peptide N-terminus
                               windows$wC - 5L)     # expected peptide C-terminus

  starts_by_gene <- split(as.integer(uniprot_peps$start), uniprot_peps$gene)  # peptide N-ends
  ends_by_gene   <- split(as.integer(uniprot_peps$end),   uniprot_peps$gene)  # peptide C-ends

  anchor_hits <- function(gene, model, anchor, tol) {
    ter <- if (model == "N") starts_by_gene[[gene]] else ends_by_gene[[gene]]
    if (is.null(ter) || is.na(anchor)) return(FALSE)
    any(abs(ter - anchor) <= tol)
  }

  windows$pep_terminus_hit <- FALSE
  idx <- which(windows$gene %in% names(starts_by_gene))   # only genes with uniprot peptides
  if (length(idx) > 0) {
    windows$pep_terminus_hit[idx] <- mapply(
      anchor_hits,
      windows$gene[idx], windows$model[idx], windows$anchor_res[idx],
      MoreArgs = list(tol = anchor_tol))
  }
  windows
}
