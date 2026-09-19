## Protein-language-model channels for the window models.
##
## `python -m lf_plm reduce` writes a long table (accession, index, AA,
## plm_01..plm_k) on [0, 1] -- see inst/python/lf_plm. These helpers slot those
## columns in beside the 26 hand-built channels of an already-built `nn_input`
## (the per-window `data` tibbles 9.2_add_contact_data.R produces), so a model
## can be benchmarked with and without them on identical windows.

#' Read a parquet file without loading R's arrow package
#'
#' Once R's `arrow` has loaded libarrow, python's `pyarrow` in the same process
#' (which `lf_dcnn_export()` uses for `meta.parquet`) refuses to initialise:
#' "Attempted to register factory for scheme 'file' but that scheme is already
#' registered". So when reticulate is up, read through pandas; fall back to R
#' arrow only when python is not in play.
#' @keywords internal
lf_read_parquet <- function(path) {
  path <- path.expand(path)
  if (requireNamespace("reticulate", quietly = TRUE) &&
      (reticulate::py_available(initialize = FALSE) || !requireNamespace("arrow", quietly = TRUE))) {
    pd <- reticulate::import("pandas", convert = FALSE)
    return(tibble::as_tibble(reticulate::py_to_r(pd$read_parquet(path))))
  }
  arrow::read_parquet(path)
}

#' Read the reduced PLM residue table
#'
#' @param path the parquet `lf_plm reduce` wrote.
#' @return a list: `by_acc`, one numeric matrix per accession with the residue
#'   `index` as rownames and the AA letters as `attr(, "AA")`; `channels`, the
#'   channel names in column order.
#' @export
lf_plm_read <- function(path) {
  tbl <- lf_read_parquet(path)
  cols <- grep("^plm_", names(tbl), value = TRUE)
  if (!length(cols)) stop("no plm_* columns in ", path)
  idx <- split(seq_len(nrow(tbl)), tbl$accession)
  by_acc <- lapply(idx, function(i) {
    m <- as.matrix(tbl[i, cols])
    rownames(m) <- tbl$index[i]
    attr(m, "AA") <- stats::setNames(tbl$AA[i], tbl$index[i])
    m
  })
  list(by_acc = by_acc, channels = cols)
}

#' Append PLM channels to every window of an nn_input
#'
#' Each window's `meta_data` carries the residue `index` (NA at padding rows)
#' and `AA`; the block for the window is the PLM rows at those indices, 0 at
#' padding -- the convention every other channel follows, with the `padding`
#' channel flagging those rows. Residue letters are checked against the table:
#' a mismatch means the window was cut from a different sequence version than
#' the embedding, and is an error rather than a silent misalignment.
#'
#' @param nn_input the `list[term][split]` window sets.
#' @param plm from [lf_plm_read()].
#' @param sequences data frame with `gene` and `accession` (what the embedding
#'   ran on); a gene with several accessions resolves to the one whose residues
#'   match the window.
#' @param channels the current channel names (`all_params3`).
#' @return `list(nn_input = <with the columns appended>, channels = <extended>)`.
#' @export
lf_plm_attach <- function(nn_input, plm, sequences, channels) {
  acc_of <- split(as.character(sequences$accession), as.character(sequences$gene))
  k <- length(plm$channels)
  zero_block <- matrix(0, nrow = 0, ncol = k)

  block_for <- function(gene, meta) {
    n <- nrow(meta)
    out <- matrix(0, nrow = n, ncol = k, dimnames = list(NULL, plm$channels))
    real <- which(!is.na(meta$index))
    if (!length(real)) return(out)
    accs <- acc_of[[gene]]
    if (is.null(accs)) stop("no accession for gene ", gene, " in `sequences`")
    idx <- as.character(meta$index[real])
    for (acc in accs) {
      m <- plm$by_acc[[acc]]
      if (is.null(m)) next
      have <- idx %in% rownames(m)
      if (!all(have)) next
      aa_ok <- attr(m, "AA")[idx] == as.character(meta$AA[real])
      if (!all(aa_ok, na.rm = TRUE)) next
      out[real, ] <- m[idx, , drop = FALSE]
      return(out)
    }
    stop(sprintf("gene %s: no accession (%s) matches the window residues %s..%s",
                 gene, paste(accs, collapse = "/"), idx[[1]], idx[[length(idx)]]))
  }

  n_done <- 0L
  for (term in names(nn_input)) {
    for (split in names(nn_input[[term]])) {
      s <- nn_input[[term]][[split]]
      if (!nrow(s)) next
      s$data <- Map(function(d, gene, meta) {
        b <- block_for(gene, meta)
        d <- as.data.frame(d)
        stopifnot(nrow(d) == nrow(b))
        cbind(d, as.data.frame(b))
      }, s$data, s$gene, s$meta_data)
      nn_input[[term]][[split]] <- s
      n_done <- n_done + nrow(s)
    }
  }
  message(sprintf("lf_plm_attach: %d windows x %d PLM channels appended", n_done, k))
  list(nn_input = nn_input, channels = c(channels, plm$channels))
}
