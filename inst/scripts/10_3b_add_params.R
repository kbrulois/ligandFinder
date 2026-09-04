## =============================================================================
## Add per-position parameters + per-window summaries to the UMAP CSV
## AFTER THE FACT -- i.e. WITHOUT re-running UMAP (the slow part).
##
## Reads the coords from a prior 10_3 run (out of `in_csv`, if present), rebuilds
## the labels + per-window summaries + raw per-position values straight from
## `nn_input` (vectorised, fast), joins on peps+terminus, and writes `out_csv`.
##
## If `in_csv` doesn't exist yet, it still writes everything except UMAP1/UMAP2,
## so you can get the parameter table immediately and merge the coords in later.
##
## Needs `nn_input`, `all_params3`, `seq_len` in the session (from 9.2).
## =============================================================================

library(tidyverse)

## ---- options ---------------------------------------------------------------
in_csv        <- path.expand("~/AF2_analysis/peptide_umap.csv")       # prior UMAP output (coords)
out_csv       <- path.expand("~/AF2_analysis/peptide_umap_full.csv")  # enriched output
load_csv      <- path.expand("~/AF2_analysis/peptide_pca_loadings.csv")  # PCA loadings (feature x PC)
feature_set   <- all_params3        # 36 per-residue features (drops end_type_ch)
top_pct       <- 0.80               # per-window "top" summary = mean of residues >= this pctile
excl_padding  <- TRUE               # summaries over real residues only (padding excluded)
add_positions <- TRUE               # append the raw per-position feature values
add_pca       <- TRUE               # add PC1..PC<pca_dims> scores + write loadings
pca_dims      <- 20                 # number of principal components
subset_peps   <- NULL               # NULL = all windows; or a character vector of peps to keep
knowns_only   <- FALSE              # TRUE = only known peptides (fast sanity run)

stopifnot(exists("nn_input"), exists("all_params3"))
if (!exists("seq_len")) seq_len <- nrow(nn_input[[1]]$all$data[[1]])
col_or_na  <- function(d, nm) if (nm %in% names(d)) d[[nm]] else NA
have_mstat <- requireNamespace("matrixStats", quietly = TRUE)   # fast row quantiles
top_tag    <- paste0("top", round(top_pct * 100))
if (!have_mstat) warning("matrixStats not installed -- the top", round(top_pct*100),
  " percentile falls back to a slow per-row apply(). install.packages('matrixStats') to make this fast.",
  immediate. = TRUE)
message("windows per terminus: ",
        paste(sprintf("%s=%d", names(nn_input),
                      sapply(nn_input, function(x) length(x$all$data))), collapse = "  "))

## ---- feature columns + wide names ------------------------------------------
first    <- nn_input[[1]]$all$data[[1]]
ch_names <- colnames(first); if (is.null(ch_names)) ch_names <- names(first)
sel      <- match(feature_set, ch_names)
stopifnot(!anyNA(sel))
P        <- length(sel)
pad_j    <- match("padding", feature_set)                 # NA if user dropped padding
resid_j  <- setdiff(seq_len(P), pad_j)
resid_nm <- feature_set[resid_j]
wide_names <- as.vector(outer(feature_set, 1:seq_len,     # pos-major, feature inner
                              function(f, p) paste0("pos", sprintf("%02d", p), "_", f)))

row_quant <- function(M) {
  if (have_mstat) matrixStats::rowQuantiles(M, probs = top_pct, na.rm = TRUE)
  else apply(M, 1, stats::quantile, probs = top_pct, names = FALSE, na.rm = TRUE)
}

## ---- build one terminus: meta + summaries + wide positions (vectorised) ----
build_term <- function(term) {
  d <- nn_input[[term]]$all
  ## optional row subsetting (fast validation runs) ---------------------------
  if (isTRUE(knowns_only))        d <- d[which(d$known == 1), ]
  if (!is.null(subset_peps))      d <- d[which(d$peps %in% subset_peps), ]
  n <- length(d$data)
  message(sprintf("  [%s] building %d windows ...", term, n)); flush.console()
  if (n == 0) return(NULL)
  ## (n, seq_len, P) array of the 36 selected features
  A <- aperm(array(unlist(lapply(d$data, function(m) as.matrix(m[, sel, drop = FALSE]))),
                   dim = c(seq_len, P, n)), c(3, 1, 2))

  pad_mask <- if (is.na(pad_j)) matrix(FALSE, n, seq_len) else A[, , pad_j] == 1
  n_resid  <- seq_len - rowSums(pad_mask); n_resid[n_resid == 0] <- NA_real_

  means <- matrix(NA_real_, n, length(resid_j), dimnames = list(NULL, paste0("mean_", resid_nm)))
  tops  <- matrix(NA_real_, n, length(resid_j), dimnames = list(NULL, paste0(top_tag, "_", resid_nm)))
  for (k in seq_along(resid_j)) {
    M  <- A[, , resid_j[k]]                       # n x seq_len
    Mr <- M; Mr[pad_mask] <- NA                   # residues only
    means[, k] <- rowSums(Mr, na.rm = TRUE) / n_resid
    thr        <- row_quant(Mr)                   # per-window top_pct percentile
    tmask      <- (M >= thr) & !pad_mask
    tops[, k]  <- rowSums(M * tmask) / rowSums(tmask)
  }

  meta <- tibble(
    peps     = col_or_na(d, "peps"),
    gene     = as.character(col_or_na(d, "gene")),
    terminus = term,
    target   = as.character(col_or_na(d, "target")),
    win_type = as.character(col_or_na(d, "win_type")),
    end_type = as.character(col_or_na(d, "end_type")),
    known    = col_or_na(d, "known"),
    set = dplyr::case_when(
      d$peps %in% nn_input[[term]]$train$peps ~ "train",
      d$peps %in% nn_input[[term]]$val$peps   ~ "val",
      TRUE                                    ~ "none"),
    pad_frac = rowSums(pad_mask) / seq_len)

  res <- bind_cols(meta, as_tibble(means), as_tibble(tops))
  W   <- matrix(aperm(A, c(1, 3, 2)), nrow = n); colnames(W) <- wide_names  # pos-major wide
  if (add_positions) res <- bind_cols(res, as_tibble(W))
  list(res = res, W = W)   # W always returned so PCA can use it even if add_positions=FALSE
}

built  <- lapply(names(nn_input), build_term)
params <- dplyr::bind_rows(lapply(built, `[[`, "res"))
Xfull  <- do.call(rbind, lapply(built, `[[`, "W"))   # n x (P*seq_len) raw per-position matrix
message(sprintf("built %d windows x %d cols (summaries%s)",
                nrow(params), ncol(params), if (add_positions) " + positions" else ""))

## ---- 20-dim PCA on the standardised per-position matrix --------------------
## Same input space as the UMAP (raw per-position features, standardised, with
## zero-variance columns dropped). Scores go in the table; loadings (feature x PC)
## are written to `load_csv` -- they're per-feature, so they can't live in the
## per-window table.
if (add_pca) {
  mu   <- colMeans(Xfull)
  sdv  <- sqrt(pmax(colMeans(Xfull * Xfull) - mu^2, 0))
  keep <- sdv > .Machine$double.eps
  Xk   <- scale(Xfull[, keep, drop = FALSE], center = mu[keep], scale = sdv[keep])
  k    <- min(pca_dims, ncol(Xk), nrow(Xk) - 1L)

  if (requireNamespace("irlba", quietly = TRUE)) {        # fast truncated PCA
    pca   <- irlba::prcomp_irlba(Xk, n = k, center = FALSE, scale. = FALSE)
    scores <- pca$x; rot <- pca$rotation; sdev <- pca$sdev
  } else {
    pca   <- prcomp(Xk, center = FALSE, scale. = FALSE)    # full PCA, then take first k
    scores <- pca$x[, 1:k, drop = FALSE]; rot <- pca$rotation[, 1:k, drop = FALSE]; sdev <- pca$sdev[1:k]
  }
  colnames(scores) <- paste0("PC", 1:k)
  params <- bind_cols(params, as_tibble(scores))

  ## loadings table: one row per kept per-position feature, PC1..PCk columns
  ve  <- sdev^2 / ncol(Xk)                                  # var explained per PC (std'd data: total var = #features)
  colnames(rot) <- paste0("PC", 1:k)
  loadings <- tibble::tibble(feature = colnames(Xk)) %>%
    dplyr::mutate(                                         # split pos12_cons_rs -> 12 + cons_rs
      position     = as.integer(stringr::str_match(feature, "^pos(\\d+)_")[, 2]),
      feature_type = stringr::str_remove(feature, "^pos\\d+_"),
      .after = feature) %>%
    bind_cols(as_tibble(rot))
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(loadings, load_csv)
  else readr::write_csv(loadings, load_csv)
  message(sprintf("PCA: %d PCs, cum. var explained = %.1f%%; loadings -> %s",
                  k, 100 * sum(ve), load_csv))
}

## ---- join the UMAP coords from a prior run (if we have them) ---------------
if (file.exists(in_csv)) {
  coord_cols <- c("peps", "terminus", "UMAP1", "UMAP2")
  um <- if (requireNamespace("data.table", quietly = TRUE))
          data.table::fread(in_csv, select = coord_cols)
        else readr::read_csv(in_csv, col_select = dplyr::any_of(coord_cols), show_col_types = FALSE)
  params <- dplyr::left_join(params, as_tibble(um), by = c("peps", "terminus"))
  message("merged UMAP1/UMAP2 from ", in_csv)
} else {
  message("no ", in_csv, " yet -- writing params without UMAP coords (merge later)")
}

## put the reduced-dim coords (UMAP + PCs) up front, positions stay at the end
params <- dplyr::relocate(params,
                          dplyr::any_of(c("UMAP1", "UMAP2")),
                          dplyr::matches("^PC[0-9]+$"),
                          .after = terminus)

## ---- write (fwrite is much faster than write_csv for the wide table) -------
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fwrite(params, out_csv)
} else {
  readr::write_csv(params, out_csv)
}
message("wrote ", out_csv, "  (", nrow(params), " x ", ncol(params), ")")
