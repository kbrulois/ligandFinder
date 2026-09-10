## ---- R <-> Python bridge for the 1D-CNN window models -----------------------
##
## The modeling lives in `inst/python/lf_dcnn` (see its README). This file is
## the whole boundary: it turns `nn_input` (as 9.2_add_contact_data.R builds it)
## into the array contract, hands the arrays to Python in memory via reticulate,
## and hands the scores back.
##
## All class bookkeeping (which column is `none`, which is `padding`, the order
## of the per-index softmax) is read FROM the Python Config rather than
## restated here, so the 1-based/0-based translation happens in exactly one
## place and cannot drift.

.lf_dcnn_env <- new.env(parent = emptyenv())

#' Locate the bundled lf_dcnn Python package
#' @keywords internal
lf_dcnn_path <- function(path = NULL) {
  if (!is.null(path)) return(normalizePath(path, mustWork = TRUE))
  candidates <- c(
    system.file("python", package = "ligandFinder"),
    file.path(getwd(), "inst", "python"),
    file.path(getwd(), "..", "inst", "python")
  )
  hit <- candidates[nzchar(candidates) &
                      file.exists(file.path(candidates, "lf_dcnn", "__init__.py"))]
  if (!length(hit))
    stop("could not find inst/python/lf_dcnn; pass path = explicitly")
  normalizePath(hit[[1]])
}

#' Import the lf_dcnn Python module (cached per session)
#'
#' @param venv virtualenv name or path. Defaults to the one keras3 already uses,
#'   which is what lets the arrays pass in memory with no conversion step.
#' @param path optional path to `inst/python`.
#' @export
lf_dcnn_python <- function(venv = "r-tensorflow", path = NULL, refresh = FALSE) {
  if (!refresh && !is.null(.lf_dcnn_env$mod)) return(.lf_dcnn_env$mod)
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("reticulate is required for the lf_dcnn bridge")
  ## Only bind the venv if python has not already been initialised (keras3 may
  ## have done it first, in which case it is already the right one).
  if (!reticulate::py_available(initialize = FALSE) && !is.null(venv))
    reticulate::use_virtualenv(venv, required = TRUE)
  .lf_dcnn_env$mod <- reticulate::import_from_path("lf_dcnn", path = lf_dcnn_path(path))
  lf_dcnn_utf8_streams()
  .lf_dcnn_env$mod
}

#' Make python's stdout/stderr UTF-8 tolerant
#'
#' Keras' progress bar draws box-drawing characters. Under `Rscript` python's
#' stdout can come up ASCII, and `fit(verbose = 1)` then dies with
#' `UnicodeEncodeError: 'ascii' codec can't encode characters`. Reconfiguring
#' the streams is a no-op wherever they are already UTF-8.
#' @keywords internal
lf_dcnn_utf8_streams <- function() {
  sys <- try(reticulate::import("sys", convert = FALSE), silent = TRUE)
  if (inherits(sys, "try-error")) return(invisible(FALSE))
  for (nm in c("stdout", "stderr"))
    try(sys[[nm]]$reconfigure(encoding = "utf-8", errors = "replace"), silent = TRUE)
  invisible(TRUE)
}

#' Build an lf_dcnn Config
#'
#' @param channel_names input channel order, i.e. `all_params3`.
#' @param r_exact reproduce 10_1dcnn_new6.R exactly, including its DB-mask
#'   overwrite and its non-functional per-epoch resample. See the package README.
#' @param ... any other Config field (`epochs`, `seed`, `noise_frac`, ...).
#' @export
lf_dcnn_config <- function(channel_names, r_exact = FALSE, ..., mod = NULL) {
  mod <- mod %||% lf_dcnn_python()
  args <- list(channel_names = as.character(channel_names),
               n_channels    = length(channel_names), ...)
  ctor <- if (isTRUE(r_exact)) mod$Config$r_exact else mod$Config
  do.call(ctor, args)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

## ---- the array contract -----------------------------------------------------

#' Turn one nn_input split into the three contract arrays
#'
#' Replaces `generate_keras_input()` in 10_1dcnn_new6.R. The `aperm` MUST happen
#' here: R is column-major and numpy is row-major, so handing over the
#' pre-`aperm` (seq_len, n_channels, n) array transposes silently and still
#' trains, just badly.
#' @keywords internal
lf_dcnn_split_arrays <- function(input, bk) {
  n <- length(input[["data"]])
  if (!n) stop("split has no windows")

  x <- array(unlist(input[["data"]]),
             dim = c(bk$seq_len, bk$n_channels, n))
  x <- aperm(x, c(3, 1, 2))                                   # (n, seq_len, C)

  ## one-hot over all classes, (seq_len, n, K) then permuted to (n, seq_len, K).
  ## unlist() walks positions fastest, which is exactly R's column-major fill
  ## order for the first axis, so the flat class vector drops straight in.
  cls <- bk$class_id[unlist(input[["known_idx"]], use.names = FALSE)]
  if (anyNA(cls))
    stop("unknown class label(s) in known_idx: ",
         paste(unique(unlist(input[["known_idx"]])[is.na(cls)]), collapse = ", "))
  onehot <- array(0, dim = c(bk$seq_len, n, bk$K))
  onehot[seq_along(cls) + cls * (bk$seq_len * n)] <- 1        # cls is 0-based
  onehot <- aperm(onehot, c(2, 1, 3))

  ## cols 1:K_cat = the [6 real, none] one-hot; col K_pi = the padding flag
  y_per_index_cat <- array(0, dim = c(n, bk$seq_len, bk$K_pi))
  y_per_index_cat[, , seq_len(bk$K_cat)] <- onehot[, , bk$cat_cols, drop = FALSE]
  y_per_index_cat[, , bk$K_pi]           <- onehot[, , bk$padding_col]

  list(x               = x,
       y_global        = matrix(as.numeric(input[["known"]]), ncol = 1),
       y_per_index_cat = y_per_index_cat)
}

#' Read the class bookkeeping out of a Python Config
#' @keywords internal
lf_dcnn_bookkeeping <- function(cfg) {
  class_names <- as.character(cfg$class_names)
  list(
    seq_len     = as.integer(cfg$seq_len),
    n_channels  = as.integer(cfg$n_channels),
    K           = length(class_names),
    K_cat       = as.integer(cfg$K_cat),
    K_pi        = as.integer(cfg$K_pi),
    class_names = class_names,
    ## 0-based ids, as the one-hot fill above expects
    class_id    = stats::setNames(seq_along(class_names) - 1L, class_names),
    ## Python indices are 0-based; R needs them 1-based
    cat_cols    = as.integer(cfg$cat_cols) + 1L,
    padding_col = as.integer(cfg$padding_col) + 1L,
    pi_names    = as.character(cfg$pi_names)
  )
}

#' Build the full nested array contract from nn_input
#'
#' @param nn_input the list of per-terminus `list(train, val, all)` sets.
#' @param cfg a Python Config (see [lf_dcnn_config()]).
#' @return `list[term][split]` of `list(x, y_global, y_per_index_cat)`.
#' @export
lf_dcnn_arrays <- function(nn_input, cfg, splits = c("train", "val", "all")) {
  bk <- lf_dcnn_bookkeeping(cfg)
  lapply(nn_input, function(term_sets) {
    use <- intersect(splits, names(term_sets))
    stats::setNames(lapply(use, function(s) lf_dcnn_split_arrays(term_sets[[s]], bk)), use)
  })
}

## ---- entry point 1: in-process ---------------------------------------------

#' Train both terminus models and score every window
#'
#' Rows of every returned array follow `names(nn_input)`, so they line up with
#' `bind_rows(lapply(nn_input, function(x) x$all))`.
#'
#' @param nn_input from 9.2_add_contact_data.R.
#' @param channel_names the input channel order, i.e. `all_params3`.
#' @param r_exact reproduce 10_1dcnn_new6.R exactly (see the package README).
#' @param verbose passed to keras `fit()`; 0 silences the per-epoch output.
#' @param n_seeds ensemble members per terminus. Each is a full retrain, and the
#'   global score and per-residue softmax come back as the mean across members
#'   with `pred_sd` / `per_index_sd` giving the spread. 1 reproduces the old
#'   single-model behaviour (sd all zero).
#' @param cache path to an .rds. If it exists and its fingerprint matches
#'   (config, n_seeds, term sizes), the whole run is loaded instead of retrained
#'   -- model weights included, so `models` still works downstream. Otherwise the
#'   run is trained and written there.
#' @param refresh TRUE to retrain and overwrite an existing cache.
#' @param keep_arrays return the built contract arrays as `$arrays`, which
#'   10_1dcnn_new6.R re-exposes as `nn_in_all` for the plotting scripts. Set
#'   FALSE to drop them once scoring is done (they are the bulk of the memory).
#' @param ... further Config fields (`epochs`, `seed`, ...).
#' @return list with `pred`, `pred_raw`, `pred_sd`, `per_index`, `per_index_sd`,
#'   `per_index_tbl`, `per_index_sd_tbl`, `emb`,
#'   `val_scores`, `val_labels`, `calibrator`, `histories`, `term_order`,
#'   `models`, `config` and (unless `keep_arrays = FALSE`) `arrays`.
#' @export
lf_dcnn_run <- function(nn_input, channel_names, r_exact = FALSE,
                        verbose = 1L, keep_arrays = TRUE, n_seeds = 5L,
                        cache = NULL, refresh = FALSE, ...,
                        venv = "r-tensorflow", path = NULL) {
  mod <- lf_dcnn_python(venv = venv, path = path)
  cfg <- lf_dcnn_config(channel_names, r_exact = r_exact, ..., mod = mod)
  cfg <- lf_dcnn_align_terms(cfg, names(nn_input))

  ## Cache fingerprint: anything that would change the answer. A cache whose
  ## fingerprint differs is ignored rather than silently reused.
  .fp <- list(cfg = cfg$to_dict(), n_seeds = as.integer(n_seeds),
              terms = names(nn_input),
              n = vapply(nn_input, \(x) nrow(x$all), integer(1)),
              npos = vapply(nn_input, \(x) sum(x$train$known), numeric(1)))
  .cp <- if (!is.null(cache)) path.expand(cache) else NULL

  if (!is.null(.cp) && file.exists(.cp) && !isTRUE(refresh)) {
    .cc <- readRDS(.cp)
    if (identical(.cc$fingerprint, .fp)) {
      message("lf_dcnn_run: reusing cache ", cache, " (no training)")
      out <- .cc$out
      ## rebuild the keras models from saved weights -- downstream scripts
      ## (10_2_per_ind_profiles.R, 10_3d_embed_umap.R) expect `models`
      wd <- paste0(tools::file_path_sans_ext(.cp), "_weights")
      if (dir.exists(wd)) {
        out$models <- stats::setNames(lapply(names(nn_input), function(tm) {
          m <- mod$build_model(cfg); m$load_weights(file.path(wd, paste0(tm, ".weights.h5"))); m
        }), names(nn_input))
      }
      if (isTRUE(keep_arrays) && is.null(out$arrays))
        out$arrays <- lf_dcnn_arrays(nn_input, cfg)
      return(out)
    }
    message("lf_dcnn_run: cache fingerprint differs -- retraining")
  }

  data <- lf_dcnn_arrays(nn_input, cfg)
  res  <- mod$run(data, cfg, verbose = as.integer(verbose),
                  n_seeds = as.integer(n_seeds))
  out  <- res$to_dict()

  out$per_index_tbl <- lf_dcnn_per_index_tibbles(out$per_index, out$pi_names)
  ## sd across ensemble members, same shape and column names as the mean, so the
  ## two line up position-for-position in the plots
  out$per_index_sd_tbl <- lf_dcnn_per_index_tibbles(out$per_index_sd, out$pi_names)
  out$pred_sd <- as.numeric(out$pred_sd)
  out$pred     <- as.numeric(out$pred)
  out$pred_raw <- as.numeric(out$pred_raw)
  ## The models come back as python keras objects. R no longer drives keras, but
  ## the downstream plotting scripts still call `predict(models[[term]], x)` on
  ## them -- and that S3 method only exists once keras3's namespace is loaded.
  ## Load it (without attaching) so the dispatch works for whoever gets these.
  if (!requireNamespace("keras3", quietly = TRUE))
    warning("keras3 is not installed: predict() on the returned models will not dispatch",
            immediate. = TRUE)
  out$models   <- res$models
  out$config   <- cfg                 # the single source of truth for the indices
  ## the downstream plotting scripts expect the built arrays as `nn_in_all`;
  ## hand back the ones we already built rather than paying for them twice.
  if (isTRUE(keep_arrays)) out$arrays <- data

  if (!is.null(.cp)) {
    dir.create(dirname(.cp), showWarnings = FALSE, recursive = TRUE)
    wd <- paste0(tools::file_path_sans_ext(.cp), "_weights")
    dir.create(wd, showWarnings = FALSE, recursive = TRUE)
    for (tm in names(out$models))
      out$models[[tm]]$save_weights(file.path(normalizePath(wd), paste0(tm, ".weights.h5")))
    ## `models` are python objects and cannot be serialised; the weights above
    ## restore them. `arrays` are rebuilt on load, so they are not stored either.
    saveRDS(list(fingerprint = .fp,
                 out = out[setdiff(names(out), c("models", "arrays"))]),
            .cp)
    message("lf_dcnn_run: cached to ", cache, "  (weights in ", basename(wd), ")")
  }
  out
}

#' Keep Config$term_order in step with names(nn_input)
#' @keywords internal
lf_dcnn_align_terms <- function(cfg, terms) {
  if (identical(as.character(cfg$term_order), as.character(terms))) return(cfg)
  cfg$evolve(term_order = as.character(terms))
}

#' Per-window per-index tibbles, one per window
#'
#' The R-side half of the old `extract_per_index()`: turns the (n, seq_len, K_pi)
#' array into a list of tidy tibbles for the downstream plots.
#' @export
lf_dcnn_per_index_tibbles <- function(arr, class_names) {
  d <- dim(arr); n <- d[[1L]]; seq_n <- d[[2L]]; K <- d[[3L]]
  if (length(class_names) != K) class_names <- class_names[seq_len(K)]
  idx <- seq_len(seq_n)
  col_names <- c("index", class_names)
  lapply(seq_len(n), function(i) {
    m <- arr[i, , , drop = TRUE]
    cols <- c(list(idx), lapply(seq_len(K), function(k) m[, k]))
    names(cols) <- col_names
    tibble::new_tibble(cols, nrow = seq_n)     # skip as_tibble/dplyr validation
  })
}

## ---- entry point 2: standalone, via disk ------------------------------------

#' Write the array contract to a directory the standalone CLI can read
#'
#' @param dir destination; gets `arrays.npz`, `config.json` and, if `meta` is
#'   supplied, `meta.parquet`.
#' @param meta optional one-row-per-window data frame for the concatenated
#'   `all` splits, in `names(nn_input)` order; carried into `predictions.parquet`.
#' @export
lf_dcnn_export <- function(nn_input, channel_names, dir, meta = NULL,
                           r_exact = FALSE, ..., venv = "r-tensorflow", path = NULL) {
  mod <- lf_dcnn_python(venv = venv, path = path)
  cfg <- lf_dcnn_config(channel_names, r_exact = r_exact, ..., mod = mod)
  cfg <- lf_dcnn_align_terms(cfg, names(nn_input))
  data <- lf_dcnn_arrays(nn_input, cfg)

  if (!is.null(meta)) {
    n_all <- sum(vapply(nn_input, function(x) length(x$all$data), integer(1)))
    if (nrow(meta) != n_all)
      stop(sprintf("meta has %d rows but the `all` splits hold %d windows",
                   nrow(meta), n_all))
    meta <- reticulate::r_to_py(as.data.frame(meta))
  }
  invisible(as.character(mod$io$save_inputs(dir, data, cfg, meta = meta)))
}

#' Read back what the standalone CLI wrote
#' @export
lf_dcnn_import <- function(dir, venv = "r-tensorflow", path = NULL) {
  mod <- lf_dcnn_python(venv = venv, path = path)
  out <- mod$io$load_outputs(dir)
  out$pred     <- as.numeric(out$pred)
  out$pred_raw <- as.numeric(out$pred_raw)
  out
}
