## ---- peptide-end tracks for the protein page ----------------------------------
## The tracks the window models contribute to a protein page, assembled here
## rather than inside make_protein_plot_win() so they can be built and tested on
## their own. Nine tracks at most, all on the residue axis:
##
##   N / C  per-residue class probability   pooled over every window covering r
##   N / C  CNN attention head              the global head, on the anchor axis
##   N / C  xgb window head                 the same quantity, trees over the surface
##   MLP / XGB per residue                  window-free, one residue's own features
##
## EVERY INPUT IS OPTIONAL. A page is rendered for all 5,616 precursors and the
## per-index surface costs ~40 s per gene per terminus of model inference, so a
## bulk run cannot pay for it: a missing input drops its track rather than
## failing the page. `lf_pepend_track_data()` returns NULL for anything it
## cannot build, and the caller skips it.
## ------------------------------------------------------------------------------

#' Where the peptide-end inputs live
#' @keywords internal
LF_TRACK_PATHS <- list(
  run   = c(N = "~/AF2_analysis/lf_pepend_run_N_near.rds",
            C = "~/AF2_analysis/lf_pepend_run_C_near.rds"),
  scan  = c(N = "~/AF2_analysis/lf_pepend_scan_N_near.rds",
            C = "~/AF2_analysis/lf_pepend_scan_C_near.rds"),
  winxgb = c(N = "~/AF2_analysis/lf_winxgb_scan_N_near.npz",
             C = "~/AF2_analysis/lf_winxgb_scan_C_near.npz"),
  cache = "~/AF2_analysis/lf_pepend_residue_cache.rds",
  resid = "~/AF2_analysis/lf_resid_preds.npz"
)

## loaded-once holders: a bulk run would otherwise re-read a 50 MB scan and
## rebuild a 20-member ensemble for every gene
.lf_tracks_env <- new.env(parent = emptyenv())

## Where THIS file was sourced from, captured at source time. `lf_dcnn_python()`
## with no `path` resolves the INSTALLED ligandFinder's inst/python first, which
## silently shadows a worktree and fails far from the cause -- so the bridge is
## always given the inst/python that sits beside this file.
.LF_TRACKS_DIR <- local({
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f, mustWork = FALSE)))
  }
  NA_character_
})

## LF_PEPEND / lf_pepend_slice come from pepend_windows.R, which a page driver
## does not necessarily source (plot_only.R sources the plotting code only), and
## which an installed ligandFinder may predate. Pull it in from beside this file
## when it is missing rather than failing deep inside the inference.
if (!exists("LF_PEPEND") && !is.na(.LF_TRACKS_DIR)) {
  .pw <- file.path(.LF_TRACKS_DIR, "pepend_windows.R")
  if (file.exists(.pw)) source(.pw)
}

#' One input path, with an option override
#'
#' Lets a caller point the tracks at a different trunk per TERMINUS -- which is
#' what reading a cross-validated page needs: a gene is held out in one fold at
#' its C end and in the other at its N end, so the honest page takes C from one
#' fold and N from the other.
#'
#'   options(lf.track_run_C = "~/AF2_analysis/lf_pepend_run_C_foldA.rds",
#'           lf.track_run_N = "~/AF2_analysis/lf_pepend_run_N_foldB.rds")
#' @keywords internal
lf_track_path <- function(kind, term) {
  o <- getOption(sprintf("lf.track_%s_%s", kind, term), NULL)
  if (!is.null(o)) return(path.expand(o))
  path.expand(LF_TRACK_PATHS[[kind]][[term]])
}

#' @keywords internal
lf_tracks_python_path <- function() {
  p <- getOption("lf.pepend_python_path", NULL)
  if (!is.null(p)) return(p)
  if (!is.na(.LF_TRACKS_DIR)) {
    cand <- file.path(dirname(.LF_TRACKS_DIR), "inst", "python")
    if (dir.exists(cand)) return(cand)
  }
  NULL
}

#' The 20 member models of one run, built once per process
#'
#' A bulk render calls this once per gene per series; rebuilding 20 Keras models
#' each time is the bulk of the per-page cost, and they are only ever used for
#' prediction, so they are cached on (weights dir, terminus).
#' @keywords internal
lf_track_models <- function(mod, iso, term) {
  cfg <- .lf_once(paste0("cfg_", iso), function() mod$Config$from_json(file.path(iso, "in", "config.json")))
  list(cfg = cfg,
       models = .lf_once(paste0("members_", iso, "_", term),
                         function() mod$scan$load_member_models(file.path(iso, "out"), cfg, term)))
}

#' @keywords internal
.lf_once <- function(key, f) {
  if (!is.null(.lf_tracks_env[[key]])) return(.lf_tracks_env[[key]])
  .lf_tracks_env[[key]] <- tryCatch(f(), error = function(e) NA)
  .lf_tracks_env[[key]]
}

#' The CNN window head for one gene, straight from the step-1 scan
#'
#' Cheap and complete: the scan already holds every mature anchor of every
#' precursor, so this track needs no inference and is available on every page.
#' @keywords internal
lf_track_window_head <- function(gene, term, path = lf_track_path("scan", term)) {
  p <- path.expand(path)
  if (!file.exists(p)) return(NULL)
  s <- .lf_once(paste0("scan_", term, "_", basename(p)), function() readRDS(p))
  if (!is.list(s) || is.null(s$scan)) return(NULL)
  d <- s$scan[s$scan$gene == gene, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  tibble::tibble(residue = as.integer(d$anchor), mean = d$score, sd = d$sd,
                 term = term, head = "CNN attention head",
                 head_label = getOption("lf.track_cnn_label", "CNN attention head"))
}

#' The xgboost window head, from the per-gene scan npz
#'
#' Only the genes that scan covers; absent gene -> NULL -> no track.
#' @keywords internal
lf_track_xgb_head <- function(gene, term, path = lf_track_path("winxgb", term)) {
  p <- path.expand(path)
  if (!file.exists(p) || !requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  z <- .lf_once(paste0("wx_", term, "_", basename(p)), function() {
    np <- reticulate::import("numpy", convert = FALSE)
    np$load(p, allow_pickle = TRUE)
  })
  if (!inherits(z, "python.builtin.object")) return(NULL)
  f <- tryCatch(as.character(reticulate::py_to_r(z$files)), error = function(e) character(0))
  k <- paste0(gene, "__anchor")
  if (!k %in% f) return(NULL)
  tibble::tibble(
    residue = as.integer(reticulate::py_to_r(z[[k]])),
    mean    = as.numeric(reticulate::py_to_r(z[[paste0(gene, "__xgb_pm")]])),
    sd      = as.numeric(reticulate::py_to_r(z[[paste0(gene, "__xgb_pm__sd")]])),
    term = term, head = "xgb window head")
}

#' The window-free per-residue models
#'
#' Trunk-independent: the same four binary targets serve both termini.
#' @keywords internal
lf_track_resid_models <- function(accession, n_prot, c_prot,
                                  path = getOption("lf.track_resid", LF_TRACK_PATHS$resid)) {
  p <- path.expand(path)
  if (!file.exists(p) || !requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  z <- .lf_once(paste0("resid_", basename(p)), function() {
    np <- reticulate::import("numpy", convert = FALSE)
    np$load(p, allow_pickle = TRUE)
  })
  if (!inherits(z, "python.builtin.object")) return(NULL)
  have <- tryCatch(as.character(reticulate::py_to_r(z$files)), error = function(e) character(0))
  cols <- grep("^(pep_pocket|pep_other|ct_context|nt_context)__(mlp|xgb)$", have, value = TRUE)
  if (!length(cols)) return(NULL)
  meta <- .lf_once(paste0("resid_meta_", basename(p)), function() list(
    acc  = as.character(reticulate::py_to_r(z[["accession"]])),
    prot = as.integer(reticulate::py_to_r(z[["prot_idx"]])),
    res  = as.integer(reticulate::py_to_r(z[["resno"]]))))
  i <- match(accession, meta$acc)
  if (is.na(i)) return(NULL)
  sel <- which(meta$prot == (i - 1L))
  if (!length(sel)) return(NULL)
  out <- lapply(cols, function(cc) {
    v <- .lf_once(paste0("resid_", cc), function() as.numeric(reticulate::py_to_r(z[[cc]])))
    sdc <- paste0(cc, "__sd")
    sv <- if (sdc %in% have)
      .lf_once(paste0("resid_", sdc), function() as.numeric(reticulate::py_to_r(z[[sdc]]))) else NULL
    tibble::tibble(residue = meta$res[sel], mean = v[sel],
                   sd = if (is.null(sv)) NA_real_ else sv[sel],
                   target = sub("__.*$", "", cc), model = toupper(sub("^.*__", "", cc)))
  })
  d <- dplyr::bind_rows(out)
  d <- d[d$residue >= n_prot & d$residue <= c_prot, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  ## lf_resid's target names ARE the per-index vocabulary under another name
  d$class <- factor(c(pep_pocket = "pep_pocket", pep_other = "pep_other",
                      ct_context = "CT_cleavage_context",
                      nt_context = "NT_cleavage_context")[d$target],
                    levels = names(nn_class_cols))
  d
}

#' Score every window of one precursor, both termini
#'
#' The expensive one: 20 members x one window per mature residue, per terminus.
#' Returns the per-index surface (for the popups), the pooled per-residue tracks
#' and the global head. NULL when the bridge or the weights are unavailable,
#' which is the normal case in a bulk run.
#' @keywords internal
lf_track_per_index <- function(accession, n_prot, c_prot, terms = c("N", "C")) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  cache_p <- path.expand(LF_TRACK_PATHS$cache)
  if (!file.exists(cache_p)) return(NULL)
  fc <- .lf_once("feat_cache", function() readRDS(cache_p))
  if (!is.list(fc)) return(NULL)
  i <- match(accession, fc$prec$accession)
  if (is.na(i)) return(NULL)

  .pp <- lf_tracks_python_path()
  mod <- tryCatch(if (is.null(.pp)) lf_dcnn_python() else lf_dcnn_python(path = .pp),
                  error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  np <- reticulate::import("numpy", convert = FALSE)
  feat <- fc$feats[[i]]
  anchors <- seq.int(n_prot, c_prot)

  out <- list()
  for (term in terms) {
    run_p <- lf_track_path("run", term)
    iso <- paste0(tools::file_path_sans_ext(run_p), "_isolated")
    if (!dir.exists(iso)) next
    got <- tryCatch({
      .m <- lf_track_models(mod, iso, term)
      cfg <- .m$cfg; models <- .m$models
      pin <- as.character(cfg$pi_names)
      A <- LF_PEPEND$anchor[[term]]
      w_start <- anchors - A + 1L
      X <- array(0, c(length(anchors), LF_PEPEND$seq_len, ncol(feat)))
      for (k in seq_along(anchors))
        X[k, , ] <- as.matrix(lf_pepend_slice(feat, w_start[k], n_prot, c_prot)[, fc$all_params3])
      xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))
      PI <- array(0, c(length(models), length(anchors), LF_PEPEND$seq_len, length(pin)))
      G  <- matrix(0, length(models), length(anchors))
      for (m in seq_along(models)) {
        o <- mod$pipeline$predict_all(models[[m]], xp)
        PI[m, , , ] <- reticulate::py_to_r(o[["per_index_cat"]])
        G[m, ] <- as.numeric(reticulate::py_to_r(o[["global"]]))
      }
      list(pi = apply(PI, c(2, 3, 4), mean),
           ## the band in the opened window: real disagreement between the 20
           ## independently trained members, which a single softmax cannot
           ## express -- it renders "confidently 0.5" and "the members split"
           ## identically
           pi_sd = apply(PI, c(2, 3, 4), stats::sd),
           g = colMeans(G),
           g_sd = apply(G, 2, stats::sd), pi_names = pin, pos = A, w_start = w_start)
    }, error = function(e) { message("  per-index ", term, " failed: ",
                                     conditionMessage(e)); NULL })
    if (!is.null(got)) out[[term]] <- got
  }
  if (!length(out)) return(NULL)
  out$anchors <- anchors
  out
}

#' An extra window-head series from another trunk, scored on the fly
#'
#' For putting a second model's window score beside the page's own -- e.g. the
#' CT-padded C model -- without a precomputed scan. Every mature residue is an
#' anchor; each of the run's members scores every window and the series is the
#' member mean with the member sd as its band, as on the scan-backed series.
#'
#' Set through an option, so the page code needs no new argument:
#'
#' Extra series are drawn on a track of their own beside each head -- one per
#' `group` (default "extra") -- coloured by terminus like the page's own.
#'
#'   options(lf.track_extra_win = list(list(
#'     name = "C terminus, CT-padded", term = "C", ct_pad = TRUE, group = "padded",
#'     run = "~/AF2_analysis/lf_pepend_run_C_near_ctpad_slack8.rds",
#'     winxgb = "~/AF2_analysis/lf_winxgb_scan_C_near_ctpad_slack8.npz")))
#'
#' `ct_pad = TRUE` cuts every C window at its own anchor, as the run was trained
#' (10_7 --ct-pad), and `nt_pad = TRUE` starts every N window at its anchor
#' (10_7 --nt-pad); scoring a padded model on unpadded windows would show it
#' inputs it has never seen. `per_residue = TRUE` also pools that run's
#' per-index class probabilities onto residues, for a class-probability track of
#' its own; each window counts only for residues it saw, so a padded window says
#' nothing about the residues under its padding. `popup = TRUE` puts its windows
#' in the click popup: a click on that series opens that run's window. `winxgb` (optional) adds the same series to the xgb
#' window-head track from a `lf_winxgb.scan_genes` npz -- which has to be built
#' from the same run, with `--ct-pad` for a padded one.
#' @keywords internal
lf_track_extra_head <- function(spec, accession, n_prot, c_prot) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  fc <- .lf_once("feat_cache", function() readRDS(path.expand(LF_TRACK_PATHS$cache)))
  if (!is.list(fc)) return(NULL)
  i <- match(accession, fc$prec$accession)
  if (is.na(i)) return(NULL)
  ## The mature range the models were TRAINED with, not the caller's: the page
  ## passes n_prot = 1, and a window over the signal peptide would show the
  ## model residues it only ever saw as padding. The scan-backed series start
  ## at this same n_prot, so the extra series lines up with them.
  n_prot <- as.integer(fc$prec$n_prot[i]); c_prot <- as.integer(fc$prec$c_prot[i])
  .pp <- lf_tracks_python_path()
  mod <- tryCatch(if (is.null(.pp)) lf_dcnn_python() else lf_dcnn_python(path = .pp),
                  error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  np <- reticulate::import("numpy", convert = FALSE)
  term <- spec$term
  iso <- paste0(tools::file_path_sans_ext(path.expand(spec$run)), "_isolated")
  if (!dir.exists(iso)) { message("  extra series '", spec$name, "': no ", iso); return(NULL) }
  tryCatch({
    .m <- lf_track_models(mod, iso, term)
    cfg <- .m$cfg; models <- .m$models
    ## options(lf.track_member = k): draw the page from ONE ensemble member,
    ## 0-based as the weight files are named, instead of averaging all 20
    .one <- getOption(sprintf("lf.track_member_%s", term), getOption("lf.track_member", NULL))
    if (!is.null(.one)) {
      k <- as.integer(.one) + 1L
      if (k < 1L || k > length(models))
        stop("lf.track_member ", .one, " is outside 0..", length(models) - 1L, call. = FALSE)
      models <- models[k]
    }
    feat <- fc$feats[[i]]
    anchors <- seq.int(n_prot, c_prot)
    w_start <- anchors - LF_PEPEND$anchor[[term]] + 1L
    ## the window's own end and start: the precursor's, or under ct_pad the C
    ## anchor / under nt_pad the N anchor
    cut <- if (isTRUE(spec$ct_pad) && term == "C") pmin(c_prot, anchors) else rep(c_prot, length(anchors))
    cut_lo <- if (isTRUE(spec$nt_pad) && term == "N") pmax(n_prot, anchors) else rep(n_prot, length(anchors))
    X <- array(0, c(length(anchors), LF_PEPEND$seq_len, ncol(feat)))
    for (k in seq_along(anchors))
      X[k, , ] <- as.matrix(lf_pepend_slice(feat, w_start[k], cut_lo[k], cut[k])[, fc$all_params3])
    xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))
    G <- matrix(0, length(models), length(anchors))
    pin <- as.character(cfg$pi_names)
    want_pi <- isTRUE(spec$per_residue) || isTRUE(spec$popup)
    PI <- PI2 <- if (want_pi)
      array(0, c(length(anchors), LF_PEPEND$seq_len, length(pin))) else NULL
    for (m in seq_along(models)) {
      o <- mod$pipeline$predict_all(models[[m]], xp)
      G[m, ] <- as.numeric(reticulate::py_to_r(o[["global"]]))
      if (want_pi) {
        p <- reticulate::py_to_r(o[["per_index_cat"]])
        PI <- PI + p; PI2 <- PI2 + p^2
      }
    }
    PSD <- NULL
    if (want_pi) {
      ## members averaged first, as on the page's own class track; the sd is
      ## across members, as in the page's own popup
      nm <- length(models)
      PSD <- sqrt(pmax(PI2 - PI^2 / nm, 0) / max(nm - 1, 1))
      PI <- PI / nm
    }
    resid <- if (isTRUE(spec$per_residue))
      lf_track_pool_residues(list(pi = PI, pi_names = pin, w_start = w_start),
                             anchors, n_prot, c_prot, term, cut = cut, cut_lo = cut_lo) else NULL
    ## for the click popup: the same shape as lf_track_per_index() returns, plus
    ## `vis`, the window positions the model actually saw -- 1-28 for a
    ## CT-padded window, 8-36 for an NT-padded one
    pix <- if (isTRUE(spec$popup))
      list(pi = PI, pi_sd = PSD, g = colMeans(G), pi_names = pin,
           pos = LF_PEPEND$anchor[[term]], w_start = w_start, anchors = anchors, term = term,
           vis = c(if (isTRUE(spec$nt_pad) && term == "N") LF_PEPEND$anchor[["N"]] else 1L,
                   if (isTRUE(spec$ct_pad) && term == "C") LF_PEPEND$anchor[["C"]]
                   else LF_PEPEND$seq_len),
           cnn_label = lf_extra_cnn_label(spec)) else NULL
    ## `head` stays "CNN attention head" whatever the run's head is: it is the
    ## KEY that puts the series on the CNN track. What the reader sees -- axis,
    ## tooltip, popup -- comes from `head_label` (spec$cnn_label).
    list(win = tibble::tibble(residue = anchors, mean = colMeans(G), sd = apply(G, 2, stats::sd),
                              term = spec$name, head = "CNN attention head",
                              head_label = lf_extra_cnn_label(spec),
                              group_ens = isTRUE(spec$group_ensemble),
                              ## spec$term, not `term`: inside tibble() that
                              ## name is the column just defined above
                              terminus = spec$term, group = lf_extra_group(spec)),
         resid = resid, pix = pix)
  }, error = function(e) { message("  extra series '", spec$name, "' failed: ",
                                   conditionMessage(e)); NULL })
}

#' The track group an extra series is drawn in
#' @keywords internal
lf_extra_group <- function(spec) if (is.null(spec$group)) "extra" else as.character(spec$group)

#' What an extra run's CNN window score is called on the page
#'
#' `spec$cnn_label`, for a run whose window head is not attention (e.g.
#' "CNN flatten head"); the attention default otherwise.
#' @keywords internal
lf_extra_cnn_label <- function(spec)
  if (is.null(spec$cnn_label)) "CNN attention head" else as.character(spec$cnn_label)

#' Pool each window's opinion onto the residues it covers
#'
#' Members are averaged first -- the ensemble mean is the prediction -- then mean
#' and sd across the COVERING WINDOWS, so the band says how much the call depends
#' on which window you look through, not how much the seeds disagree.
#' @keywords internal
lf_track_pool_residues <- function(S, anchors, n_prot, c_prot, term,
                                   classes = c("pep_pocket", "pep_other",
                                               "CT_cleavage_context", "NT_cleavage_context"),
                                   cut = NULL, cut_lo = NULL) {
  idx <- which(S$pi_names %in% classes)
  if (!length(idx)) return(NULL)
  res <- lapply(idx, function(k) {
    acc <- lapply(seq_along(anchors), function(i) {
      r <- S$w_start[i] + seq_len(LF_PEPEND$seq_len) - 1L
      ## `cut`: where THIS window's input ends (its anchor, for a CT-padded
      ## window) -- it has no opinion about residues under its padding
      ok <- r >= (if (is.null(cut_lo)) n_prot else max(n_prot, cut_lo[i])) &
            r <= (if (is.null(cut)) c_prot else cut[i])
      tibble::tibble(residue = r[ok], v = S$pi[i, ok, k])
    })
    d <- dplyr::bind_rows(acc)
    d <- dplyr::summarise(dplyr::group_by(d, residue),
                          mean = mean(v), sd = stats::sd(v), n_win = dplyr::n(),
                          .groups = "drop")
    d$class <- S$pi_names[k]; d$term <- term
    d
  })
  out <- dplyr::bind_rows(res)
  out$class <- factor(out$class, levels = names(nn_class_cols))
  out
}

#' Everything the page needs, or NULL for whatever cannot be built
#' Which fold held each known end out, for this precursor
#'
#' Under 5-fold CV a known end is TRAINING in four folds and held out in exactly
#' one, so on a per-fold track there is a single series whose value at that
#' residue is honest -- the other four trained on it. That one point gets a mark
#' of its own. Reads `lf.cv_knowns`, a list(N =, C =) of the CV knowns csvs.
#' @keywords internal
lf_track_heldout_marks <- function(accession) {
  paths <- getOption("lf.cv_knowns", NULL)
  if (is.null(paths) || !length(paths)) return(NULL)
  out <- lapply(names(paths), function(tm) {
    p <- path.expand(paths[[tm]])
    if (!file.exists(p)) return(NULL)
    d <- utils::read.csv(p, stringsAsFactors = FALSE)
    if (!all(c("accession", "anchor", "fold") %in% names(d))) return(NULL)
    d <- d[d$accession == accession & !is.na(d$fold), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    data.frame(term = sprintf("%s fold %d", tm, as.integer(d$fold)),
               terminus = tm, residue = as.integer(d$anchor),
               fold = as.integer(d$fold), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out) || !nrow(out)) NULL else out
}

#' Held-out marks for extra runs that are not the 5-fold CV
#'
#' `options(lf.heldout_extra = <data.frame>)` with columns `term` (the extra
#' series' name, e.g. "C flatten fold A"), `terminus`, `accession`, `residue`:
#' the series whose fold held that known end out, drawn as triangles exactly
#' like the 5-fold ones.
#' @keywords internal
lf_track_heldout_extra <- function(accession) {
  x <- getOption("lf.heldout_extra", NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  x <- x[x$accession == accession, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  data.frame(term = as.character(x$term), terminus = as.character(x$terminus),
             residue = as.integer(x$residue), fold = NA_integer_, stringsAsFactors = FALSE)
}

#' The marker a tooltip carries when its series has no window popup
#' @keywords internal
## (no apostrophe: it lands in an svg attribute, and in the JS that reads it)
LF_NO_POPUP <- "no window popup (scores from the full-proteome scan of this arm)"

#' An arm's two window heads for one precursor, from its tracks rds
#'
#' `options(lf.track_arm_tracks = list(list(label = "ESM-C only",
#'     N = <lf_pepend_cv5_tracks_N_flatten_esmonly36.rds>, C = <..._C_...>), ...))`
#'
#' The tracks rds is what an arm's full-proteome scan assembles
#' (tracks_esmonly36.R): every anchor of every precursor, the CNN linear head
#' and the per-member xgb head per fold (`cnn__f1..5`, `xgb_pm__f1..5`) and
#' their mean over all 50 members (`*__all50`). That covers ANY gene, so an arm
#' can be drawn on a page without the per-gene scan and PLM table that
#' `lf.track_extra_win` needs -- which exist only for the genes an arm's own
#' gallery was built for. The values are the same scan's: the curated arm's
#' all-50 rds matches its per-gene page scan to 5e-7.
#'
#' What this source does not have: the per-seed spread within a fold (no
#' ribbon) and the per-index surface (no click-to-open window). Rows come back
#' in lf.track_extra_win's shape -- one group per arm, named by `label`, five
#' thin fold series and the all-50 as the group's own heavy line -- with
#' `no_popup = TRUE`.
#' @keywords internal
lf_track_arm_heads <- function(arm, accession) {
  out <- list()
  for (tm in c("N", "C")) {
    p <- arm[[tm]]
    if (is.null(p) || !file.exists(path.expand(p))) next
    b <- .lf_once(paste0("armtracks_", p), function() {
      x <- data.table::as.data.table(readRDS(path.expand(p)))
      data.table::setkeyv(x, "accession")
      x
    })
    if (!data.table::is.data.table(b)) next
    ## NOT b[list(accession)]: inside [ that is b's own column, a self-join
    acc <- accession
    d <- b[list(acc), nomatch = NULL]
    if (!nrow(d)) next
    folds <- sub("^cnn__f", "", grep("^cnn__f[0-9]+$", names(d), value = TRUE))
    for (hd in c("cnn", "xgb_pm")) {
      key <- if (hd == "cnn") "CNN attention head" else "xgb window head"
      hl  <- if (hd == "cnn") sprintf("CNN linear head (%s)", arm$label) else NA_character_
      one <- function(col, term, ens)
        tibble::tibble(residue = as.integer(d$anchor), mean = as.numeric(d[[col]]),
                       sd = NA_real_, term = term, terminus = tm, head = key,
                       head_label = hl, group = arm$label, group_ens = ens, no_popup = TRUE)
      for (k in folds)
        out[[length(out) + 1L]] <- one(sprintf("%s__f%s", hd, k),
                                       sprintf("%s %s fold %s", arm$label, tm, k), FALSE)
      if (paste0(hd, "__all50") %in% names(d))
        out[[length(out) + 1L]] <- one(paste0(hd, "__all50"),
                                       sprintf("%s all 50, %s terminus", arm$label, tm), TRUE)
    }
  }
  if (!length(out)) NULL else dplyr::bind_rows(out)
}

#' MultiPep's class columns, in the candidate table's own `mp_*` names
#'
#' `mp_signalling_max` is a derived column (max of three of these), not a class,
#' so it never competes for "best".
#' @keywords internal
lf_mp_classes <- function(nms) setdiff(grep("^mp_", nms, value = TRUE), "mp_signalling_max")

#' Candidate sets to draw: one spec, or a list of them
#'
#' `options(lf.track_candidates = <spec> | list(<spec>, <spec>, ...))`, each spec
#'
#'   list(label = "xgb all 50 >= 0.5",   # names the set on both of its tracks
#'        segments = <end_segments_*.csv>, candidates = <candidate_peptides_*.csv>,
#'        tol = 4, thresh = 0.5)
#'
#' One spec per ARM or per WAY OF DEFINING END SEGMENTS; every set draws its own
#' pair of tracks, in the order given. The page knows nothing about how a set's
#' segments were defined -- it reads the segment file the set names, which needs
#' term, accession, start, end, centre, centre_score, n_anchors, motif_centre,
#' n_candidates -- so a new definition is a new file, not a page change. The
#' label must not contain "xgb window": the popup wires clicks on any track
#' whose text matches that.
#' @keywords internal
lf_track_candidate_specs <- function() {
  o <- getOption("lf.track_candidates", NULL)
  if (is.null(o)) return(list())
  if (!is.null(o$segments)) o <- list(o)               # a bare spec
  o
}

#' End segments and the candidate peptides built from them, for one precursor
#'
#' For one spec (see lf_track_candidate_specs): the table `candidate_pipeline.sh`
#' writes and the `end_segments_*.csv` its stage 1 writes beside it. The
#' segments are READ, not re-derived from the page's own xgb series: they are
#' then by construction the ones the candidates came from. A missing file or a
#' precursor with no segment -> NULL, and neither track is drawn.
#'
#' Each candidate gets `best_type` / `best_score`: the MultiPep class with the
#' highest score and that score. A table that has not been through stage 4 has
#' no `mp_*` columns; its candidates are then all "no MultiPep score".
#' @keywords internal
lf_track_candidates <- function(o, accession) {
  if (is.null(o$segments) || !file.exists(path.expand(o$segments))) return(NULL)
  seg_all <- .lf_once(paste0("cand_seg_", o$segments), function()
    as.data.frame(data.table::fread(path.expand(o$segments))))
  if (!is.data.frame(seg_all)) return(NULL)
  seg <- seg_all[seg_all$accession == accession, , drop = FALSE]
  if (!nrow(seg)) return(NULL)
  cand <- NULL
  if (!is.null(o$candidates) && file.exists(path.expand(o$candidates))) {
    cand_all <- .lf_once(paste0("cand_tab_", o$candidates), function() {
      d <- data.table::fread(path.expand(o$candidates))
      mp <- lf_mp_classes(names(d))
      if (length(mp)) {
        m <- as.matrix(d[, mp, with = FALSE])
        k <- max.col(replace(m, is.na(m), -Inf), ties.method = "first")
        d[, best_type := sub("^mp_", "", mp[k])]
        d[, best_score := m[cbind(seq_len(nrow(m)), k)]]
      } else {
        d[, `:=`(best_type = "no MultiPep score", best_score = NA_real_)]
      }
      as.data.frame(d)
    })
    if (is.data.frame(cand_all)) {
      cand <- cand_all[cand_all$accession == accession, , drop = FALSE]
      if (!nrow(cand)) cand <- NULL
    }
  }
  list(seg = seg, cand = cand,
       label = if (is.null(o$label)) basename(o$segments) else as.character(o$label),
       tol = if (is.null(o$tol)) NA_integer_ else as.integer(o$tol),
       thresh = if (is.null(o$thresh)) NA_real_ else as.numeric(o$thresh))
}

#' @export
lf_pepend_track_data <- function(gene, accession, n_prot, c_prot,
                                 want_per_index = TRUE) {
  ## ORDER MATTERS. `lf_dcnn_python()` binds the r-tensorflow venv only if
  ## reticulate has not already been initialised -- and reading any npz below
  ## initialises it against whatever interpreter reticulate picks by default,
  ## which has no tensorflow. Loading the members would then fail, be caught,
  ## and the class-probability tracks would silently disappear while the npz
  ## tracks rendered fine. So the bridge goes first, every time.
  pix <- if (isTRUE(want_per_index))
    lf_track_per_index(accession, n_prot, c_prot) else NULL

  win <- dplyr::bind_rows(
    lf_track_window_head(gene, "N"), lf_track_window_head(gene, "C"),
    lf_track_xgb_head(gene, "N"),    lf_track_xgb_head(gene, "C"))
  ## extra series (another trunk's window head), after the bridge is up
  extra_resid <- list(); extra_pix <- list()
  for (spec in getOption("lf.track_extra_win", list())) {
    ex <- lf_track_extra_head(spec, accession, n_prot, c_prot)
    win <- dplyr::bind_rows(win, ex$win)
    if (!is.null(ex$resid) && nrow(ex$resid)) extra_resid[[spec$name]] <- ex$resid
    if (!is.null(ex$pix)) extra_pix[[spec$name]] <- ex$pix
    if (!is.null(spec$winxgb)) {
      xh <- lf_track_xgb_head(gene, spec$term, spec$winxgb)
      if (!is.null(xh)) {
        xh$term <- spec$name; xh$terminus <- spec$term; xh$group <- lf_extra_group(spec)
        xh$group_ens <- isTRUE(spec$group_ensemble)
        win <- dplyr::bind_rows(win, xh)
      }
    }
  }
  ## other arms, from their full-proteome tracks: after the page's own groups,
  ## in the order given
  for (arm in getOption("lf.track_arm_tracks", list()))
    win <- dplyr::bind_rows(win, lf_track_arm_heads(arm, accession))
  if (!nrow(win)) win <- NULL
  rp <- lf_track_resid_models(accession, n_prot, c_prot)
  resid <- NULL
  if (!is.null(pix)) {
    resid <- dplyr::bind_rows(lapply(c("N", "C"), function(t) {
      if (is.null(pix[[t]])) return(NULL)
      lf_track_pool_residues(pix[[t]], pix$anchors, n_prot, c_prot, t)
    }))
    if (!nrow(resid)) resid <- NULL
  }
  cand <- Filter(Negate(is.null), lapply(lf_track_candidate_specs(), lf_track_candidates,
                                         accession = accession))
  if (!length(cand)) cand <- NULL
  if (is.null(win) && is.null(rp) && is.null(resid) && is.null(cand)) return(NULL)
  if (!is.null(pix) && length(extra_pix)) pix$extra <- extra_pix
  list(win = win, rp = rp, resid = resid, extra_resid = extra_resid, pix = pix,
       cand = cand,
       heldout = rbind(lf_track_heldout_marks(accession), lf_track_heldout_extra(accession)),
       anchors = if (is.null(pix)) NULL else pix$anchors)
}

## ---- the plots -----------------------------------------------------------------
## Each track is its own ggplot and, at the call site, its own girafe widget:
## ggplot2 4.x + patchwork 1.3.1 + ggiraph 0.8.13 truncate a patchwork of more
## than two panels inside dsvg, silently. `data_id = residue` on every point is
## what links hover across all of them, so pointing at a residue on one track
## marks it on every other.

#' The two termini, where a track draws both
#'
#' Deliberately outside `nn_class_cols`: on a combined head track the colour
#' means WHICH TERMINUS, not which class, and reusing a class colour there would
#' say the opposite. Green/orange: a blue N sat next to `pep_pocket`'s #197EC0
#' on the same page and read as the same thing. The green is a mid-saturation
#' one, well clear of `pep_other`'s pale #D5E4A2 without going neon.
#' @export
LF_TERM_COLS <- c(N = "#2BA84A", C = "#D95F0E")

#' The same two termini, darkened, for the all-50 ensemble drawn over a fold track
#'
#' Same hue as the fold series so the terminus still reads at a glance, but
#' clearly heavier -- the ensemble is the reference line, the folds are the
#' spread around it.
#' @export
LF_TERM_COLS_ENS <- c(N = "#13632A", C = "#8A3C08")



#' @keywords internal
lf_track_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 7.5, lineheight = 0.95),
      axis.text.y = ggplot2::element_text(size = 6.5),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 1, r = 10, b = 1, l = 80))
}

#' The page's shared residue axis
#'
#' seq_p, the detail panels and the strip all use these three settings together
#' (limits 0..max_index+1, `expansion(add = c(seq_offset, -0.4))`, oob_keep).
#' A track that sets only `limits` lines up at the N-terminus and drifts by a
#' residue or two by the C-terminus, because the expansion differs -- which is
#' the misalignment this page has had before.
#' @keywords internal
lf_track_x <- function(x_range, seq_offset = -0.45) {
  ggplot2::scale_x_continuous(
    limits = x_range,
    expand = ggplot2::expansion(add = c(seq_offset, -0.4)),
    oob = scales::oob_keep)
}

#' A class-probability track, averaged over the two termini
#'
#' An N window and a C window both have an opinion about every residue, and they
#' mostly agree about what the residue IS -- it is the ends they disagree about.
#' Averaging the two gives one statement per residue per class. The band is the
#' mean of the two across-window sds, not a pooled sd: it is the typical spread
#' within a terminus, not the spread of the average.
#' @keywords internal
lf_plot_class_track <- function(d, x_range, aa, gene, legend = FALSE,
                                seq_offset = -0.45, ylab = NULL,
                                note = "mean of N and C") {
  if (is.null(d) || !nrow(d)) return(NULL)
  ## FORCE the lazily-used arguments. `note` is referenced only inside the
  ## aes() sprintf below, which ggplot evaluates when the plot is BUILT -- so a
  ## caller looping `for (k in seq_along(terms))` and passing
  ## `sprintf(..., terms[k])` had every track's tooltip read the LAST k. The N
  ## class track rendered its points correctly but labelled them "C terminus
  ## windows only". `ylab` escaped it only because `lab <-` forces it here.
  force(note); force(gene); force(aa)
  n_terms <- length(unique(d$term))
  d <- dplyr::summarise(dplyr::group_by(d, residue, class),
                        mean = mean(mean), sd = mean(sd), n_win = sum(n_win),
                        .groups = "drop")
  d <- droplevels(d)
  lab <- if (!is.null(ylab)) ylab else
    if (n_terms > 1) "N+C per-residue\nclass probability" else "per-residue\nclass probability"
  ggplot2::ggplot(d, ggplot2::aes(residue, mean, colour = class)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(mean - sd, 0),
                                      ymax = pmin(mean + sd, 1), fill = class),
                         alpha = 0.16, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(data_id = residue,
                   tooltip = sprintf("%s  residue %d (%s)\n%s %.3f +/- %.3f  (%s)",
                                     gene, residue, aa[residue],
                                     nn_class_labels[as.character(class)], mean, sd, note)),
      ## hollow, the same as the head tracks
      shape = 21, fill = "white", size = 1.7, stroke = 0.6, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = nn_class_cols, labels = nn_class_label_fn,
                                 name = NULL) +
    ggplot2::scale_fill_manual(values = nn_class_cols, guide = "none") +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(y = lab) +
    lf_track_theme() +
    ggplot2::theme(legend.position = if (legend) "top" else "none",
                   legend.direction = "horizontal", legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = grid::unit(8, "pt"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
}

#' One window-level head, both termini on the same track
#'
#' The N and C heads answer different questions -- "does a peptide START here",
#' "does one END here" -- so they are two series, not an average. Colour says
#' which.
#' @keywords internal
lf_plot_win_track <- function(d, head_nm, x_range, aa, gene, legend = FALSE,
                              seq_offset = -0.45, group = NA_character_, ylab = NULL,
                              overlay = NULL, held = NULL, overlay_own = FALSE) {
  ## the page's own series have no group; an extra run's carry one, and each
  ## group gets a track of its own
  if (is.null(d$group)) d$group <- NA_character_
  if (is.null(d$terminus)) d$terminus <- d$term
  d$terminus[is.na(d$terminus)] <- d$term[is.na(d$terminus)]
  keep <- d$head == head_nm & (if (is.na(group)) is.na(d$group) else d$group %in% group)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  ## `overlay`: the all-50 ensemble for this head, drawn INTO this track rather
  ## than as a track of its own -- a heavier line with square points, against
  ## which the fold lines read as spread. Relabelled so it is its own series.
  ## `overlay_own = TRUE`: the overlay is instead THIS group's own ensemble (its
  ## folds pooled), which keeps its series name -- the name is what the popup
  ## resolves a click by, and it must open that ensemble, not the page's own.
  d$ens <- FALSE
  if (!is.null(overlay) && nrow(overlay)) {
    overlay <- overlay[overlay$head == head_nm, , drop = FALSE]
    if (nrow(overlay)) {
      if (is.null(overlay$terminus)) overlay$terminus <- overlay$term
      overlay$terminus[is.na(overlay$terminus)] <- overlay$term[is.na(overlay$terminus)]
      if (!overlay_own) overlay$term <- paste0("all 50 (", overlay$terminus, ")")
      overlay$group <- NA_character_
      overlay$ens   <- TRUE
      d <- dplyr::bind_rows(d, overlay[, intersect(names(overlay), names(d)), drop = FALSE])
    }
  }
  ## colour says WHICH TERMINUS, on every head track alike
  tc  <- getOption("lf.term_cols", LF_TERM_COLS)
  tce <- getOption("lf.term_cols_ens", LF_TERM_COLS_ENS)
  series  <- unique(as.character(d$term))
  term_of <- d$terminus[match(series, d$term)]
  is_ens  <- d$ens[match(series, d$term)]
  ## within a terminus the ensemble sorts first, so it heads its half of the key
  ord <- order(match(term_of, names(tc)), !is_ens, series)
  series <- series[ord]; term_of <- term_of[ord]; is_ens <- is_ens[ord]
  cols <- stats::setNames(unname(ifelse(is_ens, tce[term_of], tc[term_of])), series)
  d$term <- factor(d$term, levels = series)
  ## "N terminus" / "C terminus" for the page's own series; an extra one by name
  lab <- function(x) ifelse(x %in% c("N", "C"), paste0(x, " terminus"), x)
  ## The tooltip is built as a COLUMN, not inside aes(): an aes() expression is
  ## evaluated when the plot is BUILT, so anything it closes over must still be
  ## what it was here -- the same trap that made every class track read the last
  ## loop index (see lf_plot_class_track).
  ## The popup resolves a clicked point to a series by looking for the series
  ## LABEL inside the tooltip text (seriesOf), so the ensemble's tooltip has to
  ## name its terminus: the all-50 per-index surface IS the page's own N/C
  ## series. Without this the square is clickable but inert.
  ser <- lab(as.character(d$term))
  if (any(d$ens) && !overlay_own)
    ser[d$ens] <- sprintf("all 50 ensemble, %s terminus", d$terminus[d$ens])
  ## the head's display name: an extra run's own (`head_label`), else the key
  hl <- if (is.null(d$head_label)) rep(head_nm, nrow(d)) else
    ifelse(is.na(d$head_label), head_nm, d$head_label)
  d$tip <- sprintf("%s  anchor %d (%s)\n%s, %s: %.3f%s",
                   gene, d$residue, aa[d$residue], ser, hl, d$mean,
                   ifelse(is.finite(d$sd), sprintf(" +/- %.3f", d$sd), ""))
  ## A series read from an arm's tracks rds has no per-index surface on the
  ## page, so a click has no window to open. Say so -- and the popup handler
  ## skips any point whose tooltip carries LF_NO_POPUP, which it must: it
  ## matches a point to a series by finding the series label INSIDE the
  ## tooltip, and "ESM-C only C fold 1" contains the curated "C fold 1".
  if (!is.null(d$no_popup))
    d$tip[d$no_popup %in% TRUE] <- paste0(d$tip[d$no_popup %in% TRUE], "\n", LF_NO_POPUP)
  ## A known end of THIS precursor is honest on exactly one fold. That whole
  ## series is then the honest one here, so every one of its residue positions
  ## is drawn as a triangle -- the fold to read is identifiable anywhere along
  ## the track, not only at the end itself. The line and ribbon keep EVERY
  ## non-ensemble row, or the series would break where a triangle replaces a
  ## circle.
  d$held <- FALSE
  if (!is.null(held) && nrow(held)) {
    d$held <- !d$ens & d$term %in% unique(held$term)
    d$tip[d$held] <- paste0(d$tip[d$held], "\nthis fold held out a known end here")
  }
  dline <- d[!d$ens, , drop = FALSE]
  dfold <- d[!d$ens & !d$held, , drop = FALSE]
  dheld <- d[ d$held, , drop = FALSE]
  dens  <- d[ d$ens, , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(residue, mean, colour = term, fill = term)) +
    ggplot2::geom_ribbon(data = dline,
                         ggplot2::aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                         alpha = 0.14, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(data = dline, linewidth = 0.4, alpha = 0.75, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      data = dfold,
      ggplot2::aes(data_id = residue, tooltip = tip),
      ## hollow (a ring in the series colour, white inside) and big enough to
      ## click; the popup fills the clicked one while its window is open
      shape = 21, fill = "white", size = 1.7, stroke = 0.6, na.rm = TRUE) +
    ## the fold that held this end out: a DOWNWARD TRIANGLE, so the one honest
    ## point on the track is findable without reading tooltips
    ggiraph::geom_point_interactive(
      data = dheld,
      ggplot2::aes(data_id = residue, tooltip = tip),
      ## show.legend = FALSE or this layer's shape becomes the glyph for EVERY
      ## key in the colour legend -- all ten fold entries drawn as triangles
      shape = 25, fill = "white", size = 1.8, stroke = 0.7, na.rm = TRUE,
      show.legend = FALSE) +
    ## the ensemble, on top: heavier line and SQUARE points, so it never reads as
    ## just one more fold. White inside like every other point on the page -- the
    ## shape carries the distinction, not the fill.
    ggplot2::geom_line(data = dens, linewidth = 0.95, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      data = dens,
      ggplot2::aes(data_id = residue, tooltip = tip),
      shape = 22, fill = "white", size = 1.55, stroke = 0.7, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, labels = lab) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(y = if (!is.null(ylab)) ylab else sub(" head", "\nhead", head_nm)) +
    lf_track_theme() +
    ggplot2::theme(legend.position = if (legend) "top" else "none",
                   legend.direction = "horizontal",
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = grid::unit(8, "pt"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
}

#' A window-free per-residue model
#' @keywords internal
lf_plot_resid_track <- function(d, mdl, x_range, aa, gene, seq_offset = -0.45) {
  d <- droplevels(d[d$model == mdl, , drop = FALSE])
  if (!nrow(d)) return(NULL)
  ggplot2::ggplot(d, ggplot2::aes(residue, mean, colour = class)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(mean - sd, 0),
                                      ymax = pmin(mean + sd, 1), fill = class),
                         alpha = 0.16, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(data_id = residue,
                   tooltip = sprintf("%s  residue %d (%s)\n%s %s %.3f",
                                     gene, residue, aa[residue], mdl,
                                     nn_class_labels[as.character(class)], mean)),
      ## hollow, the same as the other tracks
      shape = 21, fill = "white", size = 1.7, stroke = 0.6, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = nn_class_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = nn_class_cols, guide = "none") +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(y = sprintf("%s per residue\n(no window)", mdl)) +
    lf_track_theme()
}

#' MultiPep class -> symbol on the candidate track, the same on every page
#'
#' Colour carries the SCORE (turbo, continuous); the symbol carries which class
#' won. MultiPep scores 20 classes; a symbol for each would be unreadable, so
#' five carry the classes that matter and the rest share an "other" cross: the
#' three signalling classes (the ones `mp_signalling_max` pools), growth
#' factors, and antimicrobial + antibacterial sharing one -- both say "looks
#' like an AMP". The tooltip always names the exact class. The five are the
#' FILLABLE shapes, drawn white with a black outline, so they read on both the
#' dark ends and the bright middle of turbo; the line-only cross gets a halo.
#' @export
LF_MP_SHAPES <- c(peptidehormone = 21, neuropeptide = 24, cellcellsignaling = 23,
                  cytokines_growthfactors = 22, antimicrobial = 25, other = 4)

#' MultiPep classes that share another class's symbol
#' @keywords internal
LF_MP_SHARE <- c(antibacterial = "antimicrobial")

#' @keywords internal
LF_MP_LABELS <- c(peptidehormone = "peptide hormone", neuropeptide = "neuropeptide",
                  cellcellsignaling = "cell-cell signalling",
                  cytokines_growthfactors = "cytokine / growth factor",
                  antimicrobial = "antimicrobial", antibacterial = "antibacterial",
                  other = "other class")
## the legend's names, one per symbol
LF_MP_KEY_LABELS <- c(LF_MP_LABELS[setdiff(names(LF_MP_SHAPES), "antimicrobial")],
                      antimicrobial = "antimicrobial / antibacterial")

## arm tokens (as candidate_peptides.R's ARMS_MAX writes them) -> short names
LF_ARM_SHORT <- c(curated = "curated", esmonly36 = "ESM-C", esmonly36x = "ESM-C (covered)",
                  curated_wAAid = "wAAid", curated_AAswap = "AAswap", motifind = "motif-ind")
.lf_arm_names <- function(x) vapply(strsplit(as.character(x), "+", fixed = TRUE), function(a) {
  a <- ifelse(a %in% names(LF_ARM_SHORT), LF_ARM_SHORT[a], a); paste(a, collapse = " + ") }, "")

## tooltip text: ggiraph escapes it itself (escaping here too shows a literal
## "&gt;"), but a single quote inside an attribute is an error rather than an
## escape (see the onclick trap), so a hit name like "5'-nucleotidase" is
## given a typographic apostrophe
.lf_tip_txt <- function(x) gsub("'", "’", as.character(x), fixed = TRUE)

#' The end segments a candidate table was cut from
#'
#' One row per terminus. The bar is the segment's run of anchors at or above
#' the threshold, the open circle its centre (the run's maximum -- the anchor
#' the candidates are built around), the thin line the +/-TOL fan-out: every
#' residue that segment offers as a peptide's first (N) or last (C) residue. A
#' segment with no partner in the length range is drawn pale: it is real, but
#' it is in no candidate.
#' @keywords internal
lf_plot_segment_track <- function(cd, x_range, aa, gene, seq_offset = -0.45) {
  s <- cd$seg
  if (is.null(s) || !nrow(s)) return(NULL)
  force(gene); force(aa)
  tc  <- getOption("lf.term_cols", LF_TERM_COLS)
  tol <- cd$tol; thr <- cd$thresh
  s$y <- ifelse(s$term == "N", 2, 1)
  s$paired <- factor(ifelse(s$n_candidates > 0, "paired", "unpaired"),
                     levels = c("paired", "unpaired"))
  has_fan <- is.finite(tol)
  s$fan_lo <- if (has_fan) pmax(1, s$centre - tol) else s$centre
  s$fan_hi <- if (has_fan) pmin(x_range[2] - 1, s$centre + tol) else s$centre
  role <- ifelse(s$term == "N", "first", "last")
  s$tip <- paste0(
    sprintf("%s  %s end segment %d-%d (%s): %d anchor%s%s\n",
            gene, s$term, s$start, s$end, cd$label, s$n_anchors,
            ifelse(s$n_anchors == 1, "", "s"),
            if (is.finite(thr)) paste(" at or above", format(thr)) else ""),
    sprintf("centre %d (%s), score %.3f, %s motif\n",
            s$centre, aa[s$centre], s$centre_score, s$motif_centre),
    if (has_fan) sprintf("offers residues %d-%d as a peptide's %s residue (+/-%d)\n",
                         s$fan_lo, s$fan_hi, role, tol) else "",
    ifelse(s$n_candidates > 0, sprintf("in %d candidate peptide%s", s$n_candidates,
                                       ifelse(s$n_candidates == 1, "", "s")),
           "no partner segment in the length range: in no candidate"))
  ## a segment cut from several arms' max (candidate_peptides.R ARMS_MAX) says
  ## which of them reach the threshold inside it
  has_arms <- !is.null(s$n_arms)
  if (has_arms)
    s$tip <- paste0(s$tip, sprintf("\nsupported by %d arm%s: %s; max at the centre from %s",
                                   s$n_arms, ifelse(s$n_arms == 1, "", "s"),
                                   .lf_arm_names(s$arms), .lf_arm_names(s$centre_arm)))
  s$tip <- .lf_tip_txt(s$tip)
  ## the bar and its fan-out light up together; the centre keeps the bare
  ## residue number so hovering it marks that anchor on every other track
  s$sid <- sprintf("seg_%s_%d", s$term, s$centre)
  ggplot2::ggplot(s) +
    ggiraph::geom_segment_interactive(
      ggplot2::aes(x = fan_lo, xend = fan_hi, y = y, yend = y, colour = term,
                   data_id = sid, tooltip = tip),
      linewidth = 0.5) +
    ggiraph::geom_rect_interactive(
      ggplot2::aes(xmin = start - 0.45, xmax = end + 0.45, ymin = y - 0.26, ymax = y + 0.26,
                   fill = term, colour = term, alpha = paired, data_id = sid, tooltip = tip),
      linewidth = 0.4) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(centre, y, colour = term, data_id = centre, tooltip = tip),
      shape = 21, fill = "white", size = 1.7, stroke = 0.6) +
    ## how many arms the segment rests on, over its bar
    (if (has_arms) ggplot2::geom_text(ggplot2::aes(x = (start + end) / 2, y = y + 0.42, label = n_arms),
                                      size = 2.3, colour = "#0b0b0b", fontface = "bold")) +
    ggplot2::scale_colour_manual(values = tc, guide = "none") +
    ggplot2::scale_fill_manual(values = tc, guide = "none") +
    ggplot2::scale_alpha_manual(
      values = c(paired = 0.85, unpaired = 0.15), drop = FALSE, name = cd$label,
      labels = c(paired = "segment in a candidate",
                 unpaired = "segment with no partner in the length range"),
      ## outlined keys: the pale one is otherwise white on white
      guide = ggplot2::guide_legend(override.aes = list(colour = "grey35"))) +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("C", "N"),
                                expand = c(0, 0)) +
    ## the set's label rides in the legend title: the axis title of a track
    ## this short has room for two short words, and is clipped past that
    ggplot2::labs(y = if (has_arms) "end\nsegments\n(digit = arms)" else "end\nsegments") +
    lf_track_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   legend.position = "top", legend.direction = "horizontal",
                   legend.justification.top = "left",
                   legend.title = ggplot2::element_text(size = 7, face = "bold"),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = grid::unit(8, "pt"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
}

#' The candidate peptides, one bar each, coloured by their best MultiPep call
#'
#' Colour = the best MultiPep score (turbo, 0..1); the symbol near the bar's N
#' end = which class it was (LF_MP_SHAPES). Every candidate is its own row -- a group's
#' candidates all overlap, since they are one pair of segments cut at different
#' offsets -- packed first-fit so that groups that do not overlap share rows.
#' Groups go in order of their centres and, within a group, the strongest
#' MultiPep call first, so the top of each group's block is its best guess.
#'
#' Returns list(p, h): the track height grows with the number of rows,
#' `options(lf.cand_row_height)` inches each -- 0.065 by default, the least at
#' which a row is tall enough to carry its symbol.
#' @keywords internal
lf_plot_candidate_track <- function(cd, x_range, aa, gene, seq_offset = -0.45,
                                    row_in = getOption("lf.cand_row_height", 0.065)) {
  d <- cd$cand
  if (is.null(d) || !nrow(d)) return(NULL)
  force(gene); force(aa)
  d <- d[order(d$n_centre, d$c_centre, -d$best_score, d$start, d$end), , drop = FALSE]
  ## first-fit packing; one clear residue between neighbours on a row
  last <- integer(0); d$row <- 0L
  for (i in seq_len(nrow(d))) {
    j <- which(last < d$start[i] - 1L)[1]
    if (is.na(j)) { last <- c(last, d$end[i]); j <- length(last) } else last[j] <- d$end[i]
    d$row[i] <- j
  }
  n_rows <- length(last)
  key <- ifelse(d$best_type %in% names(LF_MP_SHARE), LF_MP_SHARE[d$best_type], d$best_type)
  d$cls <- factor(ifelse(key %in% names(LF_MP_SHAPES), key, "other"), levels = names(LF_MP_SHAPES))
  ## the symbol sits inside the bar, one residue in from its N end
  ## under ARMS_MAX each end carries how many arms its segment rests on: drawn
  ## as a cap on that end of the bar, so the class symbol steps one residue in
  has_arms <- all(c("n_es_arms", "c_es_arms") %in% names(d))
  d$sym_x <- d$start + if (has_arms) 1.6 else 1

  ## ---- tooltip ----------------------------------------------------------------
  col <- function(n) if (n %in% names(d)) d[[n]] else rep(NA, nrow(d))
  has <- function(x) !is.na(x) & nzchar(as.character(x)) & as.character(x) != "none"
  src <- c(ligand_list = "curated ligand", peptide = "UniProt peptide",
           propeptide = "UniProt propeptide", none = "none")
  lbl <- function(t) ifelse(t %in% names(LF_MP_LABELS), LF_MP_LABELS[t], t)
  mp <- lf_mp_classes(names(d))
  runner_up <- if (length(mp) >= 3) {
    m <- as.matrix(d[, mp, drop = FALSE])
    apply(m, 1, function(v) {
      o <- order(v, decreasing = TRUE)[2:3]
      paste(sprintf("%s %.2f", lbl(sub("^mp_", "", mp[o])), v[o]), collapse = ", ")
    })
  } else rep(NA_character_, nrow(d))
  gsz <- stats::ave(d$group_rank, d$group_id, FUN = length)
  ll  <- col("ll_hit_name"); ov <- col("overlap_name"); lr <- col("lr_prob")
  bh  <- col("mmseqs_best_hit_name")
  d$tip <- paste0(
    sprintf("%s  residues %d-%d (%d aa)\n%s\n", d$peptide_id, d$start, d$end, d$length, d$sequence),
    ifelse(is.na(d$best_score), "no MultiPep score\n",
           sprintf("best MultiPep: %s %.3f%s\n", lbl(d$best_type), d$best_score,
                   ifelse(is.na(runner_up), "", paste0("  (then ", runner_up, ")")))),
    sprintf("N end %d (score %.3f, %s, %+d from centre) | C end %d (%.3f, %s, %+d)\n",
            d$start, d$n_score_end, col("n_motif_end"), d$n_dist,
            d$end, d$c_score_end, col("c_motif_end"), d$c_dist),
    if (has_arms) sprintf("N end segment: %d arm%s (%s) | C end segment: %d arm%s (%s)\n",
                          d$n_es_arms, ifelse(d$n_es_arms == 1, "", "s"), .lf_arm_names(d$n_es_arm_list),
                          d$c_es_arms, ifelse(d$c_es_arms == 1, "", "s"), .lf_arm_names(d$c_es_arm_list))
    else "",
    sprintf("group %s: rank %d of %d by end-score product %.3f",
            d$group_id, d$group_rank, gsz, d$pair_score),
    ifelse(is.na(lr), "", sprintf("\nranker lr_prob %.3f", lr)),
    ## the closest known peptide of any source (the target set is cut to
    ## peptide lengths in stage 2a), then the closest CURATED one when that is a
    ## different peptide -- HCRT's curated list has orexin-B but not orexin-A
    ifelse(has(bh), sprintf("\nbest known-peptide match: %s (%s, %s%% identity, covers %s%%, e-value %s)",
                            bh, src[as.character(col("mmseqs_hit_type"))], round(col("mmseqs_pident")),
                            round(100 * col("mmseqs_tcov")), signif(col("mmseqs_evalue"), 2)),
           "\nno known-peptide match"),
    ifelse(has(ll) & col("mmseqs_hit_type") != "ligand_list",
           sprintf("\nclosest curated ligand: %s (%s%% identity, covers %s%%)", ll,
                   round(col("ll_pident")), round(100 * col("ll_tcov"))), ""),
    ifelse(has(ov), sprintf("\noverlaps annotated: %s%s", ov,
                            ifelse(col("overlap_exact") %in% TRUE, " (exact range)", "")), ""))
  d$tip <- .lf_tip_txt(d$tip)

  present <- levels(droplevels(d$cls))
  n_grp <- length(unique(d$group_id))
  ## at least an inch, or the axis title outruns a short pile and is clipped
  panel_h <- max(n_rows * row_in, 1)
  oth <- d[d$cls == "other", , drop = FALSE]
  caps <- NULL
  if (has_arms) {
    K <- max(3L, d$n_es_arms, d$c_es_arms, na.rm = TRUE)
    ramp <- grDevices::colorRampPalette(c("#ffffff", "#0b0b0b"))(K)
    lv <- sprintf("%d arm%s", seq_len(K), ifelse(seq_len(K) == 1, "", "s"))
    caps <- rbind(data.frame(x = d$start, row = d$row, k = d$n_es_arms, id = d$peptide_id, tip = d$tip),
                  data.frame(x = d$end,   row = d$row, k = d$c_es_arms, id = d$peptide_id, tip = d$tip))
    caps$lab <- factor(lv[caps$k], levels = lv)
  }
  p <- ggplot2::ggplot(d) +
    ggiraph::geom_rect_interactive(
      ggplot2::aes(xmin = start - 0.45, xmax = end + 0.45, ymin = row - 0.42, ymax = row + 0.42,
                   fill = best_score, data_id = peptide_id, tooltip = tip),
      colour = NA) +
    ## arm-support caps, one per end. Fill is a per-row constant (the track's
    ## fill scale is the turbo score), and the legend comes from `colour`
    ## with the ramp put back into its keys.
    (if (has_arms) ggiraph::geom_rect_interactive(
      data = caps, ggplot2::aes(xmin = x - 0.45, xmax = x + 0.45, ymin = row - 0.42, ymax = row + 0.42,
                                colour = lab, data_id = id, tooltip = tip),
      fill = ramp[caps$k], linewidth = 0.2)) +
    (if (has_arms) ggplot2::scale_colour_manual(
      values = stats::setNames(rep("#0b0b0b", K), lv), drop = FALSE, name = "end segment rests on:",
      guide = ggplot2::guide_legend(order = 3, override.aes = list(fill = ramp)))) +
    ## the cross is line-only: a white stroke under a black one is its halo
    ggplot2::geom_point(data = oth, ggplot2::aes(sym_x, row), shape = 4, size = 1.5,
                        stroke = 1.3, colour = "white") +
    ggiraph::geom_point_interactive(
      ggplot2::aes(sym_x, row, shape = cls, data_id = peptide_id, tooltip = tip),
      fill = "white", colour = "black", size = 1.5, stroke = 0.35) +
    ggplot2::scale_fill_viridis_c(option = "turbo", limits = c(0, 1),
                                  breaks = c(0, 0.25, 0.5, 0.75, 1), na.value = "grey70",
                                  name = paste0(cd$label, "   best MultiPep score:")) +
    ggplot2::scale_shape_manual(values = LF_MP_SHAPES, breaks = present,
                                labels = LF_MP_KEY_LABELS[present], name = "its class:") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(order = 1, barwidth = grid::unit(110, "pt"),
                                                    barheight = grid::unit(7, "pt")),
                    shape = ggplot2::guide_legend(order = 2,
                                                  override.aes = list(size = 2.2, stroke = 0.5))) +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_reverse(limits = c(n_rows + 0.6, 0.4), expand = c(0, 0)) +
    ggplot2::labs(y = sprintf("candidates\n%d in %d group%s", nrow(d), n_grp,
                              if (n_grp == 1) "" else "s")) +
    lf_track_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   legend.position = "top", legend.direction = "horizontal",
                   legend.box = "horizontal", legend.justification.top = "left",
                   legend.title = ggplot2::element_text(size = 7, face = "bold"),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = grid::unit(8, "pt"),
                   legend.margin = ggplot2::margin(0, 12, 0, 0))
  list(p = p, h = panel_h + getOption("lf.track_legend_height", 0.27) + 0.05)
}

#' Every track that could be built, top to bottom, with its height in inches
#' @export
#' One height for every line track
#'
#' `options(lf.track_height)` (inches, default 0.77 -- the plotting height the
#' CNN head track had when a two-row legend squeezed it, which is the size Kevin
#' asked every track to match) is a track with no legend; a track that carries
#' its legend INSIDE gets `lf.track_legend_height` (default 0.27) on top, so its
#' plotting area matches every other track's.
#' @keywords internal
lf_track_h <- function(legend = FALSE)
  getOption("lf.track_height", 0.77) +
    if (isTRUE(legend)) getOption("lf.track_legend_height", 0.27) else 0

#' The key for the window-head tracks, drawn once above them
#'
#' Replaces the per-series legend (N fold 1, C fold 1, ... all 50 (N)) with what
#' the marks MEAN: colour = terminus, open dot = a fold that trained on this
#' protein's ends, downward triangle = the fold that held one out, square = the
#' mean of all folds. Built from off-panel keys so only the legend shows, on the
#' tracks' own residue axis so the page aligns it like any other track.
#' @keywords internal
## `validation`: does this precursor carry a known end that some fold held out?
## Only then do the triangles exist, and only then does "train" mean anything --
## on a precursor no fold ever trained on, every fold series is out-of-sample, so
## the open dot is just that fold's own ensemble and is named as such.
lf_plot_win_legend <- function(x_range, seq_offset = -0.45, validation = TRUE) {
  tc <- getOption("lf.term_cols", LF_TERM_COLS)
  far <- x_range[1] - 1e6                    # off the panel: keys only
  ln <- data.frame(x = far, y = c(0.4, 0.6, 0.4, 0.6),
                   key = factor(rep(c("N terminus", "C terminus"), each = 2),
                                levels = c("N terminus", "C terminus")))
  keys <- if (validation) c("train", "validation", "mean of all folds")
          else c("individual fold (mean of 10 seeds)", "mean of all folds")
  pt <- data.frame(x = far, y = 0.5, key = factor(keys, levels = keys))
  ggplot2::ggplot() +
    ggplot2::geom_line(data = ln, ggplot2::aes(x, y, colour = key, group = key), linewidth = 0.7) +
    ggplot2::geom_point(data = pt, ggplot2::aes(x, y, shape = key),
                        fill = "white", colour = "grey20", size = 1.9, stroke = 0.6) +
    ggplot2::scale_colour_manual(values = c(`N terminus` = tc[["N"]], `C terminus` = tc[["C"]]),
                                 name = NULL) +
    ggplot2::scale_shape_manual(
      values = stats::setNames(if (validation) c(21, 25, 22) else c(21, 22), keys),
      name = NULL) +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    lf_track_theme() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "inside",
                   legend.position.inside = c(0, 0.5),
                   legend.justification = c(0, 0.5),
                   legend.direction = "horizontal", legend.box = "horizontal",
                   legend.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 7.5),
                   legend.key.size = grid::unit(10, "pt"),
                   legend.spacing.x = grid::unit(6, "pt"),
                   legend.margin = ggplot2::margin(0, 12, 0, 0))
}

lf_pepend_track_plots <- function(dat, x_range, aa, gene, seq_offset = -0.45) {
  mk <- list()
  add <- function(p, h) if (!is.null(p)) mk[[length(mk) + 1L]] <<- list(p = p, h = h)
  ## options(lf.class_tracks_by_term = TRUE): one class track per terminus
  ## instead of their average -- the N and C windows mostly agree about what a
  ## residue IS, but where they disagree is exactly what the average hides
  terms <- if (is.null(dat$resid)) character(0) else
    intersect(c("N", "C"), unique(as.character(dat$resid$term)))
  if (isTRUE(getOption("lf.class_tracks_by_term", FALSE)) && length(terms)) {
    for (k in seq_along(terms))
      add(lf_plot_class_track(dat$resid[dat$resid$term == terms[k], , drop = FALSE],
                              x_range, aa, gene, k == 1L, seq_offset,
                              ylab = sprintf("%s terminus windows\nper-residue\nclass probability", terms[k]),
                              note = sprintf("%s terminus windows only", terms[k])),
          lf_track_h(legend = k == 1L))
  } else {
    add(lf_plot_class_track(dat$resid, x_range, aa, gene, TRUE, seq_offset), lf_track_h(TRUE))
  }
  for (nm in names(dat$extra_resid))
    add(lf_plot_class_track(dat$extra_resid[[nm]], x_range, aa, gene, FALSE, seq_offset,
                            ylab = paste0(nm, "\nper-residue\nclass probability"),
                            note = paste(nm, "windows only")), lf_track_h())
  ## the window-free per-residue models first: they answer a different question
  ## from the window heads and read as the page's baseline above them
  if (!is.null(dat$rp)) {
    add(lf_plot_resid_track(dat$rp, "MLP", x_range, aa, gene, seq_offset), lf_track_h())
    add(lf_plot_resid_track(dat$rp, "XGB", x_range, aa, gene, seq_offset), lf_track_h())
  }
  if (!is.null(dat$win)) {
    ## each head as the page's own track, then one per extra group (e.g. padded)
    grp <- if (is.null(dat$win$group)) character(0) else unique(stats::na.omit(dat$win$group))
    ## options(lf.win_overlay_all50 = TRUE): the page's own series (group NA --
    ## the all-50 ensemble) stops being a track of its own and is drawn INTO each
    ## grouped track instead. A fold is then read against the ensemble on one
    ## pair of axes rather than by eye across two tracks.
    ov <- isTRUE(getOption("lf.win_overlay_all50", FALSE)) && length(grp)
    own <- dat$win[is.na(dat$win$group), , drop = FALSE]
    ## Per-fold tracks get ONE key above them (lf_plot_win_legend) instead of a
    ## per-series legend inside the first; options(lf.win_custom_legend = FALSE)
    ## restores the per-series one.
    custom <- length(grp) > 0 && isTRUE(getOption("lf.win_custom_legend", TRUE))
    has_val <- !is.null(dat$heldout) && nrow(dat$heldout) > 0
    if (custom) add(lf_plot_win_legend(x_range, seq_offset, validation = has_val),
                    getOption("lf.track_legend_height", 0.27))
    for (hd in c("CNN attention head", "xgb window head"))
      if (!ov)
        add(lf_plot_win_track(dat$win, hd, x_range, aa, gene, hd == "CNN attention head", seq_offset),
            lf_track_h(hd == "CNN attention head"))
    ## GROUP-major: each arm's CNN and xgb tracks sit together, in the order the
    ## arms were given (curated first). Head-major interleaves them, which put
    ## one arm's tracks between another's two.
    for (g in grp) {
      for (hd in c("CNN attention head", "xgb window head")) {
        ## the axis names the head the group's runs actually have, when they
        ## agree on one (a flatten trunk's CNN track is not attention)
        hl <- if (is.null(dat$win$head_label)) character(0) else
          unique(stats::na.omit(dat$win$head_label[dat$win$head == hd & dat$win$group %in% g]))
        lbl <- if (length(hl) == 1L) hl else hd
        ## An arm's head_label already names the arm ("CNN linear head (ESM-C
        ## only)") because the TOOLTIP has no group to carry it; the axis does,
        ## so strip it here rather than say it twice. Fixed-string, not regex: an
        ## arm label may contain "+" or "(".
        .sfx <- paste0(" (", g, ")")
        if (nzchar(g) && endsWith(lbl, .sfx)) lbl <- substr(lbl, 1L, nchar(lbl) - nchar(.sfx))
        ## A group that carries its OWN ensemble (spec$group_ensemble: its folds
        ## pooled) overlays that, not the page's all-50: the heavy line on a
        ## track must be the average of the folds drawn on it. Only a group
        ## without one -- the 5-fold per-fold group -- gets the page's own.
        ge <- if (is.null(dat$win$group_ens)) rep(FALSE, nrow(dat$win)) else
          dat$win$group %in% g & dat$win$group_ens %in% TRUE
        own_g <- any(ge)
        leg <- !custom && hd == "CNN attention head"
        add(lf_plot_win_track(dat$win[!ge, , drop = FALSE], hd, x_range, aa, gene,
                              leg, seq_offset,
                              group = g, ylab = sprintf("%s\n(%s)", sub(" head", "\nhead", lbl), g),
                              overlay = if (own_g) dat$win[ge, , drop = FALSE] else if (ov) own else NULL,
                              overlay_own = own_g,
                              held = dat$heldout), lf_track_h(leg))
      }
    }
  }
  ## per candidate set: the end segments it was cut from, then its candidates --
  ## last, so they sit directly above the sequence they span
  for (cd in dat$cand) {
    add(lf_plot_segment_track(cd, x_range, aa, gene, seq_offset),
        getOption("lf.seg_track_height", 0.5) + getOption("lf.track_legend_height", 0.27))
    cp <- lf_plot_candidate_track(cd, x_range, aa, gene, seq_offset)
    if (!is.null(cp)) add(cp$p, cp$h)
  }
  mk
}

## ---- the click popup -------------------------------------------------------------
## Every mature residue is a window anchor, so instead of pre-rendering a panel
## per window the page carries all of their per-index surfaces as data and draws
## the one you ask for. Clicking a residue on ANY track opens the N and C windows
## anchored there, side by side.
##
## The click is bound to `data-id`, which is the bare residue number -- it has to
## be, because that is what links hover across tracks -- so a click cannot say
## which track it came from. Showing both termini is therefore not a compromise:
## it is the only well-defined answer, and the more useful one.

#' The popup payload and handler, as a <script> tag
#'
#' @param pix the `pix` element of lf_pepend_track_data()
#' @param max_anchors keep only this many anchors, highest window score first.
#'   36 x K floats per anchor per terminus is ~3 KB; a 1,600-residue precursor
#'   would otherwise add 10 MB to the page.
#' @export
lf_pepend_popup_tag <- function(pix, aa, gene, accession,
                                max_anchors = getOption("lf.pepend_popup_max", 600L)) {
  if (is.null(pix) || is.null(pix$anchors)) return(NULL)
  terms <- intersect(c("N", "C"), names(pix))
  if (!length(terms)) return(NULL)

  ## `label` is how the series is named in the head tracks' tooltips, which is
  ## how a click finds it; `title` heads the opened panel
  pack <- function(S, anchors, label, title, term) {
    keep <- seq_along(anchors)
    if (length(keep) > max_anchors)
      keep <- sort(utils::head(order(S$g, decreasing = TRUE), max_anchors))
    list(label = label, title = title, term = term,
         ## an extra run's own label, else the page's (options(lf.track_cnn_label):
         ## what the page's own trunk's window head is -- "CNN flatten head" for
         ## the flatten models)
         cnn_label = if (is.null(S$cnn_label))
           getOption("lf.track_cnn_label", "CNN attention head") else S$cnn_label,
         anchor = as.integer(anchors[keep]),
         pos = as.integer(S$pos),
         vis = as.integer(if (is.null(S$vis)) c(1L, LF_PEPEND$seq_len) else S$vis),
         g = round(S$g[keep], 4),
         pi = lapply(keep, function(i) round(S$pi[i, , ], 3)),
         sd = lapply(keep, function(i) round(S$pi_sd[i, , ], 3)))
  }
  payload <- list(gene = gene, accession = accession, aa = aa,
                  seq_len = as.integer(LF_PEPEND$seq_len),
                  classes = pix[[terms[1]]]$pi_names,
                  labels = unname(nn_class_labels[pix[[terms[1]]]$pi_names]),
                  cols = unname(nn_class_cols[pix[[terms[1]]]$pi_names]))
  series <- list()
  for (t in terms)
    series[[t]] <- pack(pix[[t]], pix$anchors, paste0(t, " terminus"), paste0(t, " window"), t)
  for (nm in names(pix$extra)) {
    S <- pix$extra[[nm]]
    if (!identical(as.character(S$pi_names), as.character(pix[[terms[1]]]$pi_names))) next
    series[[nm]] <- pack(S, S$anchors, nm, sub("^([NC]) terminus", "\\1 window", nm), S$term)
  }
  payload$series <- series

  htmltools::tags$script(htmltools::HTML(sprintf(
    "window.LF_PEP = %s;\n%s",
    jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 4, null = "null"),
    lf_pepend_popup_js())))
}

#' @keywords internal
lf_pepend_popup_js <- function() '
(function () {
  var d = window.LF_PEP;
  if (!d) return;
  /* one open window per track per terminus, keyed track id + "|" + N or C */
  var PANEL_H = 150, slots = {};

  /* The clickable track is the xgb window head, and only that one: it is the
     head that actually localises a boundary, so it is the one worth opening a
     window from. Found by its y-axis label rather than by index, since a track
     is dropped when its input is missing. */
  /* every xgb track -- the page own one and any extra run (padded) -- since
     each opens the windows of the series drawn on it */
  function xgbSvgs() {
    var all = [].slice.call(document.querySelectorAll("svg.ggiraph-svg")), out = [];
    for (var i = 0; i < all.length; i++) {
      var t = all[i].querySelectorAll("text");
      for (var j = 0; j < t.length; j++)
        if (/xgb window/.test(t[j].textContent)) { out.push(all[i]); break; }
    }
    return out;
  }

  /* residue -> client x, read off the track itself. Two points far apart give
     the linear map, so the opened window lands on the SAME residue positions as
     every other track rather than on its own 36-wide axis. */
  function residueScale(svg) {
    var pts = [].slice.call(svg.querySelectorAll("circle[data-id]"))
      .map(function (c) { var v = c.getAttribute("data-id");
                          return /^[0-9]+$/.test(v) ? { r: +v, el: c } : null; })
      .filter(Boolean);
    if (pts.length < 2) return null;
    pts.sort(function (a, b) { return a.r - b.r; });
    var lo = pts[0], hi = pts[pts.length - 1];
    var xa = lo.el.getBoundingClientRect(), xb = hi.el.getBoundingClientRect();
    var ax = xa.left + xa.width / 2, bx = xb.left + xb.width / 2;
    if (hi.r === lo.r) return null;
    var b = (bx - ax) / (hi.r - lo.r);
    return { x: function (r) { return ax + (r - lo.r) * b; }, step: b };
  }

  function slotKey(svg, term) { return (svg.id || "svg") + "|" + term; }

  /* scale an SVG radius, keeping its unit -- ggiraph writes r as e.g. 0.99pt,
     and a bare number would be read as px */
  function scaleR(r, f) {
    var m = /^([0-9.]+)([a-z%]*)$/.exec(String(r || "2"));
    return m ? (parseFloat(m[1]) * f).toFixed(2) + m[2] : r;
  }

  /* the clicked dot stays filled and enlarged while its window is open */
  function selectDot(el) {
    var keep = { el: el, fill: el.getAttribute("fill"), r: el.getAttribute("r"),
                 sw: el.getAttribute("stroke-width") };
    el.setAttribute("fill", el.getAttribute("stroke") || "#000");
    el.setAttribute("r", scaleR(keep.r, 1.6));
    el.setAttribute("stroke-width", "1.4");
    if (el.parentNode) el.parentNode.appendChild(el);   /* draw it on top */
    return keep;
  }
  function unselectDot(keep) {
    if (!keep) return;
    ["fill", "r", "stroke-width"].forEach(function (a, i) {
      var v = [keep.fill, keep.r, keep.sw][i];
      if (v === null) keep.el.removeAttribute(a); else keep.el.setAttribute(a, v);
    });
  }

  function closeSlot(k) {
    var o = slots[k];
    if (!o) return;
    var wrap = o.host.parentNode;
    unselectDot(o.dot);
    o.host.remove(); delete slots[k];
    if (wrap && !wrap.children.length) wrap.remove();
    /* a highlight left on another track would survive the panel */
    document.querySelectorAll("svg.ggiraph-svg [data-id]").forEach(function (n) {
      if (n.__lfhl) { n.removeAttribute("stroke"); delete n.__lfhl; }
    });
  }
  function closeAll() { Object.keys(slots).forEach(closeSlot); }

  /* Every window opened from one track sits in one stack directly under it.
     The panel must live in the SAME coordinate space as the track, not the
     page. These svgs are drawn at natural size -- NTN1 is about 86 inches --
     so the girafe container scrolls horizontally. A width:100% div outside it
     takes the PAGE width: every x computed from the track then lands past its
     right edge and is clipped. So the stack goes INSIDE the container,
     directly after the svg, at exactly the width of that svg. And ggiraph
     CENTRES a chart narrower than the page (.girafe_container_std is
     text-align:center) while a block div sits at the left, so on a short
     protein the panel started left of its track and was cut off on the right;
     the stack is lined up with the svg itself.
     NOTE: no apostrophes anywhere in this handler -- it is an R single-quoted
     string, and one would end it. */
  function stackFor(svg) {
    var nx = svg.nextSibling;
    if (nx && nx.classList && nx.classList.contains("lf-pep-stack")) return nx;
    var svgRect = svg.getBoundingClientRect();
    var w = document.createElement("div");
    w.className = "lf-pep-stack";
    w.style.cssText = "position:relative;margin:0;text-align:left;width:" +
      svgRect.width.toFixed(0) + "px;";
    svg.parentNode.insertBefore(w, svg.nextSibling);
    var dx = svgRect.left - w.getBoundingClientRect().left;
    if (Math.abs(dx) > 0.5) w.style.marginLeft = dx.toFixed(1) + "px";
    return w;
  }

  /* the head and value the clicked dot carries, from its own tooltip: a click
     comes from an xgb track, while the payload score is the CNN global head,
     so a header naming only one of them reads as wrong on the other track */
  function clickedScore(dot) {
    var t = dot && dot.getAttribute("title");
    var m = t && /([A-Za-z]+ [A-Za-z]+ head): ([0-9.]+)/.exec(t.replace(/&[a-z]+;/g, " "));
    return m ? { head: m[1], val: m[2] } : null;
  }

  function open(key, res, svg, dot) {
    var S = d.series && d.series[key];
    if (!S) return;
    var i = S.anchor.indexOf(res);
    if (i < 0) return;
    var sc = residueScale(svg);
    if (!sc) return;
    var term = S.term || key, k = slotKey(svg, term);
    closeSlot(k);

    var svgRect = svg.getBoundingClientRect();
    var stack = stackFor(svg);
    var host = document.createElement("div");
    host.className = "lf-pep-window";
    host.setAttribute("data-term", term);
    host.style.cssText = "position:relative;margin:0;width:100%;";
    /* N above C, whichever was opened first */
    var before = null;
    if (term === "N")
      [].slice.call(stack.children).forEach(function (c) {
        if (!before && c.getAttribute("data-term") === "C") before = c;
      });
    stack.insertBefore(host, before);
    slots[k] = { host: host, key: key, res: res, dot: dot ? selectDot(dot) : null };
    var hostRect = host.getBoundingClientRect();
    var N = d.seq_len, w0 = res - S.pos + 1, K = d.classes.length;
    /* positions the model saw (1-based, inclusive): all 36, 1-28 for a
       CT-padded window, 8-36 for an NT-padded one */
    var V0 = (S.vis && S.vis[0]) || 1, V1 = (S.vis && S.vis[1]) || N;
    var X = function (r) { return sc.x(r) - hostRect.left; };
    var Y = function (v) { return PANEL_H - 26 - v * (PANEL_H - 44); };

    /* overflow visible: a window anchored near either end of the protein runs
       past the track axis, and those positions should show, not be clipped */
    var g = "<svg width=\\"" + svgRect.width.toFixed(0) + "\\" height=\\"" + PANEL_H +
            "\\" style=\\"display:block;overflow:visible\\">";
    /* the window frame, on the residue axis */
    g += "<rect x=" + X(w0 - 0.5).toFixed(1) + " y=6 width=" +
         Math.max(1, (X(w0 + N - 0.5) - X(w0 - 0.5))).toFixed(1) +
         " height=" + (PANEL_H - 32) + " fill=#ffffff stroke=#d9d8d3 />";
    /* a padded window never saw positions outside V0-V1 -- say so, rather
       than draw predictions about residues it was not shown */
    function padBlock(p0, p1) {
      var px0 = X(w0 + p0 - 1.5), px1 = X(w0 + p1 - 0.5);
      g += "<rect x=" + px0.toFixed(1) + " y=6 width=" + Math.max(1, px1 - px0).toFixed(1) +
           " height=" + (PANEL_H - 32) + " fill=#ecebe6 />";
      g += "<text x=" + ((px0 + px1) / 2).toFixed(1) + " y=" + (PANEL_H / 2).toFixed(0) +
           " text-anchor=middle font-size=10 fill=#8a8a8a>padding</text>";
    }
    if (V0 > 1) padBlock(1, V0 - 1);
    if (V1 < N) padBlock(V1 + 1, N);
    /* the anchor: the residue this window is ABOUT */
    g += "<rect x=" + (X(res) - 1).toFixed(1) + " y=6 width=2 height=" +
         (PANEL_H - 32) + " fill=#000 opacity=0.35 />";
    for (var k = 0; k < K; k++) {
      if (d.classes[k] === "none" || d.classes[k] === "padding") continue;
      var pts = [];
      for (var j = V0 - 1; j < V1; j++)
        pts.push(X(w0 + j).toFixed(1) + "," + Y(S.pi[i][j][k]).toFixed(1));
      /* +/- 1 sd across the members, as a band under the line */
      if (S.sd) {
        var up = [], dn = [];
        for (var j3 = V0 - 1; j3 < V1; j3++) {
          var m3 = S.pi[i][j3][k], e3 = S.sd[i][j3][k];
          up.push(X(w0 + j3).toFixed(1) + "," + Y(Math.min(1, m3 + e3)).toFixed(1));
          dn.push(X(w0 + j3).toFixed(1) + "," + Y(Math.max(0, m3 - e3)).toFixed(1));
        }
        dn.reverse();
        g += "<polygon points=\\"" + up.concat(dn).join(" ") + "\\" fill=\\"" +
             d.cols[k] + "\\" fill-opacity=0.16 stroke=none />";
      }
      g += "<polyline points=\\"" + pts.join(" ") + "\\" fill=none stroke=\\"" +
           d.cols[k] + "\\" stroke-width=1.6 />";
      /* a point per residue per class, carrying the residue it belongs to so a
         hover can light that column up on every other track */
      for (var j2 = V0 - 1; j2 < V1; j2++) {
        var rr = w0 + j2, vv = S.pi[i][j2][k];
        g += "<circle cx=" + X(rr).toFixed(1) + " cy=" + Y(vv).toFixed(1) +
             " r=2 fill=\\"" + d.cols[k] + "\\" data-res=" + rr +
             " style=\\"cursor:crosshair\\"><title>residue " + rr + " (" +
             (d.aa[rr - 1] || ".") + ")  " + d.labels[k] + " " + vv.toFixed(3) +
             (S.sd ? " +/- " + S.sd[i][j2][k].toFixed(3) : "") + "</title></circle>";
      }
    }
    g += "<line x1=0 y1=" + (PANEL_H - 26) + " x2=" + svgRect.width.toFixed(0) + " y2=" + (PANEL_H - 26) +
         " stroke=#ccc />";
    g += "</svg>";

    /* centred over the window it describes, not pinned to the left margin --
       the title belongs to the 36-mer, which can sit anywhere along the protein */
    var cx = (X(w0) + X(w0 + N - 1)) / 2;
    var cs = clickedScore(dot);
    var hdr = "<div class=\\"lf-pep-hdr\\" style=\\"position:absolute;top:4px;left:" +
      cx.toFixed(0) + "px;transform:translateX(-50%);white-space:nowrap;text-align:center;" +
      "font:11px ui-sans-serif,system-ui;color:#52514e\\"><b>" + S.title +
      "</b> anchored at residue " + res +
      " (" + (d.aa[res - 1] || "?") + ") &middot; covers " + w0 + "&ndash;" + (w0 + N - 1) +
      (cs && cs.head !== (S.cnn_label || "CNN attention head") ? " &middot; " + cs.head + " " + cs.val : "") +
      " &middot; " + (S.cnn_label || "CNN attention head") + " " + S.g[i].toFixed(3) + "</div>";
    var btn = "<div class=\\"lf-pep-close\\" title=\\"close\\" style=\\"position:absolute;right:10px;top:2px;" +
      "cursor:pointer;font:16px ui-sans-serif,system-ui;color:#8a8a8a;padding:0 5px\\">&times;</div>";
    g += "<line class=\\"lf-pep-guide\\" x1=0 y1=6 x2=0 y2=" + (PANEL_H - 26) +
         " stroke=#000 stroke-width=1 opacity=0 />";
    host.innerHTML = g + hdr + btn;
    host.querySelector(".lf-pep-close").addEventListener("click", function () { closeSlot(k); });

    /* a window near either terminus would push the centred title off the panel,
       so nudge it back inside once its width is known */
    var hd = host.querySelector(".lf-pep-hdr");
    if (hd) {
      var hb = hd.getBoundingClientRect(), pb = host.getBoundingClientRect();
      var shift = 0;
      if (hb.left < pb.left + 6) shift = (pb.left + 6) - hb.left;
      else if (hb.right > pb.right - 34) shift = (pb.right - 34) - hb.right;
      if (shift !== 0) hd.style.left = (cx + shift).toFixed(0) + "px";
    }

    /* Hovering a point in the opened window marks that residue on every track
       on the page -- the same data_id linking the tracks use between
       themselves, done by hand because this panel is not a girafe widget. */
    var guide = host.querySelector(".lf-pep-guide"), held = [];
    function clearHl() {
      held.forEach(function (o) {
        if (o[1] === null) o[0].removeAttribute("stroke"); else o[0].setAttribute("stroke", o[1]);
        if (o[2] === null) o[0].removeAttribute("stroke-width"); else o[0].setAttribute("stroke-width", o[2]);
        if (o[3] !== null) o[0].setAttribute("r", o[3]);
      });
      held = [];
      if (guide) guide.setAttribute("opacity", 0);
    }
    function markResidue(r) {
      clearHl();
      document.querySelectorAll("svg.ggiraph-svg [data-id=\\"" + r + "\\"]").forEach(function (n) {
        held.push([n, n.getAttribute("stroke"), n.getAttribute("stroke-width"),
                   n.tagName.toLowerCase() === "circle" ? n.getAttribute("r") : null]);
        n.setAttribute("stroke", "#000");
        n.setAttribute("stroke-width", "1.4");
        if (n.tagName.toLowerCase() === "circle")
          n.setAttribute("r", scaleR(n.getAttribute("r"), 1.5));
      });
      if (guide) {
        guide.setAttribute("x1", X(r).toFixed(1));
        guide.setAttribute("x2", X(r).toFixed(1));
        guide.setAttribute("opacity", 0.3);
      }
    }
    host.addEventListener("mouseover", function (ev) {
      var r = ev.target && ev.target.getAttribute && ev.target.getAttribute("data-res");
      if (r) markResidue(parseInt(r, 10));
    });
    host.addEventListener("mouseleave", clearHl);
  }

  /* which series a dot belongs to -- its tooltip names it. Longest matching
     label wins, so a dot of "C terminus, CT-padded" is that run and not the
     plain C series its label also contains */
  /* a series with no surface on the page (an arm drawn from its tracks rds,
     tooltip marked by LF_NO_POPUP) has nothing to open -- and its label can
     CONTAIN a curated one ("ESM-C only C fold 1" holds "C fold 1"), so without
     this it would open the curated window */
  function noPopup(n) { return (n.getAttribute("title") || "").indexOf("no window popup") >= 0; }
  function seriesOf(n) {
    if (noPopup(n)) return null;
    var t = n.getAttribute("title") || "", key = null, best = -1;
    Object.keys(d.series || {}).forEach(function (k2) {
      var lab = d.series[k2].label;
      if (lab && t.indexOf(lab) >= 0 && lab.length > best) { key = k2; best = lab.length; }
    });
    return key;
  }

  function wire() {
    xgbSvgs().forEach(wireOne);
  }

  /* ---- a peptide in the URL fragment ---------------------------------------
     #CREG2x90x116 (or #90x116) opens the N window at 90 and the C window at
     116 on every xgb track carrying those series, and centres the view on the
     residue nearest their average; #103 just centres on residue 103. */
  function openAt(svg, term, r) {
    ptsWithId(svg, r).forEach(function (n) {
      var key = seriesOf(n);
      if (!key || (d.series[key].term || key) !== term) return;
      var cur = slots[slotKey(svg, term)];
      if (cur && cur.key === key && cur.res === r) return;
      open(key, r, svg, n);
    });
  }
  function centreOn(r, svg) {
    var sc = residueScale(svg);
    if (!sc) return;
    /* the visible area starts right of the fixed annotation box, when it is open */
    var left = 0, box = document.getElementById("annotation-box");
    if (box) {
      var b = box.getBoundingClientRect();
      if (b.right > 0 && b.right < window.innerWidth / 2) left = b.right;
    }
    var top = 0, bar = document.getElementById("top-meta-fixed");
    if (bar) top = bar.getBoundingClientRect().bottom;
    var dx = sc.x(r) - (left + window.innerWidth) / 2;
    var dy = svg.getBoundingClientRect().top - top - 12;
    window.scrollTo({ left: window.scrollX + dx, top: window.scrollY + dy });
  }
  function fromHash(tries) {
    var h = decodeURIComponent((location.hash || "").replace(/^#/, "")).trim();
    if (!h) return;
    var svgs = xgbSvgs();
    /* the widgets draw asynchronously; wait for a track with its dots */
    if (!svgs.length || !residueScale(svgs[0])) {
      if ((tries || 0) < 20) setTimeout(function () { fromHash((tries || 0) + 1); }, 500);
      return;
    }
    var m = /^(?:(.+?)x)?([0-9]+)x([0-9]+)$/.exec(h);
    if (m) {
      if (m[1] && d.gene && m[1].toUpperCase() !== String(d.gene).toUpperCase()) {
        console.warn("peptide " + h + " names " + m[1] + ", this page is " + d.gene);
        return;
      }
      var nr = parseInt(m[2], 10), cr = parseInt(m[3], 10);
      svgs.forEach(function (svg) { openAt(svg, "N", nr); openAt(svg, "C", cr); });
      centreOn(Math.round((nr + cr) / 2), svgs[0]);
      return;
    }
    var one = /^([0-9]+)$/.exec(h);
    if (one) centreOn(parseInt(one[1], 10), svgs[0]);
  }
  /* a small box to type one in: GENExNxC, NxC or a residue */
  function inputBox() {
    if (document.getElementById("lf-pep-goto")) return;
    var w = document.createElement("div");
    w.id = "lf-pep-goto";
    w.style.cssText = "position:fixed;top:6px;right:10px;z-index:10001;background:#fff;" +
      "border:1px solid #d9d8d3;border-radius:5px;padding:3px 6px;font:11px ui-sans-serif,system-ui;";
    w.innerHTML = "<input type=\\"text\\" placeholder=\\"" + (d.gene || "GENE") +
      "x90x116 or 103\\" style=\\"width:150px;font:11px ui-monospace,monospace;border:0;outline:0\\">";
    document.body.appendChild(w);
    var inp = w.querySelector("input");
    inp.addEventListener("keydown", function (ev) {
      if (ev.key !== "Enter") return;
      var v = inp.value.trim();
      if (!v) return;
      if (location.hash === "#" + v) fromHash(0); else location.hash = v;
    });
  }
  /* ggiraph draws shape 21 as <circle>, but shape 22 (the all-50 squares) as
     <rect> and shape 25 (the held-out fold triangles) as <polygon>. Binding to
     "circle[data-id]" alone left both of those inert -- clickable-looking and
     doing nothing. Ribbons carry no data-id, so they are never matched, and
     this only ever runs on the xgb head tracks. */
  var PT_SEL = "circle[data-id], polygon[data-id], rect[data-id]";
  function ptsWithId(svg, r) {
    return [].slice.call(svg.querySelectorAll(
      \'circle[data-id="\' + r + \'"], polygon[data-id="\' + r + \'"], rect[data-id="\' + r + \'"]\'));
  }

  function wireOne(svg) {
    [].slice.call(svg.querySelectorAll(PT_SEL)).forEach(function (n) {
      if (n.__lfpep) return;
      var v = n.getAttribute("data-id");
      if (!/^[0-9]+$/.test(v) || noPopup(n)) return;
      n.__lfpep = true;
      n.style.cursor = "pointer";
      n.addEventListener("click", function (ev) {
        ev.stopPropagation();
        var key = seriesOf(n);
        if (!key) return;
        /* the same point again closes its window; another point of that
           terminus replaces it; the other terminus opens beside it */
        var r = parseInt(v, 10), sk = slotKey(svg, d.series[key].term || key), cur = slots[sk];
        if (cur && cur.key === key && cur.res === r) closeSlot(sk);
        else open(key, r, svg, n);
      });
    });
  }
  function start() { setTimeout(wire, 400); setTimeout(function () { wire(); inputBox(); fromHash(0); }, 1500); }
  if (document.readyState === "complete") start();
  else window.addEventListener("load", start);
  setTimeout(wire, 4000);
  window.addEventListener("hashchange", function () { wire(); fromHash(0); });
  window.addEventListener("resize", closeAll);
})();
'
