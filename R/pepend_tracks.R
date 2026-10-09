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
                 term = term, head = "CNN attention head")
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
      cfg <- mod$Config$from_json(file.path(iso, "in", "config.json"))
      models <- mod$scan$load_member_models(file.path(iso, "out"), cfg, term)
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

#' Pool each window's opinion onto the residues it covers
#'
#' Members are averaged first -- the ensemble mean is the prediction -- then mean
#' and sd across the COVERING WINDOWS, so the band says how much the call depends
#' on which window you look through, not how much the seeds disagree.
#' @keywords internal
lf_track_pool_residues <- function(S, anchors, n_prot, c_prot, term,
                                   classes = c("pep_pocket", "pep_other",
                                               "CT_cleavage_context", "NT_cleavage_context")) {
  idx <- which(S$pi_names %in% classes)
  if (!length(idx)) return(NULL)
  res <- lapply(idx, function(k) {
    acc <- lapply(seq_along(anchors), function(i) {
      r <- S$w_start[i] + seq_len(LF_PEPEND$seq_len) - 1L
      ok <- r >= n_prot & r <= c_prot
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
  if (is.null(win) && is.null(rp) && is.null(resid)) return(NULL)
  list(win = win, rp = rp, resid = resid, pix = pix,
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
                                seq_offset = -0.45) {
  if (is.null(d) || !nrow(d)) return(NULL)
  n_terms <- length(unique(d$term))
  d <- dplyr::summarise(dplyr::group_by(d, residue, class),
                        mean = mean(mean), sd = mean(sd), n_win = sum(n_win),
                        .groups = "drop")
  d <- droplevels(d)
  lab <- if (n_terms > 1) "N+C per-residue\nclass probability" else "per-residue\nclass probability"
  ggplot2::ggplot(d, ggplot2::aes(residue, mean, colour = class)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(mean - sd, 0),
                                      ymax = pmin(mean + sd, 1), fill = class),
                         alpha = 0.16, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(data_id = residue,
                   tooltip = sprintf("%s  residue %d (%s)\n%s %.3f +/- %.3f  (mean of N and C)",
                                     gene, residue, aa[residue],
                                     nn_class_labels[as.character(class)], mean, sd)),
      size = 0.9, na.rm = TRUE) +
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
                              seq_offset = -0.45) {
  d <- d[d$head == head_nm, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$term <- factor(d$term, levels = names(LF_TERM_COLS))
  ggplot2::ggplot(d, ggplot2::aes(residue, mean, colour = term, fill = term)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                         alpha = 0.14, colour = NA, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = 0.45, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(data_id = residue,
                   tooltip = sprintf("%s  anchor %d (%s)\n%s terminus, %s: %.3f +/- %.3f",
                                     gene, residue, aa[residue], term, head_nm, mean, sd)),
      size = 0.9, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = LF_TERM_COLS, name = NULL,
                                 labels = function(x) paste0(x, " terminus")) +
    ggplot2::scale_fill_manual(values = LF_TERM_COLS, guide = "none") +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(y = sub(" head", "\nhead", head_nm)) +
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
      size = 0.9, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = nn_class_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = nn_class_cols, guide = "none") +
    lf_track_x(x_range, seq_offset) +
    ggplot2::scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(y = sprintf("%s per residue\n(no window)", mdl)) +
    lf_track_theme()
}

#' Every track that could be built, top to bottom, with its height in inches
#' @export
lf_pepend_track_plots <- function(dat, x_range, aa, gene, seq_offset = -0.45) {
  mk <- list()
  add <- function(p, h) if (!is.null(p)) mk[[length(mk) + 1L]] <<- list(p = p, h = h)
  add(lf_plot_class_track(dat$resid, x_range, aa, gene, TRUE, seq_offset), 1.7)
  if (!is.null(dat$win)) {
    add(lf_plot_win_track(dat$win, "CNN attention head", x_range, aa, gene, TRUE, seq_offset), 1.15)
    add(lf_plot_win_track(dat$win, "xgb window head",    x_range, aa, gene, FALSE, seq_offset), 1.15)
  }
  if (!is.null(dat$rp)) {
    add(lf_plot_resid_track(dat$rp, "MLP", x_range, aa, gene, seq_offset), 1.1)
    add(lf_plot_resid_track(dat$rp, "XGB", x_range, aa, gene, seq_offset), 1.1)
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

  pack <- function(S) {
    keep <- seq_along(pix$anchors)
    if (length(keep) > max_anchors)
      keep <- sort(utils::head(order(S$g, decreasing = TRUE), max_anchors))
    list(anchor = as.integer(pix$anchors[keep]),
         pos = as.integer(S$pos),
         g = round(S$g[keep], 4),
         pi = lapply(keep, function(i) round(S$pi[i, , ], 3)),
         sd = lapply(keep, function(i) round(S$pi_sd[i, , ], 3)))
  }
  payload <- list(gene = gene, accession = accession, aa = aa,
                  seq_len = as.integer(LF_PEPEND$seq_len),
                  classes = pix[[terms[1]]]$pi_names,
                  labels = unname(nn_class_labels[pix[[terms[1]]]$pi_names]),
                  cols = unname(nn_class_cols[pix[[terms[1]]]$pi_names]))
  for (t in terms) payload[[t]] <- pack(pix[[t]])

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
  var PANEL_H = 150, host = null, track = null;

  /* The clickable track is the xgb window head, and only that one: it is the
     head that actually localises a boundary, so it is the one worth opening a
     window from. Found by its y-axis label rather than by index, since a track
     is dropped when its input is missing. */
  function xgbSvg() {
    var all = [].slice.call(document.querySelectorAll("svg.ggiraph-svg"));
    for (var i = 0; i < all.length; i++) {
      var t = all[i].querySelectorAll("text");
      for (var j = 0; j < t.length; j++)
        if (/xgb window/.test(t[j].textContent)) return all[i];
    }
    return null;
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

  function close() {
    if (!host) return;
    /* a highlight left on another track would survive the panel */
    document.querySelectorAll("svg.ggiraph-svg [data-id]").forEach(function (n) {
      if (n.__lfhl) { n.removeAttribute("stroke"); delete n.__lfhl; }
    });
    host.remove(); host = null;
  }

  function open(term, res, svg) {
    var S = d[term];
    if (!S) return;
    var i = S.anchor.indexOf(res);
    if (i < 0) return;
    var sc = residueScale(svg);
    if (!sc) return;
    close();

    /* The panel must live in the SAME coordinate space as the track, not the
       page. These svgs are drawn at natural size -- NTN1 is about 86 inches --
       so the girafe container scrolls horizontally. A width:100% div outside it
       takes the PAGE width: every x computed from the track then lands past its
       right edge and is clipped, which shows up as a window stuck at the left
       with the points missing. So it goes INSIDE the container, directly after
       the svg, at exactly the width of that svg.
       NOTE: no apostrophes anywhere in this block -- the whole handler is an R
       single-quoted string, and one would end it. */
    var svgRect = svg.getBoundingClientRect();
    var cont = svg.closest(".girafe_container_std") || svg.parentElement;
    host = document.createElement("div");
    host.id = "lf-pep-window";
    host.style.cssText = "position:relative;margin:0;width:" +
      svgRect.width.toFixed(0) + "px;";
    svg.parentNode.insertBefore(host, svg.nextSibling);

    var hostRect = host.getBoundingClientRect();
    var N = d.seq_len, w0 = res - S.pos + 1, K = d.classes.length;
    var X = function (r) { return sc.x(r) - hostRect.left; };
    var Y = function (v) { return PANEL_H - 26 - v * (PANEL_H - 44); };

    var g = "<svg width=\\"" + svgRect.width.toFixed(0) + "\\" height=\\"" + PANEL_H +
            "\\" style=\\"display:block\\">";
    /* the window frame, on the residue axis */
    g += "<rect x=" + X(w0 - 0.5).toFixed(1) + " y=6 width=" +
         Math.max(1, (X(w0 + N - 0.5) - X(w0 - 0.5))).toFixed(1) +
         " height=" + (PANEL_H - 32) + " fill=#ffffff stroke=#d9d8d3 />";
    /* the anchor: the residue this window is ABOUT */
    g += "<rect x=" + (X(res) - 1).toFixed(1) + " y=6 width=2 height=" +
         (PANEL_H - 32) + " fill=#000 opacity=0.35 />";
    for (var k = 0; k < K; k++) {
      if (d.classes[k] === "none" || d.classes[k] === "padding") continue;
      var pts = [];
      for (var j = 0; j < N; j++)
        pts.push(X(w0 + j).toFixed(1) + "," + Y(S.pi[i][j][k]).toFixed(1));
      /* +/- 1 sd across the members, as a band under the line */
      if (S.sd) {
        var up = [], dn = [];
        for (var j3 = 0; j3 < N; j3++) {
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
      for (var j2 = 0; j2 < N; j2++) {
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
    var hdr = "<div id=\\"lf-pep-hdr\\" style=\\"position:absolute;top:4px;left:" +
      cx.toFixed(0) + "px;transform:translateX(-50%);white-space:nowrap;text-align:center;" +
      "font:11px ui-sans-serif,system-ui;color:#52514e\\"><b>" + term +
      " window</b> anchored at residue " + res +
      " (" + (d.aa[res - 1] || "?") + ") &middot; covers " + w0 + "&ndash;" + (w0 + N - 1) +
      " &middot; score " + S.g[i].toFixed(3) + "</div>";
    var btn = "<div id=\\"lf-pep-close\\" title=\\"close\\" style=\\"position:absolute;right:10px;top:2px;" +
      "cursor:pointer;font:16px ui-sans-serif,system-ui;color:#8a8a8a;padding:0 5px\\">&times;</div>";
    g += "<line id=\\"lf-pep-guide\\" x1=0 y1=6 x2=0 y2=" + (PANEL_H - 26) +
         " stroke=#000 stroke-width=1 opacity=0 />";
    host.innerHTML = g + hdr + btn;
    host.querySelector("#lf-pep-close").addEventListener("click", close);

    /* a window near either terminus would push the centred title off the panel,
       so nudge it back inside once its width is known */
    var hd = host.querySelector("#lf-pep-hdr");
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
    var guide = host.querySelector("#lf-pep-guide"), held = [];
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
          n.setAttribute("r", (parseFloat(n.getAttribute("r") || 2) + 1.2).toFixed(1));
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

  function wire() {
    var svg = xgbSvg();
    if (!svg) return;
    [].slice.call(svg.querySelectorAll("circle[data-id]")).forEach(function (n) {
      if (n.__lfpep) return;
      var v = n.getAttribute("data-id");
      if (!/^[0-9]+$/.test(v)) return;
      n.__lfpep = true;
      n.style.cursor = "pointer";
      n.addEventListener("click", function (ev) {
        ev.stopPropagation();
        /* which series was clicked -- the tooltip says which terminus */
        var t = n.getAttribute("title") || "";
        var term = /N terminus/.test(t) ? "N" : (/C terminus/.test(t) ? "C" : null);
        if (term) open(term, parseInt(v, 10), svg);
      });
    });
  }
  if (document.readyState === "complete") setTimeout(wire, 400);
  else window.addEventListener("load", function () { setTimeout(wire, 400); });
  setTimeout(wire, 1500); setTimeout(wire, 4000);
  window.addEventListener("resize", close);
})();
'
