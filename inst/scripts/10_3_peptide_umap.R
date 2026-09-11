## =============================================================================
## Peptide-window UMAP + parameters + PCA + model scores  -- ONE script.
##
## Consolidates what used to be three passes over the same data:
##   10_3_peptide_umap.R   flatten windows -> UMAP -> peptide_umap.csv
##   10_3b_add_params.R    rebuild params + PCA, re-join the coords by CSV
##   10_3c_add_scores.R    re-read that CSV just to attach model scores
## The wide matrix and the per-window summaries were being built twice, by two
## different implementations, and the coords/scores were round-tripped through
## CSV to get back to objects that were in the session the whole time. Now it is
## built once, in memory, and written once.
##
## Flattens every `seq_len`-position window into ONE row -- one column per
## (position x feature) -- and emits, per window:
##   labels        peps / gene / terminus / target / win_type / end_type /
##                 known / set / pad_frac / uniprot_known
##   coords        UMAP1_<metric>, UMAP2_<metric> for each metric in
##                 `umap_metrics`, plus PC1..PC<pca_dims>
##   summaries     mean_<feature> and top<NN>_<feature> per window
##   positions     the raw per-position values used as the UMAP input
##   scores        pred_nn -- the current session's model scores
##
## Outputs per UMAP metric: a rasterised static svg, and an all-vector html in
## which every window is hoverable and links to its per-gene page under
## `link_base`.
##
## "one parameter/position per column, one row per peptide" == a pivot_wider of
## the per-residue feature matrix. We flatten the arrays directly (vectorised),
## which gives the same wide table as pivot_longer() |> pivot_wider() but is far
## faster for many thousands of windows.
##
## Run AFTER 9.2 (needs `nn_input`, `all_params3`, `seq_len` in the session).
## Scores come only from the current session: `add_scores` additionally needs
## `nn_input_comb` from 10_1dcnn_new6.R, and errors if that object isn't there.
## =============================================================================

library(tidyverse)

## ---- options ---------------------------------------------------------------
out_csv       <- path.expand("~/AF2_analysis/peptide_umap_full.csv")     # the one output table
## Cold-start inputs. Building `params` needs a training session (nn_input +
## all_params3); re-PLOTTING it does not. These are the saved equivalents, used
## when the session is empty -- see the fast path under "setup / guards".
comb_rds      <- path.expand("~/AF2_analysis/nn_input_comb_ensemble.rds")
load_csv      <- path.expand("~/AF2_analysis/peptide_pca_loadings.csv")  # PCA loadings (feature x PC)
svg_stem      <- path.expand("~/AF2_analysis/peptide_umap")              # <stem>_<metric>.svg
html_name     <- "ligandFinder_v5.html"  # fixed filename for the interactive page, written
                                         # beside svg_stem. Only used when a single layout is
                                         # plotted -- with more than one it would collide, so
                                         # the <stem>_<layout>.html form takes over. NULL = always
                                         # use the stem form. The svg is unaffected.
## ---- linked side panels ----------------------------------------------------
## Extra panels drawn in the SAME girafe as the UMAP and keyed on the same window
## id, so hovering a window in one panel highlights it in all of them. Two
## separate girafe objects would not link, however they are laid out.
panel_score   <- TRUE                    # panel 1: score by stratum (violin + beeswarm)
umap_rel_w    <- 1.15                    # UMAP width relative to one side panel
panel_rel_w   <- 1                       # ... shrinking the UMAP is just this ratio
panel_all     <- TRUE                    # panels cover EVERY window, not only the
                                         # gated ones the UMAP draws. The score
                                         # distribution is the point of panel 1, and
                                         # the gate removes 56% of the windows.
set_cols      <- c(train = "#1b7837", val = "#2166ac", none = "grey75")
## Hover styling, shared by every interactive layer in the girafe. The hovered
## window (in BOTH panels, since they share data_id) gets the first rule; every
## other interactive mark gets the second, which is what makes the highlight
## read against 85k neighbours -- a thicker stroke alone disappears in the bulk.
hover_css     <- "stroke:#FF6600;stroke-width:3.5;fill:#FF6600;fill-opacity:1;opacity:1;cursor:pointer;"
hover_inv_css <- "opacity:0.3;"
## Only windows scoring ABOVE this are interactive -- hover text, click-through
## to the per-gene page, and the gene-search highlight. Everything at or below
## it is still drawn, in the UMAP and in the score panel alike, but as plain
## static marks: 58k of the 59k windows sit at ~0, and giving each of them a
## tooltip, an id and an onclick made the page heavy and put the whole cloud
## under the hover dimming. NULL = every window interactive.
interactive_min <- 0.1
panel_unk_split <- 0.1                   # the unknown column splits at this score:
                                         # nearly all 59k unknowns sit at ~0 and
                                         # would otherwise bury the few that score
panel_pt_size <- 2.2                     # panel marks are the whole panel now (no
                                         # violin behind them), so they carry the
                                         # distribution and need to be readable
panel_swarm_w <- 0.42                    # swarm half-width; wider than it could be
                                         # with a violin outline to respect

point_size    <- 1.2                     # data-mark size, shared by the static and the
                                         # interactive layer so the two cannot drift apart
term_shapes   <- c(C = 16, N = 17)       # C-terminal windows draw as points, N-terminal as
                                         # triangles: terminus is categorical and the colour
                                         # aesthetic is already spent on the score

## The per-residue features -- only the UMAP/PCA build reads this, so on the
## replot path (no training session) it is simply absent. The build branch keeps
## its own stopifnot(exists("all_params3")), so a genuinely missing feature set
## still stops the run rather than producing an empty layout.
feature_set   <- if (exists("all_params3")) all_params3 else NULL
standardize   <- TRUE               # center/scale each column so continuous feats don't dominate
seed          <- 42

## One UMAP per metric; each contributes its own UMAP1_<metric>/UMAP2_<metric>
## pair. euclidean sees absolute magnitude differences between windows, cosine
## only the shape/direction of the profile -- windows that differ mainly in
## overall conservation or burial level separate under euclidean and can collapse
## together under cosine, so the two views are worth having side by side.
umap_metrics  <- "euclidean"        # final figure is the euclidean layout only
umap_args     <- list(n_neighbors = 15, min_dist = 0.1)   # metric is supplied per run

## A THIRD layout, for comparison: UMAP of the model's own 16-d representation
## (the dense layer named "embed" in 10_1dcnn_new6.R) instead of the input
## features. The two answer different questions -- the input-space layouts show
## which windows are biophysically alike, so structure there is NOT evidence the
## network learned anything; the embedding layout shows how the trained global
## ranker organises the same windows. Needs `emb_all` (and `nn_input_comb`, to
## align its rows) in the session; skipped with a message otherwise.
## The embedding is L2-normalised at source, so cosine is the metric that
## matches it -- euclidean on unit vectors is a monotone transform of the same
## thing and would just reproduce the layout.
embed_umap    <- FALSE              # off for the final figure: the gate below is in
                                    # euclidean-layout coordinates and would cut the
                                    # embedding layout at meaningless places
embed_metric  <- "cosine"

make_plot     <- TRUE               # write a coloured UMAP svg per metric
## Which score the UMAP is coloured by and the score panel plots. "pred_nn" is
## nn_input_comb$pred -- the calibrated ENSEMBLE mean over the 5 seeds, i.e. the
## same number the per-gene pages show; "pred_raw" is the uncalibrated sigmoid.
## They rank similarly but are not on the same footing, so the viewer and the
## gene pages should quote one of them, not one each.
score_col     <- "pred_nn"
## Axis/legend text, so the reader is not left to interpret a column name.
score_label   <- c(pred_nn  = "ensemble score (mean of 5 seeds)",
                   pred_raw = "raw score (uncalibrated)")
density_bins  <- 8                  # number of 2-D density contour levels (NULL = no contours)
contour_col   <- "black"             # contour line colour
contour_lw    <- 0.35               # contour line width
contour_glow  <- TRUE               # halo the contours via ggfx so they read over dense points
glow_col      <- "white"            # halo colour -- should contrast with contour_col
## The contour-with-halo look comes from scTools::plot_subsets() (R/plots.R);
## that version uses sigma = 2, expand = 3 for a softer, thicker halo, which
## separates the line from a dense point field more strongly than these do.
glow_sigma    <- 1                  # halo softness
glow_expand   <- 2                  # halo thickness (px)
contour_dpi   <- 300                # resolution the contour layer is rasterised at (via ragg)
raster_points <- TRUE               # rasterise the point layer of the static svg (needs ggrastr);
                                    # the html is all-vector by necessity -- see note below
known_col     <- "red"              # known uniprot peptides: ring colour
known_rings   <- FALSE              # draw those rings at all. Off: the known ends
                                    # are now a column of their own in the score
                                    # panel, and the rings were punching holes in
                                    # the densest part of the cloud
                                    # Amidation is NOT a mark -- it reads out in the hover
                                    # text instead, so it can coexist with the known ring on
                                    # the same window without one obscuring the other.
## Split the figure into motif / no-motif panels. Scales stay fixed, so the two
## panels are directly comparable -- the question is WHERE in the layout the
## motif windows sit, which a free scale would destroy. Note the contours become
## a density per panel rather than of everything, which is the point: it shows
## whether motif windows concentrate somewhere or track the bulk.
facet_amidation <- FALSE            # one panel; the motif still reads out on hover
## No printed labels: which peptide a known ring belongs to reads out in the
## hover text instead. Labels could only ever name the top few dozen without
## colliding, and they punched opaque holes in the densest part of the cloud --
## exactly where the knowns are.
## Gate: keep only windows inside this box in the plotted layout's coordinates.
## It restricts the figure entirely -- points, density contours and labels are
## all computed from the gated subset, not merely zoomed to it. NULL disables.
gate_umap1    <- c(-11, -1)         # UMAP1 range to keep
gate_umap2    <- c(-11, 4)          # UMAP2 range to keep. Extended from -7 to take in the
                                    # N-terminal satellite at UMAP2 ~ -9 -- a small, sparse,
                                    # spatially separate cluster that is heavily enriched for
                                    # known N-terminal peptide ends (apelin, kisspeptin,
                                    # somatostatin, dynorphin, neurotensin, insulin). At -7 the
                                    # gate cut straight through the gap between it and the main
                                    # blob and excluded 19 known ends.
add_positions <- TRUE               # append the raw per-position feature values
top_pct       <- 0.80               # per-window "top" summary = mean of residues >= this pctile
excl_padding  <- TRUE               # summaries over real residues only (padding excluded)
add_pca       <- TRUE               # add PC1..PC<pca_dims> scores + write loadings
pca_dims      <- 20                 # number of principal components
add_scores    <- TRUE               # join pred_nn from the current session
link_base     <- "https://d2v3leolhhovg9.cloudfront.net"
                                    # base URL of the per-gene plot pages; NULL = no HTML output.
                                    # Points link to <base>/<gene>.html#<peps>, which the page
                                    # resolves via its data-peps deep link.
                                    # EVERY window is hoverable and clickable, so the html has one
                                    # vector mark per window and nothing rasterised: it is a big
                                    # file. The static svg is the small one.
## Annotate every window whose dibasic anchor sits at a KNOWN peptide boundary,
## independent of whether that peptide is in the training set. The training
## positives are deliberately narrow -- they need a cognate receptor-ligand model
## and a terminal (not middle) insertion -- so most real peptide ends are absent
## from `known`. This reads the boundaries straight off ligand_list.rds and needs
## nothing from the pipeline session.
known_end_ref <- "/Users/kbrulois/R_projects/ligandFinder/inst/extdata/ligand_list.rds"
known_end_len <- 100                # max reference peptide length (see the audit:
                                    # ligand_type is NA for most rows and unusable)
known_end_gap <- 5                  # anchor-to-boundary distance to accept; matches
                                    # 9.2's own db_spacer < 6 rule

## Gene search box in the html: type a symbol, its best windows light up in the
## layout and are listed underneath. Needs jsonlite.
gene_search   <- TRUE
gene_search_n <- 5                  # how many windows per gene to rank and light up
gene_hit_col  <- "#00E5FF"          # highlight ring -- cyan reads against magma
                                    # (purple/yellow) and against the red known rings
subset_peps   <- NULL               # NULL = all windows; or a character vector of peps to keep
knowns_only   <- FALSE              # TRUE = only known peptides (fast sanity run)

## ---- setup / guards --------------------------------------------------------
## Everything from here to the PCA block rebuilds `params` from the session --
## the expensive half of this script, and the only half that needs nn_input +
## all_params3. The annotation and plotting sections below only ever READ
## `params`, so when the training session is gone but the saved table is there,
## reload it and skip straight to them. Set reuse_csv explicitly to force a path.
.from_csv <- if (exists("reuse_csv")) isTRUE(reuse_csv) else
             (!exists("nn_input") && file.exists(out_csv))

if (.from_csv) {
  message("reuse: loading params from ", out_csv, " (skipping the UMAP/PCA rebuild)")
  ## Only the label/coord columns -- the ~1,000 per-position columns are input to
  ## the UMAP, which is not being recomputed, and reading them costs minutes.
  params <- readr::read_csv(
    out_csv,
    col_select = c(peps, gene, terminus, target, win_type, end_type, known, set,
                   pad_frac, dplyr::starts_with("UMAP"), dplyr::num_range("PC", 1:5)),
    show_col_types = FALSE, progress = FALSE)
  ## Drop any stale scores/annotations so the sections below re-derive them from
  ## the CURRENT nn_input_comb, instead of silently keeping older values under
  ## the same column names.
  params <- dplyr::select(params, -dplyr::any_of(
    c("pred_nn", "pred_raw", "pred_april17", "amidation",
      "known_end", "known_end_name", "uniprot_known")))
  message(sprintf("  %d windows x %d columns", nrow(params), ncol(params)))
} else {
stopifnot(exists("nn_input"), exists("all_params3"))
if (!exists("seq_len")) seq_len <- nrow(nn_input[[1]]$all$data[[1]])
if (!requireNamespace("uwot", quietly = TRUE))
  stop("install.packages('uwot') for UMAP")

col_or_na  <- function(d, nm) if (nm %in% names(d)) d[[nm]] else NA
have_mstat <- requireNamespace("matrixStats", quietly = TRUE)   # fast row quantiles
top_tag    <- paste0("top", round(top_pct * 100))
if (!have_mstat)
  warning("matrixStats not installed -- the top", round(top_pct * 100),
          " percentile falls back to a slow per-row apply(). ",
          "install.packages('matrixStats') to make this fast.", immediate. = TRUE)

write_csv_fast <- function(d, p) if (requireNamespace("data.table", quietly = TRUE))
                                   data.table::fwrite(d, p) else readr::write_csv(d, p)

message("windows per terminus: ",
        paste(sprintf("%s=%d", names(nn_input),
                      sapply(nn_input, function(x) length(x$all$data))), collapse = "  "))

## ---- feature columns + wide names ------------------------------------------
## wide column names: pos01_cons_rs, pos01_cons_rs_n, ..., pos36_padding
## (position is the outer / slower-varying index, feature the inner one -- this
##  matches the aperm below, which reads feature-fastest within a position)
first    <- nn_input[[1]]$all$data[[1]]
ch_names <- colnames(first); if (is.null(ch_names)) ch_names <- names(first)
sel      <- match(feature_set, ch_names)
stopifnot(!anyNA(sel), length(sel) == length(feature_set))
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
## For each window and feature: mean over residues, and mean of the top (>=
## top_pct percentile) residues within that window. `excl_padding` restricts both
## to real residues (padding channel == 0); pad_frac is reported separately.
build_term <- function(term) {
  d <- nn_input[[term]]$all
  if (isTRUE(knowns_only))   d <- d[which(d$known == 1), ]
  if (!is.null(subset_peps)) d <- d[which(d$peps %in% subset_peps), ]
  n <- length(d$data)
  message(sprintf("  [%s] building %d windows ...", term, n)); flush.console()
  if (n == 0) return(NULL)

  ## (n, seq_len, P) array of the selected features
  A <- aperm(array(unlist(lapply(d$data, function(m) as.matrix(m[, sel, drop = FALSE]))),
                   dim = c(seq_len, P, n)), c(3, 1, 2))

  ## matrix(A[, , j], nrow = n) rather than A[, , j]: the latter drops to a plain
  ## vector when n == 1 (reachable via knowns_only / subset_peps), which breaks
  ## every rowSums below.
  pad_mask <- if (is.na(pad_j) || !excl_padding) matrix(FALSE, n, seq_len)
              else matrix(A[, , pad_j], nrow = n) == 1
  n_resid  <- seq_len - rowSums(pad_mask); n_resid[n_resid == 0] <- NA_real_

  means <- matrix(NA_real_, n, length(resid_j), dimnames = list(NULL, paste0("mean_", resid_nm)))
  tops  <- matrix(NA_real_, n, length(resid_j), dimnames = list(NULL, paste0(top_tag, "_", resid_nm)))
  for (k in seq_along(resid_j)) {
    M  <- matrix(A[, , resid_j[k]], nrow = n)     # n x seq_len (see pad_mask note)
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
      d$peps %in% nn_input[[term]]$train$peps ~ "train",   # known/control used for training
      d$peps %in% nn_input[[term]]$val$peps   ~ "val",     # known/control held out
      TRUE                                    ~ "none"),   # candidate (not in the split)
    pad_frac = rowSums(pad_mask) / seq_len)

  res <- bind_cols(meta, as_tibble(means), as_tibble(tops))
  W   <- matrix(aperm(A, c(1, 3, 2)), nrow = n); colnames(W) <- wide_names  # pos-major wide
  if (add_positions) res <- bind_cols(res, as_tibble(W))
  list(res = res, W = W)   # W always returned so UMAP/PCA can use it even if add_positions=FALSE
}

built  <- lapply(names(nn_input), build_term)
params <- dplyr::bind_rows(lapply(built, `[[`, "res"))
Xfull  <- do.call(rbind, lapply(built, `[[`, "W"))   # n x (P*seq_len) raw per-position matrix
rm(built); gc()
stopifnot(nrow(Xfull) == nrow(params))
message(sprintf("flattened %d windows x %d features (%d pos x %d params)",
                nrow(Xfull), ncol(Xfull), seq_len, P))

## ---- drop zero-variance columns, standardise (once, shared by UMAP + PCA) ---
mu   <- colMeans(Xfull)
sdv  <- sqrt(pmax(colMeans(Xfull * Xfull) - mu^2, 0))
keep <- sdv > .Machine$double.eps
Xk   <- Xfull[, keep, drop = FALSE]
rm(Xfull); gc()
Xs   <- scale(Xk, center = mu[keep], scale = sdv[keep])   # PCA always uses the standardised matrix
X_umap <- if (standardize) Xs else Xk
message(sprintf("reduced-dim input: %d features (%d zero-variance dropped)",
                ncol(Xk), sum(!keep)))

## ---- UMAP, one run per metric ----------------------------------------------
## Plain UMAP1/UMAP2 are aliased to the FIRST metric so anything downstream that
## expects the old column names keeps working.
for (mt in umap_metrics) {
  message("  UMAP [", mt, "] ..."); flush.console()
  set.seed(seed)
  um <- do.call(uwot::umap, c(list(X = X_umap, metric = mt), umap_args))
  params[[paste0("UMAP1_", mt)]] <- um[, 1]
  params[[paste0("UMAP2_", mt)]] <- um[, 2]
}
params$UMAP1 <- params[[paste0("UMAP1_", umap_metrics[1])]]
params$UMAP2 <- params[[paste0("UMAP2_", umap_metrics[1])]]

## ---- PCA on the standardised per-position matrix ---------------------------
## Same input space as the UMAP. Scores go in the table; loadings (feature x PC)
## are per-feature so they can't live in the per-window table -- they go to
## `load_csv`.
if (add_pca) {
  k <- min(pca_dims, ncol(Xs), nrow(Xs) - 1L)
  if (requireNamespace("irlba", quietly = TRUE)) {         # fast truncated PCA
    pca <- irlba::prcomp_irlba(Xs, n = k, center = FALSE, scale. = FALSE)
    scores <- pca$x; rot <- pca$rotation; sdev <- pca$sdev
  } else {
    pca <- prcomp(Xs, center = FALSE, scale. = FALSE)      # full PCA, then take first k
    scores <- pca$x[, 1:k, drop = FALSE]; rot <- pca$rotation[, 1:k, drop = FALSE]; sdev <- pca$sdev[1:k]
  }
  colnames(scores) <- paste0("PC", 1:k)
  params <- bind_cols(params, as_tibble(scores))

  ## loadings table: one row per kept per-position feature, PC1..PCk columns
  ve <- sdev^2 / ncol(Xs)                     # var explained per PC (std'd data: total var = #features)
  colnames(rot) <- paste0("PC", 1:k)
  loadings <- tibble::tibble(feature = colnames(Xs)) %>%
    dplyr::mutate(                            # split pos12_cons_rs -> 12 + cons_rs
      position     = as.integer(stringr::str_match(feature, "^pos(\\d+)_")[, 2]),
      feature_type = stringr::str_remove(feature, "^pos\\d+_"),
      .after = feature) %>%
    bind_cols(as_tibble(rot))
  write_csv_fast(loadings, load_csv)
  message(sprintf("PCA: %d PCs, cum. var explained = %.1f%%; loadings -> %s",
                  k, 100 * sum(ve), load_csv))
}

}  # end of the session-only build (see .from_csv above)

## ---- embedding UMAP --------------------------------------------------------
## Same points, laid out by the model's 16-d "embed" representation rather than
## by the input features. Lands in UMAP1_embed / UMAP2_embed, so the plot loop
## below picks it up as just another layout.
if (embed_umap) {
  ekey <- intersect(c("peps", "target"), names(params))
  if (!exists("emb_all")) {
    message("emb_all not in session -- skipping the embedding UMAP ",
            "(it is built near the end of 10_1dcnn_new6.R)")
  } else if (!exists("nn_input_comb")) {
    message("nn_input_comb not in session -- cannot align emb_all rows; skipping embedding UMAP")
  } else if (nrow(emb_all) != nrow(nn_input_comb)) {
    message(sprintf("emb_all (%d rows) and nn_input_comb (%d) disagree -- skipping embedding UMAP",
                    nrow(emb_all), nrow(nn_input_comb)))
  } else {
    ## emb_all rows align with nn_input_comb, NOT with params: params is rebuilt
    ## by build_term() and can be filtered (knowns_only / subset_peps) or ordered
    ## differently. Match on the same peps+target key the score join uses instead
    ## of trusting the two tables to line up. Duplicate keys resolve to the first
    ## occurrence, as match() does -- the score join collapses them by max, so
    ## duplicates do exist.
    kk  <- function(d) do.call(paste, c(unname(as.list(d[ekey])), sep = "\r"))
    idx <- match(kk(params), kk(nn_input_comb))
    ok  <- !is.na(idx)
    message(sprintf("embedding UMAP: %d/%d windows matched to emb_all (%d-d), metric %s",
                    sum(ok), nrow(params), ncol(emb_all), embed_metric))

    if (any(ok)) {
      set.seed(seed)
      um <- do.call(uwot::umap, c(list(X = emb_all[idx[ok], , drop = FALSE],
                                       metric = embed_metric), umap_args))
      ## unmatched windows stay NA and are dropped when that layout is drawn,
      ## rather than being silently collapsed onto a shared coordinate
      params$UMAP1_embed <- NA_real_
      params$UMAP2_embed <- NA_real_
      params$UMAP1_embed[ok] <- um[, 1]
      params$UMAP2_embed[ok] <- um[, 2]
    }
  }
}

## ---- model scores ----------------------------------------------------------
## Joined in-session on peps+target (C / N / loop_C / loop_N). Duplicate keys are
## collapsed by taking the max score.
if (add_scores) {
  key <- intersect(c("peps", "target"), names(params))
  stopifnot("peps" %in% key)
  mx  <- function(x) { x <- x[!is.na(x)]; if (length(x)) max(x) else NA_real_ }

  score_by_key <- function(tbl, score_col, new) {
    tbl %>%
      dplyr::select(dplyr::all_of(c(key, score_col))) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
      dplyr::summarise("{new}" := mx(.data[[score_col]]), .groups = "drop")
  }

  ## nn_input_comb is the only score source, so a missing object stops the run
  ## rather than being skipped with a message: an unscored table would otherwise
  ## only announce itself much later, as a UMAP silently coloured by win_type.
  ## Cold start: the saved ensemble table IS the object 10_1dcnn_new6.R leaves
  ## in the session, so load it rather than demanding a training session.
  if (!exists("nn_input_comb") && file.exists(comb_rds)) {
    message("loading nn_input_comb from ", comb_rds)
    nn_input_comb <- readRDS(comb_rds)
  }
  if (!exists("nn_input_comb"))
    stop("nn_input_comb not in session and ", comb_rds, " does not exist -- ",
         "run 10_1dcnn_new6.R first, or set add_scores <- FALSE")
  sc <- intersect(c("pred", "pred_raw", "pred_cal"), names(nn_input_comb))[1]
  if (is.na(sc))
    stop("nn_input_comb has no pred/pred_raw/pred_cal column")

  params <- dplyr::left_join(params, score_by_key(nn_input_comb, sc, "pred_nn"), by = key)
  message(sprintf("pred_nn (from nn_input_comb$%s): matched %d / %d windows",
                  sc, sum(!is.na(params$pred_nn)), nrow(params)))

  ## The figure colours by the RAW (uncalibrated) score, so carry it as its own
  ## column rather than redefining pred_nn -- the calibrated score stays
  ## available, and the legend title then says which one is on screen.
  if ("pred_raw" %in% names(nn_input_comb) && !identical(sc, "pred_raw")) {
    params <- dplyr::left_join(params, score_by_key(nn_input_comb, "pred_raw", "pred_raw"), by = key)
    message(sprintf("pred_raw: matched %d / %d windows",
                    sum(!is.na(params$pred_raw)), nrow(params)))
  } else if (!"pred_raw" %in% names(nn_input_comb)) {
    message("nn_input_comb has no pred_raw column -- the figure will fall back to pred_nn")
  }

  ## Amidation-motif flag, computed in 10_1dcnn_new6.R. Collapsed with any(): a
  ## peps+target key can cover more than one row and the motif is a property of
  ## the window, so one flagged row flags the key.
  if ("amidation" %in% names(nn_input_comb)) {
    params <- dplyr::left_join(
      params,
      nn_input_comb %>%
        dplyr::select(dplyr::all_of(c(key, "amidation"))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
        dplyr::summarise(amidation = any(amidation, na.rm = TRUE), .groups = "drop"),
      by = key)
    params$amidation[is.na(params$amidation)] <- FALSE
    message(sprintf("amidation: %d of %d windows carry the motif",
                    sum(params$amidation), nrow(params)))
  } else {
    ## Recompute rather than skip. The motif is a fixed-position lookup, not a
    ## search: 9.2 anchors every db window on db_ind, so the dibasic pair always
    ## lands at the same local index -- 31 (2nd basic) for C-target windows, 6
    ## (1st basic) for N-target ones -- with the glycine immediately 5' of it.
    ## Same rule as the amidation block in 10_1dcnn_new6.R; kept here so a cold
    ## start off the saved table is not silently missing the flag.
    .ws <- if (exists("win_size")) c(N = win_size$N$start, C = win_size$C$start)
           else c(N = -5L, C = -30L)                          # 9.2 defaults
    .mpos <- function(tg) if (tg %in% c("C", "loop_C")) {
        a <- 1L - .ws[["C"]]; c(g = a - 2L, b1 = a - 1L, b2 = a)
      } else { a <- 1L - .ws[["N"]]; c(g = a - 1L, b1 = a, b2 = a + 1L) }
    .aa <- function(md, i) { v <- as.character(md[["AA"]])
                             if (i < 1L || i > length(v)) NA_character_ else v[[i]] }
    .amid <- vapply(seq_len(nrow(nn_input_comb)), function(i) {
      tg <- as.character(nn_input_comb$target[[i]])
      if (!identical(as.character(nn_input_comb$win_type[[i]]), "db")) return(FALSE)
      if (!tg %in% c("N", "loop_N", "C", "loop_C")) return(FALSE)
      q <- .mpos(tg); md <- nn_input_comb$meta_data[[i]]
      isTRUE(identical(.aa(md, q[["g"]]), "G") &&
             .aa(md, q[["b1"]]) %in% c("K", "R") && .aa(md, q[["b2"]]) %in% c("K", "R"))
    }, logical(1))
    ## Guard: a wrong offset makes every window come back FALSE, which reads as
    ## "no motifs found" rather than as a bug. Confirm the anchor really is basic.
    .chk <- vapply(utils::head(which(nn_input_comb$win_type == "db"), 2000), function(i) {
      q <- .mpos(as.character(nn_input_comb$target[[i]])); md <- nn_input_comb$meta_data[[i]]
      isTRUE(.aa(md, q[["b1"]]) %in% c("K", "R") && .aa(md, q[["b2"]]) %in% c("K", "R"))
    }, logical(1))
    message(sprintf("amidation: recomputed; dibasic anchor confirmed in %.1f%% of sampled db windows",
                    100 * mean(.chk)))
    if (mean(.chk) < 0.9)
      warning("amidation: the dibasic anchor is often NOT at the expected local ",
              "position -- check win_size against 9.2", immediate. = TRUE)
    nn_input_comb$amidation <- .amid
    params <- dplyr::left_join(
      params,
      nn_input_comb %>%
        dplyr::select(dplyr::all_of(c(key, "amidation"))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
        dplyr::summarise(amidation = any(amidation, na.rm = TRUE), .groups = "drop"),
      by = key)
    params$amidation[is.na(params$amidation)] <- FALSE
    message(sprintf("amidation: %d of %d windows carry the motif",
                    sum(params$amidation), nrow(params)))
    rm(.ws, .mpos, .aa, .amid, .chk)
  }
}

## ---- known peptide ends ----------------------------------------------------
## A db window is built with window_origin = db_ind, the dibasic site 9.2 found
## within db_spacer < 6 of a peptide boundary, and wN/wC are fixed offsets from
## it (C: -30/+5, N: -5/+30). So the anchor is recoverable from `peps` alone --
## and recovered from the NEAR edge in each case (wC for C-targets, wN for N),
## since that edge is only 5 residues from the anchor and is therefore the one
## least likely to have been clamped at a protein terminus.
##
## A window is then a known peptide end when some reference peptide of the same
## gene has its boundary 1..known_end_gap residues from that anchor -- the C
## branch measures to the peptide's END, the N branch to its START, mirroring
## how 9.2 computed db_spacer in each direction.
if (file.exists(known_end_ref)) {
  .im <- readRDS("/Users/kbrulois/R_projects/ligandFinder/data/id_mapping.rds")
  .a2s <- setNames(.im[["Gene Names (primary)"]], .im[["Entry"]])
  .e2s <- setNames(.im[["Gene Names (primary)"]], .im[["Entry Name"]])

  ref <- readRDS(known_end_ref) %>%
    dplyr::mutate(gene  = dplyr::coalesce(.a2s[accession],
                                          .e2s[paste0(uniprot_name, "_HUMAN")]),
                  start = as.integer(start), end = as.integer(end),
                  name  = dplyr::coalesce(as.character(final_name), "unnamed")) %>%
    dplyr::filter(!is.na(gene), !is.na(start), !is.na(end),
                  end - start + 1L <= known_end_len)

  .co <- stringr::str_match(params$peps, "_w(\\d+)-(\\d+)$")
  .wN <- as.integer(.co[, 2]); .wC <- as.integer(.co[, 3])
  .isC <- as.character(params$target) %in% c("C", "loop_C")

  ## TWO candidate anchors, one recovered from each edge. Either edge can have
  ## been clamped at a protein terminus, and which one cannot be told from the
  ## coordinates alone -- a clamped wC reads low, a clamped wN reads high, and
  ## both look like a plain disagreement. So derive both and accept a boundary
  ## that matches either. (QRFP and NPFF are the concrete case: 136 aa and 113 aa
  ## precursors whose C-terminal windows are clamped, which the near-edge-only
  ## rule missed.) They agree whenever the window is a full 36 residues.
  anchor_a <- ifelse(.isC, .wC -  5L, .wN +  5L)
  anchor_b <- ifelse(.isC, .wN + 30L, .wC - 30L)

  ## boundary the anchor is measured against, per direction
  bnd <- dplyr::bind_rows(
    ref %>% dplyr::transmute(gene, name, side = "C", pos = end),
    ref %>% dplyr::transmute(gene, name, side = "N", pos = start))

  key <- paste0(params$gene, "|", ifelse(.isC, "C", "N"))
  bnd$key <- paste0(bnd$gene, "|", bnd$side)
  bl <- split(bnd[, c("pos", "name")], bnd$key)

  hit <- vapply(seq_along(key), function(i) {
    b <- bl[[key[i]]]
    if (is.null(b)) return(NA_character_)
    best <- NA_character_; bestd <- Inf
    for (a in c(anchor_a[i], anchor_b[i])) {
      if (is.na(a)) next
      d  <- if (.isC[i]) a - b$pos else b$pos - a
      ok <- which(d >= 1L & d <= known_end_gap)
      if (length(ok) && min(d[ok]) < bestd) {
        bestd <- min(d[ok]); best <- b$name[ok[which.min(d[ok])]]
      }
    }
    best
  }, character(1))

  params$known_end      <- !is.na(hit)
  params$known_end_name <- hit

  ## Three strata for the score panel. "GPCR peptide" is a matched reference
  ## ligand that has a receptor annotated (GtoPdb / GPCRdb); "peptide" is a
  ## matched ligand with none. The distinction is a property of the LIGAND, so
  ## it is joined on known_end_name rather than recomputed per window.
  .rcpt <- readRDS(known_end_ref) %>%
    dplyr::transmute(name = dplyr::coalesce(as.character(final_name), "unnamed"),
                     has_receptor = !is.na(receptor)) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(has_receptor = any(has_receptor), .groups = "drop")
  params <- dplyr::left_join(params, .rcpt, by = c("known_end_name" = "name"))
  params$stratum <- factor(
    dplyr::case_when(!params$known_end              ~ "unknown",
                     params$has_receptor %in% TRUE  ~ "GPCR peptide",
                     TRUE                           ~ "peptide"),
    levels = c("unknown", "peptide", "GPCR peptide"))
  message("stratum: ", paste(sprintf("%s=%d", names(table(params$stratum)),
                                     table(params$stratum)), collapse = "  "))
  rm(.rcpt)
  message(sprintf("known peptide ends: %d windows across %d genes (%d also in the training set)",
                  sum(params$known_end), dplyr::n_distinct(params$gene[params$known_end]),
                  sum(params$known_end & params$known == 1)))
  ## sanity: the training knowns are a subset of this rule by construction, so a
  ## low recovery rate means the anchor arithmetic is off, not that data is missing
  .rec <- if (sum(params$known == 1)) mean(params$known_end[params$known == 1]) else NA_real_
  message(sprintf("  recovers %.0f%% of the %d training knowns",
                  100 * .rec, sum(params$known == 1)))
  if (!is.na(.rec) && .rec < 0.8)
    warning("known-end rule recovers <80% of the training knowns -- check the anchor offsets",
            immediate. = TRUE)
  rm(.im, .a2s, .e2s, .co, .wN, .wC, .isC, .rec)
} else {
  message("no ", known_end_ref, " -- skipping the known-peptide-end annotation")
}

## ---- final column order + write --------------------------------------------
## labels, then scores, then reduced-dim coords, then summaries, positions last.
params <- params %>%
  mutate(uniprot_known = ifelse(known == 1, gene, NA_character_)) %>%   # label knowns only
  dplyr::relocate(dplyr::any_of(c("uniprot_known", "known_end", "known_end_name",
                                  "amidation", "pred_nn", "pred_raw")),
                  dplyr::any_of(c("UMAP1", "UMAP2")),
                  dplyr::matches("^UMAP[12]_"),
                  dplyr::matches("^PC[0-9]+$"),
                  .after = terminus)

## Only on the full build. On the fast path `params` is a column subset of this
## very file, so writing it back would silently discard the ~1,000 per-position
## columns that the UMAP is built from.
if (!.from_csv) {
  write_csv_fast(params, out_csv)
  message("wrote ", out_csv, "  (", nrow(params), " x ", ncol(params), ")")
} else {
  message("replot: leaving ", out_csv, " alone")
}

## ---- UMAP plot per metric: 2-D density + model score ------------------------
## Points are coloured by the model score, with 2-D density contours over the top
## so the shape of the window population stays readable where the points saturate.
## Knowns are drawn as open rings rather than filled dots -- a solid marker would
## hide the score colour underneath, which is the thing being looked at.

## ---- reader explanation ----------------------------------------------------
## Placed ABOVE the search box and open by default, inside a <details> the reader
## can collapse. Position matters here: the caveats (input-space layout, the
## known-peptide floor) change how the picture should be read, so they cannot sit
## below it -- but a permanently-expanded seven-paragraph block would push the
## plot itself off the first screen. <details open> gives the first-time reader
## everything and lets a returning one fold it away in one click.
gene_intro_html <- '
<details open style="font:13px/1.55 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;
         max-width:1100px;margin:14px auto 0;padding:10px 14px;
         border:1px solid #d8d8d8;border-radius:6px;background:#fff">
  <summary style="cursor:pointer;font-weight:700;font-size:15px;outline:none">
    Peptide-window UMAP &mdash; how to read it
  </summary>
  <div style="margin-top:10px;color:#222">
    <p><b>Each point is one candidate cleavage window:</b> a 36-residue slice of a
    secreted protein, anchored on a dibasic site (KK/KR/RK/RR). 59,184 windows
    were scored.</p>

    <p><b>Position</b> comes from a UMAP of the input features: each window is
    flattened to one row of position &times; feature values &mdash; conservation,
    AlphaMissense, relative solvent accessibility, DSSP secondary structure and
    H-bond energies, backbone angles, amino-acid properties &mdash; then z-scored.
    So two points sit close together when their residue-level biophysical
    profiles are alike.</p>

    <p><b>Colour</b> is the model&rsquo;s window-level prediction score &mdash; the
    1D-CNN ensemble&rsquo;s global ranking output, averaged over 5 seeds. Dark =
    high. Circles are C-terminal windows, triangles N-terminal.</p>

    <p><b>Score panel</b> (right) shows the same score for <i>every</i> scored
    window, including those outside the gated UMAP, split by what is known about
    the window: unknown (divided at 0.1 so the few that score are not buried
    under the 58k that do not), a known peptide end, or a known GPCR-peptide end.
    Known ends come from UniProt, GPCRdb and Guide to Pharmacology; read them as
    a floor, not a census &mdash; the scored set is restricted to dibasic-anchored
    windows, so the chemokine family, ADM, AVP, APLN and other peptides lacking
    dibasic sites are absent. Colour there marks the training/validation set.</p>

    <p><b>Interaction.</b> Windows scoring above 0.1 are live, in both panels:
    hover one for its window id, score, peptide name if known, and amidation
    status (a glycine immediately 5&prime; of the dibasic pair, &hellip;X-G |
    K/R-K/R, the signal for a C-terminally amidated peptide); the UMAP and the
    score panel are linked, so hovering a window in one highlights it in the
    other; click to open that protein&rsquo;s per-residue page at the window.
    Windows at or below 0.1 &mdash; the great majority &mdash; are drawn but
    static. The search box takes a gene symbol and highlights its 5
    best-scoring windows in cyan, listing them below; windows that rank in the top
    5 but fall outside the gate are listed greyed rather than dropped. Once in the
    per-window page, click &ldquo;Visualize in ChimeraX&rdquo; to view predictions
    on the AlphaFold-DB structure of the precursor protein.</p>
  </div>
</details>'

## ---- gene search widget template -------------------------------------------
## sprintf slots, in order: 1 datalist <option>s, 2 lookup JSON, 3 link base,
## 4 score column name, 5 top-N, 6 highlight colour, 7 interactive threshold
## (text only -- it labels the static entries in the result list).
##
## Highlighting draws an OVERLAY circle into the same parent <g> as the matched
## point rather than restyling the point itself. Same parent means the same
## transform, so cx/cy copy across untouched; and clearing is one removal pass
## with nothing to restore and no z-order left mutated. (A resize probe confirmed
## girafe scales by viewBox and never re-creates these nodes, so overlays stick.)
gene_search_template <- '
<div style="font:13px/1.4 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;
            max-width:1100px;margin:10px auto 0;padding:8px 10px;
            border:1px solid #d8d8d8;border-radius:6px;background:#fafafa">
  <label for="gene-search" style="font-weight:600;margin-right:6px">Gene</label>
  <input id="gene-search" list="gene-list" placeholder="e.g. NPY" autocomplete="off"
         style="padding:4px 7px;border:1px solid #bbb;border-radius:4px;width:190px">
  <button id="gene-clear" type="button"
          style="margin-left:6px;padding:4px 9px;border:1px solid #bbb;
                 border-radius:4px;background:#fff;cursor:pointer">clear</button>
  <span id="gene-msg" style="margin-left:10px;color:#555"></span>
  <datalist id="gene-list">%s</datalist>
  <ol id="gene-results" style="margin:8px 0 0;padding-left:22px"></ol>
</div>
<script>
(function(){
  var LUT = %s, BASE = "%s", SCORE = "%s", TOPN = %d, HIT = "%s", MIN = "%s";
  var IDX = null;

  // Built lazily, not at parse time: this script runs before girafe has
  // rendered the svg, so there would be nothing to index yet.
  function index() {
    if (IDX && IDX.size) return IDX;
    IDX = new Map();
    // BOTH shapes: C-terminal windows render as <circle>, N-terminal ones as
    // <polygon> triangles. Indexing only circles silently made every N-terminal
    // window unsearchable.
    document.querySelectorAll("circle[data-id],polygon[data-id]").forEach(function(c){
      var k = c.getAttribute("data-id");
      if (!IDX.has(k)) IDX.set(k, []);
      IDX.get(k).push(c);
    });
    return IDX;
  }
  // A Map lookup rather than an attribute selector: ids carry "/" and "_"
  // (UNQ6494/PRO21346, HERVK_113) and would need escaping in a selector.

  function clearHits(){
    document.querySelectorAll("circle.gene-hit").forEach(function(e){ e.remove(); });
  }

  function light(peps){
    var els = index().get(peps) || [], n = 0;
    els.forEach(function(c){
      var o = document.createElementNS("http://www.w3.org/2000/svg","circle");
      o.setAttribute("class","gene-hit");
      // centre via getBBox rather than cx/cy: a triangle is a <polygon> and has
      // no cx/cy at all. getBBox is in the element own user space, and the
      // overlay goes into the same parent, so the transform already matches.
      var bb = c.getBBox();
      o.setAttribute("cx", bb.x + bb.width  / 2);
      o.setAttribute("cy", bb.y + bb.height / 2);
      o.setAttribute("r", 5);
      o.setAttribute("fill","none");
      o.setAttribute("stroke", HIT);
      o.setAttribute("stroke-width", 1.8);
      c.parentNode.appendChild(o);   // same parent => same transform, drawn last
      n++;
    });
    return n;
  }

  var inp = document.getElementById("gene-search"),
      out = document.getElementById("gene-results"),
      msg = document.getElementById("gene-msg");

  function run(){
    clearHits(); out.innerHTML = ""; msg.textContent = "";
    var g = (inp.value || "").trim();
    if (!g) return;
    var rows = LUT[g] || LUT[g.toUpperCase()];
    if (!rows){ msg.textContent = "no windows for " + g; return; }
    if (!Array.isArray(rows)) rows = [rows];

    var lit = 0, gated = 0, stat = 0;
    rows.forEach(function(r){
      if (r.o === 1) { gated++; } else if (r.o === 2) { stat++; } else { lit += light(r.p); }
      var li = document.createElement("li");
      var a  = document.createElement("a");
      // gene taken from the window id, not from what was typed: a lowercase or
      // otherwise off-case query still has to build the real page URL.
      // lastIndexOf rather than a regex -- this lives in an R string literal,
      // where a regex escape is itself a parse error and not worth the trouble.
      a.href = BASE + "/" + r.p.slice(0, r.p.lastIndexOf("_w")) + ".html#" + r.p;
      a.target = "_blank";
      a.textContent = r.p;
      li.appendChild(a);
      var tail = "  " + SCORE + " " + r.s + (r.f ? "  [" + r.f + "]" : "");
      li.appendChild(document.createTextNode(tail));
      if (r.o){
        li.style.color = "#999";
        li.appendChild(document.createTextNode(
          r.o === 2 ? "  (score <= " + MIN + " - drawn static, not highlightable)"
                    : "  (outside gate - not drawn)"));
      }
      out.appendChild(li);
    });
    msg.textContent = "top " + rows.length + " of " + g +
                      " by " + SCORE + " - " + lit + " highlighted" +
                      (gated ? ", " + gated + " outside the gate" : "") +
                      (stat  ? ", " + stat  + " static (score <= " + MIN + ")" : "");
  }

  inp.addEventListener("input", run);
  inp.addEventListener("change", run);
  document.getElementById("gene-clear").addEventListener("click", function(){
    inp.value = ""; run();
  });
})();
</script>'

if (make_plot) {
  sc <- intersect(c(score_col, "pred_nn", "pred_raw"), names(params))[1]
  ## readable name for whichever column won
  .sc_lab <- function(x) if (!is.na(x) && x %in% names(score_label))
                           unname(score_label[[x]]) else x

  if (is.na(sc) || all(is.na(params[[sc]]))) {
    message("no usable score column (", score_col, ") -- colouring by win_type instead")
    sc <- NA
  } else {
    message("colouring UMAP by ", sc, " (", sum(!is.na(params[[sc]])), " scored windows)")
  }

  ## resolved once, not per metric
  use_glow <- isTRUE(contour_glow) && !is.null(density_bins) &&
              requireNamespace("ggfx", quietly = TRUE)
  ## The same halo goes behind the known-peptide rings. That one has nothing to
  ## do with whether contours are drawn, so it must not inherit the density_bins
  ## condition -- setting density_bins <- NULL should drop the contours, not
  ## silently strip the halo off the knowns as well.
  use_glow_marks <- isTRUE(contour_glow) && requireNamespace("ggfx", quietly = TRUE)
  raster_ok      <- requireNamespace("ggrastr", quietly = TRUE) &&
                    requireNamespace("ragg", quietly = TRUE)
  if (isTRUE(contour_glow) && !use_glow && !is.null(density_bins))
    message("ggfx not installed -- plain contours. install.packages('ggfx') for the halo")

  ## Which layouts to draw: the input-space UMAPs, plus the embedding one when
  ## its columns are there. Driven by what is actually in `params` rather than by
  ## the options alone, so this section can be re-run on its own -- or off the
  ## written CSV -- and still find every layout that exists.
  plot_layouts <- intersect(c(umap_metrics, "embed"),
                            sub("^UMAP1_", "", grep("^UMAP1_", names(params), value = TRUE)))

  for (mt in plot_layouts) {
    u1 <- paste0("UMAP1_", mt); u2 <- paste0("UMAP2_", mt)

    ## what this layout is a picture OF -- the two spaces are not comparable and
    ## the panel should not let you forget which one you are looking at
    space_lab <- if (mt == "embed") paste0("model embedding, ", embed_metric)
                 else               paste0("input features, ", mt)

    ## Draw order IS the z-order in ggplot: rows plotted last land on top. Sort
    ## ascending by score so the top-scoring windows survive overplotting instead
    ## of being buried under the low-scoring bulk. NAs first (bottom of the pile).
    pdat <- if (is.na(sc)) params else params[order(params[[sc]], na.last = FALSE), ]
    ## windows with no coordinate in THIS layout (unmatched embedding rows)
    pdat <- pdat[!is.na(pdat[[u1]]) & !is.na(pdat[[u2]]), , drop = FALSE]

    ## Gate. Subsetting rather than coord_cartesian(): the contours are a density
    ## OF the plotted windows, so zooming would leave them describing the whole
    ## cloud while showing only part of it, and the labelled "top" knowns would
    ## still be picked from windows outside the box.
    if (!is.null(gate_umap1) || !is.null(gate_umap2)) {
      n_before <- nrow(pdat)
      if (!is.null(gate_umap1))
        pdat <- pdat[pdat[[u1]] > min(gate_umap1) & pdat[[u1]] < max(gate_umap1), , drop = FALSE]
      if (!is.null(gate_umap2))
        pdat <- pdat[pdat[[u2]] > min(gate_umap2) & pdat[[u2]] < max(gate_umap2), , drop = FALSE]
      message(sprintf("  [%s] gate: %d of %d windows kept", mt, nrow(pdat), n_before))
      if (nrow(pdat) == 0) { message("  [", mt, "] gate is empty -- skipping"); next }
    }

    ## Every layer above the points is the same in the static svg and in the
    ## interactive html; only the point layer differs, and it has to sit at the
    ## bottom of the stack. So build the plot FROM a point layer instead of
    ## adding interactive points on top of an already-finished plot.
    ## facet column first, so knowns and their labels inherit it and land in the
    ## right panel rather than being drawn into both
    faceted <- isTRUE(facet_amidation) && "amidation" %in% names(pdat)
    if (faceted) {
      lv <- c("no amidation motif", "amidation motif (G | dibasic)")
      pdat$.facet <- factor(lv[1L + as.integer(pdat$amidation)], levels = lv)
      message(sprintf("  [%s] panels: %s", mt,
                      paste(sprintf("%s=%d", lv, tabulate(pdat$.facet, 2L)), collapse = "  ")))
    }

    ## Ring every window at a known peptide end, not just the training positives:
    ## the positives are a deliberately narrow selection (cognate receptor model,
    ## terminal insertion, dibasic anchor), so ringing only those understates
    ## where real peptide boundaries fall. Training membership is still readable
    ## -- it goes in the hover text.
    knowns <- if ("known_end" %in% names(pdat)) pdat[which(pdat$known_end), , drop = FALSE]
              else                              dplyr::filter(pdat, known == 1)

    decorate <- function(point_layer) {
      p <- ggplot(pdat, aes(.data[[u1]], .data[[u2]])) + point_layer

      ## 2-D density contours. A thin dark line disappears over the saturated point
      ## field, so ggfx paints a contrasting halo around it (with_outer_glow): the
      ## contour then reads over both the pale low-score bulk and the dark
      ## high-score points, without having to thicken the line until it obscures
      ## the data.
      ##
      ## The layer is then rasterised explicitly through ggrastr with dev="ragg"
      ## at contour_dpi, the way scTools::plot_subsets() does it. ggfx already
      ## rasterises whatever it filters, but it picks the device and resolution
      ## itself; going through ragg pins both, so the halo is resampled by a
      ## known AGG renderer at a resolution we choose rather than at whatever
      ## the active graphics device happens to be.
      if (!is.null(density_bins)) {
        dens <- stat_density_2d(color = contour_col, linewidth = contour_lw,
                                alpha = 0.9, bins = density_bins)
        if (use_glow)
          dens <- ggfx::with_outer_glow(dens, colour = glow_col,
                                        sigma = glow_sigma, expand = glow_expand)
        if (raster_ok)
          dens <- ggrastr::rasterise(dens, dev = "ragg", dpi = contour_dpi)
        else
          message("ggrastr/ragg not both installed -- marks left to ggfx's own rasteriser")
        p <- p + dens
      }

      if (!is.na(sc))
        p <- p + scale_color_viridis_c(option = "magma", direction = -1,
                                       na.value = "grey88", name = .sc_lab(sc))

      ## Annotation rings stay a fixed open circle rather than following this
      ## scale: one shape scale cannot serve both solid data marks (16/17) and
      ## open rings (21/24), and a ring reads as an annotation either way.
      p <- p + scale_shape_manual(values = term_shapes, name = "terminus",
                                  na.value = 16)

      ## Known-peptide rings get the contour treatment: a halo so a thin ring
      ## stays visible against the dark high-score points it sits on, then the
      ## same ragg rasterisation so every mark layer is resampled by one renderer
      ## at one resolution. They are not the interactive layer, so rasterising
      ## them costs no hover or click.
      ring_layer <- function(d, col, size) {
        g <- geom_point(data = d, shape = 21, fill = NA,
                        color = col, stroke = 0.7, size = size)
        if (use_glow_marks)
          g <- ggfx::with_outer_glow(g, colour = glow_col,
                                     sigma = glow_sigma, expand = glow_expand)
        if (raster_ok) g <- ggrastr::rasterise(g, dev = "ragg", dpi = contour_dpi)
        g
      }

      if (isTRUE(known_rings)) p <- p + ring_layer(knowns, known_col, 2.4)
      p <- p +
        theme_bw() +
        labs(title = paste0("Peptide-window UMAP -- ", space_lab),
             subtitle = paste0(if (is.na(sc)) "coloured by win_type" else paste0("coloured by ", .sc_lab(sc)),
                               "; contours = 2-D density of the plotted windows",
                               if (isTRUE(known_rings))
                                 paste0("; ", known_col, " rings = known peptide ends (named on hover)") else "",
                               "; points = C-terminal, triangles = N-terminal",
                               if ("amidation" %in% names(pdat))
                                 "; amidation motif shown on hover" else "",
                               if (faceted) "; panels split by motif" else ""),
             x = "UMAP1", y = "UMAP2")

      if (faceted) p <- p + facet_wrap(ggplot2::vars(.data[[".facet"]]))

      p
    }

    ## ---- static svg --------------------------------------------------------
    ## An SVG with tens of thousands of vector points is huge and slow to open;
    ## rasterise the point layer. Axes and labels stay vector (so does the contour
    ## line, unless contour_glow is on -- ggfx rasterises what it filters).
    pt <- if (is.na(sc))
            geom_point(aes(color = win_type, shape = terminus), size = point_size, alpha = 0.55)
          else
            geom_point(aes(color = .data[[sc]], shape = terminus), size = point_size, alpha = 0.75)
    if (raster_points && requireNamespace("ggrastr", quietly = TRUE))
      pt <- ggrastr::rasterise(pt, dpi = 300)

    ## a 12x10 canvas split in two gives tall narrow panels; shorten it so each
    ## panel stays roughly square and the layout is still readable
    fig_h    <- if (faceted) 7 else 10
    svg_path <- paste0(svg_stem, "_", mt, ".svg")
    ggsave(svg_path, decorate(pt), width = 12, height = fig_h)
    message("wrote ", svg_path)

    ## ---- clickable version -------------------------------------------------
    ## EVERY window is hoverable and clickable, so the interactive layer IS the
    ## point layer here -- not a top-N overlay on a rasterised bulk. Nothing is
    ## rasterised in this file, so it carries one vector mark per window: expect
    ## it to be a lot bigger and slower to open than the svg.
    ## Each point opens <link_base>/<gene>.html#<peps> -- the per-gene page matches
    ## the fragment against its panels' data-peps and centres that window.
    if (!is.null(link_base) && requireNamespace("ggiraph", quietly = TRUE)) {
      linkable  <- !is.na(pdat$peps) & !is.na(pdat$gene)
      pdat$.url <- ifelse(linkable,
                          sprintf("%s/%s.html#%s", sub("/+$", "", link_base),
                                  pdat$gene, pdat$peps),
                          NA_character_)
      ## ggiraph renders an NA aesthetic as the literal string "NA", so every
      ## interactive attribute needs a real value: a tooltip reading "NA", a
      ## data-id that collides across rows, or onclick='NA' (a ReferenceError
      ## when clicked). Windows without a peps id are unusual but must degrade
      ## to inert points, not broken ones.
      pdat$.lab <- ifelse(is.na(pdat$peps), "(unlabelled window)", pdat$peps)
      pdat$.did <- ifelse(is.na(pdat$peps),
                          paste0("row", seq_len(nrow(pdat))), pdat$peps)
      pdat$.tip <- if (is.na(sc)) pdat$.lab else
        paste0(pdat$.lab, "\n", .sc_lab(sc), ": ", round(pdat[[sc]], 3))
      if ("terminus" %in% names(pdat))
        pdat$.tip <- paste0(pdat$.tip, "\nterminus: ", pdat$terminus)
      ## Both of these are called out only where they apply. A "no" on each of
      ## the other ~24k windows would be noise to read past on hover, and about a
      ## megabyte of tooltip text in the file for no information.
      if ("known_end" %in% names(pdat))
        pdat$.tip <- paste0(pdat$.tip,
                            ifelse(pdat$known_end,
                                   paste0("\nknown peptide end: ",
                                          dplyr::coalesce(pdat$known_end_name, pdat$gene),
                                          ifelse(pdat$known == 1, " (training positive)", "")),
                                   ""))
      else if ("known" %in% names(pdat))
        pdat$.tip <- paste0(pdat$.tip,
                            ifelse(pdat$known == 1,
                                   paste0("\nknown peptide: ", pdat$gene), ""))
      if ("amidation" %in% names(pdat))
        pdat$.tip <- paste0(pdat$.tip,
                            ifelse(pdat$amidation, "\namidation motif (G | dibasic)", ""))

      ## Plain double quotes, NOT &quot;. ggiraph HTML-escapes the attribute
      ## value itself, so a pre-escaped &quot; comes out as &amp;quot; and the
      ## browser hands JavaScript the literal text &quot; -- a syntax error, and
      ## the click silently does nothing. The per-gene pages appear to use
      ## &quot; successfully only because they never rely on ggiraph's onclick:
      ## they have their own delegated document click handler keyed on data-id.
      ## Windows with no gene/peps get no onclick at all rather than a dead link.
      pdat$.click <- ifelse(linkable,
                            sprintf('window.open("%s","_blank")', pdat$.url), "")

      ## Two point layers: the static bulk at or below interactive_min, then the
      ## interactive windows above it. Both map the same colour and shape, so
      ## they share one scale and look identical; the static layer goes first
      ## so the interactive (higher-scoring) marks land on top of it, which is
      ## the order the ascending sort already gave the single layer.
      is_int <- if (is.na(sc) || is.null(interactive_min)) rep(TRUE, nrow(pdat))
                else !is.na(pdat[[sc]]) & pdat[[sc]] > interactive_min
      pdat_int <- pdat[is_int, , drop = FALSE]
      pdat_stc <- pdat[!is_int, , drop = FALSE]
      if (!is.na(sc))
        message(sprintf("  [%s] interactive: %d of %d plotted windows (%s > %s); %d static",
                        mt, nrow(pdat_int), nrow(pdat), sc, format(interactive_min),
                        nrow(pdat_stc)))

      pt_int <- if (is.na(sc))
                  ggiraph::geom_point_interactive(
                    data = pdat_int,
                    aes(color = win_type, shape = terminus, tooltip = .data[[".tip"]],
                        data_id = .data[[".did"]], onclick = .data[[".click"]]),
                    size = point_size, alpha = 0.55)
                else
                  ggiraph::geom_point_interactive(
                    data = pdat_int,
                    aes(color = .data[[sc]], shape = terminus, tooltip = .data[[".tip"]],
                        data_id = .data[[".did"]], onclick = .data[[".click"]]),
                    size = point_size, alpha = 0.75)
      if (nrow(pdat_stc)) {
        pt_stc <- if (is.na(sc))
                    geom_point(data = pdat_stc, aes(color = win_type, shape = terminus),
                               size = point_size, alpha = 0.55)
                  else
                    geom_point(data = pdat_stc, aes(color = .data[[sc]], shape = terminus),
                               size = point_size, alpha = 0.75)
        pt_int <- list(pt_stc, pt_int)
      }

      ## ---- panel 1: score by stratum ---------------------------------------
      ## data_id is the window id, identical to the UMAP layer's, and both go into
      ## ONE girafe below -- that pairing is what links the panels. Drawn over
      ## every scored window rather than the gated subset (see panel_all).
      p_panels <- list()
      if (isTRUE(panel_score) && !is.na(sc) && "stratum" %in% names(params)) {
        sdat <- if (isTRUE(panel_all)) params else pdat
        sdat <- sdat[!is.na(sdat[[sc]]) & !is.na(sdat$stratum), , drop = FALSE]
        ## none first so the 59k grey bulk cannot bury the 1.6k train/val marks
        sdat <- sdat[order(!is.na(sdat$set) & sdat$set != "none"), , drop = FALSE]
        sdat$.did <- ifelse(is.na(sdat$peps),
                            paste0("row", seq_len(nrow(sdat))), sdat$peps)
        sdat$.set <- factor(ifelse(is.na(sdat$set), "none", as.character(sdat$set)),
                            levels = c("train", "val", "none"))
        ## panel-local strata: the unknowns split at panel_unk_split so the handful
        ## that score are not lost in the 59k that do not. `stratum` itself stays
        ## the three-level biological label.
        .lv <- c(paste0("unknown <", panel_unk_split), paste0("unknown >=", panel_unk_split),
                 "peptide", "GPCR peptide")
        sdat$.strat <- factor(dplyr::case_when(
          sdat$stratum == "unknown" & sdat[[sc]] <  panel_unk_split ~ .lv[1],
          sdat$stratum == "unknown"                                  ~ .lv[2],
          TRUE                                                       ~ as.character(sdat$stratum)),
          levels = .lv)
        sdat$.tip <- paste0(
          ifelse(is.na(sdat$peps), "(unlabelled window)", sdat$peps),
          "\n", .sc_lab(sc), ": ", round(sdat[[sc]], 3),
          "\n", ifelse(sdat$terminus %in% "C", "C-terminal window",
                                                "N-terminal window"),
          "\namidation motif: ",
          if ("amidation" %in% names(sdat)) ifelse(sdat$amidation, "yes", "no") else "unknown",
          "\nset: ", as.character(sdat$.set),
          "\nstratum: ", as.character(sdat$stratum))

        ## Same deep link as the UMAP points: <link_base>/<gene>.html#<peps>,
        ## which the per-gene page matches against its panels' data-peps and
        ## centres. Plain double quotes, NOT &quot; -- ggiraph escapes the
        ## attribute itself, so a pre-escaped entity reaches JavaScript as
        ## literal text and the click dies with a syntax error.
        .plink <- !is.na(sdat$peps) & !is.na(sdat$gene)
        sdat$.click <- ifelse(
          .plink,
          sprintf('window.open("%s/%s.html#%s","_blank")',
                  sub("/+$", "", link_base), sdat$gene, sdat$peps),
          "")

        ## quasirandom gives the beeswarm its width-proportional spread; without
        ## ggbeeswarm fall back to jitter rather than dropping the panel.
        swarm <- if (requireNamespace("ggbeeswarm", quietly = TRUE))
                   ggbeeswarm::position_quasirandom(width = panel_swarm_w)
                 else position_jitter(width = panel_swarm_w * 0.8, height = 0)

        ## Points only -- no violin outline. data_id is the window id, shared with
        ## the UMAP layer, so ggiraph highlights the same window in both panels on
        ## hover; that only works because both go into ONE girafe.
        ## Same static/interactive split as the UMAP, at the same threshold. The
        ## two layers are swarmed separately, which matters only where a column
        ## holds both (peptide, GPCR peptide) and only right at the boundary.
        s_int <- if (is.null(interactive_min)) rep(TRUE, nrow(sdat))
                 else !is.na(sdat[[sc]]) & sdat[[sc]] > interactive_min
        p_panels$score <- ggplot(sdat, aes(x = .strat, y = .data[[sc]])) +
          geom_point(data = sdat[!s_int, , drop = FALSE], aes(colour = .set),
                     position = swarm, size = panel_pt_size, alpha = 0.7) +
          ggiraph::geom_point_interactive(
            data = sdat[s_int, , drop = FALSE],
            aes(colour = .set, tooltip = .tip, data_id = .did, onclick = .click),
            position = swarm, size = panel_pt_size, alpha = 0.7) +
          scale_colour_manual(values = set_cols, name = "set", drop = FALSE) +
          ## limits pinned: with two layers the scale is trained on each layer's
          ## USED levels in turn, and a level present only in the second layer
          ## (unknown >= split, all interactive) would drop to the end
          scale_x_discrete(limits = .lv, labels = function(x)
            paste0(x, "\n(n=", tabulate(sdat$.strat, nlevels(sdat$.strat))[match(x, levels(sdat$.strat))], ")")) +
          labs(x = NULL, y = .sc_lab(sc)) +
          theme_bw() +
          theme(legend.position = "bottom",
                panel.grid.minor = element_blank()) +
          guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
        message(sprintf("  [%s] score panel: %d windows  (%s)", mt, nrow(sdat),
                        paste(sprintf("%s=%d", levels(sdat$.strat),
                                      tabulate(sdat$.strat, nlevels(sdat$.strat))),
                              collapse = "  ")))
      }

      ## One girafe over the whole patchwork. width_svg grows with the panel count
      ## so adding a panel makes the UMAP smaller rather than squeezing everything.
      gob <- decorate(pt_int)
      if (length(p_panels)) {
        ## `ggplot + ggplot` only composes once patchwork's namespace is loaded
        ## (it registers the ggplot_add method). The `::` on plot_layout() below
        ## does load it -- but R evaluates the Reduce() on its left first, so in
        ## a session where nothing else has loaded patchwork yet the Reduce()
        ## fails with "Can't add `x[[i]]` to a <ggplot> object".
        if (!requireNamespace("patchwork", quietly = TRUE))
          stop("install.packages('patchwork') to draw the side panels")
        gob <- Reduce(`+`, p_panels, init = gob) +
               patchwork::plot_layout(widths = c(umap_rel_w,
                                                 rep(panel_rel_w, length(p_panels))))
      }
      gir <- ggiraph::girafe(
        ggobj = gob,
        width_svg = 12 + 4 * length(p_panels), height_svg = fig_h,
        options = list(
          ggiraph::opts_sizing(rescale = TRUE),
          ggiraph::opts_selection(type = "none"),
          ggiraph::opts_hover(css = hover_css),
          ggiraph::opts_hover_inv(css = hover_inv_css)
        ))

      ## ---- gene search UI --------------------------------------------------
      ## Ships a PRECOMPUTED top-N per gene rather than every window: the page
      ## only ever needs the head of each gene's list, which turns a ~59k-row
      ## payload into ~3k x N and makes the lookup a plain object access.
      ## Ranking is on the colour column, so the highlight agrees with what is
      ## on screen. Windows that rank in a gene's top N but fall outside the gate
      ## (or the plotted layout) are still listed, greyed -- dropping them would
      ## quietly under-report a gene.
      search_ui <- NULL
      if (isTRUE(gene_search) && !is.na(sc) &&
          requireNamespace("jsonlite", quietly = TRUE)) {

        pool <- params[!is.na(params[[u1]]) & !is.na(params[[u2]]) &
                       !is.na(params[[sc]]), , drop = FALSE]
        pool <- pool[!duplicated(pool$peps), , drop = FALSE]   # peps repeat across targets
        pool <- pool[order(-pool[[sc]]), , drop = FALSE]
        ## why a listed window cannot be lit: 0 = drawn and interactive,
        ## 1 = outside the gate (not drawn), 2 = drawn but static (no data-id
        ## to find, so nothing to overlay the highlight on)
        pool$.out <- dplyr::case_when(!pool$peps %in% pdat$peps     ~ 1L,
                                      !pool$peps %in% pdat_int$peps ~ 2L,
                                      TRUE                          ~ 0L)
        pool$.pan <- if (faceted && "amidation" %in% names(pool))
                       ifelse(pool$amidation, "motif", "no motif") else ""

        top <- pool %>% dplyr::group_by(gene) %>%
          dplyr::slice_head(n = gene_search_n) %>% dplyr::ungroup()

        lut <- lapply(split(seq_len(nrow(top)), top$gene), function(ix)
          data.frame(p = top$peps[ix], s = round(top[[sc]][ix], 3),
                     f = top$.pan[ix], o = top$.out[ix], stringsAsFactors = FALSE))

        genes_sorted <- sort(unique(top$gene))
        message(sprintf("  [%s] gene search: %d genes indexed, %d windows",
                        mt, length(lut), nrow(top)))

        search_ui <- htmltools::HTML(sprintf(
          gene_search_template,
          paste0(sprintf("<option value=\"%s\"></option>",
                         htmltools::htmlEscape(genes_sorted, attribute = TRUE)),
                 collapse = ""),
          as.character(jsonlite::toJSON(lut, dataframe = "rows", auto_unbox = TRUE)),
          sub("/+$", "", link_base), sc, gene_search_n, gene_hit_col,
          if (is.null(interactive_min)) "" else format(interactive_min)))
      }

      html_path <- if (!is.null(html_name) && length(plot_layouts) == 1) {
                     file.path(dirname(svg_stem), html_name)
                   } else {
                     paste0(svg_stem, "_", mt, ".html")
                   }
      htmltools::save_html(
        htmltools::tagList(htmltools::HTML(gene_intro_html), search_ui, gir),
        file = html_path, libdir = "dependency_files")
      ## same inlining the per-gene pages use, so the file stands alone
      if (exists("nn_inline_deps", mode = "function"))
        nn_inline_deps(html_path, file.path(dirname(html_path), "dependency_files"))
      message("wrote ", html_path, "  (", nrow(pdat), " windows drawn, ",
              nrow(pdat_int), " interactive, ", sum(linkable & is_int), " clickable)")
    }
  }
}
