# R translation of ~/Python_scripts/cap_peptides.py plus a converter that turns
# a peptide window (e.g. "CAMP_w125-160" from nn_input$peps) into a mature
# peptide: remove the 7-residue inserting-end padding/cleavage/DB block, then
# score the non-inserting end with the cap logic to drop excess residues while
# keeping dibasic sites (KK|KR|RK|RR) out of the mature peptide.
#
# Scanner / native-cap functions are generic over the terminus being capped:
#   scan_terminal_caps(peptide, terminus = "N" | "C", ...)
#   compute_native_cap(peptide, terminus = "N" | "C")
# with back-compat wrappers scan_n_term_caps() / compute_native_ncap().
#
# The Jupyter-only helpers in the Python original (read_peptides_from_clipboard,
# print_alignment_block_from_cterm, analyze_clipboard_peptides) are not ported.


# ---- Scoring constants (match cap_peptides.py) ----------------------------

KD_HYDRO <- c(
  I = 4.5, V = 4.2, L = 3.8, F = 2.8, C = 2.5, M = 1.9,
  A = 1.8, G = -0.4, T = -0.7, S = -0.8, W = -0.9, Y = -1.3,
  P = -1.6, H = -3.2, E = -3.5, Q = -3.5, D = -3.5, N = -3.5,
  K = -3.9, R = -4.5
)

AROMATIC <- c("F", "W", "Y")
ALIPHATIC <- c("L", "I", "V", "M")
HYDROPHOBIC <- union(AROMATIC, ALIPHATIC)
POLAR_RES <- c("D", "E", "S", "T", "N", "Q", "G")

DIST_MIN_DEFAULT <- 15L
DIST_MAX_DEFAULT <- 25L
MIN_PEPTIDE_LEN <- 7L
TOP_K_CAPS <- 3L

DIBASIC_REGEX <- "KK|KR|RK|RR"


# ---- Low-level helpers ----------------------------------------------------

#' Clean a peptide string (uppercase, drop non-letters).
#' @keywords internal
clean_peptide <- function(seq) {
  seq <- toupper(as.character(seq))
  gsub("[^A-Z]", "", seq)
}

#' Average Kyte-Doolittle hydrophobicity.
#' @keywords internal
avg_hydrophobicity <- function(seq) {
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return(0)
  vals <- KD_HYDRO[chars]
  vals[is.na(vals)] <- 0
  mean(vals)
}

#' Count polar residues (D,E,S,T,N,Q,G).
#' @keywords internal
count_polar <- function(seq) {
  sum(strsplit(seq, "", fixed = TRUE)[[1]] %in% POLAR_RES)
}

#' Count strongly hydrophobic residues (F,W,Y,L,I,V,M).
#' @keywords internal
count_strong_hydrophobic <- function(seq) {
  sum(strsplit(seq, "", fixed = TRUE)[[1]] %in% HYDROPHOBIC)
}

#' Net charge at ~pH 7 (K,R = +1; D,E = -1; H ignored).
#' @keywords internal
net_charge <- function(seq) {
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  sum(chars %in% c("K", "R")) - sum(chars %in% c("D", "E"))
}

#' Longest consecutive hydrophobic run.
#' @keywords internal
longest_hydrophobic_run <- function(seq) {
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  if (length(chars) == 0) return(0L)
  r <- rle(chars %in% HYDROPHOBIC)
  runs <- r$lengths[r$values]
  if (length(runs) == 0) 0L else as.integer(max(runs))
}

#' Cap length rule (3 for short peptides, else 4).
#' @keywords internal
choose_cap_len <- function(peptide_len) {
  if (peptide_len <= 15) 3L else 4L
}

#' Reverse a string.
#' @keywords internal
str_reverse <- function(x) {
  vapply(strsplit(x, "", fixed = TRUE),
         function(chars) paste(rev(chars), collapse = ""),
         character(1))
}


# ---- Core scoring ---------------------------------------------------------

#' Compute cap features and cap_score for a peptide cap sequence.
#'
#' Lower cap_score indicates a better N-terminal blocking cap (less
#' insertion-prone). Returns a one-row tibble so callers can `bind_rows()`
#' across multiple caps.
#'
#' @param cap_seq Character. A short peptide cap sequence.
#' @return A one-row tibble with columns `cap_seq, n_aromatic, n_W,
#'   n_aliphatic, n_hydrophobic, hydro_run, n_polar, n_neg, n_pos,
#'   net_charge, has_pro, avgKD, cap_score`.
#' @export
compute_cap_features <- function(cap_seq) {
  cap_seq <- toupper(as.character(cap_seq))
  chars <- strsplit(cap_seq, "", fixed = TRUE)[[1]]

  n_aromatic    <- sum(chars %in% AROMATIC)
  n_W           <- sum(chars == "W")
  n_aliphatic   <- sum(chars %in% ALIPHATIC)
  n_hydrophobic <- sum(chars %in% HYDROPHOBIC)
  hydro_run     <- longest_hydrophobic_run(cap_seq)
  n_polar       <- sum(chars %in% POLAR_RES)
  n_neg         <- sum(chars %in% c("D", "E"))
  n_pos         <- sum(chars %in% c("K", "R"))
  net_ch        <- n_pos - n_neg
  has_pro       <- as.integer("P" %in% chars)
  avgKD         <- avg_hydrophobicity(cap_seq)

  cap_score <-
      5 * n_aromatic    +
      3 * n_W           +
      2 * n_aliphatic   +
      1 * hydro_run     +
      1 * max(net_ch, 0) -
      2 * n_neg         -
      1 * n_polar       -
      2 * has_pro

  tibble::tibble(
    cap_seq       = cap_seq,
    n_aromatic    = n_aromatic,
    n_W           = n_W,
    n_aliphatic   = n_aliphatic,
    n_hydrophobic = n_hydrophobic,
    hydro_run     = hydro_run,
    n_polar       = n_polar,
    n_neg         = n_neg,
    n_pos         = n_pos,
    net_charge    = net_ch,
    has_pro       = has_pro,
    avgKD         = avgKD,
    cap_score     = cap_score
  )
}

#' Map cap_score to a qualitative label, clamping W-containing caps.
#'
#' @param cap_score Numeric.
#' @param cap_seq Character.
#' @return One of "great", "good", "fair", "poor".
#' @export
qualitative_from_score <- function(cap_score, cap_seq) {
  label <- if (cap_score <= -2) "great"
           else if (cap_score <= 1)  "good"
           else if (cap_score <= 4)  "fair"
           else                      "poor"

  if (grepl("W", cap_seq, fixed = TRUE)) {
    if (label %in% c("great", "good")) label <- "fair"
  }
  label
}

#' Integer tier from qualitative label (0=great ... 3=poor).
#' @param label Character.
#' @export
tier_from_label <- function(label) {
  unname(c(great = 0L, good = 1L, fair = 2L, poor = 3L)[label])
}


# ---- Rank utilities -------------------------------------------------------

#' Build the lexicographic rank key used to sort caps.
#' @keywords internal
make_rank_key <- function(tier, cap_score, n_aromatic, n_W, n_aliphatic,
                          n_polar, net_charge) {
  c(tier,
    cap_score,
    n_aromatic + n_W,
    n_aliphatic,
    -n_polar,
    max(net_charge, 0))
}

#' 1-based rank where a new_key would be inserted into a sorted (best→worst)
#' list of rank keys.
#' @param rank_keys List of numeric vectors, already sorted best→worst.
#' @param new_key Numeric vector of the same length.
#' @return Integer.
#' @export
rank_position <- function(rank_keys, new_key) {
  for (i in seq_along(rank_keys)) {
    k <- rank_keys[[i]]
    cmp <- which(new_key != k)
    if (length(cmp) > 0 && new_key[cmp[1]] < k[cmp[1]]) return(i)
  }
  length(rank_keys) + 1L
}


# ---- Scanner --------------------------------------------------------------

#' Core N-terminal scan. Operates only on the given orientation; callers
#' reverse first if they want C-terminal scanning.
#' @keywords internal
scan_n_caps_core <- function(peptide, dist_min, dist_max) {
  n <- nchar(peptide)
  if (n < MIN_PEPTIDE_LEN) return(tibble::tibble())

  cap_len <- choose_cap_len(n)

  if (n < dist_min) {
    start_indices <- 0L
  } else {
    start_min <- max(0L, n - dist_max)
    start_max <- max(0L, n - dist_min)
    if (start_min > start_max) return(tibble::tibble())
    start_indices <- seq.int(start_min, start_max)
  }

  records <- vector("list", length(start_indices))
  for (i in seq_along(start_indices)) {
    start <- start_indices[i]
    subpep <- substr(peptide, start + 1L, n)
    if (nchar(subpep) < cap_len) next

    cap <- substr(subpep, 1L, cap_len)
    feats <- compute_cap_features(cap)
    label <- qualitative_from_score(feats$cap_score, feats$cap_seq)
    tier  <- tier_from_label(label)

    rk <- make_rank_key(tier, feats$cap_score, feats$n_aromatic, feats$n_W,
                        feats$n_aliphatic, feats$n_polar, feats$net_charge)

    records[[i]] <- tibble::tibble(
      peptide             = peptide,
      peptide_length      = n,
      cap_len             = cap_len,
      start_index_0based  = start,
      start_index_1based  = start + 1L,
      end_index_1based    = n,
      modeled_subpeptide  = subpep,
      modeled_length      = nchar(subpep),
      cap_seq             = feats$cap_seq,
      cap_score           = feats$cap_score,
      avgKD               = feats$avgKD,
      n_aromatic          = feats$n_aromatic,
      n_W                 = feats$n_W,
      n_aliphatic         = feats$n_aliphatic,
      n_hydrophobic       = feats$n_hydrophobic,
      hydro_run           = feats$hydro_run,
      n_polar             = feats$n_polar,
      n_neg               = feats$n_neg,
      n_pos               = feats$n_pos,
      net_charge          = feats$net_charge,
      has_pro             = feats$has_pro,
      qualitative         = label,
      tier                = tier,
      rank_key            = list(rk)
    )
  }

  records <- records[!vapply(records, is.null, logical(1))]
  if (length(records) == 0) return(tibble::tibble())

  df <- dplyr::bind_rows(records)
  ord <- order(vapply(df$rank_key, `[`, numeric(1), 1),
               vapply(df$rank_key, `[`, numeric(1), 2),
               vapply(df$rank_key, `[`, numeric(1), 3),
               vapply(df$rank_key, `[`, numeric(1), 4),
               vapply(df$rank_key, `[`, numeric(1), 5),
               vapply(df$rank_key, `[`, numeric(1), 6))
  df <- df[ord, , drop = FALSE]
  df <- tibble::add_column(df, cap_rank = seq_len(nrow(df)), .before = 1L)
  df
}

#' Scan a peptide for terminal blocking caps (N- or C-terminal).
#'
#' For `terminus = "N"`, iterates over candidate start positions such that the
#' modeled length (start → C-terminus) lies in `[dist_min, dist_max]`, and
#' scores the first `cap_len` residues as the cap.
#'
#' For `terminus = "C"`, iterates over candidate end positions such that the
#' modeled length (N-terminus → end) lies in `[dist_min, dist_max]`, and scores
#' the last `cap_len` residues as the cap. Internally this reverses the
#' peptide, runs the N-terminal scan (scores are orientation-invariant), and
#' re-orients the output: `modeled_subpeptide` and `cap_seq` are returned in
#' natural N→C orientation, and `start_index_*` / `end_index_1based` are
#' reported in the original peptide's coordinates.
#'
#' The cap-scoring features (aromatic/aliphatic/polar counts, net charge,
#' longest hydrophobic run, proline presence) are all symmetric under
#' reversal, so scores match what a direct C-terminal scan would produce.
#'
#' @param peptide Character. Peptide sequence (single string).
#' @param terminus `"N"` or `"C"`. Which terminus to cap.
#' @param dist_min,dist_max Integer. Modeled-length bounds.
#' @return A tibble of ranked caps (best→worst), or an empty tibble if no
#'   valid cap. Includes a `terminus` column.
#' @export
scan_terminal_caps <- function(peptide,
                               terminus = c("N", "C"),
                               dist_min = DIST_MIN_DEFAULT,
                               dist_max = DIST_MAX_DEFAULT) {
  terminus <- match.arg(terminus)
  pep <- clean_peptide(peptide)
  n <- nchar(pep)

  scan_pep <- if (terminus == "C") str_reverse(pep) else pep
  df <- scan_n_caps_core(scan_pep, dist_min = dist_min, dist_max = dist_max)
  if (nrow(df) == 0) return(df)

  if (terminus == "C") {
    # Re-orient: cap/subpep are in reversed order; flip them back.
    df$modeled_subpeptide <- str_reverse(df$modeled_subpeptide)
    df$cap_seq            <- str_reverse(df$cap_seq)
    # For C-term scans the modeled_subpeptide always starts at position 1 of
    # the original peptide and ends at `end_index_1based`.
    df$end_index_1based   <- n - df$start_index_0based
    df$start_index_1based <- 1L
    df$start_index_0based <- 0L
    df$peptide            <- pep
  }
  df$terminus <- terminus
  df
}

#' Scan a C-hot peptide for N-terminal blocking caps.
#'
#' Thin back-compat wrapper around [scan_terminal_caps()] with
#' `terminus = "N"`.
#' @inheritParams scan_terminal_caps
#' @return A tibble of ranked N-terminal caps, or an empty tibble.
#' @export
scan_n_term_caps <- function(peptide,
                             dist_min = DIST_MIN_DEFAULT,
                             dist_max = DIST_MAX_DEFAULT) {
  scan_terminal_caps(peptide, terminus = "N",
                     dist_min = dist_min, dist_max = dist_max)
}

#' Compute features for the native terminal cap of a peptide.
#'
#' Uses the same cap_len rule as designed caps. `context` is the 6 AA
#' immediately inward of the cap (toward the peptide interior), reported in
#' natural N→C orientation: for an N-terminal cap that is the 6 AA after the
#' cap; for a C-terminal cap that is the 6 AA before the cap.
#'
#' @param peptide Character.
#' @param terminus `"N"` or `"C"`.
#' @return One-row tibble, or empty tibble if peptide is empty.
#' @export
compute_native_cap <- function(peptide, terminus = c("N", "C")) {
  terminus <- match.arg(terminus)
  peptide <- clean_peptide(peptide)
  n <- nchar(peptide)
  if (n == 0) return(tibble::tibble())

  cap_len <- choose_cap_len(n)

  if (terminus == "N") {
    cap_seq <- substr(peptide, 1L, cap_len)
    context <- if (cap_len + 1L <= n)
                 substr(peptide, cap_len + 1L, min(n, cap_len + 6L))
               else ""
  } else {
    cap_start <- n - cap_len + 1L
    cap_seq   <- substr(peptide, cap_start, n)
    ctx_end   <- cap_start - 1L
    ctx_start <- max(1L, ctx_end - 5L)
    context   <- if (ctx_end >= 1L) substr(peptide, ctx_start, ctx_end) else ""
  }

  feats <- compute_cap_features(cap_seq)
  label <- qualitative_from_score(feats$cap_score, feats$cap_seq)
  tier  <- tier_from_label(label)

  tibble::tibble(
    terminus      = terminus,
    cap_seq       = feats$cap_seq,
    context       = context,
    cap_len       = cap_len,
    cap_score     = feats$cap_score,
    n_aromatic    = feats$n_aromatic,
    n_W           = feats$n_W,
    n_aliphatic   = feats$n_aliphatic,
    n_hydrophobic = feats$n_hydrophobic,
    hydro_run     = feats$hydro_run,
    n_polar       = feats$n_polar,
    n_neg         = feats$n_neg,
    n_pos         = feats$n_pos,
    net_charge    = feats$net_charge,
    has_pro       = feats$has_pro,
    avgKD         = feats$avgKD,
    qualitative   = label,
    tier          = tier
  )
}

#' Compute features for the native N-terminal cap of a peptide.
#'
#' Thin back-compat wrapper around [compute_native_cap()] with
#' `terminus = "N"`.
#' @param peptide Character.
#' @return One-row tibble, or empty tibble if peptide is empty.
#' @export
compute_native_ncap <- function(peptide) {
  compute_native_cap(peptide, terminus = "N")
}


# ---- Window → mature-peptide pipeline -------------------------------------

#' Parse `nn_input$peps` strings into gene + coordinates.
#'
#' Matches the naming established by `get_pep_data()` in
#' `9.2_add_contact_data.R` (`paste0(gene, "_w", wN, "-", wC)`).
#'
#' @param peps Character vector, e.g. `"CAMP_w125-160"`.
#' @return A tibble with columns `peps, gene, wN, wC`.
#' @export
parse_peps_window <- function(peps) {
  m <- stringr::str_match(peps, "^(.+)_w(-?\\d+)-(-?\\d+)$")
  tibble::tibble(
    peps = peps,
    gene = m[, 2],
    wN   = as.integer(m[, 3]),
    wC   = as.integer(m[, 4])
  )
}

#' Pull the AA+index pair out of a per-window meta_data data.frame.
#'
#' Takes a 36-row data.frame (one element of `nn_input$N$all$meta_data` or
#' `nn_input$C$all$meta_data`) and returns parallel character/integer vectors
#' of length 36, with `NA` at padding rows.
#'
#' @param md A data.frame with `AA` and `index` columns.
#' @return A list with `aa` (character) and `index` (integer).
#' @export
meta_window_vectors <- function(md) {
  stopifnot(all(c("AA", "index") %in% names(md)))
  list(aa    = as.character(md$AA),
       index = as.integer(md$index))
}

#' Trim the fixed 7-residue padding/cleavage-context/DB block from the
#' inserting end of a 36-position window, carrying indices along.
#'
#' Mirrors the layout assigned by `make_known_idx()` in
#' `9.2_add_contact_data.R` (positions 1-7 for N-inserting, positions 30-36
#' for C-inserting; 5 cleavage_context + 2 DB residues). Operates on the
#' raw 36-slot AA / index vectors from `meta_data` so that padding is stripped
#' correctly even when the window extends past the protein terminus. Remaining
#' NA residues (from the non-inserting side extending past a protein end) are
#' then dropped, keeping AA and index aligned.
#'
#' @param aa Character vector of length-36 residue codes (NA = padding).
#' @param index Integer vector of length-36 uniprot positions (NA = padding).
#' @param target `"N"`, `"C"`, `"loop_N"`, or `"loop_C"`.
#' @return A list with `seq` (concatenated AA string, NAs removed) and
#'   `index` (matching integer vector).
#' @export
trim_inserting_end <- function(aa, index, target) {
  if (length(target) != 1L) {
    stop("`target` must be scalar; vectorise over windows with purrr/map.")
  }
  if (is.na(target)) return(list(seq = NA_character_, index = integer()))

  n <- length(aa)
  if (length(index) != n) stop("`aa` and `index` must be the same length.")
  if (n <= 7L) return(list(seq = "", index = integer()))

  keep <- if (target %in% c("N", "loop_N")) {
    seq.int(8L, n)
  } else if (target %in% c("C", "loop_C")) {
    seq.int(1L, n - 7L)
  } else {
    stop("Unrecognised target: ", target)
  }

  aa_k  <- aa[keep]
  idx_k <- index[keep]
  # Drop any residual padding rows (NA AA) that were on the non-inserting side.
  ok <- !is.na(aa_k) & nzchar(aa_k) & aa_k != "-"
  list(seq   = paste(aa_k[ok], collapse = ""),
       index = idx_k[ok])
}

#' Trim the non-inserting end of a peptide using cap scoring, carrying
#' indices along.
#'
#' Delegates to [scan_terminal_caps()] with `terminus` set to the
#' non-inserting end: `"C"` when `target` is N-inserting, `"N"` when `target`
#' is C-inserting. Candidates containing an internal dibasic site
#' (`KK|KR|RK|RR`) are dropped when `forbid_internal_dibasic = TRUE`. The
#' corresponding `index` vector is sliced with the same cut so that AA and
#' uniprot positions stay aligned.
#'
#' @param seq Character. A single peptide sequence (inserting-end already
#'   removed).
#' @param index Integer vector of uniprot positions, same length as `seq`.
#' @param target Character. `"N"`, `"C"`, `"loop_N"`, or `"loop_C"`.
#' @param dist_min,dist_max Integer. Modeled-length bounds for scanning.
#' @param forbid_internal_dibasic Logical. Drop candidates with internal
#'   dibasic sites.
#' @return One-row tibble with `mature_peptide, mature_index` (list-column
#'   of integer vectors), `cap_seq, cap_score, qualitative, modeled_length`
#'   (all NA / `integer(0)` if no valid candidate).
#' @export
trim_non_inserting_end <- function(seq, index, target,
                                   dist_min = DIST_MIN_DEFAULT,
                                   dist_max = DIST_MAX_DEFAULT,
                                   forbid_internal_dibasic = TRUE) {
  empty <- tibble::tibble(
    mature_peptide = NA_character_,
    mature_index   = list(integer()),
    cap_seq        = NA_character_,
    cap_score      = NA_real_,
    qualitative    = NA_character_,
    modeled_length = NA_integer_
  )
  if (is.na(seq) || is.na(target) || !nzchar(seq)) return(empty)
  if (length(index) != nchar(seq)) {
    stop("`index` length (", length(index), ") must match `seq` length (",
         nchar(seq), ").")
  }

  non_inserting <- if (target %in% c("N", "loop_N")) "C" else "N"

  scans <- scan_terminal_caps(seq, terminus = non_inserting,
                              dist_min = dist_min, dist_max = dist_max)
  if (nrow(scans) == 0) return(empty)

  if (forbid_internal_dibasic) {
    scans <- scans[!grepl(DIBASIC_REGEX, scans$modeled_subpeptide), , drop = FALSE]
    if (nrow(scans) == 0) return(empty)
  }

  top <- scans[1L, , drop = FALSE]

  # Slice the index vector to match the kept sub-peptide:
  #   N-term cap: keep positions start_index_1based .. peptide_length
  #   C-term cap: keep positions 1 .. end_index_1based
  if (non_inserting == "N") {
    idx_keep <- index[top$start_index_1based:top$peptide_length]
  } else {
    idx_keep <- index[1L:top$end_index_1based]
  }

  tibble::tibble(
    mature_peptide = top$modeled_subpeptide,
    mature_index   = list(idx_keep),
    cap_seq        = top$cap_seq,
    cap_score      = top$cap_score,
    qualitative    = top$qualitative,
    modeled_length = top$modeled_length
  )
}

#' Convert `nn_input` meta_data windows into mature peptides + uniprot indices.
#'
#' Pipeline: [meta_window_vectors()] → [trim_inserting_end()] →
#' [trim_non_inserting_end()], carrying the uniprot index vector through each
#' trim so AA positions stay aligned.
#'
#' Vectorised: `meta_data` is a list of per-window data.frames (e.g.
#' `nn_input$N$all$meta_data`) and `target` is a same-length character vector
#' (e.g. `nn_input$N$all$target`). A single data.frame + scalar target is
#' also accepted.
#'
#' @param meta_data A list of per-window data.frames with `AA` and `index`
#'   columns (36 rows each), or a single such data.frame.
#' @param target Character vector matching `length(meta_data)`.
#' @param peps Optional character vector of window names for labelling.
#' @param gene Optional character vector of gene names for labelling.
#' @param dist_min,dist_max Integer. Modeled-length bounds.
#' @param forbid_internal_dibasic Logical.
#' @return A tibble with one row per input containing `peps, gene, target,
#'   window_seq` (residues with NA padding dropped), `window_index`
#'   (list-column), `post_inserting_trim, post_inserting_index` (list-column),
#'   `mature_peptide, mature_index` (list-column of uniprot positions),
#'   `cap_seq, cap_score, qualitative, modeled_length`.
#' @examples
#' \dontrun{
#' nn <- readRDS("~/AF2_analysis/nn_input.rds")
#' dat <- nn$N$all
#' out <- window_to_mature_peptide(
#'   meta_data = dat$meta_data,
#'   target    = dat$target,
#'   peps      = dat$peps,
#'   gene      = dat$gene
#' )
#' out[, c("peps", "mature_peptide", "mature_index")]
#' }
#' @export
window_to_mature_peptide <- function(meta_data,
                                     target,
                                     peps = NULL,
                                     gene = NULL,
                                     dist_min = DIST_MIN_DEFAULT,
                                     dist_max = DIST_MAX_DEFAULT,
                                     forbid_internal_dibasic = TRUE) {
  if (inherits(meta_data, "data.frame")) meta_data <- list(meta_data)
  N <- length(meta_data)
  if (length(target) != N) {
    stop("`target` length (", length(target), ") must match ",
         "`meta_data` length (", N, ").")
  }
  if (!is.null(peps) && length(peps) != N) stop("`peps` length must match meta_data length.")
  if (!is.null(gene) && length(gene) != N) stop("`gene` length must match meta_data length.")

  pieces <- vector("list", N)
  for (i in seq_len(N)) {
    vecs <- meta_window_vectors(meta_data[[i]])
    trimmed <- trim_inserting_end(vecs$aa, vecs$index, target[i])

    # Window sequence/index with padding NAs removed (for reporting).
    ok_w <- !is.na(vecs$aa) & nzchar(vecs$aa) & vecs$aa != "-"
    window_seq <- paste(vecs$aa[ok_w], collapse = "")
    window_idx <- vecs$index[ok_w]

    mature <- trim_non_inserting_end(
      seq                     = trimmed$seq,
      index                   = trimmed$index,
      target                  = target[i],
      dist_min                = dist_min,
      dist_max                = dist_max,
      forbid_internal_dibasic = forbid_internal_dibasic
    )

    pieces[[i]] <- tibble::tibble(
      peps                  = if (is.null(peps)) NA_character_ else peps[i],
      gene                  = if (is.null(gene)) NA_character_ else gene[i],
      target                = target[i],
      window_seq            = window_seq,
      window_index          = list(window_idx),
      post_inserting_trim   = trimmed$seq,
      post_inserting_index  = list(trimmed$index),
      mature_peptide        = mature$mature_peptide,
      mature_index          = mature$mature_index,
      cap_seq               = mature$cap_seq,
      cap_score             = mature$cap_score,
      qualitative           = mature$qualitative,
      modeled_length        = mature$modeled_length
    )
  }
  dplyr::bind_rows(pieces)
}
