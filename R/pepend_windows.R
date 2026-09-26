## Peptide-END windows: a second, separate training set for the window models.
##
## The production set (9.2_add_contact_data.R) anchors every window on a
## dibasic site -- DB at positions 30-31 (C) / 6-7 (N) -- and keeps only the
## docked peptides' INSERTING terminus, and only where a dibasic site sits
## within 6 residues of it. This set instead anchors on the peptide's own
## terminal residue and takes BOTH ends of every docked peptide:
##
##   N window  first peptide residue at position 8   (7 residues upstream)
##   C window  last  peptide residue at position 28  (8 residues downstream)
##
## and drops the DB and gap classes: whatever precedes the peptide is
## NT_cleavage_context, whatever follows it CT_cleavage_context, dibasic or not.
## Candidate windows are anchored the same way, at the residue a cleavage there
## would leave as the peptide terminus.

#' Window geometry of the peptide-end set
#' @keywords internal
LF_PEPEND <- list(seq_len = 36L, anchor = c(N = 8L, C = 28L))

#' The class vocabulary of the peptide-end set (no DB, no gap)
#'
#' Order is the order the one-hot columns are built in; pass it to
#' [lf_dcnn_config()] as `class_names`.
#' @export
LF_PEPEND_CLASSES <- c("CT_cleavage_context", "NT_cleavage_context",
                       "pep_other", "pep_pocket", "padding", "none")

#' First precursor residue of a peptide-end window
#'
#' @param term `"N"` or `"C"` (vectorised).
#' @param anchor the peptide's first (N) or last (C) residue, in precursor
#'   coordinates.
#' @return the precursor coordinate at window position 1; the window runs
#'   through that + 35.
#' @export
lf_pepend_start <- function(term, anchor) {
  as.integer(anchor) - LF_PEPEND$anchor[as.character(term)] + 1L
}

#' Per-residue labels of one peptide-end window
#'
#' Term-agnostic: the anchoring decides where the peptide falls, so the same
#' rule labels N and C windows. Residues outside the mature precursor
#' `[n_prot, c_prot]` are padding (the signal peptide counts as outside, as in
#' 9.2); before the peptide is NT context, after it CT context, and the
#' peptide is `pep_other` except at `pocket` residues.
#'
#' @param w_start precursor coordinate of window position 1.
#' @param pep_start,pep_end the peptide, in precursor coordinates.
#' @param n_prot,c_prot first and last residue of the mature precursor.
#' @param pocket precursor coordinates the docked model puts in the pocket.
#' @return character vector, length 36.
#' @export
lf_pepend_labels <- function(w_start, pep_start, pep_end, n_prot, c_prot,
                             pocket = integer(0)) {
  r <- w_start + seq_len(LF_PEPEND$seq_len) - 1L
  lab <- ifelse(r < pep_start, "NT_cleavage_context",
         ifelse(r > pep_end,   "CT_cleavage_context",
         ifelse(r %in% pocket, "pep_pocket", "pep_other")))
  lab[r < n_prot | r > c_prot] <- "padding"
  lab
}

#' Candidate peptide-end anchors in one precursor
#'
#' Where a cleavage would leave a peptide terminus, anchored as the knowns are:
#' \describe{
#'   \item{`db`}{a dibasic pair KK/KR/RK/RR at `d, d+1` (overlapping pairs
#'     all count). Prohormone convertases cut after the pair and
#'     carboxypeptidase E trims the basics, so the upstream peptide ends at
#'     `d - 1` (C) and the downstream one starts at `d + 2` (N).}
#'   \item{`db_amid`}{`d - 2` as well when `d - 1` is G: the glycine is the
#'     amide donor and leaves, so an amidated peptide ends one residue earlier.}
#'   \item{`pep_end`}{the mature N-terminus `n_prot` (N) and the last residue
#'     (C): a peptide can end at the precursor's own terminus.}
#' }
#' @param seq precursor sequence.
#' @param n_prot first residue of the mature precursor (after the signal peptide).
#' @return tibble `term`, `anchor`, `win_type`, one row per distinct
#'   (term, anchor); a `db` reading wins over `db_amid` and `pep_end`.
#' @export
lf_pepend_candidates <- function(seq, n_prot) {
  c_prot <- nchar(seq)
  d <- stringr::str_locate_all(seq, "(?=(KK|KR|RK|RR))")[[1]][, "start"]
  aa <- strsplit(seq, "")[[1]]
  amid <- d[d > 2L & aa[pmax(d - 1L, 1L)] == "G"]
  out <- dplyr::bind_rows(
    tibble::tibble(term = "C", anchor = d - 1L,    win_type = "db"),
    tibble::tibble(term = "N", anchor = d + 2L,    win_type = "db"),
    tibble::tibble(term = "C", anchor = amid - 2L, win_type = "db_amid"),
    tibble::tibble(term = c("N", "C"), anchor = c(n_prot, c_prot), win_type = "pep_end"))
  out %>%
    dplyr::filter(anchor >= n_prot, anchor <= c_prot) %>%
    dplyr::distinct(term, anchor, .keep_all = TRUE)
}

#' What sits just outside a peptide terminus
#'
#' A coarse label for the plots and diagnostics: the two residues after a C
#' terminus (three when the first is an amidation G), the two before an N
#' terminus.
#' @return one of `"dibasic"`, `"G + dibasic"`, `"monobasic"`, `"terminus"`,
#'   `"other"` per row.
#' @export
lf_pepend_motif <- function(seq, term, anchor, n_prot) {
  vapply(seq_along(seq), function(i) {
    s <- seq[[i]]; a <- anchor[[i]]
    if (term[[i]] == "C") {
      if (a >= nchar(s)) return("terminus")
      nxt <- substring(s, a + 1L, a + 3L)
      if (grepl("^(KK|KR|RK|RR)", nxt)) return("dibasic")
      if (grepl("^G(KK|KR|RK|RR)", nxt)) return("G + dibasic")
      if (grepl("^G?[KR]", nxt)) return("monobasic")
    } else {
      if (a <= n_prot[[i]]) return("terminus")
      prv <- substring(s, max(a - 2L, 1L), a - 1L)
      if (grepl("(KK|KR|RK|RR)$", prv)) return("dibasic")
      if (grepl("[KR]$", prv)) return("monobasic")
    }
    "other"
  }, character(1))
}

#' Slice one window out of a precursor's residue-feature matrix
#'
#' @param feat `(nchar(seq), n_channels)` matrix from 9.2's channel recipe,
#'   with a `padding` column.
#' @param w_start precursor coordinate of window position 1.
#' @param n_prot,c_prot residues outside this range are padding: every channel
#'   0 and `padding` = 1, as 9.2 pads.
#' @return `(36, n_channels)` tibble, the `data` element of the contract.
#' @export
lf_pepend_slice <- function(feat, w_start, n_prot, c_prot) {
  r  <- w_start + seq_len(LF_PEPEND$seq_len) - 1L
  ok <- r >= n_prot & r <= c_prot
  x  <- matrix(0, LF_PEPEND$seq_len, ncol(feat), dimnames = list(NULL, colnames(feat)))
  x[ok, ] <- feat[r[ok], , drop = FALSE]
  x[!ok, "padding"] <- 1
  tibble::as_tibble(x)
}
