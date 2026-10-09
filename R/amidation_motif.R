## The G|K/R-K/R amidation motif at a db window's anchor.
##
## Lifted from 10_1dcnn_new6.R so the benchmark CSVs and the pipeline agree.

#' Flag db windows whose dibasic anchor carries the amidation motif
#'
#' An amidated peptide is cut at a dibasic site with a glycine immediately 5'
#' of it: ...X-G | K/R-K/R. The G is the amide donor, so a db window carrying
#' one is a candidate amidation site.
#'
#' The anchor sits at a FIXED position in every db window, which is what makes
#' this a lookup rather than a search: 9.2 sets window_origin = db_ind and then
#' wN = db_ind + win_size[[t]]$start, so db_ind always lands at local position
#' `1 - start` (31 for C windows, 6 for N). Clamped windows are padded back out
#' to seq_len at the front, so the offset holds there too. Two asymmetries:
#' db_ind is the SECOND basic residue for C-target windows (the C branch of 9.2
#' adds a full lookahead offset) but the FIRST for N-target ones; and BOTH
#' termini are eligible -- the motif is a property of the dibasic SITE, not of
#' the peptide you approach it from, so an N-anchored window reads its G at
#' local position 5, a C-anchored one at 29.
#'
#' @param windows data frame with `win_type`, `target` and the per-window
#'   `meta_data` tibbles (with an `AA` column at local positions).
#' @param win_start the `start` of `win_size` per terminus as 9.2 built the
#'   windows; the 9.2 defaults.
#' @return `windows` with a logical `amidation` column.
#' @export
lf_amidation_motif <- function(windows, win_start = c(N = -5L, C = -30L)) {
  amid_targets <- c("N", "loop_N", "C", "loop_C")

  ## local positions of the glycine and the two basic residues, per target
  motif_pos <- function(tg) {
    if (tg %in% c("C", "loop_C")) {
      a <- 1L - win_start[["C"]]                  # db_ind = 2nd basic
      c(g = a - 2L, b1 = a - 1L, b2 = a)
    } else {
      a <- 1L - win_start[["N"]]                  # db_ind = 1st basic
      c(g = a - 1L, b1 = a, b2 = a + 1L)
    }
  }
  aa_at <- function(md, i) {
    aa <- as.character(md[["AA"]])
    if (i < 1L || i > length(aa)) NA_character_ else aa[[i]]
  }

  tg <- as.character(windows$target)
  wt <- as.character(windows$win_type)

  ## Sanity check FIRST: if the anchor offset were wrong, every window would
  ## quietly come back FALSE and look like "no amidation motifs found". Confirm
  ## the two anchor positions really are basic residues before trusting the G.
  db_i <- which(wt == "db")
  chk <- vapply(db_i, function(i) {
    p <- motif_pos(tg[[i]]); md <- windows$meta_data[[i]]
    isTRUE(aa_at(md, p[["b1"]]) %in% c("K", "R") && aa_at(md, p[["b2"]]) %in% c("K", "R"))
  }, logical(1))
  message(sprintf("amidation: dibasic anchor confirmed at the expected offset in %d/%d db windows (%.1f%%)",
                  sum(chk), length(chk), 100 * mean(chk)))
  if (length(chk) && mean(chk) < 0.9)
    warning("amidation: the dibasic anchor is often NOT at the expected local position -- ",
            "check win_start against 9.2 before using the `amidation` column", immediate. = TRUE)

  windows$amidation <- vapply(seq_len(nrow(windows)), function(i) {
    if (!identical(wt[[i]], "db") || !tg[[i]] %in% amid_targets) return(FALSE)
    p <- motif_pos(tg[[i]]); md <- windows$meta_data[[i]]
    isTRUE(identical(aa_at(md, p[["g"]]), "G") &&
           aa_at(md, p[["b1"]]) %in% c("K", "R") &&
           aa_at(md, p[["b2"]]) %in% c("K", "R"))
  }, logical(1))
  message(sprintf("amidation: %d windows carry the G|dibasic motif (%d of them known peptides)",
                  sum(windows$amidation), sum(windows$amidation & windows$known == 1)))
  windows
}
