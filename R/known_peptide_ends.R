## Which windows sit at a KNOWN peptide boundary, and what kind of peptide.
##
## Lifted from the viewer (10_3_peptide_umap.R) so the benchmark CSVs and the
## viewer agree on the labels. The training positives are deliberately narrow --
## they need a cognate receptor-ligand model and a terminal (not middle)
## insertion -- so most real peptide ends are absent from `known`; this reads
## the boundaries straight off ligand_list.rds instead.

#' Annotate windows at known peptide boundaries, with the peptide's class
#'
#' A window is a known peptide end when some reference peptide of the same gene
#' has its boundary 1..`known_end_gap` residues from the window's anchor -- the
#' C branch measures to the peptide's END, the N branch to its START, mirroring
#' how 9.2 computed db_spacer in each direction. Two candidate anchors are
#' derived, one from each window edge, because either edge can have been
#' clamped at a protein terminus.
#'
#' Matched windows get a `stratum`: `"GPCR peptide"` when the reference ligand
#' has a receptor annotated (GtoPdb / GPCRdb), `"peptide"` otherwise. GPCR
#' peptides are further split by how THIS terminus sits in the docked
#' receptor-ligand models (`docked_ref`, 9.2's own rule on `lig1_end`):
#' `"end insertion"` (pocket-inserting residue < 3 from the terminus),
#' `"loop insertion"` (3-15), `"non-inserting end"` (the other terminus or the
#' middle inserts), `"no model"`. The class is per WINDOW: both ends of a
#' peptide are known ends but only one inserts.
#'
#' @param windows data frame with `peps` (`<gene>_w<N>-<C>`), `gene`, `target`
#'   and `known`.
#' @param known_end_ref path to ligand_list.rds.
#' @param docked_ref path to knowns.rds (the docked models); NULL skips the
#'   insertion split.
#' @param id_map the id_mapping data frame (`Entry`, `Entry Name`,
#'   `Gene Names (primary)`).
#' @param known_end_len max reference peptide length to consider.
#' @param known_end_gap anchor-to-boundary distance to accept.
#' @param docked_tol residues by which a docked peptide's end may miss the
#'   reference boundary and still count as the same peptide.
#' @return `windows` with `known_end`, `known_end_name`, `stratum`, `insertion`
#'   and one combined `peptide_class` label: `"unknown"`,
#'   `"peptide (non-GPCR ligand)"`, or `"GPCR peptide: <insertion>"`.
#' @export
lf_known_peptide_ends <- function(windows, known_end_ref, docked_ref = NULL, id_map,
                                  known_end_len = 100L, known_end_gap = 5L, docked_tol = 2L) {
  a2s <- stats::setNames(id_map[["Gene Names (primary)"]], id_map[["Entry"]])
  e2s <- stats::setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])

  ref <- readRDS(path.expand(known_end_ref)) %>%
    dplyr::mutate(gene  = dplyr::coalesce(a2s[accession], e2s[paste0(uniprot_name, "_HUMAN")]),
                  start = as.integer(start), end = as.integer(end),
                  name  = dplyr::coalesce(as.character(final_name), "unnamed")) %>%
    dplyr::filter(!is.na(gene), !is.na(start), !is.na(end),
                  end - start + 1L <= known_end_len)

  co  <- stringr::str_match(windows$peps, "_w(\\d+)-(\\d+)$")
  wN  <- as.integer(co[, 2]); wC <- as.integer(co[, 3])
  isC <- as.character(windows$target) %in% c("C", "loop_C")
  side <- ifelse(isC, "C", "N")
  ## TWO candidate anchors, one recovered from each edge (either may be clamped
  ## at a protein terminus); they agree whenever the window is a full 36 residues
  anchor_a <- ifelse(isC, wC -  5L, wN +  5L)
  anchor_b <- ifelse(isC, wN + 30L, wC - 30L)

  bnd <- dplyr::bind_rows(
    ref %>% dplyr::transmute(gene, name, side = "C", pos = end),
    ref %>% dplyr::transmute(gene, name, side = "N", pos = start))
  bl  <- split(bnd[, c("pos", "name")], paste0(bnd$gene, "|", bnd$side))
  key <- paste0(windows$gene, "|", side)

  hit <- vapply(seq_along(key), function(i) {
    b <- bl[[key[i]]]
    if (is.null(b)) return(NA_character_)
    best <- NA_character_; bestd <- Inf
    for (a in c(anchor_a[i], anchor_b[i])) {
      if (is.na(a)) next
      d  <- if (isC[i]) a - b$pos else b$pos - a
      ok <- which(d >= 1L & d <= known_end_gap)
      if (length(ok) && min(d[ok]) < bestd) {
        bestd <- min(d[ok]); best <- b$name[ok[which.min(d[ok])]]
      }
    }
    best
  }, character(1))
  windows$known_end      <- !is.na(hit)
  windows$known_end_name <- hit

  ## stratum is a property of the LIGAND: joined on the matched name
  rcpt <- readRDS(path.expand(known_end_ref)) %>%
    dplyr::transmute(name = dplyr::coalesce(as.character(final_name), "unnamed"),
                     has_receptor = !is.na(receptor)) %>%
    dplyr::group_by(name) %>% dplyr::summarise(has_receptor = any(has_receptor), .groups = "drop")
  has_rcpt <- rcpt$has_receptor[match(windows$known_end_name, rcpt$name)]
  windows$stratum <- factor(
    dplyr::case_when(!windows$known_end    ~ "unknown",
                     has_rcpt %in% TRUE    ~ "GPCR peptide",
                     TRUE                  ~ "peptide"),
    levels = c("unknown", "peptide", "GPCR peptide"))

  ## insertion type of the GPCR-peptide windows, from the docked models: same
  ## funnel 9.2 applied before it labelled targets (relevant site, rank-1 model,
  ## per docked peptide the model whose inserting residue sits closest to a
  ## terminus, ties by iptm); several docked forms can share a boundary, so the
  ## best class wins (end > loop > middle > other end)
  windows$insertion <- NA_character_
  gi <- which(windows$stratum == "GPCR peptide")
  if (!is.null(docked_ref) && file.exists(path.expand(docked_ref)) && length(gi)) {
    dk <- readRDS(path.expand(docked_ref)) %>%
      dplyr::filter(location == "relevant", rank == 1) %>%
      dplyr::transmute(gene = e2s[p2_name],
                       ps   = as.integer(stringr::str_extract(p2_range, "^\\d+")),
                       pe   = as.integer(stringr::str_extract(p2_range, "\\d+$")),
                       ind  = as.integer(stringr::str_extract(lig1_end, "\\d+")),
                       tt   = stringr::str_remove(lig1_end, "\\d+"),
                       iptm) %>%
      dplyr::filter(!is.na(gene), !is.na(ind)) %>%
      dplyr::group_by(gene, ps, pe) %>%
      dplyr::arrange(ind, dplyr::desc(iptm), .by_group = TRUE) %>%
      dplyr::slice(1) %>% dplyr::ungroup()
    dkl <- split(dk, dk$gene)
    lb  <- dplyr::bind_rows(ref %>% dplyr::transmute(name, gene, side = "C", pos = end),
                            ref %>% dplyr::transmute(name, gene, side = "N", pos = start)) %>%
      dplyr::distinct(name, gene, side, .keep_all = TRUE)
    pos <- lb$pos[match(paste(windows$known_end_name, windows$gene, side),
                        paste(lb$name, lb$gene, lb$side))]
    cls <- c(end = "end insertion", loop = "loop insertion",
             middle = "non-inserting end", other = "non-inserting end")
    gi2 <- gi[!is.na(pos[gi])]
    windows$insertion[gi2] <- vapply(gi2, function(i) {
      d <- dkl[[windows$gene[i]]]
      if (is.null(d)) return("no model")
      d <- d[abs((if (side[i] == "C") d$pe else d$ps) - pos[i]) <= docked_tol, , drop = FALSE]
      if (!nrow(d)) return("no model")
      k <- ifelse(d$tt != side[i], "other",
                  ifelse(d$ind < 3, "end", ifelse(d$ind < 16, "loop", "middle")))
      unname(cls[[names(cls)[min(match(k, names(cls)))]]])
    }, character(1))
    windows$insertion[gi][is.na(windows$insertion[gi])] <- "no model"
  }

  windows$peptide_class <- factor(
    dplyr::case_when(windows$stratum == "unknown"      ~ "unknown",
                     windows$stratum == "peptide"      ~ "peptide (non-GPCR ligand)",
                     !is.na(windows$insertion)         ~ paste0("GPCR peptide: ", windows$insertion),
                     TRUE                              ~ "GPCR peptide"),
    levels = c("unknown", "peptide (non-GPCR ligand)", "GPCR peptide",
               "GPCR peptide: end insertion", "GPCR peptide: loop insertion",
               "GPCR peptide: non-inserting end", "GPCR peptide: no model"))

  message("peptide_class: ",
          paste(sprintf("%s=%d", names(table(windows$peptide_class)), table(windows$peptide_class)),
                collapse = "  "))
  rec <- if (any(windows$known == 1)) mean(windows$known_end[windows$known == 1]) else NA_real_
  message(sprintf("  known-end rule recovers %.0f%% of the %d training knowns", 100 * rec, sum(windows$known == 1)))
  if (!is.na(rec) && rec < 0.8)
    warning("known-end rule recovers <80% of the training knowns -- check the anchor offsets", immediate. = TRUE)
  windows
}
