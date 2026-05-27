## ---------------------------------------------------------------------------
## Cleave a peptide range out of an AlphaFold-DB structure.
##
## AF-DB residue numbering matches UniProt sequence position (1-based), so a
## (uniprot, start, end) triple maps directly onto atoms in the cached PDB.
##
## The cached PDB is expected at:
##   <cache_dir>/pdb/AF-<UNIPROT>-F1-model_<af_version>.pdb
## (the layout written by inst/scripts/af_cache_download.R).  If the file is
## absent and `download_if_missing=TRUE`, it is pulled from EBI on demand.
##
## Designed to produce ligand PDBs consumable by
## `position_ligand_initial_guess()`.
## ---------------------------------------------------------------------------

#' Cleave a peptide range from an AlphaFold-DB structure.
#'
#' @param uniprot      UniProt accession of the parent protein.
#' @param start,end    1-based residue range in the parent's UniProt sequence.
#' @param cache_dir    Cache root containing
#'   \code{pdb/AF-<UNIPROT>-F1-model_v4.pdb}.  Created if missing when
#'   downloading.
#' @param output_pdb   Path to write the cleaved PDB.  If \code{NULL}, a
#'   tempfile is used.
#' @param expected_seq Optional 1-letter amino acid string; if supplied, the
#'   extracted residue sequence must match exactly (a sanity check on
#'   UniProt/AF alignment).
#' @param renumber     If \code{TRUE} (default), output residues are
#'   renumbered \code{1..N}.  If \code{FALSE}, original UniProt numbering is
#'   preserved.
#' @param af_version   EBI AlphaFold DB version tag.  Default \code{"v4"}.
#' @param download_if_missing  If \code{TRUE} (default), fetch the parent PDB
#'   from EBI when not present in \code{cache_dir}.  Set \code{FALSE} to
#'   require pre-cached entries.
#'
#' @details
#' AF-DB ships single-fragment (\code{F1}) predictions for proteins up to
#' ~2700 aa.  Longer entries are split into overlapping \code{F1, F2, ...}
#' fragments — this function only consults \code{F1} and raises an error if
#' \code{[start, end]} falls outside its coverage.  Handle those rare cases
#' manually.
#'
#' The cleaved fragment is in its precursor-context conformation, not
#' necessarily the bioactive fold.  Disulfides and PTMs are absent.  Adequate
#' as an initial-guess seed; do not interpret as a refined structure.
#'
#' @return Path to the written PDB (invisible).
#' @export
cleave_peptide_from_af <- function(uniprot, start, end,
                                   cache_dir,
                                   output_pdb = NULL,
                                   expected_seq = NULL,
                                   renumber = TRUE,
                                   af_version = "v4",
                                   download_if_missing = TRUE) {

  if (!is.character(uniprot) || length(uniprot) != 1L || !nzchar(uniprot)) {
    stop("`uniprot` must be a single non-empty string")
  }
  start <- as.integer(start); end <- as.integer(end)
  if (is.na(start) || is.na(end) || start < 1L || end < start) {
    stop("Invalid residue range: start=", start, " end=", end)
  }

  pdb_dir <- file.path(cache_dir, "pdb")
  src_pdb <- file.path(
    pdb_dir,
    sprintf("AF-%s-F1-model_%s.pdb", uniprot, af_version)
  )

  if (!file.exists(src_pdb)) {
    if (!download_if_missing) {
      stop("Cache miss for ", uniprot, " and download_if_missing=FALSE: ", src_pdb)
    }
    dir.create(pdb_dir, recursive = TRUE, showWarnings = FALSE)
    url <- sprintf("https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_%s.pdb",
                   uniprot, af_version)
    message("Cache miss for ", uniprot, " — fetching from EBI")
    ok <- tryCatch({
      suppressWarnings(utils::download.file(url, src_pdb, mode = "wb",
                                            quiet = TRUE, method = "libcurl"))
      file.exists(src_pdb) && file.info(src_pdb)$size > 0
    }, error = function(e) {
      if (file.exists(src_pdb)) file.remove(src_pdb)
      FALSE
    })
    if (!ok) {
      stop("Failed to download AF prediction for ", uniprot,
           " — entry may be absent from EBI AF DB, or be a multi-fragment ",
           "structure (F1 not sufficient).")
    }
  }

  pdb <- bio3d::read.pdb(src_pdb, verbose = FALSE)

  min_resno <- min(pdb$atom$resno, na.rm = TRUE)
  max_resno <- max(pdb$atom$resno, na.rm = TRUE)
  if (start < min_resno || end > max_resno) {
    stop(sprintf(
      "Range %d-%d falls outside F1 coverage (%d-%d) for %s. Likely a multi-fragment AF entry — handle manually.",
      start, end, min_resno, max_resno, uniprot
    ))
  }

  sel <- bio3d::atom.select(pdb, resno = start:end)
  if (length(sel$atom) == 0L) {
    stop("No atoms in range ", start, "-", end, " for ", uniprot)
  }
  cleaved <- bio3d::trim.pdb(pdb, sel)

  if (!is.null(expected_seq)) {
    ca <- cleaved$atom[cleaved$atom$elety == "CA", ]
    ca <- ca[order(ca$resno), ]
    seq_extracted <- paste0(bio3d::aa321(ca$resname), collapse = "")
    if (!identical(seq_extracted, as.character(expected_seq))) {
      stop(sprintf(
        "Extracted sequence does not match expected:\n  got:      %s\n  expected: %s",
        seq_extracted, expected_seq
      ))
    }
  }

  if (renumber) {
    orig <- cleaved$atom$resno
    cleaved$atom$resno <- match(orig, sort(unique(orig)))
    cleaved$atom$eleno <- seq_len(nrow(cleaved$atom))
  }

  if (is.null(output_pdb)) output_pdb <- tempfile(fileext = ".pdb")
  bio3d::write.pdb(cleaved, file = output_pdb)
  invisible(output_pdb)
}
