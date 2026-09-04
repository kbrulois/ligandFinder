## ---------------------------------------------------------------------------
## Read structures out of an AlphaFold-DB cache.
##
## The cache is the layout produced by inst/scripts/af_db_bulk_download.sh:
##   <cache_dir>/cif/AF-<UNI>-F1-model_v6.cif    (always present after bulk download)
##   <cache_dir>/pdb/AF-<UNI>-F1-model_v6.pdb    (optional; written on demand here)
##   <cache_dir>/pae/AF-<UNI>-F1-predicted_aligned_error_v6.json
##
## EBI's per-proteome tarballs ship CIF only.  Some downstream tools want PDB
## (e.g. position_ligand_initial_guess()).  af_structure_to_pdb() converts CIF
## to PDB on first request and caches the result under pdb/ — subsequent calls
## are instant.  If neither file is present and download_if_missing=TRUE, the
## PDB is pulled from EBI's per-entry endpoint.
## ---------------------------------------------------------------------------

#' Resolve a UniProt accession to a PDB file path, converting from CIF if needed.
#'
#' @param uniprot   UniProt accession.
#' @param cache_dir Cache root (the dir passed to af_db_bulk_download.sh).
#' @param af_version  AlphaFold DB version tag.  Default \code{"v6"}.
#' @param download_if_missing  If \code{TRUE} (default), fetch from EBI if
#'   neither CIF nor PDB is in the cache.
#'
#' @details
#' Lookup order:
#' \enumerate{
#'   \item \code{<cache>/pdb/AF-<UNI>-F1-model_v6.pdb} — return as is.
#'   \item \code{<cache>/cif/AF-<UNI>-F1-model_v6.cif} — convert to PDB,
#'         cache under \code{pdb/}, return.
#'   \item Download PDB from EBI to \code{<cache>/pdb/}, return.
#' }
#' Only \code{F1} is consulted; multi-fragment entries error with a clear
#' message.
#'
#' @return Path to a PDB file (invisible).
#' @export
af_structure_to_pdb <- function(uniprot, cache_dir,
                                af_version = "v6",
                                download_if_missing = TRUE) {

  if (!is.character(uniprot) || length(uniprot) != 1L || !nzchar(uniprot)) {
    stop("`uniprot` must be a single non-empty string")
  }

  pdb_dir <- file.path(cache_dir, "pdb")
  cif_dir <- file.path(cache_dir, "cif")
  pdb_path <- file.path(pdb_dir, sprintf("AF-%s-F1-model_%s.pdb", uniprot, af_version))
  cif_path <- file.path(cif_dir, sprintf("AF-%s-F1-model_%s.cif", uniprot, af_version))

  if (file.exists(pdb_path)) return(invisible(pdb_path))

  dir.create(pdb_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(cif_path)) {
    cif_to_pdb(cif_path, pdb_path)
    return(invisible(pdb_path))
  }

  if (!download_if_missing) {
    stop("No cached structure for ", uniprot, " and download_if_missing=FALSE")
  }

  url <- sprintf("https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_%s.pdb",
                 uniprot, af_version)
  message("Cache miss for ", uniprot, " - fetching PDB from EBI")
  ok <- tryCatch({
    suppressWarnings(utils::download.file(url, pdb_path, mode = "wb",
                                          quiet = TRUE, method = "libcurl"))
    file.exists(pdb_path) && file.info(pdb_path)$size > 0
  }, error = function(e) {
    if (file.exists(pdb_path)) file.remove(pdb_path)
    FALSE
  })
  if (!ok) {
    stop("Failed to fetch AF prediction for ", uniprot,
         " - entry may be absent from EBI AF DB or be a multi-fragment structure.")
  }
  invisible(pdb_path)
}


#' Read an AlphaFold-DB structure into a bio3d object.
#'
#' Prefers PDB when present (faster parser); otherwise converts the cached CIF
#' via \code{cif_to_pdb()} and reads that.  Falls back to downloading PDB from
#' EBI if neither is in the cache.
#'
#' @inheritParams af_structure_to_pdb
#' @return A bio3d \code{"pdb"} object.
#' @export
read_af_structure <- function(uniprot, cache_dir,
                              af_version = "v6",
                              download_if_missing = TRUE) {

  pdb_path <- file.path(cache_dir, "pdb",
                        sprintf("AF-%s-F1-model_%s.pdb", uniprot, af_version))
  cif_path <- file.path(cache_dir, "cif",
                        sprintf("AF-%s-F1-model_%s.cif", uniprot, af_version))

  if (file.exists(pdb_path)) {
    return(bio3d::read.pdb(pdb_path, verbose = FALSE))
  }
  if (file.exists(cif_path)) {
    # convert (and cache) rather than read.cif directly -- see cif_to_pdb()
    return(bio3d::read.pdb(cif_to_pdb(cif_path, pdb_path), verbose = FALSE))
  }
  if (!download_if_missing) {
    stop("No cached structure for ", uniprot, " and download_if_missing=FALSE")
  }
  pdb_path <- af_structure_to_pdb(uniprot, cache_dir,
                                  af_version = af_version,
                                  download_if_missing = TRUE)
  bio3d::read.pdb(pdb_path, verbose = FALSE)
}


#' Convert an AlphaFold mmCIF to PDB.
#'
#' Parses the \code{_atom_site} loop directly instead of going through
#' \code{bio3d::read.cif()}.  That reader is beta and populates \code{elety}
#' from the element symbol rather than the atom name, so every CA-based
#' selection downstream silently matches nothing -- \code{cleave_peptide_from_af()}
#' reports an empty extracted sequence, \code{position_ligand_initial_guess()}
#' has no backbone to work with, and the pLDDT B-factor column is unreliable.
#' Writing that object out to PDB bakes the damage into the cache, so the
#' conversion has to avoid read.cif entirely.
#'
#' Atom names come from \code{label_atom_id}, residue numbering from
#' \code{auth_seq_id} (which is UniProt position for AF-DB entries), and
#' pLDDT from \code{B_iso_or_equiv}.
#'
#' Assumes whitespace-delimited values with no quoted fields containing spaces,
#' which holds for AF-DB mmCIF but not for mmCIF in general.
#'
#' @param cif_path Path to the mmCIF file.
#' @param pdb_path Path to write the PDB to.
#' @return \code{pdb_path}, invisibly.
#' @export
cif_to_pdb <- function(cif_path, pdb_path) {

  ln  <- readLines(cif_path, warn = FALSE)
  hdr <- grep("^_atom_site\\.", ln)
  if (!length(hdr)) stop("no _atom_site loop in ", cif_path)
  field <- sub("^_atom_site\\.", "", ln[hdr])

  body <- ln[(max(hdr) + 1L):length(ln)]
  body <- body[grepl("^(ATOM|HETATM)", body)]
  if (!length(body)) stop("no atom records in ", cif_path)

  f <- do.call(rbind, strsplit(trimws(body), "\\s+"))
  if (ncol(f) != length(field))
    stop("_atom_site has ", length(field), " fields but rows have ", ncol(f),
         " values: ", cif_path)
  colnames(f) <- field

  g <- function(nm, alt = NULL) {
    if (nm %in% field) return(f[, nm])
    if (!is.null(alt) && alt %in% field) return(f[, alt])
    stop("missing CIF field: ", nm)
  }

  nm   <- g("label_atom_id", "auth_atom_id")
  # PDB convention: atom names shorter than 4 characters start in column 14
  nm_f <- ifelse(nchar(nm) >= 4, nm, sprintf(" %-3s", nm))

  rec <- sprintf(
    "%-6s%5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s",
    g("group_PDB"), g("id"), nm_f, " ",
    g("label_comp_id", "auth_comp_id"),
    g("auth_asym_id", "label_asym_id"),
    g("auth_seq_id", "label_seq_id"), " ",
    as.numeric(g("Cartn_x")), as.numeric(g("Cartn_y")), as.numeric(g("Cartn_z")),
    as.numeric(g("occupancy")), as.numeric(g("B_iso_or_equiv")),
    g("type_symbol"))

  dir.create(dirname(pdb_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(rec, "END"), pdb_path)
  invisible(pdb_path)
}
