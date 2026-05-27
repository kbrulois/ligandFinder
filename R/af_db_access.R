## ---------------------------------------------------------------------------
## Read structures out of an AlphaFold-DB cache.
##
## The cache is the layout produced by inst/scripts/af_db_bulk_download.sh:
##   <cache_dir>/cif/AF-<UNI>-F1-model_v4.cif    (always present after bulk download)
##   <cache_dir>/pdb/AF-<UNI>-F1-model_v4.pdb    (optional; written on demand here)
##   <cache_dir>/pae/AF-<UNI>-F1-predicted_aligned_error_v4.json
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
#' @param af_version  AlphaFold DB version tag.  Default \code{"v4"}.
#' @param download_if_missing  If \code{TRUE} (default), fetch from EBI if
#'   neither CIF nor PDB is in the cache.
#'
#' @details
#' Lookup order:
#' \enumerate{
#'   \item \code{<cache>/pdb/AF-<UNI>-F1-model_v4.pdb} — return as is.
#'   \item \code{<cache>/cif/AF-<UNI>-F1-model_v4.cif} — convert to PDB,
#'         cache under \code{pdb/}, return.
#'   \item Download PDB from EBI to \code{<cache>/pdb/}, return.
#' }
#' Only \code{F1} is consulted; multi-fragment entries error with a clear
#' message.
#'
#' @return Path to a PDB file (invisible).
#' @export
af_structure_to_pdb <- function(uniprot, cache_dir,
                                af_version = "v4",
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
    pdb <- bio3d::read.cif(cif_path, verbose = FALSE)
    bio3d::write.pdb(pdb, file = pdb_path)
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
#' Prefers PDB when present (faster parser); otherwise reads CIF.  Falls back
#' to downloading PDB from EBI if neither is in the cache.
#'
#' @inheritParams af_structure_to_pdb
#' @return A bio3d \code{"pdb"} object.
#' @export
read_af_structure <- function(uniprot, cache_dir,
                              af_version = "v4",
                              download_if_missing = TRUE) {

  pdb_path <- file.path(cache_dir, "pdb",
                        sprintf("AF-%s-F1-model_%s.pdb", uniprot, af_version))
  cif_path <- file.path(cache_dir, "cif",
                        sprintf("AF-%s-F1-model_%s.cif", uniprot, af_version))

  if (file.exists(pdb_path)) {
    return(bio3d::read.pdb(pdb_path, verbose = FALSE))
  }
  if (file.exists(cif_path)) {
    return(bio3d::read.cif(cif_path, verbose = FALSE))
  }
  if (!download_if_missing) {
    stop("No cached structure for ", uniprot, " and download_if_missing=FALSE")
  }
  pdb_path <- af_structure_to_pdb(uniprot, cache_dir,
                                  af_version = af_version,
                                  download_if_missing = TRUE)
  bio3d::read.pdb(pdb_path, verbose = FALSE)
}
