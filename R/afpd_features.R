

## Input features for AlphaPulldown.
##
## Every job-file token -- "BKRB2", "CXL14,95-102" -- is resolved by
## AlphaPulldown to a pickle named after it in the feature directory.  A name
## with no pickle is not caught at submission time: the array task starts,
## burns its GPU allocation, and dies in parse_fold() with a FileNotFoundError
## listing the directory it searched.  The checks below are here so that
## failure happens in R, before sbatch.


afpd_features_dir <- function() {
  getOption(
    "ligandFinder.features_dir",
    "/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"
  )
}


## Names the feature directory can serve.  Accepts several directories, since
## run_structure_prediction.py takes --features_directory more than once and
## custom peptides are better kept out of the curated Homo_sapiens set.
afpd_feature_db <- function(features_dir = afpd_features_dir()) {

  files <- unlist(lapply(features_dir, list.files), use.names = FALSE)

  files <- files[grepl("\\.pkl(\\.xz)?$", files)]

  unique(stringr::str_remove(files, "\\.pkl(\\.xz)?$"))
}


## "CXL14,95-102" -> "CXL14".  A residue range is sliced out of the
## full-length pickle at prediction time, so it never needs features of its
## own -- only the protein it is a window of does.
afpd_feature_name <- function(x) {
  stringr::str_remove(stringr::str_trim(x), ",.*$")
}


## Every distinct protein a job file will ask for.
afpd_job_names <- function(lines, delim_proteins = ";") {

  lines <- lines[!grepl("^[[:space:]]*(#|$)", lines)]

  parts <- unlist(strsplit(lines, delim_proteins, fixed = TRUE), use.names = FALSE)

  unique(afpd_feature_name(parts))
}


## Which names the feature directory can serve, and which cannot.  Pass job
## file lines through afpd_job_names() first, or hand it names directly.
afpd_check_features <- function(x, features_dir = afpd_features_dir()) {

  db <- afpd_feature_db(features_dir)

  out <- tibble::tibble(input = x,
                        name = afpd_feature_name(x),
                        present = afpd_feature_name(x) %in% db)

  n_missing <- sum(!out[["present"]])

  if (n_missing > 0) {
    message(n_missing, " of ", nrow(out), " name(s) have no features in:\n  ",
            paste(features_dir, collapse = "\n  "), "\n  missing: ",
            paste(unique(out[["name"]][!out[["present"]]]), collapse = ", "))
  } else {
    message("all ", nrow(out), " name(s) resolve to features")
  }

  out
}


## Write a FASTA for sequences that need features built.
##
## `seqs` is a named character vector: names become FASTA headers, and
## AlphaPulldown names each pickle after the header -- so the header is
## exactly the token the job file will use.  Follow the directory's
## convention: UniProt entry name minus the _HUMAN suffix.
afpd_write_feature_fasta <- function(seqs,
                                     file,
                                     features_dir = afpd_features_dir(),
                                     skip_existing = TRUE,
                                     out_dir = NULL) {

  if (is.null(names(seqs)) || any(!nzchar(names(seqs)))) {
    stop("`seqs` must be a named character vector -- names become the feature names")
  }

  seqs <- stats::setNames(toupper(gsub("[[:space:]]", "", seqs)), names(seqs))

  ## A name carrying a job-file delimiter would split a fold in two, or read
  ## as a residue range, long after the features are built.
  bad_name <- grepl("[;,:+[:space:]]", names(seqs))
  if (any(bad_name)) {
    stop("these names contain a job-file delimiter (; , : + or space): ",
         paste(names(seqs)[bad_name], collapse = ", "))
  }

  if (any(duplicated(names(seqs)))) {
    stop("duplicate names: ", paste(unique(names(seqs)[duplicated(names(seqs))]), collapse = ", "))
  }

  bad_seq <- grepl("[^ACDEFGHIKLMNPQRSTVWY]", seqs)
  if (any(bad_seq)) {
    warning("non-standard residues in: ", paste(names(seqs)[bad_seq], collapse = ", "))
  }

  if (skip_existing) {
    have <- names(seqs) %in% afpd_feature_db(features_dir)
    if (any(have)) {
      message("skipping ", sum(have), " name(s) that already have features: ",
              paste(names(seqs)[have], collapse = ", "))
    }
    seqs <- seqs[!have]
  }

  if (length(seqs) == 0) {
    message("nothing to build")
    return(invisible(NULL))
  }

  writeLines(as.vector(rbind(paste0(">", names(seqs)), unname(seqs))), file)

  message("wrote ", length(seqs), " record(s) to ", file)
  message("build them with:")
  message("  sbatch --array=1-", length(seqs), " ",
          "$R_LIBS/ligandFinder/scripts/afpd_create_features.sh ",
          file, if (!is.null(out_dir)) paste0(" ", out_dir) else "")

  invisible(file)
}
