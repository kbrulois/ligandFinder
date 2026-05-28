#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Generic AlphaFold-DB cache downloader.
##
## Given an RDS, TSV, CSV, or newline-delimited text file listing UniProt
## accessions, downloads per receptor:
##   <out>/pdb/AF-<UNI>-F1-model_v6.pdb           (pLDDT in B-factor)
##   <out>/pae/AF-<UNI>-F1-predicted_aligned_error_v6.json
## from the EBI AlphaFold Protein Structure Database (https://alphafold.ebi.ac.uk).
##
## Writes <out>/index.tsv with status + paths.  Idempotent — re-runs skip
## files already present unless --refresh is passed.
##
## Examples:
##   # GPCR cache  (column "uniprot_id" in gpcr_list.rds)
##   Rscript inst/scripts/af_cache_download.R \
##     --list inst/extdata/gpcr_list.rds \
##     --col  uniprot_id \
##     --out  /oak/stanford/groups/ebutcher/deorphan-AI-ze/gpcr_af_db
##
##   # Secretome cache  (column "accession")
##   Rscript inst/scripts/af_cache_download.R \
##     --list /oak/stanford/groups/ebutcher/deorphan-AI-ze/secretome_uniprot.rds \
##     --col  accession \
##     --out  /oak/stanford/groups/ebutcher/deorphan-AI-ze/secretome_af_db
##
##   # From a plain text list of accessions, one per line
##   Rscript inst/scripts/af_cache_download.R \
##     --list ids.txt \
##     --out  /oak/.../my_af_db
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "Input list: .rds (data frame), .tsv/.csv (with column), or .txt (one ID per line). REQUIRED."),
  make_option(c("-c", "--col"), type = "character", default = "uniprot_id",
              help = "Column name holding UniProt accessions (ignored for plain .txt). [default: %default]"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output cache directory. REQUIRED."),
  make_option("--refresh", action = "store_true", default = FALSE,
              help = "Re-download even if files exist locally."),
  make_option("--afver", type = "character", default = "v6",
              help = "EBI AlphaFold DB version tag [default: %default]"),
  make_option("--sleep", type = "double", default = 0.1,
              help = "Seconds to sleep between requests, to be polite to EBI [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$list) || !nzchar(opt$list)) stop("--list is required")
if (is.null(opt$out)  || !nzchar(opt$out))  stop("--out is required")
if (!file.exists(opt$list)) stop("List file not found: ", opt$list)

ext <- tolower(tools::file_ext(opt$list))
ids <- if (ext == "rds") {
  d <- readRDS(opt$list)
  if (!opt$col %in% names(d)) stop("Column '", opt$col, "' not found in ", opt$list,
                                   "\nAvailable: ", paste(names(d), collapse = ", "))
  d[[opt$col]]
} else if (ext %in% c("tsv", "csv")) {
  sep <- if (ext == "tsv") "\t" else ","
  d <- utils::read.delim(opt$list, sep = sep, stringsAsFactors = FALSE)
  if (!opt$col %in% names(d)) stop("Column '", opt$col, "' not found in ", opt$list,
                                   "\nAvailable: ", paste(names(d), collapse = ", "))
  d[[opt$col]]
} else if (ext %in% c("txt", "")) {
  readLines(opt$list)
} else {
  stop("Unrecognised list extension: ", ext, " (expected rds/tsv/csv/txt)")
}

ids <- unique(stats::na.omit(trimws(ids)))
ids <- ids[nzchar(ids)]

pdb_dir <- file.path(opt$out, "pdb")
pae_dir <- file.path(opt$out, "pae")
dir.create(pdb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pae_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://alphafold.ebi.ac.uk/files"
pdb_name <- function(u) sprintf("AF-%s-F1-model_%s.pdb", u, opt$afver)
pae_name <- function(u) sprintf("AF-%s-F1-predicted_aligned_error_%s.json", u, opt$afver)

message(sprintf("Input: %s  (column '%s')", opt$list, opt$col))
message(sprintf("IDs to fetch: %d unique", length(ids)))
message("Cache dir: ", opt$out)
message("Mode: ", if (opt$refresh) "refresh (re-download all)" else "incremental (skip existing)")

download_one <- function(url, dest) {
  tryCatch({
    suppressWarnings(
      utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
    )
    file.exists(dest) && file.info(dest)$size > 0
  }, error = function(e) {
    if (file.exists(dest)) file.remove(dest)
    FALSE
  })
}

results <- vector("list", length(ids))
t0 <- Sys.time()
for (i in seq_along(ids)) {
  u <- ids[[i]]
  pdb_path <- file.path(pdb_dir, pdb_name(u))
  pae_path <- file.path(pae_dir, pae_name(u))

  pdb_ok <- if (!opt$refresh && file.exists(pdb_path) && file.info(pdb_path)$size > 0) TRUE
            else download_one(file.path(base_url, pdb_name(u)), pdb_path)
  pae_ok <- if (!opt$refresh && file.exists(pae_path) && file.info(pae_path)$size > 0) TRUE
            else download_one(file.path(base_url, pae_name(u)), pae_path)

  results[[i]] <- data.frame(
    uniprot       = u,
    pdb_path      = if (pdb_ok) pdb_path else NA_character_,
    pae_path      = if (pae_ok) pae_path else NA_character_,
    pdb_status    = if (pdb_ok) "ok" else "missing",
    pae_status    = if (pae_ok) "ok" else "missing",
    downloaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE
  )

  if (i %% 50 == 0 || i == length(ids)) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf("  %5d / %5d  (%.1fs)  last: %s pdb=%s pae=%s",
                    i, length(ids), elapsed, u,
                    if (pdb_ok) "ok" else "MISS",
                    if (pae_ok) "ok" else "MISS"))
  }
  if (opt$sleep > 0) Sys.sleep(opt$sleep)
}

idx <- do.call(rbind, results)
n_pdb_ok <- sum(idx$pdb_status == "ok")
n_pae_ok <- sum(idx$pae_status == "ok")
message(sprintf("\nDone. PDB: %d / %d  PAE: %d / %d", n_pdb_ok, nrow(idx), n_pae_ok, nrow(idx)))

miss <- idx$uniprot[idx$pdb_status != "ok"]
if (length(miss)) {
  message("Missing in EBI (", length(miss), "): ",
          paste(utils::head(miss, 20), collapse = ", "),
          if (length(miss) > 20) ", ..." else "")
}

idx_path <- file.path(opt$out, "index.tsv")
utils::write.table(idx, idx_path, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "")
message("Wrote index: ", idx_path)
