#!/usr/bin/env Rscript
#
# Which cell type expresses each gene most: the companion to the eec_* block.
#
# eec_top_other_celltype names the strongest cell type EXCLUDING enteroendocrine
# cells, because that block asks how EECs compare with the rest of the body.
# This one asks the plainer question -- where does the gene peak, EECs included
# -- and answers per gene from the same table:
#
#   top_celltype            HPA single-cell type with the highest nTPM
#                           ("not_detected" when nothing reaches DETECT)
#   top_celltype_ntpm       that nTPM
#   top_celltype_fc_vs_2nd  fold over the second-highest cell type; ~1 means the
#                           gene has no single home, >= 4 is HPA's enrichment bar
#
# HPA v23 single-cell consensus (81 cell types), the last release with a
# standalone "Enteroendocrine cells" cluster -- the same table the eec_* columns
# were read from, so the two blocks agree on every gene. The v23 download is
# cached under <cache_dir>; nothing is fetched while it is there.
#
# Two outputs, both optional:
#   <out_dir>/hpa_top_celltype.csv   one row per HPA gene (symbol + Ensembl id)
#   --append <predictions.csv>       the three columns appended IN PLACE to a
#                                    predictions table keyed by gene symbol, a
#                                    timestamped .bak written first. Rows are
#                                    never re-serialised: the fields are pasted
#                                    onto the raw lines, so everything already
#                                    there survives byte for byte. Symbols are
#                                    matched directly, then rescued through HPA
#                                    synonyms only when the synonym is
#                                    unambiguous and the HPA gene is not already
#                                    a row of its own (the HLA-H / HFE trap).
#
# Usage:
#   Rscript inst/scripts/hpa_top_celltype.R [out_dir] [cache_dir] [--append predictions.csv]
#   defaults: out_dir ~/AF2_analysis, cache_dir <out_dir>/hpa_cache

suppressPackageStartupMessages({
  library(data.table)
})

DETECT <- 1     # nTPM >= this counts as detected (HPA convention, as in eec_*)

args      <- commandArgs(trailingOnly = TRUE)
ap        <- match("--append", args)
pred_path <- if (!is.na(ap)) path.expand(args[ap + 1]) else NA_character_
if (!is.na(ap)) args <- args[-c(ap, ap + 1)]
out_dir   <- path.expand(if (length(args) >= 1) args[1] else "~/AF2_analysis")
cache_dir <- path.expand(if (length(args) >= 2) args[2] else file.path(out_dir, "hpa_cache"))
dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

## ---- fetch (cached) ---------------------------------------------------------
fetch_hpa <- function(url, zip_name, tsv_name) {
  zip_path <- file.path(cache_dir, zip_name)
  tsv_path <- file.path(cache_dir, tsv_name)
  if (!file.exists(tsv_path)) {
    if (!file.exists(zip_path)) {
      message("downloading ", url)
      utils::download.file(url, zip_path, mode = "wb", quiet = TRUE)
    }
    utils::unzip(zip_path, files = tsv_name, exdir = cache_dir)
  }
  tsv_path
}
sc_path <- fetch_hpa("https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip",
                     "rna_single_cell_type.tsv.zip", "rna_single_cell_type.tsv")
pa_path <- fetch_hpa("https://v23.proteinatlas.org/download/proteinatlas.tsv.zip",
                     "proteinatlas.tsv.zip", "proteinatlas.tsv")

## ---- per gene: the peak cell type -------------------------------------------
sc <- data.table::fread(sc_path, showProgress = FALSE)
setnames(sc, c("Gene", "Gene name", "Cell type", "nTPM"),
             c("ensembl", "gene", "cell_type", "ntpm"))
stopifnot("Enteroendocrine cells" %in% sc$cell_type)

wide <- data.table::dcast(sc, ensembl + gene ~ cell_type, value.var = "ntpm")
key  <- wide[, .(ensembl, gene)]
mat  <- as.matrix(wide[, !c("ensembl", "gene"), with = FALSE])
stopifnot(!anyNA(mat))

## ties go to the first column alphabetically, as eec_top_other_celltype did
top_i  <- max.col(mat, ties.method = "first")
top    <- mat[cbind(seq_len(nrow(mat)), top_i)]
second <- apply(mat, 1, function(v) sort(v, decreasing = TRUE)[2])

res <- data.table(
  gene                   = key$gene,
  ensembl                = key$ensembl,
  top_celltype           = ifelse(top >= DETECT, colnames(mat)[top_i], "not_detected"),
  top_celltype_ntpm      = round(top, 1),
  top_celltype_fc_vs_2nd = ifelse(top >= DETECT, round(top / pmax(second, 0.01), 2), NA_real_)
)

out_file <- file.path(out_dir, "hpa_top_celltype.csv")
data.table::fwrite(res, out_file)

cat("HPA v23 single-cell consensus | cell types:", ncol(mat), "| genes:", nrow(res), "\n")
cat("wrote ", out_file, "\n\n", sep = "")
cat("genes by top cell type (top 15 of ", ncol(mat), ")\n", sep = "")
print(res[, .N, by = top_celltype][order(-N)][1:15])
cat("\nEEC is the top cell type for", sum(res$top_celltype == "Enteroendocrine cells"),
    "genes; of those,", sum(res$top_celltype == "Enteroendocrine cells" &
                            res$top_celltype_fc_vs_2nd >= 4), "clear a 4-fold gap\n")

## ---- append to a predictions table, in place ----------------------------------
if (!is.na(pred_path)) {
  stopifnot(file.exists(pred_path))
  add <- c("top_celltype", "top_celltype_ntpm", "top_celltype_fc_vs_2nd")

  lines <- readLines(pred_path, warn = FALSE)
  pred  <- data.table::fread(pred_path)
  if (length(lines) != nrow(pred) + 1L)
    stop("raw line count does not match parsed rows; embedded newlines?")
  stopifnot("gene" %in% names(pred))
  dup <- intersect(add, names(pred))
  if (length(dup))
    stop("already present in ", pred_path, ": ", paste(dup, collapse = ", "),
         " -- refusing to append twice")

  ## gene -> row of res: direct, then the guarded synonym rescue
  idx <- match(pred$gene, res$gene)
  if (anyNA(idx)) {
    miss <- setdiff(unique(pred$gene[is.na(idx)]), res$gene)
    pa <- data.table::fread(pa_path, quote = "\"",
                            select = c("Gene", "Gene synonym", "Ensembl"))
    setnames(pa, c("gene", "synonym", "ensembl"))
    syn <- pa[nzchar(synonym), .(s = trimws(unlist(strsplit(synonym, ",")))),
              by = .(gene, ensembl)][nzchar(s) & s %in% miss]
    syn <- unique(syn[s %in% syn[, .(n = uniqueN(ensembl)), by = s][n == 1, s]], by = "s")
    syn <- syn[!gene %in% pred$gene]
    if (nrow(syn)) {
      hit <- is.na(idx) & pred$gene %in% syn$s
      idx[hit] <- match(syn$gene[match(pred$gene[hit], syn$s)], res$gene)
    }
  }

  ## genes with no HPA record: "no_expression_data", kept apart from a measured
  ## "not_detected", exactly as eec_class does it
  fmt <- function(x) {
    if (is.character(x)) ifelse(is.na(x), '"no_expression_data"', paste0('"', x, '"'))
    else ifelse(is.na(x), "NA", as.character(x))
  }
  block <- do.call(paste, c(lapply(add, function(col) fmt(res[[col]][idx])), sep = ","))

  backup <- paste0(pred_path, ".bak-", format(Sys.time(), "%Y%m%d%H%M%S"))
  file.copy(pred_path, backup, overwrite = FALSE)
  writeLines(c(paste0(lines[1], ",", paste0('"', add, '"', collapse = ",")),
               paste0(lines[-1], ",", block)), pred_path)

  check <- data.table::fread(pred_path)
  stopifnot(nrow(check) == nrow(pred),
            identical(names(check), c(names(pred), add)),
            identical(check$peps, pred$peps))
  cat("\nbackup  ", backup, "\nmodified", pred_path,
      "\nrows", nrow(check), "| columns", ncol(pred), "->", ncol(check), "\n")
  cat("matched genes:", sum(!is.na(idx)), "/", nrow(pred), "rows (",
      length(unique(pred$gene[is.na(idx)])), "genes without an HPA record )\n")
  cat("\nwindows by top cell type (top 12)\n")
  print(check[, .(windows = .N, genes = uniqueN(gene)), by = top_celltype][order(-windows)][1:12])
}
