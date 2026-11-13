#!/usr/bin/env Rscript
#
# Append an eec_class column to a predictions CSV *in place*.
#
# The original rows are never re-serialised: the file is read as raw lines and
# the new field is appended textually, so numeric formatting, quoting and column
# order survive byte for byte. A timestamped .bak copy is written first.
#
# Refuses to run twice -- an existing eec_class column is an error, not an
# invitation to append a second one.
#
# Usage:
#   Rscript inst/scripts/add_eec_class_to_predictions.R <predictions.csv> \
#           <eec_gene_class.csv> [proteinatlas.tsv]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("need <predictions.csv> <eec_gene_class.csv>")

pred_path <- path.expand(args[1])
eec_path  <- path.expand(args[2])
pa_path   <- if (length(args) >= 3) path.expand(args[3]) else NA_character_

stopifnot(file.exists(pred_path), file.exists(eec_path))

lines <- readLines(pred_path, warn = FALSE)
pred  <- data.table::fread(pred_path)
eec   <- data.table::fread(eec_path)

if ("eec_class" %in% names(pred))
  stop("eec_class already present in ", pred_path, " -- refusing to append twice")
if (length(lines) != nrow(pred) + 1L)
  stop("raw line count does not match parsed rows; embedded newlines?")
stopifnot("gene" %in% names(pred), all(c("gene", "class") %in% names(eec)))

## ---- gene symbol -> class -------------------------------------------------

lut <- setNames(eec$class, eec$gene)
cls <- unname(lut[pred$gene])

# Rescue unmatched symbols through unambiguous HPA synonyms, but never when the
# HPA target already stands on its own in this file (e.g. HLA-H vs HFE).
if (!is.na(pa_path) && file.exists(pa_path) && anyNA(cls)) {
  missing <- setdiff(unique(pred$gene[is.na(cls)]), eec$gene)
  pa  <- data.table::fread(pa_path, quote = "\"",
                           select = c("Gene", "Gene synonym", "Ensembl"))
  setnames(pa, c("gene", "synonym", "ensembl"))
  syn <- pa[nzchar(synonym),
            .(s = trimws(unlist(strsplit(synonym, ",")))), by = .(gene, ensembl)]
  syn <- syn[nzchar(s) & s %in% missing]
  syn <- unique(syn[s %in% syn[, .(n = uniqueN(ensembl)), by = s][n == 1, s]], by = "s")
  syn <- syn[!gene %in% pred$gene]
  if (nrow(syn)) {
    rescue <- setNames(unname(lut[syn$gene]), syn$s)
    hit    <- is.na(cls) & pred$gene %in% names(rescue)
    cls[hit] <- unname(rescue[pred$gene[hit]])
    cat("rescued via synonym:", nrow(syn), "symbols\n")
  }
}

# No HPA record is not a measured zero, so it gets its own value.
cls[is.na(cls)] <- "no_expression_data"
stopifnot(length(cls) == nrow(pred), !anyNA(cls))

## ---- write ----------------------------------------------------------------

backup <- paste0(pred_path, ".bak-", format(Sys.time(), "%Y%m%d%H%M%S"))
file.copy(pred_path, backup, overwrite = FALSE)

out <- c(paste0(lines[1], ',"eec_class"'),
         paste0(lines[-1], ',"', cls, '"'))
writeLines(out, pred_path)

## ---- verify ---------------------------------------------------------------

check <- data.table::fread(pred_path)
stopifnot(nrow(check) == nrow(pred),
          identical(names(check), c(names(pred), "eec_class")),
          identical(check$peps, pred$peps),
          identical(check$score_mean, pred$score_mean),
          identical(check$eec_class, cls))

cat("backup   ", backup, "\n")
cat("modified ", pred_path, "\n")
cat("rows", nrow(check), "| columns", ncol(pred), "->", ncol(check), "\n\n")
print(check[, .N, by = eec_class][order(-N)])
