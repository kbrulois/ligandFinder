## Per-residue extra channels for the page genes, as the pages want them.
## Kevin, 2026-10-25.
##
## WHY THIS EXISTS. A page re-scores every series ON THE FLY from its members,
## so a series whose trunk was trained with `--plm` needs those channels on the
## feature table at render time. `R/pepend_tracks.R` takes them from
## `spec$plm_rds`: a named list, accession -> matrix of (residues x channels),
## one row per residue of the precursor IN ORDER, which it `cbind`s onto the
## curated 26 before slicing windows.
##
## The ESM arm's `esm_c_pca36_pagegenes.rds` was built by hand and the builder
## was never kept, so the next arm had nothing to copy. This is that builder,
## parameterised: any `plm_*` parquet in, one pagegenes rds out.
##
##   PLM=~/AF2_analysis/lf_plm/aa_onehot.parquet \
##   OUT=~/AF2_analysis/lf_plm/aa_onehot_pagegenes.rds \
##   Rscript plm_pagegenes.R
##
## Residues the parquet does not cover get ZEROS -- the same policy 10_8 applies
## at training time and `plm_residues` applies to the scan table, so the page
## scores the same windows the scan did. The count is printed, not assumed.
suppressPackageStartupMessages({library(dplyr)})
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
A2  <- path.expand("~/AF2_analysis")
S   <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
PLM <- path.expand(Sys.getenv("PLM", file.path(A2, "lf_plm", "aa_onehot.parquet")))
OUT <- path.expand(Sys.getenv("OUT", sub("\\.parquet$", "_pagegenes.rds", PLM)))
GENES_F <- Sys.getenv("GENES_FILE", file.path(S, "page_genes.txt"))

source(lf_src("R/plm_features.R"))   # lf_read_parquet

genes <- readLines(GENES_F)
genes <- genes[nzchar(genes)]
prec  <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))$prec
i     <- match(genes, prec$gene)
if (anyNA(i)) stop("genes not in the residue cache: ", paste(genes[is.na(i)], collapse = ", "))
prec  <- prec[i, , drop = FALSE]
message(sprintf("page genes: %d, %s residues", nrow(prec), format(sum(nchar(prec$seq)), big.mark = ",")))

tbl  <- lf_read_parquet(PLM)
cols <- grep("^plm_", names(tbl), value = TRUE)
if (!length(cols)) stop("no plm_* columns in ", PLM)
message(sprintf("parquet: %s rows x %d channels (%s .. %s)",
                format(nrow(tbl), big.mark = ","), length(cols), cols[1], cols[length(cols)]))

tbl <- tbl[tbl$accession %in% prec$accession, , drop = FALSE]
sp  <- split(seq_len(nrow(tbl)), tbl$accession)

zero_rows <- 0L; no_rows <- character(0)
out <- lapply(seq_len(nrow(prec)), function(k) {
  acc <- as.character(prec$accession[k]); n <- nchar(prec$seq[k])
  m <- matrix(0, nrow = n, ncol = length(cols), dimnames = list(NULL, cols))
  r <- sp[[acc]]
  if (is.null(r)) { no_rows <<- c(no_rows, acc); zero_rows <<- zero_rows + n; return(m) }
  idx <- as.integer(tbl$index[r])
  ok  <- which(idx >= 1L & idx <= n)
  m[idx[ok], ] <- as.matrix(tbl[r[ok], cols])
  zero_rows <<- zero_rows + (n - length(ok))
  m
})
names(out) <- as.character(prec$accession)

## the shape the page asserts: one matrix per accession, rows == residues
stopifnot(all(vapply(out, is.matrix, TRUE)),
          identical(unname(vapply(out, nrow, 1L)), as.integer(nchar(prec$seq))))
message(sprintf("zero-filled residues: %s of %s (%.2f%%)%s",
                format(zero_rows, big.mark = ","),
                format(sum(nchar(prec$seq)), big.mark = ","),
                100 * zero_rows / sum(nchar(prec$seq)),
                if (length(no_rows)) sprintf("; %d precursor(s) absent entirely: %s",
                                             length(no_rows),
                                             paste(utils::head(no_rows, 5), collapse = ", ")) else ""))
saveRDS(out, OUT)
message("-> ", OUT, "  (", length(out), " precursors x ", length(cols), " channels)")
