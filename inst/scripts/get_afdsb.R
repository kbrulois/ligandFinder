## =============================================================================
## AlphaFold-predicted disulfides (afdsb) for every secretome accession.
##
##   /usr/local/bin/Rscript inst/scripts/get_afdsb.R          # all accessions
##   /usr/local/bin/Rscript inst/scripts/get_afdsb.R 200      # first 200 (bench)
##
## Writes a small table to ~/AF2_analysis/afdsb.rds rather than rewriting the
## 1.8 GB secretome.rds: plot_all.R merges it into secretome[["features"]] at
## load time, so this stays reversible (delete the file and the plots go back to
## UniProt-only disulfides).
##
## Same geometric criterion as get_dsb.R: CYS SG-SG pairs with 1.6 < d < 3 A.
## Parsed straight out of the PDB text -- bio3d::read.pdb parses every atom in
## every one of the ~5,600 models, and only the SG atoms are wanted here.
## =============================================================================
suppressMessages(library(tidyverse))

pdb_dir       <- "~/peptide_alg/UP000005640_9606_HUMAN_v4"
secretome_rds <- Sys.getenv("LF_SECRETOME", "~/AF2_analysis/secretome.rds")
out_path      <- "~/AF2_analysis/afdsb.rds"

extract_afdsb <- function(acc) {
  f <- file.path(path.expand(pdb_dir), paste0("AF-", acc, "-F1-model_v4.pdb"))
  if (!file.exists(f)) return(NULL)
  ln <- readLines(f, warn = FALSE)
  ln <- ln[startsWith(ln, "ATOM") &
             substr(ln, 13, 16) == " SG " & substr(ln, 18, 20) == "CYS"]
  if (length(ln) < 2) return(NULL)
  resno <- as.integer(substr(ln, 23, 26))
  d <- as.matrix(dist(cbind(as.numeric(substr(ln, 31, 38)),
                            as.numeric(substr(ln, 39, 46)),
                            as.numeric(substr(ln, 47, 54)))))
  hits <- which(d < 3 & d > 1.6, arr.ind = TRUE)
  hits <- hits[hits[, 1] < hits[, 2], , drop = FALSE]
  if (!nrow(hits)) return(NULL)
  tibble(accession = acc, type = "afdsb", evidence = "afdb", source = "afdb",
         start = resno[hits[, 1]], end = resno[hits[, 2]],
         description = as.character(d[hits]))
}

accs <- readRDS(path.expand(secretome_rds))[["accession"]]
n <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
if (!is.na(n)) accs <- head(accs, n)

message("scanning ", length(accs), " AlphaFold models in ", pdb_dir)
t0  <- Sys.time()
res <- vector("list", length(accs))
for (i in seq_along(accs)) {
  if (i %% 500 == 0)
    message(sprintf("  %d/%d  (%.1f min)", i, length(accs),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res[[i]] <- extract_afdsb(accs[i])
}
afdsb <- bind_rows(res)

missing <- sum(purrr::map_lgl(res, is.null))
message(sprintf("done in %.1f min: %d bonds across %d proteins (%d accessions had none or no model)",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                nrow(afdsb), dplyr::n_distinct(afdsb$accession), missing))

if (is.na(n)) {
  saveRDS(afdsb, path.expand(out_path))
  message("wrote ", out_path, "  (", round(file.size(path.expand(out_path)) / 2^10), " KB)")
} else {
  message("benchmark run (", n, " accessions) -- nothing written")
  print(head(afdsb, 10))
}
