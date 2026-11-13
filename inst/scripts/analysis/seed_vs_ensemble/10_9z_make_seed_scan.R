## A scan rds carrying ONE member's score, shaped like the ensemble scans so the
## page's window-head track can read it unchanged. sd is 0: a single seed has no
## spread, and drawing the ensemble's would misrepresent it.
suppressMessages(library(dplyr))
A2 <- path.expand("~/AF2_analysis")
np <- reticulate::import("numpy", convert = FALSE)
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))$prec
for (spec in list(list(tag = "C_near", seed = 14L, out = "lf_pepend_scan_C_seed14.rds"),
                  list(tag = "N_near", seed = 7L,  out = "lf_pepend_scan_N_seed7.rds"))) {
  z <- np$load(file.path(A2, sprintf("lf_pepend_scanmem_%s.npz", spec$tag)))
  M <- reticulate::py_to_r(z[["members"]])
  prot <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
  anc  <- as.integer(reticulate::py_to_r(z[["anchor"]]))
  sc   <- as.numeric(M[spec$seed + 1L, ])
  ens  <- as.numeric(reticulate::py_to_r(z[["score"]]))
  scan <- tibble(accession = fc$accession[prot + 1L], gene = fc$gene[prot + 1L],
                 anchor = anc, score = sc, sd = 0,
                 rank_all = rank(-sc, ties.method = "min"))
  saveRDS(list(scan = scan, seed = spec$seed), file.path(A2, spec$out))
  message(sprintf("%s seed %d: %s anchors, cor with the ensemble %.3f -> %s",
                  spec$tag, spec$seed, format(nrow(scan), big.mark = ","),
                  cor(sc, ens), spec$out))
}
