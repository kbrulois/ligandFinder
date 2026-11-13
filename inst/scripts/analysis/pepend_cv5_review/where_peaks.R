## For each gene: where does the window score peak, and how far is that from the
## KNOWN peptide boundary? Peak height is not the question -- position is.
suppressMessages(library(dplyr))
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_allres.rds")
kn <- x$knowns %>% distinct(gene, term, anchor, pep_id, split, .keep_all = FALSE)
arms <- list(
  old = c(C = "~/AF2_analysis/lf_pepend_scan_C_near.rds", N = "~/AF2_analysis/lf_pepend_scan_N_near.rds"),
  new = c(C = "~/AF2_analysis/lf_pepend_scan_C_allres.rds", N = "~/AF2_analysis/lf_pepend_scan_N_allres.rds"))
genes <- c("NPY","ADM","MT-RNR2","CXCL14","ANO8","NUCB1")
out <- list()
for (arm in names(arms)) for (tm in c("N","C")) {
  p <- path.expand(arms[[arm]][[tm]]); if (!file.exists(p)) next
  s <- readRDS(p)$scan
  for (g in genes) {
    d <- s[s$gene == g, ]; if (!nrow(d)) next
    top <- d$anchor[which.max(d$score)]
    k <- kn$anchor[kn$gene == g & kn$term == tm]
    out[[length(out)+1]] <- tibble(gene = g, term = tm, arm = arm,
      peak_at = top, peak_score = round(max(d$score), 3),
      known = if (length(k)) paste(sort(k), collapse=",") else "-",
      dist_to_nearest_known = if (length(k)) min(abs(k - top)) else NA_integer_)
  }
}
r <- bind_rows(out) %>% tidyr::pivot_wider(id_cols = c(gene, term, known),
      names_from = arm, values_from = c(peak_at, peak_score, dist_to_nearest_known))
print(as.data.frame(r), row.names = FALSE)
