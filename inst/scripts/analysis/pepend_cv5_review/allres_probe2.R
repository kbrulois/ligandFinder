suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
x <- readRDS(file.path(S2, "sa1_slim.rds")); chk <- readRDS(file.path(S2, "sa1_chk.rds"))
bad <- chk$acc[!(chk$idx_eq_pos %in% TRUE)]
pat <- bind_rows(lapply(bad, function(a) {
  r <- x[x$accession == a, ]; pos <- seq_len(nrow(r)); ok <- !is.na(r$index)
  off <- r$index[ok] - pos[ok]
  data.frame(acc = a, gene = r$gene[1], n = nrow(r), n_idx = sum(ok), n_off = sum(off != 0),
             offsets = paste(head(names(sort(table(off), decreasing = TRUE)), 3), collapse = ","),
             first_bad_pos = pos[ok][which(off != 0)[1]],
             relASA_na_where_off = round(mean(is.na(r$relASA[ok][off != 0])), 2),
             dup_idx = sum(duplicated(r$index[ok])))
}))
print(pat, row.names = FALSE)
cat("\nresidues affected (index != position):", sum(pat$n_off), "\n")
