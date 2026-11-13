suppressMessages({library(dplyr)})
S2 <- commandArgs(TRUE)[1]
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec
st_all <- readRDS(file.path(S2, "sa_status.rds"))
sa <- readRDS(file.path(S2, "sa_slim.rds"))
sec <- readRDS("~/AF2_analysis/secretome_latest.rds")
cat("secretome_latest:", nrow(sec), "rows\n"); print(intersect(c("accession","gene","location","topo","has_topo","sequence_uni","dssp","af_missense","cons"), names(sec)))
cat("location in secretome_latest:\n"); print(table(sec$location, useNA="ifany"))
sec <- sec %>% distinct(accession, .keep_all = TRUE)
ix <- match(prec$accession, sec$accession)
topo <- sec$topo[ix]; has_topo <- sec$has_topo[ix]; loc <- sec$location[ix]
cat("class of topo:", class(topo), "\n"); print(head(substr(topo, 1, 80), 3))

## topology letter of every cached residue, crossed with status
tl <- unlist(lapply(seq_len(nrow(prec)), function(i) {
  n <- nchar(prec$seq[i]); t <- topo[i]
  if (is.na(t) || nchar(t) != n) rep(if (is.na(t)) "no topo string" else "topo length != seq", n)
  else strsplit(t, "")[[1]]
}))
st <- unlist(st_all, use.names = FALSE)
ht <- rep(has_topo, nchar(prec$seq))
cat("\nresidue status x topology letter (counts):\n")
print(addmargins(table(status = st, topo = tl)))
cat("\nstatus x has_topo:\n"); print(table(status = st, has_topo = ht, useNA = "ifany"))

## absent precursors: why?
ab <- prec$accession[vapply(st_all, function(s) s[1] == "precursor absent", TRUE)]
ia <- match(ab, sec$accession)
cat("\nabsent precursors by location:\n"); print(table(sec$location[ia], useNA = "ifany"))
cat("absent precursors, fraction of residues topo 'e':\n")
fe <- sapply(sec$topo[ia], function(t) if (is.na(t)) NA else mean(strsplit(t,"")[[1]] == "e"))
print(summary(fe)); cat("has_topo:\n"); print(table(sec$has_topo[ia], useNA = "ifany"))

## NA-index rows
na <- sa[is.na(sa$index), ]
cat(sprintf("\nNA-index rows: %d in %d accessions; AA NA in %d; any of 4 ch non-NA in %d\n",
    nrow(na), length(unique(na$accession)), sum(is.na(na$AA)),
    sum(rowSums(!is.na(na[, c("cons_rs_n","min_afm","mean_afm","relASA")])) > 0)))
print(head(sort(table(na$accession), decreasing = TRUE), 8))
