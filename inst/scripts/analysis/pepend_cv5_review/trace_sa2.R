suppressMessages({library(dplyr)})
S2 <- commandArgs(TRUE)[1]
slim_p <- file.path(S2, "sa_slim.rds")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")
prec <- fc$prec; ap <- fc$all_params3
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
if (!file.exists(slim_p)) {
  sa <- readRDS("~/peptide_alg/build_residue_db/processed/secretome_aa.rds")
  keep <- intersect(c("accession","gene","index","AA","topo","topo2","has_topo","location",
                      "cons_rs","species_limit", ch4, "SS","Phi","dssp_index"), names(sa))
  sa <- as.data.frame(sa)[, keep]
  saveRDS(sa, slim_p); message("slim cols: ", paste(keep, collapse=", "))
} else sa <- readRDS(slim_p)
message("secretome_aa rows: ", nrow(sa), "; NA index: ", sum(is.na(sa$index)),
        "; NA accession: ", sum(is.na(sa$accession)))
cat("\ntopo2 in secretome_aa:\n"); print(table(sa$topo2, useNA="ifany"))
if ("location" %in% names(sa)) { cat("\nlocation (rows):\n"); print(table(sa$location, useNA="ifany")) }

## per residue status in cached precursors
sa_rows <- split(seq_len(nrow(sa)), sa$accession)
st_all <- vector("list", nrow(prec))
for (i in seq_len(nrow(prec))) {
  acc <- prec$accession[i]; aa <- strsplit(prec$seq[i], "")[[1]]; n <- length(aa)
  st <- rep(if (is.null(sa_rows[[acc]])) "precursor absent" else "residue absent", n)
  r <- sa_rows[[acc]]
  if (length(r)) {
    s <- sa[r, ]; s <- s[which(!is.na(s$index) & s$index >= 1 & s$index <= n), ]
    mm <- s$AA != aa[s$index]
    st[s$index[mm]] <- "letter mismatch"
    ok <- s[!mm, ]; ok <- ok[!duplicated(ok$index), ]
    nna <- rowSums(!is.na(ok[, ch4]))
    st[ok$index] <- ifelse(nna == 0, "row, 4 ch all NA", ifelse(nna < 4, "row, some ch NA", "row, all 4 present"))
  }
  st_all[[i]] <- st
}
names(st_all) <- prec$accession
## agree with the all-four-zero marker in the cache?
z <- unlist(lapply(prec$accession, function(a) { f <- fc$feats[[a]]; rowSums(f[, ch4] == 0) == 4 }))
st <- unlist(st_all)
cat("\nresidue status, % of", length(st), "cached residues, and % all-four-zero in cache:\n")
print(cbind(pct = round(100*prop.table(table(st)), 2),
            pct_zero_in_cache = round(100*tapply(z, st, mean), 1)))
cat(sprintf("\nall-four-zero overall: %.1f%%\n", 100*mean(z)))
saveRDS(st_all, file.path(S2, "sa_status.rds"))

## gene-symbol hypothesis for the absent precursors
absent <- prec[vapply(prec$accession, function(a) is.null(sa_rows[[a]]), TRUE), ]
g_in <- absent$gene %in% sa$gene
cat(sprintf("\nabsent precursors: %d; their gene symbol IS in secretome_aa under another accession: %d\n",
            nrow(absent), sum(g_in)))
if (sum(g_in)) {
  ex <- absent[g_in, c("accession","gene")] %>% head(15)
  ex$sa_acc <- sapply(ex$gene, function(g) paste(unique(sa$accession[sa$gene == g]), collapse=","))
  print(ex)
}
saveRDS(absent, file.path(S2, "sa_absent_prec.rds"))
