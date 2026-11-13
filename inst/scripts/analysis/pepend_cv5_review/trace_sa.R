suppressMessages({library(dplyr)})
t0 <- Sys.time()
sa <- readRDS("~/peptide_alg/build_residue_db/processed/secretome_aa.rds")
message("loaded in ", round(difftime(Sys.time(), t0, units="secs")), " s; class ", paste(class(sa), collapse="/"))
message("dim ", paste(dim(sa), collapse=" x "))
print(names(sa))
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")
prec <- fc$prec; ap <- fc$all_params3
scale_cols <- setdiff(ap, c("AA_hydro","AA_charge","AA_mw","AA_pI","padding"))
message("scale cols: ", paste(scale_cols, collapse=", "))
saveRDS(list(names = names(sa), n = nrow(sa)), file.path(commandArgs(TRUE)[1], "sa_meta.rds"))

## 1. coverage: residues in prec vs rows in sa
sa <- as.data.frame(sa)
sa_acc <- unique(sa$accession)
prec$in_sa <- prec$accession %in% sa_acc
message(sprintf("precursors in cache: %d; with ANY secretome_aa row: %d (%.1f%%)",
                nrow(prec), sum(prec$in_sa), 100*mean(prec$in_sa)))
tot_res <- sum(nchar(prec$seq))
message(sprintf("residues in cache: %d; residues in precursors absent from sa: %d (%.1f%%)",
                tot_res, sum(nchar(prec$seq[!prec$in_sa])), 100*sum(nchar(prec$seq[!prec$in_sa]))/tot_res))

## 2. per residue classification, for precursors present in sa
sa_rows <- split(seq_len(nrow(sa)), sa$accession)
res <- lapply(which(prec$in_sa), function(i) {
  acc <- prec$accession[i]; aa <- strsplit(prec$seq[i], "")[[1]]; n <- length(aa)
  s <- sa[sa_rows[[acc]], c("index","AA", scale_cols)]
  s <- s[s$index >= 1 & s$index <= n, ]
  st <- rep("absent", n)
  mism <- s$index[s$AA != aa[s$index]]
  okr <- s[s$AA == aa[s$index], ]; okr <- okr[!duplicated(okr$index), ]
  st[mism] <- "letter mismatch"
  allna <- rowSums(!is.na(okr[, scale_cols])) == 0
  st[okr$index] <- ifelse(allna, "row present, all NA", "row present, has data")
  ## which channels are NA where row present
  data.frame(acc = acc, n = n, absent = sum(st=="absent"), mism = sum(st=="letter mismatch"),
             allna = sum(st=="row present, all NA"), ok = sum(st=="row present, has data"),
             sa_rows = length(sa_rows[[acc]]), sa_maxidx = max(sa$index[sa_rows[[acc]]]),
             sa_dupidx = sum(duplicated(sa$index[sa_rows[[acc]]])))
})
res <- bind_rows(res)
saveRDS(res, file.path(commandArgs(TRUE)[1], "sa_trace_res.rds"))
tot <- colSums(res[, c("absent","mism","allna","ok")])
abs_prec <- sum(nchar(prec$seq[!prec$in_sa]))
cat("\nresidue status over all", tot_res, "cached residues:\n")
print(round(100 * c(precursor_absent = abs_prec, tot) / tot_res, 2))

## 3. channel-wise NA rate where a row is present and letter matches
m <- sa[, scale_cols]
cat("\nNA rate per channel in secretome_aa rows:\n")
print(round(100*colMeans(is.na(m)), 2))
