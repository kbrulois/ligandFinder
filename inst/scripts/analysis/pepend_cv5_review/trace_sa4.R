suppressMessages({library(dplyr); library(purrr)})
S2 <- commandArgs(TRUE)[1]
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec
st_all <- readRDS(file.path(S2, "sa_status.rds"))
sa <- readRDS(file.path(S2, "sa_slim.rds"))

## ---- A. the NA-index rows: who, how long, recoverable by position?
na_acc <- unique(sa$accession[is.na(sa$index)])
L <- nchar(prec$seq[match(na_acc, prec$accession)])
cat(sprintf("NA-index accessions: %d (in cache: %d); length median %d, min %d; >2700 aa: %d\n",
    length(na_acc), sum(!is.na(L)), median(L, na.rm=TRUE), min(L, na.rm=TRUE), sum(L > 2700, na.rm=TRUE)))
cat("share of each accession's rows with NA index:\n")
print(summary(sapply(na_acc, function(a) mean(is.na(sa$index[sa$accession == a])))))
## do the NA rows sit past the last indexed residue (a fragment boundary)?
ex <- lapply(na_acc[1:min(6, length(na_acc))], function(a) {
  r <- sa[sa$accession == a, ]; w <- which(is.na(r$index))
  data.frame(acc = a, gene = r$gene[1], n_rows = nrow(r), len = nchar(prec$seq[prec$accession == a]),
             max_index = max(r$index, na.rm=TRUE), first_na_row = min(w), last_na_row = max(w),
             na_contiguous = all(diff(w) == 1))
}); print(bind_rows(ex))
## residue-absent residues in these accessions
in_na <- names(st_all) %in% na_acc
ra <- sapply(st_all, function(s) sum(s == "residue absent"))
cat(sprintf("residue-absent residues: %d total, %d (%.0f%%) inside the NA-index accessions\n",
    sum(ra), sum(ra[in_na]), 100*sum(ra[in_na])/sum(ra)))

## ---- B. what training windows look like: share of non-padding residues with no data
nn <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")$nn_input
frac0 <- function(d) {
  d <- as.data.frame(d); keep <- d$padding == 0
  if (!any(keep)) return(NA_real_)
  mean(rowSums(d[keep, ch4] == 0) == 4)
}
res <- map_dfr(c("N","C"), function(term) map_dfr(c("train","val"), function(sp) {
  t <- nn[[term]][[sp]]
  tibble(term = term, split = sp, win_type = t$win_type, f0 = map_dbl(t$data, frac0))
}))
cat("\nshare of a window's (non-padding) residues with no data, by window type:\n")
print(res %>% group_by(term, win_type) %>%
        summarise(n = n(), mean_f0 = round(mean(f0, na.rm=TRUE), 3),
                  any_gap = round(mean(f0 > 0, na.rm=TRUE), 3),
                  half_or_more = round(mean(f0 >= 0.5, na.rm=TRUE), 3),
                  all_missing = round(mean(f0 == 1, na.rm=TRUE), 3), .groups = "drop"), n = 40)
saveRDS(res, file.path(S2, "window_f0.rds"))
