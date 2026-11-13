## PCA over peptide-end windows, now in THREE sets:
##   known                        178 ends, training in 4 folds, held out in 1
##   never seen (peptide, < 60 aa)  short cleaved peptides no fold ever saw
##   large cleavage product ends    >= 60 aa products -- GDF5/GDF6, the BMPs,
##                                  TGFBs, inhibins: ends the length filter
##                                  had been excluding all along
## One PCA over all of them, so the facets share axes. Scores come from the
## exhaustive CV scans, which already cover every mature residue, so nothing
## is re-scored here.
suppressMessages({library(dplyr); library(data.table); library(tidyr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
source("R/pepend_windows.R")
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds")); prec <- fc$prec
ap <- fc$all_params3
uns <- fread(file.path(S2, "cv5_unseen_ends.csv"))
KNOWN <- "known (train in 4 folds, held out in 1)"

feat_of <- function(acc) {
  j <- match(acc, prec$accession); if (is.na(j)) return(NULL)
  list(m = fc$feats[[j]], n = prec$n_prot[j], c = prec$c_prot[j], s = prec$seq[j])
}

rows <- list()
add <- function(set, term, label, gene, fold, motif, s_cv, s_all, len, v)
  rows[[length(rows) + 1L]] <<- list(set = set, term = term, label = label, gene = gene,
    fold = fold, motif = motif, score_cv = s_cv, score_all50 = s_all, len = len, v = v)

for (tm in c("C", "N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  b[, key := paste(accession, anchor)]; setkey(b, key)
  score <- function(a, anc) { s <- b[.(paste(a, anc))]; if (!nrow(s)) NA_real_ else s$xgb_pm__all50[1] }

  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  for (i in seq_len(nrow(kn))) {
    f <- feat_of(kn$accession[i]); if (is.null(f)) next
    x <- lf_pepend_slice(f$m, lf_pepend_start(tm, kn$anchor[i]), f$n, f$c)
    add(KNOWN, tm, kn$peps[i], sub("_w[0-9]+-[0-9]+$", "", kn$peps[i]),
        kn$fold[i], kn$motif[i],
        kn$xgb_pm__cv[i], score(kn$accession[i], kn$anchor[i]), NA_real_,
        as.numeric(as.matrix(x[, ap])))
  }

  u <- uns[term == tm]
  for (i in seq_len(nrow(u))) {
    f <- feat_of(u$accession[i]); if (is.null(f)) next
    x <- lf_pepend_slice(f$m, lf_pepend_start(tm, u$anchor[i]), f$n, f$c)
    sc <- score(u$accession[i], u$anchor[i])          # honest: no fold trained on it
    add(u$set[i], tm, u$name[i], u$gene[i], NA_integer_,
        lf_pepend_motif(list(f$s), tm, list(u$anchor[i]), list(f$n)),
        sc, sc, u$len[i], as.numeric(as.matrix(x[, ap])))
  }
}

M <- do.call(rbind, lapply(rows, `[[`, "v"))
chr <- function(f) vapply(rows, function(r) { v <- r[[f]]
  if (length(v) != 1L) stop("field ", f, " is not a scalar") else as.character(v) }, "")
num <- function(f) vapply(rows, function(r) { v <- r[[f]]
  if (length(v) != 1L) stop("field ", f, " is not a scalar") else as.numeric(v) }, 0)
meta <- tibble(set = chr("set"), term = chr("term"), label = chr("label"),
               gene = chr("gene"), fold = num("fold"), motif = chr("motif"),
               score_cv = num("score_cv"), score_all50 = num("score_all50"),
               len = num("len"))
keep <- apply(M, 2, sd) > 0
pc <- prcomp(M[, keep], center = TRUE, scale. = TRUE)
ve <- round(100 * summary(pc)$importance[2, 1:4], 1)
meta <- bind_cols(meta, as_tibble(pc$x[, 1:4]))
cen <- colMeans(pc$x[, 1:10, drop = FALSE])
meta$d_cen <- sqrt(rowSums(sweep(pc$x[, 1:10, drop = FALSE], 2, cen)^2))
saveRDS(list(meta = meta, ve = ve), file.path(S2, "cv5_pca4.rds"))

cat(sprintf("\nwindows %d (%d channels used) | PC1-4 %s\n", nrow(M), sum(keep),
            paste0(ve, "%", collapse = " ")))
cat("\n=== windows per set ===\n")
print(meta %>% count(set, term) %>% pivot_wider(names_from = term, values_from = n) %>%
      as.data.frame(), row.names = FALSE)
cat("\n=== end type by set ===\n"); print(table(meta$set, meta$motif))
cat("\n=== score by set (honest score: held-out fold for knowns, all-50 otherwise) ===\n")
print(meta %>% filter(!is.na(score_cv)) %>% group_by(set, term) %>%
      summarise(n = n(), median = round(median(score_cv), 3),
                q25 = round(quantile(score_cv, .25), 3),
                q75 = round(quantile(score_cv, .75), 3), .groups = "drop") %>%
      as.data.frame(), row.names = FALSE)
cat("\n=== GDF5 / GDF6 and friends ===\n")
print(meta %>% filter(gene %in% c("GDF5","GDF6","BMP2","BMP4","TGFB1","INHBB")) %>%
      transmute(gene, term, motif, score = round(score_cv, 3), PC1 = round(PC1, 2),
                PC2 = round(PC2, 2)) %>% arrange(gene, term) %>% as.data.frame(),
      row.names = FALSE)
