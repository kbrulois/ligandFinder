## Recall vs rank for the xgb window head.
##   A. the three evaluation sets, negatives = every anchor that is not an
##      annotated end anywhere
##   B. the never-seen peptides by motif class, negatives = anchors of the SAME
##      class, so recall reaches 100% in every panel and the classes compare
##      despite very different prevalences
## A known end is ranked against the pool using the fold that HELD IT OUT; the
## never-seen sets use all 50, which is honest for them.
suppressMessages({library(dplyr); library(data.table); library(tidyr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
source("R/pepend_windows.R")
u  <- fread(file.path(S2, "cv5_unseen_ends.csv"))
ll <- as.data.frame(readRDS("inst/extdata/ligand_list.rds"))
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
seqs <- setNames(fc$prec$seq, fc$prec$accession)
NGRID <- unique(round(10^seq(0, 6.48, length.out = 150)))
curve <- function(rk, n_pool) tibble(N = NGRID, recall = sapply(NGRID, \(k) mean(rk <= k)))
A <- list(); B <- list()
for (tm in c("N","C")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  setorder(b, accession, anchor)
  b[, motif := { s <- seqs[[.BY$accession]]
      if (is.null(s)) rep(NA_character_, .N) else lf_pepend_motif_vec(s, tm, anchor, 1L) }, by = accession]
  b[, key := paste(accession, anchor)]
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  ann <- ll %>% filter(!is.na(start), !is.na(end)) %>%
    transmute(accession, anchor = if (tm=="N") start else end) %>% distinct()
  all_ann <- union(paste(ann$accession, ann$anchor), paste(kn$accession, kn$anchor))

  ## ---- A1: cross-validated, each end ranked by its own held-out fold --------
  rk <- integer(0)
  for (k in 1:5) {
    pos <- kn[fold == k]; pk <- paste(pos$accession, pos$anchor)
    neg <- b[!(key %in% all_ann)][[sprintf("xgb_pm__f%d", k)]]
    neg <- sort(neg[!is.na(neg)], decreasing = TRUE)
    sp  <- b[key %in% pk][[sprintf("xgb_pm__f%d", k)]]
    rk  <- c(rk, findInterval(-sp, -neg) + 1L)
  }
  A[[length(A)+1]] <- curve(rk) %>% mutate(term = tm, set = "cross-validated", n = length(rk))

  ## ---- A2/A3: the never-seen sets, all-50 ----------------------------------
  negall <- sort(b[!(key %in% all_ann)]$xgb_pm__all50, decreasing = TRUE)
  for (nm in c("never seen (peptide, < 60 aa)", "large cleavage product ends")) {
    p <- u[set == nm & term == tm]
    sp <- b[key %in% paste(p$accession, p$anchor)]$xgb_pm__all50
    sp <- sp[!is.na(sp)]
    A[[length(A)+1]] <- curve(findInterval(-sp, -negall) + 1L) %>%
      mutate(term = tm, set = if (grepl("^never", nm)) "never seen: peptides"
                              else "never seen: large products", n = length(sp))
  }

  ## ---- B: never-seen peptides by motif, negatives of the SAME class --------
  p <- u[set == "never seen (peptide, < 60 aa)" & term == tm]
  pm <- b[key %in% paste(p$accession, p$anchor), .(key, motif, s = xgb_pm__all50)]
  for (g in c("dibasic","monobasic","G + dibasic")) {
    ps <- pm[motif == g & !is.na(s)]$s
    if (!length(ps)) next
    ng <- sort(b[motif == g & !(key %in% all_ann)]$xgb_pm__all50, decreasing = TRUE)
    B[[length(B)+1]] <- curve(findInterval(-ps, -ng) + 1L) %>%
      mutate(term = tm, group = g, n = length(ps), pool = length(ng))
  }
  rm(b); invisible(gc()); message(tm, " done")
}
saveRDS(list(A = bind_rows(A), B = bind_rows(B)), file.path(S2, "rr_curves.rds"))
cat("\n-- A: ends per set --\n"); print(bind_rows(A) %>% distinct(term, set, n) %>% as.data.frame(), row.names=FALSE)
cat("\n-- B: never-seen peptides by motif --\n"); print(bind_rows(B) %>% distinct(term, group, n, pool) %>% as.data.frame(), row.names=FALSE)
