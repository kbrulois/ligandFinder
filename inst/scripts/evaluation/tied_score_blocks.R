## Tied-score blocks in an arm's full-proteome tracks (found 2026-11-10, while
## reading the AAswap-vs-ESM-C metrics plots).
##
## Every ESM-C (esmonly36) track holds ~60,000 anchors tied at ONE score. The
## 434 precursors esm_c_pca36.parquet does not cover get all-zero ESM-C input,
## so all their windows that hold no padding -- ~60,000 of their ~75,000
## anchors -- are one and the same input. A block that size moves any ranking
## metric wherever it sits. On N it sits high -- all-50 linear 0.2310 at ranks
## ~54.6k-114.6k, xgb 0.2567 at ~52.6k-111.6k, fold-1 xgb 0.5934 at
## 11.5k-70.5k -- which makes the flat stretches in ESM-C's N recall-vs-rank
## curves and inflates its N anchors-for-80%-recall (e.g. cross-validated N xgb
## 90,751). C's block sits deeper (~190k-360k). Compare arms on embedded
## precursors when the depth of a ranking matters.
## The block_unembedded column checks the attribution: 2026-11-12, 0.9999-1.0
## on all eight esmonly36 columns.
##
## A 6-decimal bin is not a tie by itself. Near zero it merges thousands of
## distinct scores: curated_AAswap's fold-1 C CNN bin 0.000016 holds 303,420
## anchors and 225,059 distinct values, and its most common N all-50 bin (~3,000
## anchors) holds ~3,000 values, the largest exact tie 22. n_exact and n_values
## tell the two apart: the ESM-C blocks are 1-13 values, one of them holding
## 98-100% of the bin (the rest within 5e-7 of it, float noise on the same
## input). An arm with real input on every precursor has no block where it
## matters: curated_AAswap's largest exact tie in these columns, 26,361 anchors
## on fold-1 C xgb, is that head's floor, the bottom of the ranking.
##
## Per arm, terminus and track column: the most common score rounded to 6
## decimals (modal) and how many anchors carry it (n_modal); of those, the most
## that hold one exact value (n_exact) and how many distinct values there are
## (n_values); the rank span of the bin (1 = the highest score; every anchor in
## a higher bin ranks above it, and ties are broken by nothing, so the block
## spans its whole width); the next most common bin's count for scale; the
## share of the bin on precursors with no row in EMB (block_unembedded; 1 =
## every anchor in it is on one) and how many anchors those precursors have
## (n_unembedded). Columns: cnn__f<k> / xgb_pm__f<k> are fold k's ten members,
## *__all50 the mean of the five folds (tracks_esmonly36.R / metrics_final.R
## flat_tracks).
##
##   LC_ALL=en_US.UTF-8 Rscript tied_score_blocks.R
##   ARMS=esmonly36 COLS=cnn__all50 HEADCMP=/tmp/tsb LC_ALL=en_US.UTF-8 Rscript tied_score_blocks.R
##
## Reads ~/AF2_analysis/lf_pepend_cv5_tracks_<T>_flatten_<arm>.rds and the
## accession column of EMB (default ~/AF2_analysis/lf_plm/esm_c_pca36.parquet).
## Writes HEADCMP/tied_score_blocks.csv (default ~/AF2_analysis/lf_headcmp).
## ~1.2 GB of RAM, ~15 s for the default four track files.
suppressMessages(library(data.table))
A2 <- path.expand("~/AF2_analysis")
ARMS <- strsplit(Sys.getenv("ARMS", "esmonly36,curated_AAswap"), ",")[[1]]
COLS <- strsplit(Sys.getenv("COLS", "cnn__f1,xgb_pm__f1,cnn__all50,xgb_pm__all50"), ",")[[1]]
EMB  <- path.expand(Sys.getenv("EMB", file.path(A2, "lf_plm", "esm_c_pca36.parquet")))
HC   <- path.expand(Sys.getenv("HEADCMP", file.path(A2, "lf_headcmp")))
dir.create(HC, showWarnings = FALSE, recursive = TRUE)

emb <- unique(arrow::read_parquet(EMB, col_select = "accession")$accession)
cat(sprintf("%s covers %d precursors\n", basename(EMB), length(emb)))

res <- list()
for (a in ARMS) for (tm in c("C", "N")) {
  f <- file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm,
                             if (a %in% c("", "curated", "flatten")) "" else paste0("_", a)))
  if (!file.exists(f)) { message("no tracks: ", f); next }
  b <- readRDS(f)
  unemb <- !b$accession %in% emb
  for (col in COLS) {
    if (!col %in% names(b)) { message(a, " ", tm, ": no column ", col); next }
    v <- b[[col]]; ok <- !is.na(v); v <- v[ok]; vr <- round(v, 6)
    n <- data.table(r = vr)[, .N, by = r][order(-N)]
    top <- n$r[1]; blk <- vr == top
    ## above the bin = in any higher bin, the next one up included (a score
    ## between top + 5e-7 and top + 1e-6 still ranks above every anchor here)
    above <- sum(vr > top)
    ex <- data.table(v = v[blk])[, .N, by = v]
    res[[length(res) + 1]] <- data.table(
      arm = a, term = tm, col = col, n_anchors = length(v), modal = top, n_modal = n$N[1],
      n_exact = max(ex$N), n_values = nrow(ex),
      rank_lo = above + 1L, rank_hi = above + n$N[1], n_next = n$N[2],
      block_unembedded = round(mean(unemb[ok][blk]), 4), n_unembedded = sum(unemb[ok]))
  }
  rm(b); invisible(gc())
}
res <- rbindlist(res)
print(res, row.names = FALSE)
fwrite(res, file.path(HC, "tied_score_blocks.csv"))
cat("wrote", file.path(HC, "tied_score_blocks.csv"), "\n")
