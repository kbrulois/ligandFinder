## PR curves and recall-vs-rank curves, within each motif class, for the three
## ensemble rules. Positives are the 37 never-seen peptides; negatives are the
## anchors of the SAME motif class, so recall reaches 100% inside every group
## and the groups are comparable despite very different prevalences.
suppressMessages({library(dplyr); library(data.table); library(tidyr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis"); L <- file.path(A2,"pepend_locus")
TOL <- 2L
source("R/pepend_windows.R")
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
ll <- fread(file.path(L,"ligand_list.csv")); up <- fread(file.path(L,"uniprot_peptides.csv")) %>%
  filter(type %in% c("peptide","propeptide"))
fc <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); seqs <- setNames(fc$prec$seq, fc$prec$accession)
GROUPS <- c("all","dibasic","monobasic","G + dibasic")
pr <- list(); rr <- list()
for (tm in c("C","N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  setorder(b, accession, anchor)
  b[, motif := { s <- seqs[[.BY$accession]]
    if (is.null(s)) rep(NA_character_,.N) else lf_pepend_motif_vec(s, tm, anchor, 1L) }, by=accession]
  b[, xgb_pm__maxfold := pmax(xgb_pm__f1,xgb_pm__f2,xgb_pm__f3,xgb_pm__f4,xgb_pm__f5)]
  b[, xgb_pm__minfold := pmin(xgb_pm__f1,xgb_pm__f2,xgb_pm__f3,xgb_pm__f4,xgb_pm__f5)]
  b[, key := paste(accession, anchor)]
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  .pa <- if (tm=="N") tst$start else tst$end
  allends <- unique(rbindlist(list(kn[, list(accession, anchor)],
    data.table(accession=tst$accession, anchor=as.integer(.pa)),
    data.table(accession=ll$accession, anchor=as.integer(if (tm=="N") ll$start else ll$end)),
    data.table(accession=up$accession, anchor=as.integer(if (tm=="N") up$start else up$end))),
    use.names=TRUE))[!is.na(anchor)]
  near <- CJ(d=-TOL:TOL, i=seq_len(nrow(allends)))[, .(k=paste(allends$accession[i], allends$anchor[i]+d))]
  b[, is_end := key %in% near$k]
  pos <- merge(data.table(accession=tst$accession, anchor=as.integer(.pa)),
               b[, .(accession, anchor, key, motif)], by=c("accession","anchor"))
  scor <- c(`mean of 50`="xgb_pm__all50", `max of folds`="xgb_pm__maxfold",
            `min of folds`="xgb_pm__minfold")
  for (g in GROUPS) {
    sub   <- if (g=="all") b else b[motif == g]
    pg    <- if (g=="all") pos else pos[motif == g]
    if (!nrow(pg)) next
    negs  <- sub[is_end == FALSE]
    for (nm in names(scor)) {
      col <- scor[[nm]]; nv <- negs[[col]]; nn <- length(nv)
      sc  <- b[match(pg$key, key)][[col]]
      ok  <- !is.na(sc); np <- sum(ok)
      s2 <- c(sc[ok], nv); y <- c(rep(TRUE,np), rep(FALSE,nn))
      o <- order(-s2); y <- y[o]
      tp <- cumsum(y); fp <- cumsum(!y)
      prec <- tp/(tp+fp); rec <- tp/np
      k <- unique(c(which(y), round(seq(1, length(y), length.out=400))))
      k <- sort(k[k <= length(y)])
      pr[[length(pr)+1]] <- tibble(term=tm, group=g, scorer=nm, n=np,
                                   recall=rec[k], precision=prec[k])
      rk <- sapply(sc[ok], function(v) sum(nv > v) + 1L)
      Ns <- unique(round(10^seq(0, log10(max(nn,10)), length.out=90)))
      rr[[length(rr)+1]] <- tibble(term=tm, group=g, scorer=nm, n=np,
                                   N=Ns, recall=sapply(Ns, function(N) mean(rk <= N)))
    }
  }
}
saveRDS(list(pr=bind_rows(pr), rr=bind_rows(rr)), file.path(S2,"cv5_curves.rds"))
cat("curve data written\n")
print(bind_rows(rr) %>% distinct(term, group, n) %>% as.data.frame(), row.names=FALSE)
