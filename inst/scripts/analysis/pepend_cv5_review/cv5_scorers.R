## Which scoring rule should ship: the mean of all 50 models, or the 10 seeds of
## one chosen fold?
##
## Protocol. The fold is CHOSEN on the knowns (each fold judged on the ends IT
## held out -- its own honest CV set), then both rules are EVALUATED on the 37
## never-seen peptides, which no fold ever trained on. Selection and evaluation
## therefore use disjoint peptides.
##
## For the knowns the all-50 mean is partly IN-SAMPLE (4 of 5 folds trained on
## each end), so knowns are reported only as the honest per-fold CV, never as a
## comparison between rules.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis"); L <- file.path(A2,"pepend_locus")
TOL <- 2L; PW <- 5L
source("R/pepend_windows.R")
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
ll <- fread(file.path(L,"ligand_list.csv")); up <- fread(file.path(L,"uniprot_peptides.csv")) %>%
  filter(type %in% c("peptide","propeptide"))
fc <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); seqs <- setNames(fc$prec$seq, fc$prec$accession)

res <- list(); peaks <- list(); foldsel <- list()
for (tm in c("C","N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  setorder(b, accession, anchor)
  b[, motif := { s <- seqs[[.BY$accession]]
    if (is.null(s)) rep(NA_character_,.N) else lf_pepend_motif_vec(s, tm, anchor, 1L) }, by=accession]
  ## held-out ends, per fold, from the CV knowns table
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  ## every end to hold OUT of the negative pool
  .tst <- as.data.table(tst); .pa <- if (tm=="N") .tst$start else .tst$end
  allends <- unique(rbindlist(list(
    kn[, list(accession, anchor)],
    data.table(accession = .tst$accession, anchor = as.integer(.pa)),
    data.table(accession = ll$accession, anchor = as.integer(if (tm=="N") ll$start else ll$end)),
    data.table(accession = up$accession, anchor = as.integer(if (tm=="N") up$start else up$end))),
    use.names=TRUE))[!is.na(anchor)]
  near <- CJ(d=-TOL:TOL, i=seq_len(nrow(allends)))[, .(k=paste(allends$accession[i], allends$anchor[i]+d))]
  ## A third rule: at each anchor take the MAX of the 5 fold means rather than
  ## their average. Any fold that fires high carries the anchor -- the sensitive
  ## end of the dial. Must be added BEFORE negpool is copied out of `b`, or the
  ## negative pool has no such column, nv comes back NULL and every rank is 1.
  b[, xgb_pm__maxfold := pmax(xgb_pm__f1, xgb_pm__f2, xgb_pm__f3, xgb_pm__f4, xgb_pm__f5)]
  b[, xgb_pm__minfold := pmin(xgb_pm__f1, xgb_pm__f2, xgb_pm__f3, xgb_pm__f4, xgb_pm__f5)]
  b[, key := paste(accession, anchor)]
  negpool <- b[!(key %in% near$k)]

  scorers <- c(all50 = "xgb_pm__all50", maxfold = "xgb_pm__maxfold",
               minfold = "xgb_pm__minfold",
               setNames(sprintf("xgb_pm__f%d",1:5), paste0("fold",1:5)))
  ## --- 1. choose the fold on the KNOWNS: each fold judged on its own val ends
  for (k in 1:5) {
    v <- kn[fold == k]
    if (!nrow(v)) next
    nv <- negpool[[sprintf("xgb_pm__f%d", k)]]
    sc <- b[match(paste(v$accession, v$anchor), key)][[sprintf("xgb_pm__f%d", k)]]
    rk <- sapply(sc, function(x) if (is.na(x)) NA_integer_ else sum(nv > x) + 1L)
    foldsel[[length(foldsel)+1]] <- tibble(term=tm, fold=k, n=sum(!is.na(rk)),
      median_rank=median(rk,na.rm=TRUE), top1k=mean(rk<=1000,na.rm=TRUE))
  }
  ## --- 2. evaluate every scorer on the 37 NEVER-SEEN peptides
  pos <- data.table(accession = .tst$accession, anchor = as.integer(.pa))
  pos <- merge(pos, b[, .(accession, anchor, key, motif)], by=c("accession","anchor"))
  for (nm in names(scorers)) {
    col <- scorers[[nm]]; nv <- negpool[[col]]; nn <- length(nv)
    stopifnot(!is.null(nv), nn > 1e6)   # an absent column would make every rank 1
    sc <- b[match(pos$key, key)][[col]]
    rk <- sapply(sc, function(x) if (is.na(x)) NA_integer_ else sum(nv > x) + 1L)
    ok <- !is.na(rk)
    y <- c(rep(TRUE,sum(ok)), rep(FALSE,nn)); s2 <- c(sc[ok], nv)
    o <- order(-s2); y <- y[o]; tp <- cumsum(y); fp <- cumsum(!y)
    prec_ <- tp/(tp+fp); rec <- tp/sum(ok); fpr <- fp/nn
    res[[length(res)+1]] <- tibble(term=tm, scorer=nm, group="all", n=sum(ok),
      roc=sum(diff(c(0,fpr))*(head(c(0,rec),-1)+c(0,rec)[-1])/2), ap=sum(prec_[y])/sum(ok),
      median_rank=median(rk,na.rm=TRUE), pct_rank=median(rk,na.rm=TRUE)/nn,
      top1k=mean(rk<=1000,na.rm=TRUE), top10k=mean(rk<=10000,na.rm=TRUE))
    for (g in unique(na.omit(pos$motif))) {
      i <- which(pos$motif==g & ok); if (!length(i)) next
      ng <- negpool[motif==g][[col]]
      rkg <- sapply(sc[i], function(x) sum(ng > x)+1L)
      res[[length(res)+1]] <- tibble(term=tm, scorer=nm, group=g, n=length(i),
        roc=mean((length(ng)-(rkg-1))/length(ng)), ap=NA_real_,
        median_rank=median(rkg), pct_rank=median(rkg)/length(ng),
        top1k=mean(rkg<=1000), top10k=mean(rkg<=10000))
    }
    ## --- 3. distance to the closest local peak (+/-5), same scorer
    sh <- c(lapply(1:PW,\(j) list(j,"lag")), lapply(1:PW,\(j) list(j,"lead")))
    b[, rmax := do.call(pmax, c(list(get(col)), lapply(sh, \(z) shift(get(col), z[[1]], fill=-Inf, type=z[[2]])))), by=accession]
    pk <- split(b[get(col) >= rmax]$anchor, b[get(col) >= rmax]$accession)
    off <- mapply(function(a,an){ p <- pk[[a]]; if (is.null(p)) NA_integer_ else p[which.min(abs(p-an))]-an },
                  pos$accession, pos$anchor)
    peaks[[length(peaks)+1]] <- tibble(term=tm, scorer=nm, n=sum(!is.na(off)),
      exact=mean(off==0,na.rm=TRUE), w2=mean(abs(off)<=2,na.rm=TRUE),
      w5=mean(abs(off)<=5,na.rm=TRUE), median_abs=median(abs(off),na.rm=TRUE))
    b[, rmax := NULL]
  }
}
F <- bind_rows(foldsel); R <- bind_rows(res); P <- bind_rows(peaks)
saveRDS(list(F=F,R=R,P=P), file.path(S2,"cv5_scorers.rds"))
cat("=== fold choice, made on the KNOWNS (each fold on the ends it held out) ===\n")
print(F %>% mutate(top1k=sprintf("%.0f%%",100*top1k)) %>% as.data.frame(), row.names=FALSE)
cat("\n=== scorers on the 37 NEVER-SEEN peptides, all anchors ===\n")
print(R %>% filter(group=="all") %>% transmute(term, scorer, n, roc=round(roc,4), ap=round(ap,4),
  median_rank, pct_rank=sprintf("%.2f%%",100*pct_rank), top1k=sprintf("%.0f%%",100*top1k),
  top10k=sprintf("%.0f%%",100*top10k)) %>% as.data.frame(), row.names=FALSE)
cat("\n=== distance to the closest local peak (never-seen peptides) ===\n")
print(P %>% transmute(term, scorer, n, exact=sprintf("%.0f%%",100*exact), within2=sprintf("%.0f%%",100*w2),
  within5=sprintf("%.0f%%",100*w5), median_abs) %>% as.data.frame(), row.names=FALSE)
