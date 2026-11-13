suppressMessages({library(tidyverse)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
source(file.path(W,"R/pepend_tracks.R")); A2 <- path.expand("~/AF2_analysis")
extra <- unlist(lapply(1:5, function(k) lapply(c("C","N"), function(tm)
  list(name=sprintf("%s fold %d",tm,k), term=tm, group="per-fold", popup=TRUE,
       run=file.path(A2,sprintf("lf_pepend_run_%s_cv5_f%d_slack8.rds",tm,k)),
       winxgb=file.path(A2,sprintf("lf_winxgb_scan_%s_cv5_f%d.npz",tm,k))))), recursive=FALSE)
options(lf.pepend_python_path=file.path(W,"inst","python"),
        lf.cv_knowns=list(N=file.path(A2,"lf_pepend_cv5_knowns_N.csv"),
                          C=file.path(A2,"lf_pepend_cv5_knowns_C.csv")),
        lf.track_scan_C=file.path(A2,"lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N=file.path(A2,"lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C=file.path(A2,"lf_winxgb_scan_C_cv5all50.npz"),
        lf.track_winxgb_N=file.path(A2,"lf_winxgb_scan_N_cv5all50.npz"),
        lf.track_extra_win=extra)
fc <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); i <- match("NPY", fc$prec$gene)
dat <- lf_pepend_track_data("NPY", fc$prec$accession[i], fc$prec$n_prot[i], fc$prec$c_prot[i], TRUE)
cat("\n== rows by head x group ==\n")
w <- dat$win; w$group[is.na(w$group)] <- "(all 50)"
print(as.data.frame(table(head=w$head, group=w$group)))
cat("\n== held-out marks ==\n"); print(dat$heldout)
cat("\n== do the marked points exist in the per-fold data? ==\n")
pf <- w[w$group=="per-fold", ]
for (k in seq_len(nrow(dat$heldout))) {
  h <- dat$heldout[k, ]
  hit <- pf[pf$term==h$term & pf$residue==h$residue, c("head","term","residue","mean")]
  cat(sprintf("  %s @ %d -> %d row(s): %s\n", h$term, h$residue, nrow(hit),
      paste(sprintf("%s=%.3f", hit$head, hit$mean), collapse=", ")))
}
