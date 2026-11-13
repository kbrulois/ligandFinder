suppressMessages({library(tidyverse); library(ggiraph)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
source(file.path(W,"R/pepend_tracks.R")); A2 <- path.expand("~/AF2_analysis")
extra <- unlist(lapply(1:5, function(k) lapply(c("C","N"), function(tm)
  list(name=sprintf("%s fold %d",tm,k), term=tm, group="per-fold", popup=TRUE,
       run=file.path(A2,sprintf("lf_pepend_run_%s_cv5_f%d_slack8.rds",tm,k)),
       winxgb=file.path(A2,sprintf("lf_winxgb_scan_%s_cv5_f%d.npz",tm,k))))), recursive=FALSE)
options(lf.pepend_python_path=file.path(W,"inst","python"), lf.win_overlay_all50=TRUE,
        lf.cv_knowns=list(N=file.path(A2,"lf_pepend_cv5_knowns_N.csv"),
                          C=file.path(A2,"lf_pepend_cv5_knowns_C.csv")),
        lf.track_scan_C=file.path(A2,"lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N=file.path(A2,"lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C=file.path(A2,"lf_winxgb_scan_C_cv5all50.npz"),
        lf.track_winxgb_N=file.path(A2,"lf_winxgb_scan_N_cv5all50.npz"),
        lf.track_extra_win=extra)
fc <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); i <- match("NPY", fc$prec$gene)
aa <- strsplit(fc$prec$seq[i],"")[[1]]
dat <- lf_pepend_track_data("NPY", fc$prec$accession[i], fc$prec$n_prot[i], fc$prec$c_prot[i], TRUE)
own <- dat$win[is.na(dat$win$group), , drop=FALSE]
mk <- function(lo, hi, tag, hd) {
  p <- lf_plot_win_track(dat$win, hd, c(lo,hi), aa, "NPY", legend=FALSE,
                         group="per-fold", ylab=paste0(hd,"\n(per-fold)"),
                         overlay=own, held=dat$heldout) +
       ggplot2::coord_cartesian(xlim=c(lo,hi))
  ggsave(file.path(A2, sprintf("zoom_%s.png", tag)), p, width=7, height=3, dpi=200)
}
mk(24, 38, "NPY_Nends_xgb", "xgb window head")
mk(56, 72, "NPY_Cend_xgb",  "xgb window head")
cat("marks:\n"); print(dat$heldout)
cat("wrote zoom_NPY_Nends_xgb.png, zoom_NPY_Cend_xgb.png\n")
