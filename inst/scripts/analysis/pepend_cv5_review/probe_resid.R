suppressMessages({library(tidyverse)}); suppressWarnings(suppressMessages(library(ligandFinder)))
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"; A2 <- path.expand("~/AF2_analysis")
source(file.path(W,"R/plot_proteins_win_new.R")); source(file.path(W,"R/pepend_tracks.R"))
options(lf.pepend_python_path=file.path(W,"inst","python"),
        lf.track_run_C=file.path(A2,"lf_pepend_run_C_cv5all50_slack8.rds"),
        lf.track_run_N=file.path(A2,"lf_pepend_run_N_cv5all50_slack8.rds"),
        lf.track_scan_C=file.path(A2,"lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N=file.path(A2,"lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C=file.path(A2,"lf_winxgb_scan_C_cv5all50.npz"),
        lf.track_winxgb_N=file.path(A2,"lf_winxgb_scan_N_cv5all50.npz"),
        lf.track_resid=file.path(A2,"lf_resid_preds_v2.npz"))
d <- lf_pepend_track_data(gene="CXCL14", accession="O95715", n_prot=1L, c_prot=111L)
cat("elements:", paste(names(d), collapse=", "), "\n")
if (!is.null(d$resid)) { cat("\nresid rows by term x class:\n"); print(table(d$resid$term, d$resid$class)) }
if (!is.null(d$win)) { cat("\nwin rows by term x head:\n"); print(table(d$win$term, d$win$head)) }
if (!is.null(d$pix)) cat("\npix: ", paste(names(d$pix), collapse=", "), "\n")

options(lf.class_tracks_by_term = TRUE)
pl <- lf_pepend_track_plots(d, x_range = c(0, 112), aa = strsplit(paste(rep("A",111),collapse=""),"")[[1]],
                            gene = "CXCL14")
cat("\ntracks returned:", length(pl), "\n")
for (i in seq_along(pl)) {
  p <- pl[[i]]$p
  yl <- tryCatch(p$labels$y, error=function(e) NA)
  nr <- tryCatch(nrow(p$data), error=function(e) NA)
  tip <- tryCatch({ l <- p$layers[[3]]; as.character(rlang::quo_get_expr(l$mapping$tooltip))[1] }, error=function(e) NA)
  cat(sprintf("  %d: ylab=%-42s rows=%-5s\n", i, gsub("\n"," | ", yl %||% "NULL"), nr))
}
