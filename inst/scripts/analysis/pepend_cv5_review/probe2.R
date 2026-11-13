suppressMessages({library(tidyverse)}); suppressWarnings(suppressMessages(library(ligandFinder)))
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"; A2 <- path.expand("~/AF2_analysis")
source(file.path(W,"R/plot_proteins_win_new.R")); source(file.path(W,"R/pepend_tracks.R"))
options(lf.pepend_python_path=file.path(W,"inst","python"),
        lf.track_run_C=file.path(A2,"lf_pepend_run_C_cv5all50_slack8.rds"),
        lf.track_run_N=file.path(A2,"lf_pepend_run_N_cv5all50_slack8.rds"),
        lf.track_resid=file.path(A2,"lf_resid_preds_v2.npz"))
d <- lf_pepend_track_data(gene="CXCL14", accession="O95715", n_prot=1L, c_prot=111L)
r <- d$resid
cat("resid: NA in `mean` by term:\n")
print(r %>% group_by(term) %>% summarise(n=n(), na_mean=sum(is.na(mean)), na_sd=sum(is.na(sd)),
      min=min(mean,na.rm=TRUE), max=max(mean,na.rm=TRUE), .groups="drop") %>% as.data.frame())
cat("\nresidue range covered by each term:\n")
print(r %>% group_by(term) %>% summarise(lo=min(residue), hi=max(residue), n_res=n_distinct(residue), .groups="drop") %>% as.data.frame())
cat("\nn_win by term (how many windows cover each residue):\n")
print(r %>% group_by(term) %>% summarise(min_nwin=min(n_win), max_nwin=max(n_win), .groups="drop") %>% as.data.frame())

options(lf.class_tracks_by_term = TRUE)
pl <- lf_pepend_track_plots(d, x_range=c(0,112), aa=strsplit(paste(rep("A",111),collapse=""),"")[[1]], gene="CXCL14")
for (i in 1:2) {
  g <- ggiraph::girafe(ggobj = pl[[i]]$p, width_svg = 10, height_svg = pl[[i]]$h)
  f <- file.path(Sys.getenv("S2OUT"), sprintf("probe_track%d.html", i))
  htmlwidgets::saveWidget(g, f, selfcontained = FALSE)
  h <- paste(readLines(f, warn=FALSE), collapse="")
  cat(sprintf("track %d (%s): tooltip hits N=%d C=%d | <circle> %d\n", i,
    gsub("\n"," | ", pl[[i]]$p$labels$y),
    lengths(regmatches(h, gregexpr("N terminus windows only", h))),
    lengths(regmatches(h, gregexpr("C terminus windows only", h))),
    lengths(regmatches(h, gregexpr("<circle", h)))))
}
