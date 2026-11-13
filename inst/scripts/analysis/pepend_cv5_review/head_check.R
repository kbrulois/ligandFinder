suppressMessages({library(tidyverse)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
source(file.path(W, "R/pepend_tracks.R")); A2 <- path.expand("~/AF2_analysis")
extra <- unlist(lapply(1:5, function(k) lapply(c("C","N"), function(tm)
  list(name = sprintf("%s fold %d", tm, k), term = tm, group = "per-fold", popup = FALSE,
       run = file.path(A2, sprintf("lf_pepend_run_%s_cv5_f%d_slack8.rds", tm, k)),
       winxgb = file.path(A2, sprintf("lf_winxgb_scan_%s_cv5_f%d.npz", tm, k))))), recursive = FALSE)
options(lf.pepend_python_path = file.path(W,"inst","python"),
        lf.track_scan_C = file.path(A2,"lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N = file.path(A2,"lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C = file.path(A2,"lf_winxgb_scan_C_cv5all50.npz"),
        lf.track_winxgb_N = file.path(A2,"lf_winxgb_scan_N_cv5all50.npz"),
        lf.track_extra_win = extra)
fc <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); i <- match("CXCL14", fc$prec$gene)
dat <- lf_pepend_track_data("CXCL14", fc$prec$accession[i], fc$prec$n_prot[i], fc$prec$c_prot[i], FALSE)
w <- dat$win; w$group[is.na(w$group)] <- "(page's own / all 50)"
print(as.data.frame(table(head = w$head, group = w$group)))
cat("\nper-fold series by head:\n")
print(w %>% filter(group == "per-fold") %>% count(head, term) %>% tidyr::pivot_wider(names_from=head, values_from=n) %>% as.data.frame())
