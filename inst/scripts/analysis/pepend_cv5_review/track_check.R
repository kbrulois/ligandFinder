suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))
W <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
source(file.path(W, "R/pepend_tracks.R"))
A2 <- path.expand("~/AF2_analysis")
extra <- unlist(lapply(1:5, function(k) lapply(c("C","N"), function(tm)
  list(name = sprintf("%s fold %d", tm, k), term = tm, group = "per-fold", popup = TRUE,
       run = file.path(A2, sprintf("lf_pepend_run_%s_cv5_f%d_slack8.rds", tm, k)),
       winxgb = file.path(A2, sprintf("lf_winxgb_scan_%s_cv5_f%d.npz", tm, k))))), recursive = FALSE)
options(lf.pepend_python_path = file.path(W, "inst", "python"),
        lf.win_overlay_all50 = TRUE,
        lf.cv_knowns = list(N = file.path(A2, "lf_pepend_cv5_knowns_N.csv"),
                            C = file.path(A2, "lf_pepend_cv5_knowns_C.csv")),
        lf.track_scan_C = file.path(A2, "lf_pepend_scan_C_cv5all50.rds"),
        lf.track_scan_N = file.path(A2, "lf_pepend_scan_N_cv5all50.rds"),
        lf.track_winxgb_C = file.path(A2, "lf_winxgb_scan_C_cv5all50.npz"),
        lf.track_winxgb_N = file.path(A2, "lf_winxgb_scan_N_cv5all50.npz"),
        lf.track_extra_win = extra)
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
GENE <- Sys.getenv("LF_GENE", "NPY")
i <- match(GENE, fc$prec$gene)
acc <- fc$prec$accession[i]; np <- fc$prec$n_prot[i]; cp <- fc$prec$c_prot[i]
aa <- strsplit(fc$prec$seq[i], "")[[1]]
dat <- lf_pepend_track_data(GENE, acc, np, cp, want_per_index = TRUE)
cat("gene:", GENE, "| held-out marks:", if (is.null(dat$heldout)) 0 else nrow(dat$heldout), "\n")
cat("win series:", paste(unique(dat$win$term), collapse=", "), "\n")
cat("groups:", paste(unique(dat$win$group), collapse=", "), " (NA = the page's own all-50)\n")
own <- dat$win[is.na(dat$win$group), , drop = FALSE]
xr <- c(1, length(aa))
p <- lf_plot_win_track(dat$win, "xgb window head", xr, aa, GENE, legend = TRUE,
                       group = "per-fold", ylab = "xgb window head\n(per-fold)", overlay = own)
ggsave(file.path(A2, "track_check_xgb.png"), p, width = 13, height = 2.4, dpi = 165)
p2 <- lf_plot_win_track(dat$win, "CNN attention head", xr, aa, GENE, legend = TRUE,
                        group = "per-fold", ylab = "CNN attention head\n(per-fold)", overlay = own)
ggsave(file.path(A2, "track_check_cnn.png"), p2, width = 13, height = 2.4, dpi = 165)
cat("wrote track_check_{xgb,cnn}.png\n")
