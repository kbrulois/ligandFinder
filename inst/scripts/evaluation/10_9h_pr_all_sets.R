#!/usr/bin/env Rscript
## ---- PR AUC and PR curves for every set a number has been quoted on -----------
## PR AUC's floor IS the prevalence, so a PR AUC means nothing without it. These
## five sets differ in prevalence by more than 2,000x, which is why the headline
## numbers elsewhere range from 0.17 to 1.00 for the same models. The table adds
## lift (AUC / prevalence) so the sets can actually be compared.
##
##   set                      n           positives  prevalence
##   pepend val split             429            29     6.76%
##   db val split                 407             7     1.72%
##   pepend all (candidates)   41,870            86     0.205%
##   db all (candidates)       29,599            22     0.0743%
##   full step-1 scan       2,980,852            86     0.0029%
##
## NOT EVERY SET CAN CARRY BOTH MODELS. The db arm predates per-member weight
## saving, so only its member 0 survives and it cannot be run on windows it did
## not already score. It therefore appears only on the two sets it scored during
## its own run. The peptide-end model has all 20 weight files and appears
## everywhere. Where a set carries one model, the table says so rather than
## implying a comparison.
##
## The full scan carries NO seed band: lf_dcnn.scan accumulates per-window
## summaries and never stores the 20 x 3M member matrix, so only the ensemble
## curve exists there. Getting a band would mean re-scanning (~20 min) to retain
## every member.
##
## Curves are exact: with P positives there are only P achievable recall levels
## k/P, and precision there is k over the rank of the k-th positive, so all seeds
## share a recall grid with no interpolation.
##
##   Rscript inst/scripts/10_9h_pr_all_sets.R
##
## Output: ~/AF2_analysis/lf_pr_all_sets_<term>.{csv,svg,png}
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(ggplot2); library(yardstick); library(tidyr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term   <- toupper(.opt("--term", "C"))
arm    <- path.expand(.opt("--arm", "~/AF2_analysis/lf_dcnn_bench_none20/unet_none_w1"))
in_dir <- path.expand(.opt("--in-dir", file.path(dirname(arm), "in")))
run_p  <- path.expand(.opt("--pepend-run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
pin    <- path.expand(.opt("--pepend-input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
scan_p <- path.expand(.opt("--scan", sprintf("~/AF2_analysis/lf_pepend_scan_%s.rds", term)))
cache  <- path.expand(.opt("--cache", sprintf("~/AF2_analysis/lf_pr_all_sets_scores_%s.rds", term)))
out_st <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pr_all_sets_%s", term)))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))

DB <- "db-anchored"; PE <- "peptide-end"

## ---- 1. assemble every set: labels, and a score matrix per available model ----
if (file.exists(cache) && !("--refresh" %in% .args)) {
  sets <- readRDS(cache); message("reusing ", basename(cache), " (--refresh to rebuild)")
} else {
  mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
  np  <- reticulate::import("numpy", convert = FALSE)
  npc <- reticulate::import("numpy", convert = TRUE)
  z   <- npc$load(file.path(in_dir, "arrays.npz"))
  dbm <- sort(Sys.glob(file.path(arm, "members", sprintf("member_*_%s.npz", term))))
  pem <- sort(Sys.glob(file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"),
                                 "out", "members", sprintf("member_*_%s.npz", term))))
  out_dir <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out")
  cfg <- mod$Config$from_json(file.path(dirname(out_dir), "in", "config.json"))
  models <- mod$scan$load_member_models(out_dir, cfg, term)
  run_pe <- function(X) {                      # peptide-end ensemble over windows
    xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))
    vapply(models, function(m)
      as.numeric(reticulate::py_to_r(mod$pipeline$predict_all(m, xp)[["global"]])),
      numeric(dim(X)[[1]]))
  }
  read_mem <- function(fs, key, n) vapply(fs, function(f) as.numeric(npc$load(f)[[key]]), numeric(n))

  x <- readRDS(pin)
  ## pepend val / all: the peptide-end model's own stored scores. The db model
  ## cannot be run here -- no ensemble weights -- so these sets carry one model.
  yv <- as.numeric(x$nn_input[[term]]$val$known)
  ya <- as.numeric(x$nn_input[[term]]$all$known)
  ## db val / all: db's stored scores, plus the peptide-end ensemble run here
  ydv <- as.numeric(z[[paste0(term, "__val__y_global")]])
  yda <- as.numeric(z[[paste0(term, "__all__y_global")]])
  ## full scan: ensemble mean only (no per-member matrix is retained)
  sc <- readRDS(scan_p)
  kn <- sc$knowns %>% filter(!is.na(rank_all)) %>% select(accession, anchor)
  ys <- as.integer(paste(sc$scan$accession, sc$scan$anchor) %in%
                     paste(kn$accession, kn$anchor))

  sets <- list(
    list(name = "pepend val split",        y = yv,
         s = list(read_mem(pem, paste0(term, "__val__x"), length(yv))) %>% setNames(PE)),
    list(name = "db val split",            y = ydv,
         s = setNames(list(read_mem(dbm, paste0(term, "__val__x"), length(ydv)),
                           run_pe(z[[paste0(term, "__val__x")]])), c(DB, PE))),
    list(name = "pepend all (candidates)", y = ya,
         s = list(read_mem(pem, paste0(term, "__global__x"), length(ya))) %>% setNames(PE)),
    list(name = "db all (candidates)",     y = yda,
         s = setNames(list(read_mem(dbm, paste0(term, "__global__x"), length(yda)),
                           run_pe(z[[paste0(term, "__all__x")]])), c(DB, PE))),
    list(name = "full step-1 scan",        y = ys,
         s = list(matrix(sc$scan$score, ncol = 1)) %>% setNames(PE)))
  saveRDS(sets, cache); message("scored -> ", basename(cache))
}

## ---- 2. exact curve + AUC ------------------------------------------------------
curve_one <- function(sc, y, P) {
  r <- which(y[order(-sc, seq_along(sc))] == 1)
  tibble(k = seq_along(r), recall = k / P, precision = k / r)
}
rows <- list(); curves <- list()
for (st in sets) {
  y <- st$y; P <- sum(y == 1); n <- length(y); prev <- P / n
  truth <- factor(as.integer(y), levels = c(0, 1))
  for (mdl in names(st$s)) {
    S <- st$s[[mdl]]; ns <- ncol(S)
    a <- apply(S, 2, function(v) pr_auc_vec(truth, v, event_level = "second"))
    cv <- bind_rows(lapply(seq_len(ns), function(i)
      curve_one(S[, i], y, P) %>% mutate(seed = i - 1L)))
    curves[[length(curves) + 1]] <- cv %>%
      mutate(set = st$name, model = mdl, baseline = prev)
    rows[[length(rows) + 1]] <- tibble(
      set = st$name, n = n, positives = P, prevalence = prev, model = mdl,
      seeds = ns, auc_mean = mean(a), auc_sd = if (ns > 1) sd(a) else NA_real_,
      lift = mean(a) / prev)
  }
}
tab <- bind_rows(rows)
cur <- bind_rows(curves)
ORD <- c("pepend val split", "db val split", "pepend all (candidates)",
         "db all (candidates)", "full step-1 scan")
tab$set <- factor(tab$set, levels = ORD); cur$set <- factor(cur$set, levels = ORD)
tab <- tab %>% arrange(set, model)

message("\n== PR AUC by set and model ==")
tab %>%
  transmute(set, n = format(n, big.mark = ","), positives,
            prevalence = sprintf("%.4f%%", 100 * prevalence), model, seeds,
            `PR AUC` = ifelse(is.na(auc_sd), sprintf("%.3f", auc_mean),
                              sprintf("%.3f +/- %.3f", auc_mean, auc_sd)),
            lift = sprintf("%.0fx", lift)) %>%
  as.data.frame() %>% print(row.names = FALSE)
readr::write_csv(tab, paste0(out_st, ".csv"))

## ---- 3. one curve panel per set ------------------------------------------------
band <- cur %>% group_by(set, model, baseline, k, recall) %>%
  summarise(mean = mean(precision), sd = sd(precision), .groups = "drop") %>%
  mutate(sd = ifelse(is.na(sd), 0, sd),
         lo = pmax(mean - sd, baseline), hi = mean + sd)
## The AUC values go in the STRIP, not inside the panel: the saturated sets pin
## their curve to the top of the panel, which is exactly where in-panel text sits.
auc_lines <- tab %>%
  mutate(txt = ifelse(is.na(auc_sd), sprintf("%s  %.3f  (1 seed)", model, auc_mean),
                      sprintf("%s  %.3f +/- %.3f", model, auc_mean, auc_sd))) %>%
  group_by(set) %>% summarise(auc_txt = paste(txt, collapse = "\n"), .groups = "drop")
strip <- tab %>% distinct(set, n, positives, prevalence) %>%
  left_join(auc_lines, by = "set") %>%
  mutate(lab = sprintf("%s\n%s windows, %d known (%.3f%%)\n%s", set,
                       format(n, big.mark = ","), positives, 100 * prevalence, auc_txt))
lev <- strip$lab[match(ORD, strip$set)]
band$facet <- factor(strip$lab[match(band$set, strip$set)], levels = lev)
base_df <- strip %>% mutate(facet = factor(lab, levels = lev))

FILL <- setNames(c(LF_VIZ$s1, LF_VIZ$s2), c(DB, PE))
p <- ggplot(band, aes(recall, mean, colour = model, fill = model)) +
  geom_hline(data = base_df, aes(yintercept = prevalence), inherit.aes = FALSE,
             colour = LF_VIZ$ink2, linewidth = 0.3, linetype = "22") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ facet, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = FILL, name = NULL) +
  scale_fill_manual(values = FILL, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_log10() +
  labs(title = sprintf("Precision-recall across every evaluated set (%s terminus)", term),
       subtitle = paste("Dashed line is each set's random baseline, i.e. its prevalence --",
                        "PR AUC's floor. Band is +/- 1 sd of precision across seeds.\nY is log10",
                        "and FREE per panel: the prevalences span 2,000x, so a shared axis would",
                        "flatten every curve."),
       x = "recall", y = "precision (log, free scale)",
       caption = paste("The db-anchored model appears only where it already scored the windows:",
                       "it predates per-member weight saving, so only member 0 survives and it",
                       "cannot be\nrun on new windows. The full scan has no band -- the scan keeps",
                       "per-window summaries, not the 20 x 3M member matrix.")) +
  lf_viz_theme() +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(colour = LF_VIZ$ink2),
        strip.text = element_text(colour = LF_VIZ$ink, size = 8, lineheight = 1.35),
        panel.spacing.x = unit(11, "pt"))

ggsave(paste0(out_st, ".svg"), p, width = 13.5, height = 4.9, device = svglite::svglite)
ggsave(paste0(out_st, ".png"), p, width = 13.5, height = 4.9, dpi = 200, bg = LF_VIZ$surface)
message("\nplot  -> ", out_st, ".{svg,png}")
message("table -> ", out_st, ".csv")
