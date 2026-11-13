#!/usr/bin/env Rscript
## ---- the insertion head's three classes, one panel each, mean +/- sd ----------
## 10_9i draws all three classes on one axis, which shows the crossover but hides
## the spread -- and the stored scan only carries `ins_sd_inserting`, so the
## other two sds are not in it at all. This re-runs the ensemble over one gene's
## windows (cheap: a 97-residue precursor is 69 windows) and keeps every member,
## so each class gets its own panel with the across-seed mean and sd.
##
## The window score is drawn on top for comparison, because the two are NOT the
## same curve and it is easy to assume they are: on NPY r = 0.62 and they differ
## by up to 0.51. They answer different questions. The global head asks "is this
## a real peptide end?". The insertion head asks "GIVEN a peptide ending here,
## where does the pocket sit?" -- a conditional, trained on the 57 labelled
## windows only, whose softmax must still sum to 1 at an anchor the global head
## rejects outright. Read the class panels where the score panel is high.
##
##   Rscript inst/scripts/10_9o_pepend_ins_class_tracks.R --gene NPY
##
## Output: <out-dir>/lf_ins_class_tracks_<gene>_<term>.svg (+ .png unless
##         --no-png) and the matching .csv
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
gene    <- .opt("--gene", "NPY")
term    <- toupper(.opt("--term", "C"))
run_p   <- path.expand(.opt("--run", "~/AF2_analysis/lf_pepend_run_C_ins.rds"))
cache_p <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
scan_p  <- path.expand(.opt("--scan", "~/AF2_analysis/lf_pepend_scan_C_ins.rds"))
ll_p    <- .opt("--ligand-list", "inst/extdata/ligand_list.rds")
out_dir <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
no_png  <- "--no-png" %in% .args

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)

ANCHOR  <- c(N = 8L, C = 28L)[[term]]
SEQ_LEN <- 36L
RES <- c(`dibasic pair (KK/KR/RK/RR)` = LF_VIZ$s2,
         `single K / R`               = LF_VIZ$s3,
         other                        = LF_VIZ$grid)
INS_PAL <- c(inserting = "#1B7837", loop = "#E08214", non_inserting = "#878787")

## ---- 1. the windows ----------------------------------------------------------
fc <- readRDS(cache_p)
p  <- fc$prec %>% filter(gene == !!gene) %>% slice(1)
if (!nrow(p)) stop("no precursor called ", gene, call. = FALSE)
feat <- fc$feats[[match(p$accession, fc$prec$accession)]]
L <- nchar(p$seq); aa <- strsplit(p$seq, "")[[1]]
anchors <- seq.int(p$n_prot, p$c_prot)
message(sprintf("%s (%s): %d residues, mature %d-%d, %d windows",
                gene, p$accession, L, p$n_prot, p$c_prot, length(anchors)))

X <- array(0, c(length(anchors), SEQ_LEN, ncol(feat)))
for (i in seq_along(anchors))
  X[i, , ] <- as.matrix(lf_pepend_slice(feat, anchors[i] - ANCHOR + 1L,
                                        p$n_prot, p$c_prot)[, fc$all_params3])
xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))

## ---- 2. every member, kept ---------------------------------------------------
cfg <- mod$Config$from_json(file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"),
                                      "in", "config.json"))
models <- mod$scan$load_member_models(
  file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out"), cfg, term)
ins_names <- tryCatch(as.character(cfg$ins_class_names), error = function(e) character(0))
if (!length(ins_names))
  stop(basename(run_p), " has no insertion head -- train with --ins-head", call. = FALSE)

G   <- matrix(0, length(models), length(anchors))
INS <- array(0, c(length(models), length(anchors), length(ins_names)))
for (m in seq_along(models)) {
  o <- mod$pipeline$predict_all(models[[m]], xp)
  G[m, ] <- as.numeric(reticulate::py_to_r(o[["global"]]))
  INS[m, , ] <- reticulate::py_to_r(o[["ins_class"]])
}
message(sprintf("scored %d window(s) x %d member(s)", length(anchors), length(models)))

tr <- tibble(anchor = anchors, mean = colMeans(G), sd = apply(G, 2, sd), class = "window score")
cls <- lapply(seq_along(ins_names), function(k)
  tibble(anchor = anchors, mean = colMeans(INS[, , k]), sd = apply(INS[, , k], 2, sd),
         class = ins_names[k])) %>% bind_rows()

## ---- 3. cross-check the means against the stored full scan -------------------
## This gene's windows are a subset of the 2,980,852 the scan covers, computed
## here by a different code path. The means must agree exactly; if they do not,
## the run and the scan are not the same model and nothing below is comparable.
sc <- tryCatch(readRDS(scan_p), error = function(e) NULL)
kn <- NULL
if (!is.null(sc)) {
  ref <- sc$scan %>% filter(accession == p$accession) %>% arrange(anchor)
  kn  <- sc$knowns %>% filter(accession == p$accession)
  if (nrow(ref) == length(anchors)) {
    d_g <- max(abs(ref$score - tr$mean))
    d_i <- max(vapply(ins_names, function(n) {
      col <- paste0("ins_p_", n)
      if (col %in% names(ref)) max(abs(ref[[col]] - cls$mean[cls$class == n])) else NA_real_
    }, numeric(1)), na.rm = TRUE)
    message(sprintf("cross-check vs %s: max|delta| global %.2e, insertion %.2e  [%s]",
                    basename(scan_p), d_g, d_i,
                    if (max(d_g, d_i) < 1e-5) "OK" else "MISMATCH"))
    if (max(d_g, d_i, na.rm = TRUE) >= 1e-5)
      warning("the run and the scan disagree -- are they the same model?", call. = FALSE)
  } else {
    message("cross-check skipped: scan has ", nrow(ref), " anchors here, run built ",
            length(anchors))
  }
}

## ---- 4. the figure -----------------------------------------------------------
LL <- tryCatch(readRDS(file.path(ROOT, ll_p)), error = function(e) NULL)
peps <- if (is.null(LL)) tibble(final_name = character(), start = integer(), end = integer()) else
  LL %>% filter(accession == p$accession) %>% select(final_name, start, end) %>%
  distinct() %>% arrange(start)

PLOT_W  <- max(11, min(42, 3.5 + L * 0.034))
per_res <- (PLOT_W - 2.2) * 25.4 / L
SEQ_SZ  <- max(0.85, min(1.9, per_res / 0.62))
SHOW_AA <- per_res / 0.62 >= 0.8
BRK_W   <- if (L > 800) 100 else if (L > 400) 50 else if (L > 200) 25 else 10

panel <- function(d, col, ylab) {
  ggplot(d, aes(anchor, mean)) +
    {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
    {if (!is.null(kn) && nrow(kn)) geom_vline(data = kn, inherit.aes = FALSE,
        aes(xintercept = anchor), colour = LF_VIZ$ink, linewidth = 0.4, linetype = "22") } +
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                fill = col, alpha = 0.22) +
    geom_line(colour = col, linewidth = 0.8) +
    {if (!is.null(kn) && nrow(kn)) geom_point(
        data = kn %>% left_join(d, by = "anchor"), inherit.aes = FALSE,
        aes(anchor, mean), colour = LF_VIZ$ink, size = 1.9) } +
    scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                       breaks = scales::breaks_width(BRK_W)) +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = ylab) +
    lf_viz_theme() +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
}

ttl <- ggplot() + theme_void() +
  labs(title = sprintf("%s (%s): window score and each insertion class, step-1 %s scan",
                       gene, p$accession, term),
       subtitle = sprintf(paste0(
         "One window per mature residue (%d of them), x = the window's anchor. Line is the ",
         "%d-seed mean, band +/- 1 sd across seeds.\nThe score and p(inserting) are NOT the same ",
         "curve (r = %.2f here): the global head asks whether this is a real peptide end, the ",
         "insertion\nhead asks where the pocket would sit GIVEN one ends here -- a conditional, ",
         "so read it where the score panel is high. Dashed: known end."),
         length(anchors), length(models), cor(tr$mean, cls$mean[cls$class == "inserting"]))) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = LF_VIZ$ink),
        plot.subtitle = element_text(size = 8, colour = LF_VIZ$ink2, lineheight = 1.2))

is_basic <- aa %in% c("K", "R"); in_pair <- rep(FALSE, L)
pr <- which(head(is_basic, -1) & tail(is_basic, -1)); in_pair[c(pr, pr + 1L)] <- TRUE
p_seq <- ggplot(tibble(pos = seq_len(L), aa = aa,
                       res = ifelse(in_pair, names(RES)[1],
                             ifelse(is_basic, names(RES)[2], names(RES)[3]))),
                aes(pos, 1)) +
  geom_tile(aes(fill = res), height = 0.42, colour = NA) +
  {if (SHOW_AA) geom_text(aes(label = aa), size = SEQ_SZ, colour = LF_VIZ$ink, vjust = 0.5) } +
  scale_fill_manual(values = RES, name = "residue", breaks = names(RES)) +
  scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                     breaks = scales::breaks_width(BRK_W)) +
  labs(x = "residue (the window's anchor = the peptide's last residue)", y = NULL) +
  lf_viz_theme() +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.justification = "left",
        legend.key.size = unit(9, "pt"),
        legend.text = element_text(colour = LF_VIZ$ink2, size = 8),
        legend.title = element_text(colour = LF_VIZ$ink2, size = 8))

panels <- c(list(panel(tr, LF_VIZ$s1, "window score")),
            lapply(ins_names, function(n)
              panel(cls %>% filter(class == n), INS_PAL[[n]], sprintf("p(%s)", n))))
fig <- wrap_plots(c(list(ttl), panels, list(p_seq)), ncol = 1) +
  plot_layout(heights = c(1.3, rep(3, length(panels)), 1.5)) +
  plot_annotation(theme = theme(plot.background =
                                  element_rect(fill = LF_VIZ$surface, colour = NA)))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
st <- file.path(out_dir, sprintf("lf_ins_class_tracks_%s_%s", gene, term))
H <- 1.5 + 1.7 * length(panels) + 1.3
ggsave(paste0(st, ".svg"), fig, width = PLOT_W, height = H, device = svglite::svglite)
if (!no_png) ggsave(paste0(st, ".png"), fig, width = PLOT_W, height = H, dpi = 200,
                    bg = LF_VIZ$surface)
readr::write_csv(bind_rows(tr, cls) %>% arrange(class, anchor), paste0(st, ".csv"))

message("\n", paste0(st, ".svg"), "\n", paste0(st, ".csv"))
bind_rows(tr, cls) %>% group_by(class) %>%
  summarise(peak = max(mean), at = anchor[which.max(mean)],
            mean_sd = mean(sd), max_sd = max(sd), max_sd_at = anchor[which.max(sd)],
            .groups = "drop") %>% as.data.frame() %>% print(row.names = FALSE, digits = 3)
