#!/usr/bin/env Rscript
## ---- protein page v2: both termini, every window a click away ------------------
## 10_4_plot_proteins.R draws one protein with the production window model's
## detail panels baked in as rectangles, and annotates peptides with the GPCRdb
## line segments. This replaces both ideas.
##
##   * BOTH TERMINI, as separate tracks. An N window and a C window answer
##     different questions about the same residue ("does a peptide start here",
##     "does one end here") and their anchors sit at different window positions
##     (8 vs 28), so they can never share a track.
##   * NO WINDOW RECTANGLES. Every mature residue is a window anchor, and all of
##     their per-index surfaces travel with the page as data. Clicking a residue
##     on a class-probability track pops up the 36-residue window anchored
##     THERE, drawn on demand -- so the page carries every window instead of the
##     handful someone chose in advance.
##   * PEPTIDE ANNOTATION, not GPCRdb segments: the known peptides' own spans,
##     overlapping ones on separate lines (`assign_overlap_layers`), the way the
##     resid-gallery figures draw ground truth.
##   * LINKED HOVER. Every track keys its points by residue, so hovering one
##     track marks the same residue on all of them.
##
##   Rscript inst/scripts/10_9t_protein_page_v2.R --gene VIP
##   ... --run-c ~/AF2_analysis/lf_pepend_run_C_near.rds --out page.html
##
## Output: one self-contained HTML.
## ------------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
gene_want <- .opt("--gene", "VIP")
## The near-end trunks by default: this page is for reading a boundary off a
## precursor, and those are the trunks that localise one (on-end 60% under the
## xgb head, against 28% on the plain trunks). They are NOT the better trunks
## for finding which protein carries a peptide -- see the near-end memory.
run_p  <- c(N = path.expand(.opt("--run-n",  "~/AF2_analysis/lf_pepend_run_N_near.rds")),
            C = path.expand(.opt("--run-c",  "~/AF2_analysis/lf_pepend_run_C_near.rds")))
scan_p <- c(N = path.expand(.opt("--scan-n", "~/AF2_analysis/lf_pepend_scan_N_near.rds")),
            C = path.expand(.opt("--scan-c", "~/AF2_analysis/lf_pepend_scan_C_near.rds")))
wx_p   <- c(N = path.expand(.opt("--winxgb-n", "~/AF2_analysis/lf_winxgb_scan_N_near.npz")),
            C = path.expand(.opt("--winxgb-c", "~/AF2_analysis/lf_winxgb_scan_C_near.npz")))
cache_p <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
rp_path <- path.expand(.opt("--resid-preds", "~/AF2_analysis/lf_resid_preds.npz"))
out_p   <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_protein_v2_%s.html", gene_want)))
## Popups for every anchor is ~36*K floats per residue. Fine for a 170-residue
## precursor; for a 1600-residue one it is megabytes, so a floor can drop the
## anchors nothing will ever be clicked on.
min_pop <- as.numeric(.opt("--popup-min-score", "0"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))
## nn_class_cols / nn_class_labels / assign_overlap_layers
source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)

SEQ_LEN <- LF_PEPEND$seq_len
DRAW_CLASSES <- c("pep_pocket", "pep_other", "CT_cleavage_context", "NT_cleavage_context")

fc <- readRDS(cache_p)
p  <- fc$prec %>% filter(gene == gene_want) %>% slice(1)
if (!nrow(p)) stop("no precursor called ", gene_want, call. = FALSE)
feat <- fc$feats[[match(p$accession, fc$prec$accession)]]
L  <- nchar(p$seq); aa <- strsplit(p$seq, "")[[1]]
anchors <- seq.int(p$n_prot, p$c_prot)
message(sprintf("%s (%s): %d residues, mature %d-%d, %d anchors per terminus",
                gene_want, p$accession, L, p$n_prot, p$c_prot, length(anchors)))

## ---- 1. score every anchor, both termini ---------------------------------------
## One forward pass per member per terminus carries BOTH the per-index surface
## (the popups, and the pooled per-residue track) and the global head.
score_term <- function(term) {
  iso <- paste0(tools::file_path_sans_ext(run_p[[term]]), "_isolated")
  cfg <- mod$Config$from_json(file.path(iso, "in", "config.json"))
  models <- mod$scan$load_member_models(file.path(iso, "out"), cfg, term)
  pi_names <- as.character(cfg$pi_names)
  A <- LF_PEPEND$anchor[[term]]
  w_start <- anchors - A + 1L

  X <- array(0, c(length(anchors), SEQ_LEN, ncol(feat)))
  for (i in seq_along(anchors))
    X[i, , ] <- as.matrix(lf_pepend_slice(feat, w_start[i], p$n_prot, p$c_prot)[, fc$all_params3])
  xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))

  PI <- array(0, c(length(models), length(anchors), SEQ_LEN, length(pi_names)))
  G  <- matrix(0, length(models), length(anchors))
  for (m in seq_along(models)) {
    o <- mod$pipeline$predict_all(models[[m]], xp)
    PI[m, , , ] <- reticulate::py_to_r(o[["per_index_cat"]])
    G[m, ]      <- as.numeric(reticulate::py_to_r(o[["global"]]))
  }
  message(sprintf("  %s: %d anchors x %d members", term, length(anchors), length(models)))
  list(pi_mean = apply(PI, c(2, 3, 4), mean), pi_sd = apply(PI, c(2, 3, 4), sd),
       g_mean = colMeans(G), g_sd = apply(G, 2, sd),
       pi_names = pi_names, anchor_pos = A, w_start = w_start, n_members = length(models))
}
TERM <- lapply(c(N = "N", C = "C"), score_term)

## ---- 2. pool each window's opinion onto the residues it covers ------------------
## A residue is seen by up to 36 windows, each from a different position in the
## frame. Members are averaged first (the ensemble mean IS the prediction), then
## mean and sd are taken ACROSS WINDOWS -- so the band says how much the call
## depends on which window you look through, not how much the seeds disagree.
pool_residues <- function(S, term) {
  idx <- which(S$pi_names %in% DRAW_CLASSES)
  out <- vector("list", 0)
  for (k in idx) {
    acc <- vector("list", length(anchors))
    for (i in seq_along(anchors)) {
      r <- S$w_start[i] + seq_len(SEQ_LEN) - 1L
      ok <- r >= p$n_prot & r <= p$c_prot
      acc[[i]] <- tibble(residue = r[ok], v = S$pi_mean[i, ok, k])
    }
    out[[S$pi_names[k]]] <- bind_rows(acc) %>% group_by(residue) %>%
      summarise(mean = mean(v), sd = sd(v), n_win = n(), .groups = "drop") %>%
      mutate(class = S$pi_names[k], term = term)
  }
  bind_rows(out)
}
resid_tracks <- bind_rows(pool_residues(TERM$N, "N"), pool_residues(TERM$C, "C")) %>%
  mutate(class = factor(class, levels = names(nn_class_cols)))

win_tracks <- bind_rows(
  tibble(residue = anchors, mean = TERM$N$g_mean, sd = TERM$N$g_sd, term = "N", head = "CNN attention head"),
  tibble(residue = anchors, mean = TERM$C$g_mean, sd = TERM$C$g_sd, term = "C", head = "CNN attention head"))

## ---- 3. the xgb window head, from the scans --------------------------------------
read_wx <- function(path, term) {
  if (!file.exists(path)) return(NULL)
  z <- np$load(path, allow_pickle = TRUE)
  f <- as.character(reticulate::py_to_r(z$files))
  if (!paste0(gene_want, "__anchor") %in% f) return(NULL)
  a <- as.integer(reticulate::py_to_r(z[[paste0(gene_want, "__anchor")]]))
  v <- as.numeric(reticulate::py_to_r(z[[paste0(gene_want, "__xgb_pm")]]))
  s <- as.numeric(reticulate::py_to_r(z[[paste0(gene_want, "__xgb_pm__sd")]]))
  tibble(residue = a, mean = v, sd = s, term = term, head = "xgb window head")
}
win_tracks <- bind_rows(win_tracks, read_wx(wx_p[["N"]], "N"), read_wx(wx_p[["C"]], "C"))

## ---- 4. the window-free residue models -------------------------------------------
read_rp <- function() {
  if (!file.exists(rp_path)) return(NULL)
  z <- np$load(rp_path, allow_pickle = TRUE)
  have <- as.character(reticulate::py_to_r(z$files))
  cols <- grep("^(pep_pocket|pep_other|ct_context|nt_context)__(mlp|xgb)$", have, value = TRUE)
  if (!length(cols)) return(NULL)
  accs <- as.character(reticulate::py_to_r(z[["accession"]]))
  pi_i <- match(p$accession, accs)
  if (is.na(pi_i)) return(NULL)
  prot <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
  res  <- as.integer(reticulate::py_to_r(z[["resno"]]))
  sel  <- which(prot == (pi_i - 1L))
  lapply(cols, function(cc) {
    tibble(residue = res[sel], mean = as.numeric(reticulate::py_to_r(z[[cc]]))[sel],
           sd = {sdc <- paste0(cc, "__sd")
                 if (sdc %in% have) as.numeric(reticulate::py_to_r(z[[sdc]]))[sel] else NA_real_},
           target = sub("__.*$", "", cc), model = toupper(sub("^.*__", "", cc)))
  }) %>% bind_rows() %>% filter(residue >= p$n_prot, residue <= p$c_prot)
}
## lf_resid's target names are the per-index vocabulary under another name
RP_CLASS <- c(pep_pocket = "pep_pocket", pep_other = "pep_other",
              ct_context = "CT_cleavage_context", nt_context = "NT_cleavage_context")
rp <- read_rp()
if (!is.null(rp)) rp <- rp %>% mutate(class = factor(RP_CLASS[target], levels = names(nn_class_cols)))

## ---- 5. the known peptides, overlapping ones on their own lines -------------------
ll <- tryCatch(readRDS(file.path(ROOT, "inst", "extdata", "ligand_list.rds")),
               error = function(e) NULL)
peps <- if (is.null(ll)) tibble(final_name = character(), start = integer(), end = integer()) else
  ll %>% filter(accession == p$accession) %>% select(final_name, start, end) %>%
  distinct() %>% arrange(start, end)
if (nrow(peps)) {
  peps$layer <- assign_overlap_layers(peps$start, peps$end, gap = 1)
  ## every residue of every peptide, carrying the class a window would call it:
  ## pocket residues are the inserting part, the rest is the peptide body, and
  ## the 8 / 7 residues outside each end are that end's cleavage context
  message(sprintf("  %d known peptide(s), %d layer(s)", nrow(peps), max(peps$layer)))
} else message("  no annotated peptides for this accession")

saveRDS(list(resid = resid_tracks, win = win_tracks, rp = rp, peps = peps,
             prec = p, anchors = anchors, term = TERM),
        sub("\\.html$", "_data.rds", out_p))
message("data -> ", sub("\\.html$", "_data.rds", out_p))

## ---- 6. the tracks ---------------------------------------------------------------
## One patchwork, one girafe. A single SVG means ggiraph links hover by data_id
## across every panel for free: `data_id = residue` on each point, so pointing at
## a residue on any track marks the same residue on all of them.
X_LIM <- c(p$n_prot - 0.5, p$c_prot + 0.5)
BRK_W <- if (L > 400) 50 else if (L > 200) 25 else 10
x_scale <- function() scale_x_continuous(limits = X_LIM, expand = c(0, 0),
                                         breaks = scales::breaks_width(BRK_W))
bare <- function() theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.grid.major.x = element_blank(),
                         axis.title.y = element_text(size = 7.5, lineheight = 0.95),
                         legend.position = "none")

## a class-probability track. Clicks are wired in JS against `data_id`, not with
## an `onclick` aesthetic -- see the handler below for why.
mk_class_track <- function(term, legend = FALSE) {
  d <- resid_tracks %>% filter(term == !!term, class %in% DRAW_CLASSES) %>% droplevels()
  ggplot(d, aes(residue, mean, colour = class)) +
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = class),
                alpha = 0.16, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      aes(data_id = residue,
          tooltip = sprintf("%s  residue %d (%s)\n%s %.3f +/- %.3f  (%d windows)",
                            gene_want, residue, aa[residue], nn_class_labels[as.character(class)],
                            mean, sd, n_win)),
      size = 1.1, na.rm = TRUE) +
    scale_colour_manual(values = nn_class_cols, labels = nn_class_label_fn, name = NULL) +
    scale_fill_manual(values = nn_class_cols, guide = "none") +
    x_scale() + scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = NULL, y = sprintf("%s per-residue\nclass probability", term)) +
    theme_bw() + bare() +
    theme(legend.position = if (legend) "top" else "none", legend.direction = "horizontal",
          legend.title = element_blank(), legend.text = element_text(size = 7),
          legend.key.size = unit(8, "pt"), legend.margin = margin(0, 0, 0, 0))
}

## a window-level head on the anchor axis: one value per "a peptide ends here"
mk_win_track <- function(term, head_nm) {
  d <- win_tracks %>% filter(term == !!term, head == head_nm)
  if (!nrow(d)) return(NULL)
  ggplot(d, aes(residue, mean)) +
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                fill = LF_VIZ$ink, alpha = 0.15, na.rm = TRUE) +
    geom_line(colour = LF_VIZ$ink, linewidth = 0.5, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      aes(data_id = residue,
          tooltip = sprintf("%s  anchor %d (%s)\n%s %.3f +/- %.3f",
                            gene_want, residue, aa[residue], head_nm, mean, sd)),
      size = 1.1, colour = LF_VIZ$ink, na.rm = TRUE) +
    x_scale() + scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = NULL, y = sprintf("%s %s", term, sub(" head", "\nhead", head_nm))) +
    theme_bw() + bare()
}

## the window-free models: no anchor, no frame, one residue's own features
mk_rp_track <- function(mdl) {
  d <- rp %>% filter(model == mdl) %>% droplevels()
  if (!nrow(d)) return(NULL)
  ggplot(d, aes(residue, mean, colour = class)) +
    geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = class),
                alpha = 0.16, colour = NA, na.rm = TRUE) +
    geom_line(linewidth = 0.5, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      aes(data_id = residue,
          tooltip = sprintf("%s  residue %d (%s)\n%s %s %.3f",
                            gene_want, residue, aa[residue], mdl,
                            nn_class_labels[as.character(class)], mean)),
      size = 1.1, na.rm = TRUE) +
    scale_colour_manual(values = nn_class_cols, labels = nn_class_label_fn, guide = "none") +
    scale_fill_manual(values = nn_class_cols, guide = "none") +
    x_scale() + scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = NULL, y = sprintf("%s per residue\n(no window)", mdl)) +
    theme_bw() + bare()
}

## the annotation track: each known peptide on its own line when it overlaps
mk_peps <- function() {
  if (!nrow(peps)) return(NULL)
  d <- peps %>% mutate(y = layer)
  ggplot(d) +
    ggiraph::geom_segment_interactive(
      aes(x = start, xend = end, y = y, yend = y, data_id = final_name,
          tooltip = sprintf("%s  %d-%d  (%d aa)", final_name, start, end, end - start + 1L)),
      colour = nn_class_cols[["pep_pocket"]], linewidth = 3, lineend = "butt") +
    geom_text(aes(x = (start + end) / 2, y = y, label = final_name),
              size = 2.4, colour = "white", fontface = "bold") +
    x_scale() +
    scale_y_continuous(limits = c(0.4, max(d$y) + 0.6), expand = c(0, 0), breaks = NULL) +
    labs(x = NULL, y = "known\npeptides") + theme_bw() + bare()
}

## the sequence strip, and the only x axis on the page
mk_seq <- function() {
  s <- tibble(residue = seq_len(L), AA = aa) %>% filter(residue >= p$n_prot, residue <= p$c_prot)
  s$kind <- "other"
  db <- which(grepl("KK|KR|RK|RR", substring(p$seq, seq_len(L), pmin(seq_len(L) + 1L, L))))
  s$kind[s$residue %in% c(db, db + 1L)] <- "dibasic pair (KK/KR/RK/RR)"
  s$kind[s$kind == "other" & s$AA %in% c("K", "R")] <- "single K / R"
  ggplot(s, aes(residue, 1)) +
    ggiraph::geom_tile_interactive(aes(fill = kind, data_id = residue,
        tooltip = sprintf("residue %d: %s", residue, AA)), width = 1, height = 1) +
    scale_fill_manual(values = c(`dibasic pair (KK/KR/RK/RR)` = LF_VIZ$s2,
                                 `single K / R` = LF_VIZ$s3, other = LF_VIZ$grid), name = NULL) +
    x_scale() + scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    labs(x = sprintf("residue of %s (%s)", gene_want, p$accession), y = NULL) +
    theme_bw() +
    theme(legend.position = "bottom", legend.text = element_text(size = 7),
          legend.key.size = unit(8, "pt"), panel.grid = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank())
}

## ---- 7. assemble ---------------------------------------------------------------
## ONE WIDGET PER TRACK, not one patchwork inside one girafe. ggplot2 4.x moved
## to S7 and patchwork 1.3.1 cannot be drawn through ggiraph's dsvg device on
## this combination: the figure comes out of `ggsave` complete and out of
## `girafe` truncated after the second panel, silently and with no warning. The
## tracks are therefore separate girafe widgets stacked in one page, which is
## also what the protein pages already do.
panels <- Filter(Negate(is.null), list(
  mk_peps(),
  mk_class_track("N", legend = TRUE), mk_win_track("N", "CNN attention head"),
  mk_win_track("N", "xgb window head"),
  mk_class_track("C"), mk_win_track("C", "CNN attention head"),
  mk_win_track("C", "xgb window head"),
  mk_rp_track("MLP"), mk_rp_track("XGB"),
  mk_seq()))
hts <- c(if (nrow(peps)) 0.9 else NULL, 2.6, 1.6, 1.6, 2.6, 1.6, 1.6,
         if (!is.null(rp)) c(1.8, 1.8) else NULL, 0.7)
stopifnot(length(panels) == length(hts))

## Separate widgets each size their own y-axis gutter, so a track with a long
## label would sit further right than its neighbours and the residue axes would
## drift apart. Pad every panel out to the widest gutter on both sides.
PLOT_W <- max(11, min(30, 3.5 + L * 0.055))
lefts  <- vapply(panels, panel_left_in,  numeric(1))
rights <- vapply(panels, panel_right_in, numeric(1))
panels <- lapply(seq_along(panels), function(i)
  panels[[i]] + theme(plot.margin = margin(t = 2, b = 2,
                                           l = (max(lefts)  - lefts[i])  * 72,
                                           r = (max(rights) - rights[i]) * 72, unit = "pt")))

mk_wgt <- function(p, h) {
  w <- ggiraph::girafe(ggobj = p, width_svg = PLOT_W, height_svg = h)
  ggiraph::girafe_options(w,
    ggiraph::opts_hover(css = "stroke:#000;stroke-width:1.4px;"),
    ## the same residue lights up on every OTHER track too: ggiraph keys hover
    ## by data_id across every girafe figure on the page
    ggiraph::opts_hover_inv(css = "opacity:0.3;"),
    ggiraph::opts_tooltip(css = paste0("background:#fcfcfb;border:1px solid ", LF_VIZ$grid,
                                       ";padding:5px 7px;font:11px ui-sans-serif,system-ui;",
                                       "border-radius:4px;white-space:pre;")),
    ggiraph::opts_selection(type = "none"),
    ggiraph::opts_sizing(rescale = FALSE),
    ggiraph::opts_toolbar(saveaspng = FALSE))
}
wgts <- Map(mk_wgt, panels, hts)

## ---- 8. the popup payload ----------------------------------------------------------
## Every anchor's 36 x K surface travels with the page, so a click draws the
## window with no round trip. Rounded to 3 dp: the popup cannot show more.
pack <- function(S) {
  keep <- if (min_pop > 0) which(S$g_mean >= min_pop) else seq_along(anchors)
  list(anchor = anchors[keep], pos = S$anchor_pos, g = round(S$g_mean[keep], 4),
       pi = lapply(keep, function(i) round(S$pi_mean[i, , ], 3)))
}
payload <- list(
  gene = gene_want, accession = p$accession, aa = aa,
  n_prot = p$n_prot, c_prot = p$c_prot, seq_len = SEQ_LEN,
  classes = TERM$C$pi_names, labels = unname(nn_class_labels[TERM$C$pi_names]),
  cols = unname(nn_class_cols[TERM$C$pi_names]),
  N = pack(TERM$N), C = pack(TERM$C))

## ggiraph 0.8.13 drops the `onclick` aesthetic on this ggplot2 -- the rendered
## svg carries data-id and title and no onclick at all -- so the click handler is
## attached here instead, to the data-id elements that ARE emitted. data-id is
## the bare residue (it has to be: that is what links hover across tracks), so a
## click cannot say which terminus it came from. The popup therefore shows BOTH
## frames anchored at that residue, which is the more useful answer anyway.
js <- sprintf("
window.LFV2 = %s;
(function () {
  var d = window.LFV2, box = null;
  function mkbox() {
    box = document.createElement('div'); box.id = 'lf-popup';
    box.style.cssText = 'position:fixed;display:none;z-index:9999;background:#fcfcfb;' +
      'border:1px solid #c9c8c3;border-radius:6px;box-shadow:0 8px 28px rgba(0,0,0,.18);' +
      'padding:10px 12px;font:12px ui-sans-serif,system-ui;max-height:92%%;overflow:auto;';
    document.body.appendChild(box);
    document.addEventListener('mousedown', function (ev) {
      if (box && !box.contains(ev.target) && !ev.target.closest('svg')) box.style.display = 'none';
    });
  }
  function frame(term, res) {
    var S = d[term], i = S.anchor.indexOf(res);
    if (i < 0) return '<div style=\"color:#888\">no ' + term + ' window at ' + res + '</div>';
    var K = d.classes.length, N = d.seq_len, W = 470, H = 150, PAD = 30;
    var sx = function (j) { return PAD + j * (W - PAD - 8) / (N - 1); };
    var sy = function (v) { return H - PAD - v * (H - PAD - 12); };
    var g = '<svg width=' + W + ' height=' + H + '>';
    var ax = sx(S.pos - 1);
    g += '<rect x=' + (ax - 1.5) + ' y=8 width=3 height=' + (H - PAD - 8) + ' fill=#000 opacity=.12/>';
    for (var k = 0; k < K; k++) {
      if (d.classes[k] === 'none' || d.classes[k] === 'padding') continue;
      var pts = [];
      for (var j = 0; j < N; j++) pts.push(sx(j).toFixed(1) + ',' + sy(S.pi[i][j][k]).toFixed(1));
      g += '<polyline points=\"' + pts.join(' ') + '\" fill=none stroke=\"' + d.cols[k] + '\" stroke-width=1.7/>';
    }
    g += '<line x1=' + PAD + ' y1=' + (H - PAD) + ' x2=' + (W - 8) + ' y2=' + (H - PAD) + ' stroke=#aaa/>';
    g += '<text x=2 y=' + (sy(1) + 4) + ' font-size=9 fill=#666>1.0</text>';
    g += '<text x=2 y=' + (H - PAD) + ' font-size=9 fill=#666>0.0</text>';
    var w0 = res - S.pos + 1;
    for (var j = 0; j < N; j++) {
      var r = w0 + j, ch = (r >= 1 && r <= d.aa.length) ? d.aa[r - 1] : '.';
      g += '<text x=' + sx(j).toFixed(1) + ' y=' + (H - 14) + ' font-size=7 text-anchor=middle fill=' +
           (r === res ? '#000' : '#999') + '>' + ch + '</text>';
    }
    g += '</svg>';
    return '<div style=\"font-weight:600;margin:6px 0 1px\">' + term +
      ' window &middot; covers ' + w0 + '&ndash;' + (w0 + N - 1) +
      ' &middot; anchor at frame position ' + S.pos +
      ' &middot; window score ' + S.g[i].toFixed(3) + '</div>' + g;
  }
  window.lfPopup = function (res, ev) {
    if (!box) mkbox();
    var leg = d.classes.map(function (c, k) {
      if (c === 'none' || c === 'padding') return '';
      return '<span style=\"color:' + d.cols[k] + ';font-weight:600;margin-right:10px\">&#9632; ' +
             d.labels[k] + '</span>'; }).join('');
    box.innerHTML = '<div style=\"font-weight:700\">' + d.gene + ' &mdash; residue ' + res +
      ' (' + d.aa[res - 1] + ')</div>' + frame('N', res) + frame('C', res) +
      '<div style=\"margin-top:6px\">' + leg + '</div>';
    box.style.display = 'block';
    var px = ev ? ev.clientX : 150, py = ev ? ev.clientY : 100;
    box.style.left = Math.max(8, Math.min(px + 16, window.innerWidth - 520)) + 'px';
    box.style.top  = Math.max(8, Math.min(py - 60, window.innerHeight - 420)) + 'px';
  };
  function wire() {
    document.querySelectorAll('svg [data-id]').forEach(function (n) {
      if (n.__lf) return;
      var v = n.getAttribute('data-id');
      if (!/^[0-9]+$/.test(v)) return;          // peptide segments carry a name, not a residue
      n.__lf = true; n.style.cursor = 'pointer';
      n.addEventListener('click', function (ev) { window.lfPopup(parseInt(v, 10), ev); });
    });
  }
  // the widgets render after this script runs, so wire on idle and once more later
  if (document.readyState === 'complete') setTimeout(wire, 300); else
    window.addEventListener('load', function () { setTimeout(wire, 300); });
  setTimeout(wire, 1500); setTimeout(wire, 4000);
})();
", jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 4, null = "null"))

page <- htmltools::tagList(
  htmltools::tags$style(htmltools::HTML(
    "body{margin:0;background:#fcfcfb;font:13px ui-sans-serif,system-ui,-apple-system,Arial}
     .lfhdr{padding:14px 18px 6px} .lfhdr h1{margin:0 0 3px;font-size:17px}
     .lfhdr .sub{color:#52514e;font-size:11.5px;line-height:1.45;max-width:1100px}
     .lftracks{padding:0 10px 24px}
     .girafe_container_std{margin:0 !important}")),
  htmltools::tags$div(class = "lfhdr",
    htmltools::tags$h1(sprintf("%s (%s) — per-residue and per-window predictions, both termini",
                               gene_want, p$accession)),
    htmltools::tags$div(class = "sub", htmltools::HTML(sprintf(
      "Click any point to open the 36-residue windows anchored THERE, for both termini. Hovering a residue marks it on every track.<br>N anchors sit at window position %d of %d, C at %d. Trunks: N <code>%s</code>, C <code>%s</code>, %d members each. %d anchors per terminus carry a popup.",
      LF_PEPEND$anchor[["N"]], SEQ_LEN, LF_PEPEND$anchor[["C"]],
      basename(tools::file_path_sans_ext(run_p[["N"]])),
      basename(tools::file_path_sans_ext(run_p[["C"]])),
      TERM$C$n_members, length(payload$C$anchor))))),
  htmltools::tags$div(class = "lftracks", wgts),
  htmltools::tags$script(htmltools::HTML(js)))

lib_dir <- file.path(dirname(out_p),
                     paste0(tools::file_path_sans_ext(basename(out_p)), "_files"))
htmltools::save_html(page, out_p, libdir = lib_dir)
n_inlined <- nn_inline_deps(out_p, lib_dir)
if (isTRUE(n_inlined > 0)) unlink(lib_dir, recursive = TRUE)
message(sprintf("page -> %s  (%.1f MB, %d tracks, %s deps inlined)",
                out_p, file.size(out_p) / 1048576, length(wgts),
                if (isTRUE(n_inlined > 0)) n_inlined else "NO"))
