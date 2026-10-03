#!/usr/bin/env Rscript
## ---- per-RESIDUE class tracks, pooled over every window that sees a residue ----
## 10_9i plots the window-level score on the ANCHOR axis, because that score is a
## property of the whole 36-residue window and belongs to the cleavage hypothesis
## "the peptide ends here". This plots something different in kind: the per-index
## head's prediction for a RESIDUE, which is a property of that residue.
##
## Every mature residue is seen by up to 36 different windows -- a window
## anchored at `a` covers a-27..a+8 (C terminus, anchor at window position 28),
## so residue r sits at window position r - a + 28 of every window with
## r-27 <= a <= r+8. Each of those windows gives its own opinion of what r is,
## from a different position in the window and a different 36-residue context.
##
## So the band here is NOT ensemble spread: it is how much the call for a residue
## depends on WHICH WINDOW you look at it through. Members are averaged first
## (the ensemble mean is the model's prediction), then mean and sd are taken
## across the covering windows. The across-member sd is kept in the output table
## as `sd_member` for comparison, but it is not what is drawn.
##
##   Rscript inst/scripts/10_9n_pepend_residue_tracks.R --gene NPY
##   ... --run ~/AF2_analysis/lf_pepend_run_C.rds --classes pep_pocket,pep_other
##
## Output: <out-dir>/lf_resid_tracks_<gene>_<term>.svg  (+ .png unless --no-png)
##         <out-dir>/lf_resid_tracks_<gene>_<term>.csv  (the tidy table)
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
## several genes share one load of the residue cache and one rebuild of the
## 20-member ensemble, which is the whole cost -- so --genes, not one run each
genes   <- strsplit(.opt("--genes", .opt("--gene", "NPY")), ",")[[1]] |> trimws()
term    <- toupper(.opt("--term", "C"))
run_p   <- path.expand(.opt("--run", "~/AF2_analysis/lf_pepend_run_C_ins.rds"))
cache_p <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
scan_p  <- path.expand(.opt("--scan", "~/AF2_analysis/lf_pepend_scan_C_ins.rds"))
in_path <- path.expand(.opt("--input", "~/AF2_analysis/lf_pepend_nn_input.rds"))
## Per-residue MLP / xgboost scores from `python -m lf_resid train`. These are a
## different model FAMILY from everything else on the page: no window, no
## convolution, one residue's own features (plus +/-2 neighbours) and nothing
## more. Absent file -> no panel.
rp_path <- path.expand(.opt("--resid-preds", "~/AF2_analysis/lf_resid_preds.npz"))
## The xgboost window heads, scanned over every anchor and averaged across the
## same 20 members as the CNN track -- so the three window-level panels are the
## same quantity computed three ways. From `python -m lf_winxgb.scan_genes`.
wx_path <- path.expand(.opt("--winxgb-scan", "~/AF2_analysis/lf_winxgb_scan_C.npz"))
## Which xgboost window heads to draw. Default is the per-member head alone:
## the stacked head is one tree model on 20 correlated views of the same 41
## positives and is not the one going forward.
wx_heads <- strsplit(.opt("--winxgb-heads", "xgb_pm"), ",")[[1]] |> trimws()
ll_p    <- .opt("--ligand-list", "inst/extdata/ligand_list.rds")
out_dir <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
## "all" expands to the whole per-index vocabulary once the Config is loaded,
## since only the model knows what that vocabulary is.
classes <- strsplit(.opt("--classes", "pep_pocket,CT_cleavage_context"), ",")[[1]] |> trimws()
no_png  <- "--no-png" %in% .args
meta_p  <- path.expand(.opt("--meta", file.path(out_dir, "lf_resid_tracks_meta.rds")))
## one stacked-area panel instead of one panel per class
combined <- "--combined" %in% .args
## Classes computed and written to the csv, but not DRAWN as lines. `none` is
## the vocabulary's bookkeeping column -- it is the argmax at most residues and
## dominates the panel, while the ground truth it sits above never contains it.
## --drop-classes "" keeps everything.
## `padding` joins it: a mature residue is never padding in a window that
## covers it, so the line is flat at ~0 by construction -- a correctness check
## worth running (and it is still in the csv) but dead space on the panel.
drop_cls <- strsplit(.opt("--drop-classes", "none,padding"), ",")[[1]] |> trimws()
drop_cls <- drop_cls[nzchar(drop_cls)]

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "dcnn_bridge.R"))
source(file.path(ROOT, "R", "viz_theme.R"))
## nn_class_cols: the class palette the protein page's detail panels use.
## Sourced rather than copied so the two figures cannot drift apart.
source(file.path(ROOT, "R", "plot_proteins_win_new.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)

## ---- the xgboost window-head scans, loaded once --------------------------------
WX <- NULL
if (file.exists(wx_path)) {
  .wz <- np$load(wx_path, allow_pickle = TRUE)
  .wf <- as.character(reticulate::py_to_r(.wz$files))
  .wg <- unique(sub("__.*$", "", grep("__anchor$", .wf, value = TRUE)))
  if (length(.wg)) {
    WX <- list(z = .wz, files = .wf, genes = .wg)
    message("xgb window heads: ", length(.wg), " gene(s) from ", basename(wx_path))
  }
}

## ---- the window-free residue scores, loaded once -------------------------------
RP <- NULL
if (file.exists(rp_path)) {
  .z <- np$load(rp_path, allow_pickle = TRUE)
  .have <- as.character(reticulate::py_to_r(.z$files))
  .cols <- grep("^(pep_pocket|ct_context)__(mlp|xgb)$", .have, value = TRUE)
  if (length(.cols)) {
    RP <- list(cols = .cols,
               prot = as.integer(reticulate::py_to_r(.z[["prot_idx"]])),
               res  = as.integer(reticulate::py_to_r(.z[["resno"]])),
               acc  = as.character(reticulate::py_to_r(.z[["accession"]])))
    for (cc in .cols) {
      RP[[cc]] <- as.numeric(reticulate::py_to_r(.z[[cc]]))
      sdc <- paste0(cc, "__sd")
      RP[[sdc]] <- if (sdc %in% .have) as.numeric(reticulate::py_to_r(.z[[sdc]]))
                   else rep(NA_real_, length(RP[[cc]]))
    }
    message("residue models: ", paste(.cols, collapse = ", "), " from ", basename(rp_path))
  }
}

ANCHOR <- c(N = 8L, C = 28L)[[term]]
SEQ_LEN <- 36L
RES <- c(`dibasic pair (KK/KR/RK/RR)` = LF_VIZ$s2,
         `single K / R`               = LF_VIZ$s3,
         other                        = LF_VIZ$grid)

fc <- readRDS(cache_p)
cfg <- mod$Config$from_json(file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"),
                                      "in", "config.json"))
models <- mod$scan$load_member_models(
  file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out"), cfg, term)
pi_names <- as.character(cfg$pi_names)
if (identical(tolower(classes), "all")) classes <- pi_names
if (!all(classes %in% pi_names))
  stop("unknown class(es): ", paste(setdiff(classes, pi_names), collapse = ", "),
       "\n  the vocabulary is: ", paste(pi_names, collapse = ", "), call. = FALSE)

make_one <- function(gene) {
  ## ---- 1. the precursor and all of its windows ---------------------------------
  p  <- fc$prec %>% filter(gene == !!gene) %>% slice(1)
  if (!nrow(p)) stop("no precursor called ", gene, call. = FALSE)
  feat <- fc$feats[[match(p$accession, fc$prec$accession)]]
  L <- nchar(p$seq); aa <- strsplit(p$seq, "")[[1]]

  anchors <- seq.int(p$n_prot, p$c_prot)          # every mature residue, as 10_8b scans
  w_start <- anchors - ANCHOR + 1L
  message(sprintf("%s (%s): %d residues, mature %d-%d, %d windows",
                  gene, p$accession, L, p$n_prot, p$c_prot, length(anchors)))

  X <- array(0, c(length(anchors), SEQ_LEN, ncol(feat)))
  for (i in seq_along(anchors))
    X[i, , ] <- as.matrix(lf_pepend_slice(feat, w_start[i], p$n_prot, p$c_prot)[, fc$all_params3])
  xp <- np$ascontiguousarray(np$asarray(X, dtype = "float32"))

  ## ---- 2. every member's per-index softmax -------------------------------------

  ## The same forward pass carries the window-level heads, so collect them here
  ## rather than paying for a second one: `global` always, `ins_class` when the
  ## run was trained with --ins-head.
  ins_names <- tryCatch(as.character(cfg$ins_class_names), error = function(e) character(0))
  PI  <- array(0, c(length(models), length(anchors), SEQ_LEN, length(pi_names)))
  G   <- matrix(0, length(models), length(anchors))
  INS <- NULL
  for (m in seq_along(models)) {
    o <- mod$pipeline$predict_all(models[[m]], xp)
    PI[m, , , ] <- reticulate::py_to_r(o[["per_index_cat"]])
    G[m, ] <- as.numeric(reticulate::py_to_r(o[["global"]]))
    if ("ins_class" %in% names(o)) {
      if (is.null(INS)) INS <- array(0, c(length(models), length(anchors), length(ins_names)))
      INS[m, , ] <- reticulate::py_to_r(o[["ins_class"]])
    }
  }
  message(sprintf("scored %d window(s) x %d member(s); classes: %s",
                  length(anchors), length(models), paste(classes, collapse = ", ")))

  ## ---- 3. pool every window's opinion of each residue --------------------------
  ## mean over members FIRST (that is the model's prediction), then across windows.
  pi_mean <- apply(PI, c(2, 3, 4), mean)          # (window, position, class)
  pi_sdm  <- apply(PI, c(2, 3, 4), sd)            # across-member sd, kept for the table

  long <- expand_grid(wi = seq_along(anchors), pos = seq_len(SEQ_LEN)) %>%
    mutate(anchor  = anchors[wi],
           ## the inverse of `r sits at position r - a + 28`
           residue = anchor - ANCHOR + pos) %>%
    ## a window position whose residue falls outside the mature range is padding,
    ## not an opinion about a residue -- drop it rather than average it in
    filter(residue >= p$n_prot, residue <= p$c_prot)

  tracks <- lapply(classes, function(cl) {
    k <- match(cl, pi_names)
    long %>%
      mutate(p_mean = pi_mean[cbind(wi, pos, k)],
             p_sdm  = pi_sdm[cbind(wi, pos, k)]) %>%
      group_by(residue) %>%
      summarise(class     = cl,
                n_windows = n(),
                mean      = mean(p_mean),
                sd        = sd(p_mean),          # ACROSS WINDOWS -- the band
                sd_member = mean(p_sdm),         # mean across-member sd, for reference
                min       = min(p_mean),
                max       = max(p_mean),
                .groups = "drop")
  }) %>% bind_rows() %>% mutate(class = factor(class, levels = classes))

  message(sprintf("windows per residue: %d-%d (median %d)",
                  min(tracks$n_windows), max(tracks$n_windows),
                  as.integer(median(tracks$n_windows))))

  ## ---- 3b. the window-level heads, on the ANCHOR axis --------------------------
  ## These are a different KIND of quantity from the tracks above. A window score
  ## or an insertion class belongs to the whole 36-residue window and to the
  ## hypothesis "the peptide ends at this anchor"; it is drawn at the anchor, not
  ## at a residue it describes. And its band is the spread across the 20 SEEDS,
  ## where the panel above is the spread across WINDOWS. Same picture language,
  ## two different uncertainties -- which the axis labels have to say out loud.
  WIN_PAL <- c(`window score` = LF_VIZ$ink, inserting = "#1B7837",
               loop = "#E08214", non_inserting = "#878787")
  win <- tibble(anchor = anchors, mean = colMeans(G), sd = apply(G, 2, sd),
                class = "window score")
  if (!is.null(INS))
    win <- bind_rows(win, lapply(seq_along(ins_names), function(k)
      tibble(anchor = anchors, mean = colMeans(INS[, , k]),
             sd = apply(INS[, , k], 2, sd), class = ins_names[k])) %>% bind_rows())
  win <- win %>% mutate(class = factor(class, levels = names(WIN_PAL)))

  ## ---- 3c. GROUND TRUTH for this gene's known window(s) ------------------------
  ## The labels the model was fitted against, on the same residue axis as the
  ## predictions: the peptide's span, its per-index classes, and where the pocket
  ## sits relative to each of its two ends -- which is the whole definition of the
  ## insertion class. `lf_pepend_insertion_class(detail = TRUE)` is called rather
  ## than the geometry re-derived here, so the annotation cannot disagree with the
  ## label the head was actually trained on.
  truth_seg <- truth_pep <- truth_end <- NULL
  kw <- tryCatch(readRDS(in_path)$nn_input[[term]]$all %>%
                   filter(gene == !!gene, known == 1L),
                 error = function(e) NULL)
  if (!is.null(kw) && nrow(kw)) {
    det <- lf_pepend_insertion_class(kw$known_idx, kw$term, kw$len, detail = TRUE)
    truth_seg <- lapply(seq_len(nrow(kw)), function(i) {
      ki <- as.character(kw$known_idx[[i]])
      res <- kw$anchor[i] - ANCHOR + seq_along(ki)
      r <- rle(ki); e <- cumsum(r$lengths); st <- c(1L, head(e, -1) + 1L)
      tibble(peps = kw$peps[i], class = r$values, start = res[st], end = res[e])
    }) %>% bind_rows() %>%
      ## padding rows are not labels about residues
      filter(class %in% names(nn_class_cols), class != "padding") %>%
      mutate(class = factor(class, levels = names(nn_class_cols)))
    truth_pep <- tibble(peps = kw$peps, name = kw$pep_name,
                        start = kw$pep_start, end = kw$pep_end,
                        ins_class = det$ins_class,
                        gap_end = det$gap_end, gap_far = det$gap_far)
    ## the two ends of the peptide: the one the pocket reaches (inserting) and the
    ## one it does not. `gap_*` are in residues; tol is 2 in the class definition.
    truth_end <- bind_rows(
      tibble(x = kw$pep_end, kind = "inserting end", gap = det$gap_end),
      tibble(x = kw$pep_start, kind = "non-inserting end", gap = det$gap_far)) %>%
      mutate(lab = sprintf("%s\n(pocket %s away)", kind,
                           ifelse(is.na(gap), "absent", sprintf("%.0f res", gap))))
    message(sprintf("ground truth: %d known window(s); %s",
                    nrow(kw), paste(sprintf("%s = %s (gap_end %s, gap_far %s)",
                      kw$peps, det$ins_class, det$gap_end, det$gap_far), collapse = "; ")))
  }

  ## ---- 3d. the CORRECTLY POSITIONED window --------------------------------
  ## The panel above pools every window that sees a residue, which is what you
  ## have when you do not know where the peptide ends. This is the opposite: the
  ## single window anchored exactly ON the known end, i.e. what the model says
  ## when handed the right frame. Its band is across the 20 SEEDS, not across
  ## windows -- there is only one window here.
  ##
  ## Read it against the ground-truth strip directly below: same residues, same
  ## classes, prediction over label.
  cw <- NULL
  if (!is.null(kw) && nrow(kw)) {
    kk <- kw %>% mutate(wi = match(anchor, anchors)) %>% filter(!is.na(wi))
    ## ONE window per gene. POMC has 4 known C ends, GCG and PTHLH 3; a panel
    ## each would make those figures unreadable and the rest of the gallery
    ## inconsistent. Highest-scoring, so the pick is deterministic and is the
    ## window the model itself is most confident about.
    if (nrow(kk) > 1) {
      kk <- kk %>% mutate(.s = colMeans(G)[wi]) %>% arrange(desc(.s))
      message(sprintf("  %s: %d known windows, drawing %s (score %.3f)",
                      gene, nrow(kk), kk$peps[1], kk$.s[1]))
      kk <- kk %>% slice_head(n = 1)
    }
    cw <- lapply(seq_len(nrow(kk)), function(i) {
      wi <- kk$wi[i]
      mu <- apply(PI[, wi, , , drop = FALSE], c(3, 4), mean)
      sg <- apply(PI[, wi, , , drop = FALSE], c(3, 4), sd)
      res <- kk$anchor[i] - ANCHOR + seq_len(SEQ_LEN)
      expand_grid(pos = seq_len(SEQ_LEN), k = seq_along(pi_names)) %>%
        mutate(residue = res[pos], class = pi_names[k],
               mean = mu[cbind(pos, k)], sd = sg[cbind(pos, k)]) %>%
        filter(residue >= p$n_prot, residue <= p$c_prot, class %in% classes,
               !class %in% drop_cls) %>%
        mutate(class = factor(class, levels = setdiff(intersect(names(nn_class_cols),
                                                                classes), drop_cls)),
               ttl = sprintf("%s\n@ %d (band: seeds)",
                             substr(kk$pep_name[i], 1, 22), kk$anchor[i]))
    })
  }

  ## ---- 3e. the window-free residue scores for THIS precursor ---------------
  rp <- NULL
  if (!is.null(RP)) {
    pi_i <- match(p$accession, RP$acc)
    if (!is.na(pi_i)) {
      sel <- which(RP$prot == (pi_i - 1L))
      if (length(sel)) {
        rp <- lapply(RP$cols, function(cc) {
          tibble(residue = RP$res[sel], mean = RP[[cc]][sel],
                 sd = RP[[paste0(cc, "__sd")]][sel],
                 target = sub("__.*$", "", cc), model = sub("^.*__", "", cc))
        }) %>% bind_rows() %>%
          filter(residue >= p$n_prot, residue <= p$c_prot) %>%
          ## colour by the CLASS so it reads against the per-index panel above,
          ## linetype by the model: same colour means the same thing predicted
          mutate(class = factor(ifelse(target == "pep_pocket", "pep_pocket",
                                       "CT_cleavage_context"),
                                levels = c("pep_pocket", "CT_cleavage_context")),
                 model = factor(model, levels = c("xgb", "mlp")))
      }
    }
  }

  ## ---- 3f. the xgboost window heads for THIS gene --------------------------
  wx <- NULL
  if (!is.null(WX) && gene %in% WX$genes) {
    .g <- function(k) as.numeric(reticulate::py_to_r(WX$z[[paste0(gene, "__", k)]]))
    .an <- as.integer(reticulate::py_to_r(WX$z[[paste0(gene, "__anchor")]]))
    wx <- bind_rows(lapply(intersect(c("xgb_pm", "xgb_st"), wx_heads), function(k)
      tibble(anchor = .an, mean = .g(k), sd = .g(paste0(k, "__sd")), head = k))) %>%
      mutate(head = factor(head, levels = intersect(c("xgb_pm", "xgb_st"), wx_heads)))
  }

  ## ---- 4. annotation: known ends and peptide spans -----------------------------
  sc <- tryCatch(readRDS(scan_p), error = function(e) NULL)
  kn <- if (!is.null(sc)) sc$knowns %>% filter(accession == p$accession) else NULL
  LL <- tryCatch(readRDS(file.path(ROOT, ll_p)), error = function(e) NULL)
  peps <- if (is.null(LL)) tibble(final_name = character(), start = integer(), end = integer()) else
    LL %>% filter(accession == p$accession) %>% select(final_name, start, end) %>%
    distinct() %>% arrange(start)

  PLOT_W  <- max(11, min(42, 3.5 + L * 0.034))
  per_res <- (PLOT_W - 2.2) * 25.4 / L
  SEQ_SZ  <- max(0.85, min(1.9, per_res / 0.62))
  SHOW_AA <- per_res / 0.62 >= 0.8
  BRK_W   <- if (L > 800) 100 else if (L > 400) 50 else if (L > 200) 25 else 10

  ## One colour per class so six stacked panels stay distinguishable at a glance.
  ## `none` and `padding` are deliberately grey: they are the vocabulary's
  ## bookkeeping columns, not biology.
  CLS_PAL <- c(pep_pocket          = LF_VIZ$s1,
               CT_cleavage_context = "#1B7837",
               NT_cleavage_context = "#E08214",
               pep_other           = "#8073AC",
               none                = "#878787",
               padding             = "#BDBDBD")
  mk_track <- function(cl) {
    d <- tracks %>% filter(class == cl)
    col <- if (cl %in% names(CLS_PAL)) CLS_PAL[[cl]] else LF_VIZ$s1
    ggplot(d, aes(residue, mean)) +
      {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                  ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
      {if (nrow(peps)) geom_rect(data = peps, inherit.aes = FALSE,
          aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
          fill = LF_VIZ$s1, alpha = 0.07) } +
      geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1)),
                  fill = col, alpha = 0.22) +
      geom_line(colour = col, linewidth = 0.8) +
      {if (!is.null(kn) && nrow(kn)) geom_vline(data = kn, aes(xintercept = anchor),
          colour = LF_VIZ$ink, linewidth = 0.4, linetype = "22") } +
      scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                         breaks = scales::breaks_width(BRK_W)) +
      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      labs(x = NULL, y = sprintf("p(%s)", cl)) +
      lf_viz_theme() +
      theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  }

  ## The subtitle used to say "the band is window-position sensitivity, not
  ## ensemble spread" -- true of the residue panel, and FALSE of the window panels
  ## added below it, whose bands are across seeds. With panels of both kinds on one
  ## canvas the blanket claim had to go; each panel's y axis names its own sd.
  ttl <- ggplot() + theme_void() +
    labs(title = sprintf("%s (%s): per-residue and per-window predictions, step-1 %s scan",
                         gene, p$accession, term),
         subtitle = sprintf(paste0(
           "TOP, a residue property: each mature residue is seen by up to %d windows (median %d ",
           "here) -- a window anchored at a covers\na-%d..a+%d, so residue r sits at window position ",
           "r - a + %d. Members averaged first, then mean and +/- 1 sd ACROSS WINDOWS.%s",
           "\nDashed line: a known peptide end.%s"),
           SEQ_LEN, as.integer(median(tracks$n_windows)), ANCHOR - 1L, SEQ_LEN - ANCHOR, ANCHOR,
           if (combined)
             paste0("\nBELOW, window properties on the anchor axis: each value belongs to a whole ",
                    "36-residue window and to the hypothesis\n\"the peptide ends at this anchor\", ",
                    "and its band is +/- 1 sd ACROSS THE ", length(models), " SEEDS. The two kinds ",
                    "of band are not comparable.")
           else "",
           if (!is.na(grp <- tryCatch(paste(sort(unique(as.character(kn$split))), collapse = "/"),
                                      error = function(e) NA_character_)) && nzchar(grp))
             sprintf("  This gene's known end is in the %s split.", toupper(grp)) else "")) +
    theme(plot.title = element_text(size = 12, face = "bold", colour = LF_VIZ$ink),
          plot.subtitle = element_text(size = 8, colour = LF_VIZ$ink2, lineheight = 1.2))

  is_basic <- aa %in% c("K", "R")
  in_pair  <- rep(FALSE, L)
  pr <- which(head(is_basic, -1) & tail(is_basic, -1))
  in_pair[c(pr, pr + 1L)] <- TRUE
  seqd <- tibble(pos = seq_len(L), aa = aa,
                 res = ifelse(in_pair, names(RES)[1],
                       ifelse(is_basic, names(RES)[2], names(RES)[3])))
  p_seq <- ggplot(seqd, aes(pos, 1)) +
    geom_tile(aes(fill = res), height = 0.42, colour = NA) +
    {if (SHOW_AA) geom_text(aes(label = aa), size = SEQ_SZ, colour = LF_VIZ$ink, vjust = 0.5) } +
    scale_fill_manual(values = RES, name = "residue", breaks = names(RES)) +
    scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                       breaks = scales::breaks_width(BRK_W)) +
    labs(x = "residue (a residue property this time, NOT the window anchor)", y = NULL) +
    lf_viz_theme() +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom", legend.justification = "left",
          legend.key.size = unit(9, "pt"),
          legend.text = element_text(colour = LF_VIZ$ink2, size = 8),
          legend.title = element_text(colour = LF_VIZ$ink2, size = 8))

  ## ---- one panel: the classes stacked --------------------------------------
  ## The per-index head is a softmax over the vocabulary, so at each residue the
  ## classes sum to 1 (checked: max |sum - 1| ~ 1e-8). A stacked area is therefore
  ## not a design choice laid over the numbers -- it IS the decomposition, and the
  ## band heights are the model's whole belief about that residue.
  ##
  ## The spread cannot ride along as a ribbon here (six overlapping ribbons is
  ## mud), so it goes underneath as a rug: the across-window sd of whichever class
  ## wins at that residue. A dark rug means the stack above it would look
  ## different through a different window.
  mk_combined <- function() {
    ## Same visual language as the protein page's detail panels
    ## (make_detail_panel in R/plot_proteins_win_new.R): every class overlaid as a
    ## line in the project's canonical class colours, its own +/- 1 sd as a
    ## translucent ribbon, a point per residue, legend in one row on top.
    ##
    ## NOT a stacked area. Stacking makes each band's height readable only against
    ## a moving baseline, so a class is easy to read at the bottom and impossible
    ## in the middle -- and it implies the classes compete for one budget, which
    ## is true of the softmax but not of what you want to see: whether pep_pocket
    ## is high HERE, against its own axis.
    ## canonical order and palette, minus anything --drop-classes excludes
    ord <- setdiff(intersect(names(nn_class_cols), classes), drop_cls)
    d <- tracks %>% filter(!as.character(class) %in% drop_cls) %>%
      mutate(class = factor(as.character(class), levels = ord))
    ggplot(d, aes(residue, mean, colour = class)) +
      {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                  ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
      {if (!is.null(kn) && nrow(kn)) geom_vline(data = kn, inherit.aes = FALSE,
          aes(xintercept = anchor), colour = LF_VIZ$ink, linewidth = 0.5, linetype = "22") } +
      geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = class),
                  alpha = 0.20, colour = NA, na.rm = TRUE) +
      scale_fill_manual(values = nn_class_cols, guide = "none") +
      geom_line(linewidth = 0.6, na.rm = TRUE) +
      ## the per-residue dot the detail panels carry; dropped on a long precursor,
      ## where one dot per residue is a smear rather than a mark
      {if (L <= 300) geom_point(pch = 21, stroke = 0.6, size = 1.6,
                                fill = LF_VIZ$surface, na.rm = TRUE) } +
      scale_colour_manual(values = nn_class_cols) +
      scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                         breaks = scales::breaks_width(BRK_W)) +
      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      guides(colour = guide_legend(nrow = 1)) +
      labs(x = NULL, y = "per-residue class probability\n(band: across windows)") +
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = "top",
            legend.direction = "horizontal", legend.justification = "center",
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill = "white", colour = NA),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  }

  ## The window-level panels, in the same language as the residue panel above.
  ## The score and the insertion classes get a panel each: they are on the same
  ## axis but they are not comparable quantities -- the score is an uncalibrated
  ## sigmoid that ranks windows, the classes are a softmax that partitions one
  ## window's probability mass. Sharing a y axis invited reading one off the other.
  win_panel <- function(d, pal, ylab, legend = TRUE) {
    ggplot(d, aes(anchor, mean, colour = class)) +
      {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                  ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
      {if (!is.null(kn) && nrow(kn)) geom_vline(data = kn, inherit.aes = FALSE,
          aes(xintercept = anchor), colour = LF_VIZ$ink, linewidth = 0.5, linetype = "22") } +
      geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = class),
                  alpha = 0.20, colour = NA, na.rm = TRUE) +
      scale_fill_manual(values = pal, guide = "none") +
      geom_line(linewidth = 0.6, na.rm = TRUE) +
      {if (L <= 300) geom_point(pch = 21, stroke = 0.6, size = 1.6,
                                fill = LF_VIZ$surface, na.rm = TRUE) } +
      scale_colour_manual(values = pal, guide = if (legend) "legend" else "none") +
      scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                         breaks = scales::breaks_width(BRK_W)) +
      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      guides(colour = guide_legend(nrow = 1)) +
      labs(x = NULL, y = ylab) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.position = if (legend) "top" else "none",
            legend.direction = "horizontal", legend.justification = "center",
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill = "white", colour = NA),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  }

  ## A thin band of what the labels SAY, directly under the predicted tracks.
  mk_truth <- function() {
    ggplot() +
      {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                  ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
      geom_rect(data = truth_seg,
                aes(xmin = start - 0.5, xmax = end + 0.5, ymin = 0.62, ymax = 1.18, fill = class),
                colour = "white", linewidth = 0.15) +
      scale_fill_manual(values = nn_class_cols, name = NULL,
                        breaks = levels(droplevels(truth_seg$class))) +
      geom_segment(data = truth_pep, aes(x = start, xend = end, y = 1.72, yend = 1.72),
                   colour = LF_VIZ$ink2, linewidth = 1.8, lineend = "butt") +
      geom_text(data = truth_pep, aes(x = (start + end) / 2, y = 2.12,
                                      label = sprintf("%s  [%s]", name, ins_class)),
                colour = LF_VIZ$ink2, size = 2.7) +
      geom_point(data = truth_end, aes(x = x, y = 1.72), colour = LF_VIZ$ink,
                 size = 1.9, shape = 18) +
      geom_text(data = truth_end, aes(x = x, y = 0.18, label = lab),
                colour = LF_VIZ$ink, size = 2.3, lineheight = 0.95, vjust = 0) +
      scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                         breaks = scales::breaks_width(BRK_W)) +
      scale_y_continuous(limits = c(0, 2.45), expand = c(0, 0)) +
      guides(fill = guide_legend(nrow = 1, keyheight = unit(7, "pt"))) +
      labs(x = NULL, y = "ground truth") +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), legend.position = "right",
            legend.text = element_text(size = 7), legend.key.size = unit(8, "pt"),
            legend.margin = margin(0, 0, 0, 0))
  }

  ## The correctly-framed window's own prediction, in the same language as the
  ## pooled panel above it. Legend suppressed: that panel already carries it.
  mk_correct <- function(d) {
    ggplot(d, aes(residue, mean, colour = class)) +
      {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                  ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
      {if (!is.null(kn) && nrow(kn)) geom_vline(data = kn, inherit.aes = FALSE,
          aes(xintercept = anchor), colour = LF_VIZ$ink, linewidth = 0.5, linetype = "22") } +
      geom_ribbon(aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = class),
                  alpha = 0.20, colour = NA, na.rm = TRUE) +
      scale_fill_manual(values = nn_class_cols, guide = "none") +
      geom_line(linewidth = 0.6, na.rm = TRUE) +
      {if (L <= 300) geom_point(pch = 21, stroke = 0.6, size = 1.6,
                                fill = LF_VIZ$surface, na.rm = TRUE) } +
      scale_colour_manual(values = nn_class_cols, guide = "none") +
      scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                         breaks = scales::breaks_width(BRK_W)) +
      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
      labs(x = NULL, y = d$ttl[1]) +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 7, lineheight = 0.95),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  }

  ## No window and no convolution: one residue at a time. A panel per model
  ## family rather than both in one -- the two disagree enough (xgboost leads on
  ## PR AUC, the MLP on ROC AUC) that overlaying them on a shared axis made the
  ## disagreement look like noise on one curve. Colour still matches the
  ## per-index panel's classes, so a colour means the same thing on every panel.
  mk_resid <- function(mdl, legend = TRUE) {
    d <- rp %>% filter(model == mdl)
    ggplot(d, aes(residue, mean, colour = class)) +
      {if (p$n_prot > 1) annotate("rect", xmin = 0.5, xmax = p$n_prot - 0.5,
                                  ymin = -Inf, ymax = Inf, fill = LF_VIZ$grid, alpha = 0.55) } +
      {if (!is.null(kn) && nrow(kn)) geom_vline(data = kn, inherit.aes = FALSE,
          aes(xintercept = anchor), colour = LF_VIZ$ink, linewidth = 0.5, linetype = "22") } +
      {if (any(!is.na(d$sd))) geom_ribbon(
          aes(ymin = pmax(mean - sd, 0), ymax = pmin(mean + sd, 1), fill = class),
          alpha = 0.18, colour = NA, na.rm = TRUE) } +
      scale_fill_manual(values = nn_class_cols, guide = "none") +
      geom_line(linewidth = 0.6, na.rm = TRUE) +
      scale_colour_manual(values = nn_class_cols, name = NULL,
                          guide = if (legend) "legend" else "none") +
      scale_x_continuous(limits = c(0.5, L + 0.5), expand = c(0, 0),
                         breaks = scales::breaks_width(BRK_W)) +
      scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.5), expand = c(0, 0)) +
      guides(colour = guide_legend(nrow = 1)) +
      labs(x = NULL, y = sprintf("%s, per residue\n(no window; band: seeds)",
                                 toupper(as.character(mdl)))) +
      theme_bw() +
      theme(legend.position = if (legend) "top" else "none",
            legend.direction = "horizontal", legend.justification = "center",
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill = "white", colour = NA),
            legend.text = element_text(size = 7), legend.key.size = unit(8, "pt"),
            axis.title.y = element_text(size = 7.5, lineheight = 0.95),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  }
  ## which model families actually have predictions here
  rp_models <- if (!is.null(rp) && nrow(rp)) levels(droplevels(rp$model)) else character(0)

  body <- if (!combined) lapply(classes, mk_track) else {
    ## one series needs no legend; its y axis already names it
    b <- list(mk_combined())
    if (!is.null(cw) && length(cw)) b <- c(b, lapply(cw, mk_correct))
    if (!is.null(truth_seg) && nrow(truth_seg)) b <- c(b, list(mk_truth()))
    if (length(rp_models))
      b <- c(b, lapply(seq_along(rp_models),
                       \(i) mk_resid(rp_models[i], legend = i == 1L)))
    b <- c(b, list(win_panel(win %>% filter(class == "window score"), WIN_PAL,
                             "CNN attention head\n(band: across seeds)", legend = FALSE)))
    ## One panel per xgboost head, in the CNN track's colour so the eye reads
    ## them as the same quantity -- they are, computed three ways.
    if (!is.null(wx) && nrow(wx)) {
      .lab <- c(xgb_pm = "xgb window head\n(band: across seeds)",
                xgb_st = "xgb head, stacked\n(band: across seeds)")
      for (h in levels(droplevels(wx$head)))
        b <- c(b, list(win_panel(wx %>% filter(head == h) %>%
                                   mutate(class = factor("window score",
                                                         levels = names(WIN_PAL))),
                                 WIN_PAL, .lab[[h]], legend = FALSE)))
    }
    if (!is.null(INS))
      b <- c(b, list(win_panel(win %>% filter(class != "window score"), WIN_PAL,
                               "insertion class\n(band: across seeds)")))
    b
  }
  fig <- wrap_plots(c(list(ttl), body, list(p_seq)), ncol = 1) +
    plot_layout(heights = c(1.15,
                            ## one entry per panel mk_* actually produced, in the
                            ## same order -- patchwork errors rather than guessing
                            ## if these drift apart
                            if (combined) c(6,
                                            if (!is.null(cw) && length(cw)) rep(2.6, length(cw)) else NULL,
                                            if (!is.null(truth_seg) && nrow(truth_seg)) 1.5 else NULL,
                                            if (length(rp_models)) rep(2.9, length(rp_models)) else NULL,
                                            3.2,
                                            if (!is.null(wx) && nrow(wx))
                                              rep(2.6, dplyr::n_distinct(wx$head)) else NULL,
                                            if (!is.null(INS)) 4.2 else NULL)
                            else rep(if (length(classes) > 3) 2.8 else 4, length(classes)),
                            1.5)) +
    plot_annotation(theme = theme(plot.background =
                                    element_rect(fill = LF_VIZ$surface, colour = NA)))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  st <- file.path(out_dir, sprintf("lf_resid_tracks_%s%s_%s",
                                   gene, if (combined) "_combined" else "", term))
  H <- if (combined) (if (is.null(INS)) 8.6 else 12.2) +
         (if (!is.null(truth_seg) && nrow(truth_seg)) 1.1 else 0) +
         (if (!is.null(cw)) 1.5 * length(cw) else 0) +
         (if (length(rp_models)) 1.7 * length(rp_models) else 0) +
         (if (!is.null(wx) && nrow(wx)) 1.5 * dplyr::n_distinct(wx$head) else 0) else
    1.6 + (if (length(classes) > 3) 1.45 else 2.1) * length(classes) + 1.3
  ggsave(paste0(st, ".svg"), fig, width = PLOT_W, height = H, device = svglite::svglite)
  if (!no_png) ggsave(paste0(st, ".png"), fig, width = PLOT_W, height = H, dpi = 200,
                      bg = LF_VIZ$surface)
  ## Both panels' data, with `axis` saying which coordinate the row is on and
  ## `sd_over` what its sd is across -- the two are not interchangeable and a bare
  ## merge of them would silently imply they are.
  out_tbl <- bind_rows(
    tracks %>% mutate(axis = "residue", sd_over = "windows") %>%
      rename(x = residue),
    if (combined) win %>% mutate(axis = "anchor", sd_over = "seeds",
                                 n_windows = NA_integer_, sd_member = NA_real_,
                                 min = NA_real_, max = NA_real_) %>%
      rename(x = anchor) else NULL
  ) %>% select(axis, x, class, mean, sd, sd_over, everything())
  readr::write_csv(out_tbl, paste0(st, ".csv"))


  ## one row per gene, in the shape 10_9k's gallery already reads
  bs <- win %>% filter(class == "window score")
  pks <- if (!is.null(sc)) sc$peaks %>% filter(accession == p$accession) else NULL
  knm <- if (!is.null(kn) && nrow(kn))
    kn %>% select(any_of(c("pep_name", "anchor", "score", "rank_all", "split", "motif")))
  else tibble(pep_name = character(), anchor = integer(), score = numeric(),
              rank_all = integer(), split = character(), motif = character())
  message(sprintf("  %-8s %4d res, %3d windows, best %.3f @ %d, %d known",
                  gene, L, length(anchors), max(bs$mean), bs$anchor[which.max(bs$mean)],
                  nrow(knm)))
  tibble(gene = gene, accession = p$accession,
         ## derived, never a flag: a gene with no known end is `other`, and one
         ## with knowns is named by the split(s) they are actually in
         group = if (!nrow(knm)) "other"
                 else paste(sort(unique(as.character(knm$split))), collapse = "/"),
         term = term, n_res = L, n_scanned = length(anchors),
         n_peaks = if (is.null(pks)) 0L else nrow(pks),
         best_score = max(bs$mean), best_anchor = bs$anchor[which.max(bs$mean)],
         width_in = PLOT_W, svg = paste0(st, ".svg"), n_known = nrow(knm),
         knowns = list(knm))
}

meta <- bind_rows(lapply(genes, function(g)
  tryCatch(make_one(g), error = function(e) {
    message("  ", g, ": skipped -- ", conditionMessage(e)); NULL })))
if (!nrow(meta)) stop("nothing rendered", call. = FALSE)
## merge, so one gene can be re-rendered without dropping the rest
old <- if (file.exists(meta_p)) readRDS(meta_p) else NULL
if (!is.null(old)) meta <- bind_rows(old %>% filter(!gene %in% meta$gene), meta)
saveRDS(meta, meta_p)
message(sprintf("\n%d figure(s) this call; %d in %s",
                length(genes), nrow(meta), basename(meta_p)))
