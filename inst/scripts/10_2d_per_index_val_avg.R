#!/usr/bin/env Rscript
## ---- per-residue profiles of the VALIDATION windows, member-averaged ----------
## 10_2c plots every member for one window; this plots every validation positive
## (plus the negatives the arm scores highest) with the per-residue softmax
## averaged over members, two ways:
##
##   all    the plain ensemble mean over every member
##   top<k> the mean over the k members that score THAT window highest -- "what
##          do the seeds that like it see?", per window. Selection is on the
##          window score, so on a window every member likes it barely differs
##          from `all`; it moves where the members disagree. On the windows every
##          member dismisses, the top k is just the k members with the highest
##          score floor -- nearly the same seeds for every such window.
##
## Drawn as the protein pages' detail panels are (make_detail_panel): the mean
## as a plain line with a ribbon at +/- 1 sd across the averaged members, no
## loess. `--band se` makes it sd / sqrt(n members averaged) instead.
##
## Reads the per-member arrays the isolated trainer leaves behind
## (<arm>/members/member_NNN_<term>.npz), so it needs no retraining.
##
##   /usr/local/bin/Rscript inst/scripts/10_2d_per_index_val_avg.R
##   ... --arms unet_none_w1 --top 5 --n-neg 10 --band se
##
## Output, one per arm and averaging:
##   ~/AF2_analysis/per_index_val_<arm>_<all|top<k>>_<sd|se>_<stamp>.svg
## Labels are ASCII: under Rscript's non-UTF-8 locale svglite writes "·" as "..".
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
cache_dir <- path.expand(.opt("--cache-dir", "~/AF2_analysis/lf_dcnn_bench_none20"))
arms      <- strsplit(.opt("--arms", "unet_none_w1,unet_base"), ",")[[1]]
term      <- .opt("--term", "C")
top_k     <- as.integer(.opt("--top", "5"))
n_neg     <- as.integer(.opt("--n-neg", "5"))       # highest-scoring val negatives to add
band      <- .opt("--band", "sd")                    # sd | se
stopifnot(band %in% c("sd", "se"))
nn_cache  <- path.expand(.opt("--nn-input", "~/AF2_analysis/lf_dcnn_compare_nn_input.rds"))
out_dir   <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
base_seed <- 42L                                    # member k (0-based) trained with seed 42 + k
stamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
invisible(lf_dcnn_python(path = file.path(ROOT, "inst", "python")))   # binds the venv
np <- reticulate::import("numpy", convert = FALSE)
pd <- reticulate::import("pandas", convert = FALSE)
.npz <- function(f, keys = NULL) {
  z <- np$load(f)
  keys <- keys %||% as.character(reticulate::py_to_r(z$files))
  stats::setNames(lapply(keys, function(k) reticulate::py_to_r(z$`__getitem__`(k))), keys)
}

## ---- the windows: val order is the order val_scores / val_labels are in --------
.cc  <- readRDS(nn_cache)
val  <- .cc$nn_input[[term]]$val
all_peps <- .cc$nn_input[[term]]$all$peps

class_cols <- c(CT_cleavage_context = "#FED439FF", DB = "#370335FF", gap = "#8A9197FF",
                NT_cleavage_context = "#D2AF81FF", pep_other = "#D5E4A2FF",
                pep_pocket = "#197EC0FF", padding = "grey85", none = "#075149FF")

for (arm_dir in arms) {
  arm_path <- file.path(cache_dir, arm_dir)
  out <- .npz(file.path(arm_path, "outputs.npz"))
  stopifnot(length(out$val_scores) == nrow(val), length(out$pred_raw) == length(all_peps),
            identical(as.integer(out$val_labels), as.integer(val$known)))
  pq <- file.path(arm_path, "predictions.parquet")
  if (file.exists(pq))
    stopifnot(identical(as.character(reticulate::py_to_r(pd$read_parquet(pq, columns = list("peps"))$peps$tolist())),
                        all_peps))

  ## positives by ensemble score, then the arm's n_neg most confident negatives
  v <- tibble(peps = val$peps, target = val$target, known = as.integer(val$known),
              score = as.numeric(out$val_scores)) %>%
    mutate(row = match(peps, all_peps)) %>%
    arrange(desc(known), desc(score)) %>%
    group_by(known) %>% filter(known == 1 | row_number() <= n_neg) %>% ungroup()

  scores <- out$pred_raw_members[, v$row, drop = FALSE]        # (n_members, n_windows)
  n_mem  <- nrow(scores)
  stopifnot(top_k <= n_mem)
  ## the members' val scores are the same windows scored the same way
  stopifnot(isTRUE(all.equal(scores, out$val_scores_members[, match(v$peps, val$peps), drop = FALSE],
                             tolerance = 1e-5)))

  mfiles <- sort(Sys.glob(file.path(arm_path, "members", sprintf("member_*_%s.npz", term))))
  stopifnot(length(mfiles) == n_mem)
  pi_names <- as.character(unlist(jsonlite::read_json(file.path(arm_path, "history.json"))$pi_names))
  ## (n_members, n_windows, seq_len, K_pi)
  pi_key <- paste0(term, "__per_index__x")
  pi_m <- simplify2array(lapply(mfiles, function(f)
    .npz(f, pi_key)[[pi_key]][v$row, , , drop = FALSE]))
  pi_m <- aperm(pi_m, c(4, 1, 2, 3))
  dimnames(pi_m)[[4]] <- pi_names

  ## ---- the sequence / true-label track -----------------------------------------
  trk <- bind_rows(lapply(seq_len(nrow(v)), function(i) {
    iv <- match(v$peps[[i]], val$peps)
    tibble(win = i, index = seq_along(val$known_idx[[iv]]),
           AA = as.character(val$meta_data[[iv]]$AA),
           known_idx = as.character(val$known_idx[[iv]]))
  }))

  for (mode in c("all", paste0("top", top_k))) {
    ## which members each window averages over
    pick <- lapply(seq_len(nrow(v)), function(i)
      if (mode == "all") seq_len(n_mem) else order(-scores[, i])[seq_len(top_k)])

    ## sample sd, as the pipeline's per_index_sd; se divides by sqrt(members averaged)
    .spread <- function(x) sd(x) / if (band == "se") sqrt(length(x)) else 1
    long <- bind_rows(lapply(seq_len(nrow(v)), function(i) {
      x  <- pi_m[pick[[i]], i, , , drop = FALSE]
      mu <- apply(x, c(3, 4), mean)                                   # (seq_len, K_pi)
      tibble(win = i, index = rep(seq_len(nrow(mu)), ncol(mu)),
             name = rep(pi_names, each = nrow(mu)), value = as.vector(mu),
             spread = as.vector(apply(x, c(3, 4), .spread)))
    }))
    lab <- v %>% mutate(
      win   = row_number(),
      avg   = vapply(win, function(i) mean(scores[pick[[i]], i]), numeric(1)),
      avg_s = vapply(win, function(i) .spread(scores[pick[[i]], i]), numeric(1)),
      n_hi  = colSums(scores > 0.5),
      seeds = vapply(pick, function(p) paste(base_seed + p - 1L, collapse = ","), character(1)),
      panel = if (mode == "all")
                sprintf("%s %s | %s\nensemble %.2f +/- %.2f | %d/%d members > .5",
                        peps, target, ifelse(known == 1, "pos", "neg"), avg, avg_s, n_hi, n_mem)
              else
                sprintf("%s %s | %s\ntop %d %.2f +/- %.2f (all %.2f) | seeds %s",
                        peps, target, ifelse(known == 1, "pos", "neg"), top_k, avg, avg_s, score, seeds))
    lab$panel <- factor(lab$panel, levels = lab$panel)

    long <- long %>% left_join(select(lab, win, panel), by = "win") %>%
      mutate(name = factor(name, levels = names(class_cols)))
    trk_p <- trk %>% left_join(select(lab, win, panel), by = "win") %>% mutate(value = -0.12)

    avg_desc <- if (mode == "all") sprintf("mean over all %d members", n_mem) else
      sprintf("mean over the %d members scoring each window highest", top_k)
    p <- ggplot(long, aes(index, value, colour = name)) +
      ## the class palette is the track's too, so the ribbons share its fill scale
      geom_ribbon(aes(ymin = pmax(value - spread, 0), ymax = pmin(value + spread, 1), fill = name),
                  alpha = 0.20, colour = NA, show.legend = FALSE) +
      geom_line(linewidth = 0.6) +
      geom_point(pch = 21, stroke = 0.8, size = 1.3) +
      geom_tile(data = trk_p, mapping = aes(x = index, y = value, fill = known_idx),
                height = 0.07, inherit.aes = FALSE) +
      geom_text(data = filter(trk_p, !is.na(AA)), aes(x = index, y = value, label = AA), inherit.aes = FALSE,
                size = 1.6, fontface = "bold", colour = "black") +
      scale_colour_manual(values = class_cols, name = NULL, drop = FALSE) +
      scale_fill_manual(values = class_cols, name = "true label", drop = FALSE) +
      ## one scale in every file, so top<k> and all compare panel for panel
      scale_y_continuous(breaks = seq(0, 1, 0.25)) +
      coord_cartesian(ylim = c(-0.15, 1.05)) +
      facet_wrap(vars(panel), ncol = 4) +
      labs(title = sprintf("Validation windows, %s terminus: per-residue predictions, %s", term, avg_desc),
           subtitle = sprintf("%s / %s: the %d validation positives, then the %d highest-scoring of %d negatives. Ribbon and strip +/- = %s across the averaged members. The `none` curve is dark teal.",
                              basename(cache_dir), arm_dir, sum(v$known), sum(v$known == 0),
                              sum(val$known == 0),
                              if (band == "sd") "1 sd" else "1 se (sd / sqrt(n))"),
           x = "position in window", y = "predicted probability") +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            strip.text = element_text(size = 7.5))

    f <- file.path(out_dir, sprintf("per_index_val_%s_%s_%s_%s.svg", arm_dir, mode, band, stamp))
    ggsave(f, p, width = 4 * 3.4 + 1, height = 3.1 * ceiling(nrow(v) / 4) + 1.6, limitsize = FALSE)
    message("wrote ", f)
  }
}
