#!/usr/bin/env Rscript
## ---- per-member per-residue profiles for ONE window ----------------------------
## The ensemble-mean per-index plot (10_2_per_ind_profiles.R) hides exactly what
## you want when the members disagree: a window whose ensemble score is 0.31
## because a quarter of the members love it and the rest do not looks identical
## to one every member scores 0.31. This plots each member separately, ordered
## by that member's window score, so "what do the seeds that like it actually
## see?" is answerable.
##
## Reads the per-member arrays the isolated trainer leaves behind
## (<arm>/members/member_NNN_<term>.npz, which carry the full per_index softmax),
## so it needs no retraining and no `models` in the session.
##
##   /usr/local/bin/Rscript inst/scripts/10_2c_per_index_members.R \
##      --cache-dir ~/AF2_analysis/lf_dcnn_bench_none20 --arm unet_none_w1 \
##      --peps ANO8_w12-47
##   ... --members good     # only the members scoring the window above --threshold
##
## Output: ~/AF2_analysis/<out-prefix>_<peps>_<stamp>.svg
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
cache_dir <- path.expand(.opt("--cache-dir", "~/AF2_analysis/lf_dcnn_bench_none20"))
arm_dir   <- .opt("--arm", "unet_none_w1")
term      <- .opt("--term", "C")
peps      <- .opt("--peps", "ANO8_w12-47")
which_m   <- .opt("--members", "all")            # all | good
thresh    <- as.numeric(.opt("--threshold", "0.5"))
nn_cache  <- path.expand(.opt("--nn-input", "~/AF2_analysis/lf_dcnn_compare_nn_input.rds"))
out_dir   <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
out_prefix <- .opt("--out-prefix", paste0("per_index_members_", arm_dir))
stamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)

## ---- locate the window ----------------------------------------------------------
pred <- arrow_free_parquet <- NULL
pq <- file.path(cache_dir, arm_dir, "predictions.parquet")
ref_pq <- if (file.exists(pq)) pq else file.path(cache_dir, "unet_base", "predictions.parquet")
pd   <- reticulate::import("pandas", convert = FALSE)
ref  <- reticulate::py_to_r(pd$read_parquet(ref_pq))
j    <- which(ref$peps == peps)
if (!length(j)) stop(peps, " not found in ", ref_pq)
j <- j[[1]]
message(sprintf("%s: row %d of %s", peps, j, basename(dirname(ref_pq))))

## per-member window score, to order and label the panels
out <- reticulate::py_to_r(np$load(file.path(cache_dir, arm_dir, "outputs.npz")))
scores <- as.numeric(out[["pred_raw_members"]][, j])
n_mem  <- length(scores)

## ---- per-member per-residue softmax ---------------------------------------------
mfiles <- sort(Sys.glob(file.path(cache_dir, arm_dir, "members", sprintf("member_*_%s.npz", term))))
stopifnot(length(mfiles) == n_mem)
pi_names <- as.character(unlist(jsonlite::read_json(file.path(cache_dir, arm_dir, "history.json"))$pi_names))
prof <- bind_rows(lapply(seq_along(mfiles), function(k) {
  z <- reticulate::py_to_r(np$load(mfiles[[k]]))
  m <- z[[paste0(term, "__per_index__x")]][j, , ]            # (seq_len, K_pi)
  colnames(m) <- pi_names
  as_tibble(m) %>% mutate(index = row_number(), member = k, seed = 41L + k,
                          score = scores[[k]])
}))

keep <- if (identical(which_m, "good")) which(scores > thresh) else seq_len(n_mem)
message(sprintf("%d/%d members plotted (%d score above %.2f: seeds %s)",
                length(keep), n_mem, sum(scores > thresh), thresh,
                paste(41L + which(scores > thresh), collapse = ", ")))
prof <- prof %>% filter(member %in% keep) %>%
  mutate(panel = sprintf("seed %d · %.2f%s", seed, score, ifelse(score > thresh, "  ✓", "")),
         panel = factor(panel, levels = unique(panel[order(-score)])))

## ---- the sequence / true-label track, as 10_2_per_ind_profiles.R draws it --------
.cc <- readRDS(nn_cache)
w   <- .cc$nn_input[[term]]$all
i_w <- which(w$peps == peps)
md  <- w$meta_data[[i_w]]
trk <- tibble(index = seq_len(nrow(md)), AA = as.character(md$AA),
              known_idx = as.character(w$known_idx[[i_w]]))

class_cols <- c(CT_cleavage_context = "#FED439FF", DB = "#370335FF", gap = "#8A9197FF",
                NT_cleavage_context = "#D2AF81FF", pep_other = "#D5E4A2FF",
                pep_pocket = "#197EC0FF", padding = "grey85", none = "#075149FF")

long <- prof %>% select(index, panel, all_of(pi_names)) %>%
  pivot_longer(all_of(pi_names), names_to = "name", values_to = "value") %>%
  mutate(name = factor(name, levels = names(class_cols)))
trk_p <- tidyr::crossing(trk, panel = levels(long$panel)) %>%
  mutate(panel = factor(panel, levels = levels(long$panel)), value = -0.12)

p <- ggplot(long, aes(index, value, colour = name)) +
  geom_point(pch = 21, stroke = 0.8, size = 1.3) +
  geom_smooth(method = "loess", span = 0.6, se = FALSE, linewidth = 0.7) +
  geom_tile(data = trk_p, mapping = aes(x = index, y = value, fill = known_idx),
            height = 0.07, inherit.aes = FALSE) +
  geom_text(data = trk_p, aes(x = index, y = value, label = AA), inherit.aes = FALSE,
            size = 1.6, fontface = "bold", colour = "black") +
  scale_colour_manual(values = class_cols, name = NULL, drop = FALSE) +
  scale_fill_manual(values = class_cols, name = "true label", drop = FALSE) +
  facet_wrap(vars(panel)) +
  labs(title = sprintf("%s — per-residue predictions by ensemble member", peps),
       subtitle = sprintf("%s / %s, %d members, ordered by that member's window score (✓ = above %.2f). The `none` curve is dark teal.",
                          basename(cache_dir), arm_dir, length(keep), thresh),
       x = "position in window", y = "predicted probability") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

nc <- ceiling(sqrt(length(keep)))
ggsave(file.path(out_dir, sprintf("%s_%s_%s.svg", out_prefix, gsub("[^A-Za-z0-9]", "", peps), stamp)),
       p, width = 3.4 * nc + 1, height = 2.9 * ceiling(length(keep) / nc) + 1.6, limitsize = FALSE)
message("wrote ", file.path(out_dir, sprintf("%s_%s_%s.svg", out_prefix, gsub("[^A-Za-z0-9]", "", peps), stamp)))
