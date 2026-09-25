#!/usr/bin/env Rscript
## ---- the 10_2 analyses on a benchmark arm -------------------------------------
## Rebuilds the 10_1 session contract for a model trained by
## 10_5_benchmark_window_model.R -- the arm's first member from its saved
## weights, the windows with the arm's channels attached, the ensemble scores
## and per-residue softmax from outputs.npz -- and then runs, per --what:
##
##   saliency   10_2_model_importance.R: vanilla-gradient saliency of the global
##              logit on the validation positives, input profiles, sequence track
##              + NJ tree; plm_* channels as extra groups of 8, most salient
##              first; plus a per-channel ranking of mean |saliency|
##   per_index  10_2_per_ind_profiles.R: per-residue class curves for the known
##              peptides (ensemble MEAN softmax) over the sequence track, and the
##              masked per-index accuracy / confusion on validation
##
## Default: the `t5` preset's ProtT5 arm, C terminus, both analyses.
##
##   /usr/local/bin/Rscript inst/scripts/10_2b_bench_arm_session.R
##   /usr/local/bin/Rscript inst/scripts/10_2b_bench_arm_session.R --what per_index
##   /usr/local/bin/Rscript inst/scripts/10_2b_bench_arm_session.R --arm unet_base --input base
##
## Output under ~/AF2_analysis, stamped:
##   <out-prefix>_<group>[_<term>]_<stamp>.svg          (saliency, as 10_2)
##   <out-prefix>_channel_ranking_<stamp>.{svg,png,csv}
##   <out-prefix>_per_index_knowns_<term>_<stamp>.svg   (per_index)
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
cache_dir   <- path.expand(.opt("--cache-dir", "~/AF2_analysis/lf_dcnn_bench_t5"))
arm_dir     <- .opt("--arm", "unet_t5")            # <cache-dir>/<arm>: outputs + members
input_name  <- .opt("--input", "t5")               # base | t5 | esm_c -- which channels the arm trained on
term        <- .opt("--term", "C")
member      <- as.integer(.opt("--member", "0"))   # only member 0 has saved weights
nn_cache    <- path.expand(.opt("--nn-input", "~/AF2_analysis/lf_dcnn_compare_nn_input.rds"))
plm_parquet <- path.expand(.opt("--plm", switch(input_name,
                                                t5    = "~/AF2_analysis/lf_plm/prot_t5_pca32.parquet",
                                                esm_c = "~/AF2_analysis/lf_plm/esm_c_pca32.parquet",
                                                "")))
seq_parquet <- path.expand(.opt("--sequences", "~/AF2_analysis/lf_plm/sequences.parquet"))
out_prefix  <- .opt("--out-prefix", paste0("model_importance_", arm_dir, "_", term))
what        <- strsplit(.opt("--what", "saliency,per_index"), ",")[[1]]
out_dir     <- path.expand("~/AF2_analysis")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
if (!exists("lf_dcnn_run"))   source(file.path(ROOT, "R", "dcnn_bridge.R"))
if (!exists("lf_plm_attach")) source(file.path(ROOT, "R", "plm_features.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))

## ---- windows + channels, exactly as the arm saw them ---------------------------
.cc <- readRDS(nn_cache)
nn_input    <- .cc$nn_input[term]
all_params3 <- .cc$all_params3
if (input_name != "base") {
  built <- lf_plm_attach(nn_input, lf_plm_read(plm_parquet), lf_read_parquet(seq_parquet), all_params3)
  nn_input <- built$nn_input; all_params3 <- built$channels
}
## an arm whose Config changes the LABEL arrays (e.g. include_none = FALSE)
## has its own export directory; --in-dir names it
in_dir  <- .opt("--in-dir", "")
in_dir  <- if (nzchar(in_dir)) path.expand(in_dir) else
           file.path(cache_dir, if (input_name == "base") "in" else paste0("in_", input_name))
cfg     <- mod$Config$from_json(file.path(in_dir, "config.json"))
stopifnot(identical(as.character(cfg$channel_names), all_params3),
          identical(as.character(cfg$term_order), names(nn_input)))
seq_len    <- as.integer(cfg$seq_len)
n_channels <- as.integer(cfg$n_channels)

## ---- the model: member 0 from its weights -------------------------------------
w <- file.path(cache_dir, arm_dir, "members", sprintf("member_%03d_%s.weights.h5", member, term))
if (!file.exists(w)) stop("no saved weights at ", w, " (only member 0 saves them)")
model <- mod$build_model(cfg); model$load_weights(normalizePath(w))
models <- stats::setNames(list(model), term)
message(sprintf("model: %s, %s params, %d channels", basename(w), format(model$count_params(), big.mark = ","), n_channels))

## ---- the 10_1 session contract --------------------------------------------------
out <- lf_dcnn_import(file.path(cache_dir, arm_dir))
nn_input_comb <- nn_input[[term]]$all
stopifnot(length(out$pred) == nrow(nn_input_comb))
nn_input_comb$pred     <- as.numeric(out$pred)
nn_input_comb$pred_raw <- as.numeric(out$pred_raw)
nn_input_comb$pred_sd  <- as.numeric(out$pred_sd)
## per-residue softmax: the ensemble MEAN (and sd) across the arm's members
pi_names_out <- as.character(unlist(out$pi_names))
pi_names     <- pi_names_out          # 10_2_per_ind_profiles.R reads this
nn_input_comb$per_index    <- lf_dcnn_per_index_tibbles(out$per_index,    pi_names_out)
nn_input_comb$per_index_sd <- lf_dcnn_per_index_tibbles(out$per_index_sd, pi_names_out)
known_dat <- nn_input_comb %>% filter(known == 1)     # 10_2 profiles these...
c_dat     <- nn_input_comb %>% filter(known == 0)     # ...against these
classes   <- stats::setNames(seq_along(cfg$class_names) - 1L, as.character(cfg$class_names))
real_cols <- as.integer(cfg$real_cols) + 1L
class_cols <- stats::setNames(c("#FED439FF", "#370335FF", "#8A9197FF", "#D2AF81FF",
                                "#D5E4A2FF", "#197EC0FF", "grey85", "#075149FF"),
                              c("CT_cleavage_context", "DB", "gap", "NT_cleavage_context",
                                "pep_other", "pep_pocket", "padding", "none"))
targs      <- list(c(term, paste0("loop_", term)))
pi_terms   <- term
pi_prefix  <- paste0(out_prefix, "_per_index_knowns_")
extra_peps <- character(0)
plm_cols   <- grep("^plm_", all_params3, value = TRUE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
suppressMessages({ library(keras3); library(tensorflow) })

if ("saliency" %in% what) {
## ---- per-channel ranking: mean |d logit / d x| over val positives and positions --
build_x0 <- function(dataset) {
  n <- length(dataset[["data"]])
  aperm(array(unlist(dataset[["data"]]), dim = c(seq_len, n_channels, n)), c(3, 1, 2))
}
val_x   <- build_x0(nn_input[[term]]$val)
val_pos <- which(nn_input[[term]]$val$known == 1)
g         <- keras3::get_layer(model, "global")
pen_model <- keras3::keras_model(model$input, g$input)
W <- tf$convert_to_tensor(g$kernel); b <- tf$convert_to_tensor(g$bias)
xt <- tf$convert_to_tensor(val_x, dtype = tf$float32)
with(tf$GradientTape() %as% tape, { tape$watch(xt); logit <- tf$matmul(pen_model(xt), W) + b })
grads <- as.array(tape$gradient(logit, xt))                       # (n, seq_len, C)
rank_tbl <- tibble(channel = all_params3,
                   group   = if_else(grepl("^plm_", all_params3), "PLM PCA", "hand-built"),
                   mean_abs_grad_pos = apply(abs(grads[val_pos, , , drop = FALSE]), 3, mean),
                   mean_abs_grad_all = apply(abs(grads), 3, mean)) %>%
  arrange(desc(mean_abs_grad_pos)) %>% mutate(rank = row_number())
write.csv(rank_tbl, file.path(out_dir, sprintf("%s_channel_ranking_%s.csv", out_prefix, stamp)), row.names = FALSE)
message("channel ranking (mean |d logit/dx| over the ", length(val_pos), " validation positives):")
print(as.data.frame(head(rank_tbl, 15)), row.names = FALSE, digits = 3)
message(sprintf("share of total |gradient| on PLM channels: %.1f%% (they are %.0f%% of channels)",
                100 * sum(rank_tbl$mean_abs_grad_pos[rank_tbl$group == "PLM PCA"]) / sum(rank_tbl$mean_abs_grad_pos),
                100 * mean(rank_tbl$group == "PLM PCA")))
p_rank <- ggplot(rank_tbl, aes(x = reorder(channel, mean_abs_grad_pos), y = mean_abs_grad_pos, fill = group)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("hand-built" = "grey55", "PLM PCA" = "#1B9E77"), name = NULL) +
  labs(x = NULL, y = "mean |d logit / d x|  (validation positives, all positions)",
       title = sprintf("Channel saliency, %s arm, %s terminus (member %d)", arm_dir, term, member)) +
  theme_bw(base_size = 10) + theme(legend.position = "top")
for (ext in c("svg", "png"))
  ggsave(file.path(out_dir, sprintf("%s_channel_ranking_%s.%s", out_prefix, stamp, ext)), p_rank,
         width = 7, height = 2 + 0.17 * n_channels, dpi = 150)

## PLM channels as extra 10_2 groups: 8 per figure, most salient first, so the
## components that move the logit are in plm_a and the tail in plm_d
if (length(plm_cols)) {
  ordered <- rank_tbl$channel[rank_tbl$channel %in% plm_cols]
  met_it_extra <- split(ordered, ceiling(seq_along(ordered) / 8)) %>%
    stats::setNames(paste0("plm_", letters[seq_along(.)]))
}

## ---- the 10_2 figures ------------------------------------------------------------
source(file.path(ROOT, "inst/scripts/10_2_model_importance.R"))
}

if ("per_index" %in% what) {
  ## per-residue class curves for the knowns + val accuracy / confusion (uses
  ## `models` for the diagnostic and nn_input_comb$per_index for the curves)
  source(file.path(ROOT, "inst/scripts/10_2_per_ind_profiles.R"))
}
message("done: ~/AF2_analysis/", out_prefix, "_*_", stamp, "*")
