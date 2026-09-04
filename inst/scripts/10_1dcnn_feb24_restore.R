## =============================================================================
## Restore the Feb-24 models (known-good) and score the CURRENT candidates.
##
## Why loading, not hand-rebuilding: the .keras files hold the exact architecture
## AND weights. The model is a functional net with 4 outputs (global, per_index,
## per_index2, per_index_logits); reconstructing that by hand is error-prone, so
## we load it verbatim and (optionally) clone_model() it to retrain.
##
## Data fix: the Feb24 nets expect 36 input channels = the all_params3 features
## only (NO end_type / position channels). The current nn_input carries 37
## channels (all_params3 + end_type_ch), so we drop end_type_ch to match.
##
## Run AFTER 9.2 (needs `nn_input`, `all_params3`, and `classes` in the session).
## =============================================================================

library(keras3)
library(tidyverse)

feb_dir <- path.expand("~/AF2_analysis")

## --- restore the models (architecture + weights) ---------------------------
## compile = FALSE so we don't need the original custom losses/metrics to load.
models_feb <- list(
  N = load_model(file.path(feb_dir, "model_N_Feb24.keras"), compile = FALSE),
  C = load_model(file.path(feb_dir, "model_C_Feb24.keras"), compile = FALSE)
)
message("restored Feb24 architecture (N):")
print(models_feb$N)                        # <- the restored 'code' (architecture summary)

## To RETRAIN this architecture on current data instead of using the old weights,
## clone it (fresh weights, same graph) -- but its 4 outputs need their losses
## wired up first; ask and I'll add that:
# feb_arch_N <- clone_model(models_feb$N)

## --- adapt current data to the Feb24 36-channel input ----------------------
if (!exists("seq_len")) seq_len <- nrow(nn_input[[1]]$all$data[[1]])

## indices of the 36 all_params3 feature columns within the current 37-col data
## (drops end_type_ch). match() is robust to column reordering.
feb_channels <- match(all_params3, names(nn_input[[1]]$all$data[[1]]))
stopifnot(!anyNA(feb_channels), length(feb_channels) == 36L)

build_x36 <- function(dataset) {
  nch <- ncol(dataset$data[[1]])
  arr <- aperm(array(unlist(dataset$data),
                     dim = c(seq_len, nch, length(dataset$data))), c(3, 1, 2))  # (n, seq, 37)
  arr[, , feb_channels, drop = FALSE]                                            # (n, seq, 36)
}

## --- score current candidates with the restored models ---------------------
feb_preds <- setNames(
  lapply(names(nn_input), function(term)
    predict(models_feb[[term]], build_x36(nn_input[[term]]$all), verbose = 0)),
  names(nn_input))

## global (ranking) score. Feb24 predict returns a named list incl 'global'.
glob <- function(p) as.numeric((if (is.list(p)) p[["global"]] else p)[, 1])

nn_input_comb_feb <- bind_rows(lapply(nn_input, function(x) x$all))
nn_input_comb_feb$pred <- do.call(c, lapply(feb_preds, glob)) %>% unname
nn_input_comb_feb <- nn_input_comb_feb %>%
  mutate(model = stringr::str_remove(as.character(target), "^loop_")) %>%
  arrange(desc(pred)) %>%
  mutate(rank = row_number(), .by = known)

## OPTIONAL: per-index tracks for plots. Feb24's `per_index` head is a 7-class
## softmax; this ASSUMES its class order matches names(classes)[c(real_cols, none)].
## Verify before trusting the labels.
## show every Feb24 output and its shape, so we pick the right one, not guess
out_dims <- sapply(feb_preds[[1]], function(x) if (is.array(x)) paste(dim(x), collapse = "x") else "?")
message("Feb24 output shapes: ",
        paste(sprintf("%s=%s", names(out_dims), out_dims), collapse = "   "))

## a per-POSITION head is 3-D (n, seq, K) with the seq axis == seq_len (36)
pos_outs <- names(feb_preds[[1]])[vapply(
  feb_preds[[1]], function(x) is.array(x) && length(dim(x)) == 3L && dim(x)[2] == seq_len, logical(1))]

if (length(pos_outs) > 0) {
  pi_out <- pos_outs[1]                                   # e.g. "per_index"
  K_feb  <- dim(feb_preds[[1]][[pi_out]])[3]
  pi_cn  <- names(classes)[c(which(!names(classes) %in% c("none", "padding")),
                             which(names(classes) == "none"))][seq_len(K_feb)]  # 6 real + none, trimmed to K
  bad <- is.na(pi_cn) | pi_cn == ""; pi_cn[bad] <- paste0("class", which(bad))
  message(sprintf("using '%s' (%d classes) as per_index -> columns: %s",
                  pi_out, K_feb, paste(pi_cn, collapse = ", ")))

  extract_pi <- function(arr) {
    d <- dim(arr)
    lapply(seq_len(d[1]), function(i) {
      m <- matrix(arr[i, , ], nrow = d[2], ncol = d[3]); colnames(m) <- pi_cn
      tibble::as_tibble(m) %>% dplyr::mutate(index = dplyr::row_number(), .before = 1)
    })
  }
  nn_input_comb_feb$per_index <- do.call(c, lapply(feb_preds, function(p) extract_pi(p[[pi_out]])))
} else {
  message("no per-position (seq=", seq_len, ") output found among {",
          paste(names(feb_preds[[1]]), collapse = ", "), "}; skipping per_index")
}

message(sprintf("scored %d windows with the Feb24 models", nrow(nn_input_comb_feb)))
nn_input_comb_feb %>%
  filter(known == 1) %>% arrange(desc(pred)) %>%
  select(peps, gene, model, win_type, pred, rank) %>%
  print(n = 25)

## nn_input_comb_feb is now a scored table (pred + rank, and per_index if present)
## you can hand to the plotting / analysis scripts in place of nn_input_comb.
