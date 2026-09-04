## =============================================================================
## UMAP over the LEARNED 16-d embedding (the "embed" dense layer).
##
## Sibling of 10_3_peptide_umap.R, which UMAPs the raw flattened per-position
## features (36 params x 36 positions). This one instead uses the 16-d
## representation the global-ranking head already learned to sit on -- the same
## tap 10_1dcnn_new6.R uses for nearest-known retrieval -- so proximity here
## means "the model views these windows alike", not "their inputs look alike".
##
## Run AFTER 10_1dcnn_new6.R, in the same session: needs `models`, `nn_in_all`
## and `nn_input_comb`. It reuses `emb_all` if that object is already present.
##
## NOTE: the embedding CANNOT be recovered from a saved .keras file. Every R
## layer_lambda serializes under the one name "anonymous_R_function" (there are
## three: the positional ramp, the op_take, and the mask multiply that closes
## over mask_matrix_cat), so load_model() cannot tell them apart -- passing one
## replacement via custom_objects silently mis-wires the graph rather than
## erroring. Hence: same-session only.
## =============================================================================

library(tidyverse)

## ---- options ---------------------------------------------------------------
out_csv   <- path.expand("~/AF2_analysis/peptide_embed16_umap.csv")
out_svg   <- path.expand("~/AF2_analysis/peptide_embed16_umap.svg")
seed      <- 42
umap_args <- list(n_neighbors = 15, min_dist = 0.1, metric = "cosine")
make_plot <- TRUE
per_terminus <- TRUE   # ALSO write one UMAP per terminus (see caveat below)

stopifnot(exists("models"), exists("nn_in_all"), exists("nn_input_comb"))
if (!requireNamespace("uwot", quietly = TRUE))
  stop("install.packages('uwot') for UMAP")

l2norm <- function(m) m / sqrt(pmax(rowSums(m^2), 1e-12))

## ---- 16-d embedding, rows aligned with nn_input_comb -----------------------
## Identical to the retrieval block in 10_1dcnn_new6.R; reuse it if it survived.
if (exists("emb_all") && nrow(emb_all) == nrow(nn_input_comb)) {
  message("10_3d: reusing existing emb_all")
  emb <- emb_all
} else {
  emb <- do.call(rbind, lapply(names(nn_in_all), function(term) {
    embed_model <- keras3::keras_model(models[[term]]$input,
                                       keras3::get_layer(models[[term]], "embed")$output)
    l2norm(predict(embed_model, nn_in_all[[term]][["all"]][["x"]]))
  }))
}
stopifnot(nrow(emb) == nrow(nn_input_comb))
emb <- l2norm(emb)                       # idempotent; guards a non-normalised emb_all
message(sprintf("10_3d: embedding %d windows x %d dims", nrow(emb), ncol(emb)))

## ---- UMAP ------------------------------------------------------------------
set.seed(seed)
um <- do.call(uwot::umap, c(list(X = emb), umap_args))

out <- nn_input_comb %>%
  transmute(peps, gene,
            terminus = stringr::str_remove(as.character(target), "^loop_"),
            target   = as.character(target),
            win_type, known,
            uni_pep  = if ("uni_pep" %in% names(.)) uni_pep else NA,
            location = if ("location" %in% names(.)) location else NA,
            pred, rank,
            nn_closest_peptide = if ("nn_closest_peptide" %in% names(.)) nn_closest_peptide else NA,
            nn_closest_sim     = if ("nn_closest_sim" %in% names(.)) nn_closest_sim else NA) %>%
  mutate(UMAP1 = um[, 1], UMAP2 = um[, 2],
         uniprot_known = if_else(known == 1, gene, NA_character_), .after = 1) %>%
  bind_cols(as_tibble(emb, .name_repair = ~ paste0("emb", sprintf("%02d", seq_len(ncol(emb))))))

## ---- per-terminus UMAPs ----------------------------------------------------
## CAVEAT: N and C are SEPARATE networks, so their 16-d spaces are not aligned --
## a joint UMAP can split on terminus for that reason alone rather than for
## anything biological. These per-terminus embeddings are each internally valid.
if (per_terminus) {
  for (tm in unique(out$terminus)) {
    i <- which(out$terminus == tm)
    set.seed(seed)
    u2 <- do.call(uwot::umap, c(list(X = emb[i, , drop = FALSE]), umap_args))
    out[[paste0("UMAP1_", tm)]] <- NA_real_; out[[paste0("UMAP1_", tm)]][i] <- u2[, 1]
    out[[paste0("UMAP2_", tm)]] <- NA_real_; out[[paste0("UMAP2_", tm)]][i] <- u2[, 2]
  }
}

readr::write_csv(out, out_csv)
message("wrote ", out_csv, "  (", nrow(out), " windows x ", ncol(out), " cols)")

## ---- quick coloured UMAP (knowns highlighted + labelled) -------------------
if (make_plot) {
  p <- ggplot(out, aes(UMAP1, UMAP2)) +
    geom_point(aes(color = win_type, shape = terminus), size = 1, alpha = 0.45) +
    geom_point(data = filter(out, known == 1), color = "black", size = 2.2) +
    theme_bw() +
    labs(title = "Peptide-window UMAP -- learned 16-d embedding",
         subtitle = "black = known uniprot peptides")
  if (requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(data = filter(out, known == 1),
                                      aes(label = gene), size = 2, max.overlaps = 20)
  ggsave(out_svg, p, width = 12, height = 10)
  message("wrote ", out_svg)
}
