# make_training_sets_v2.R
# Drop-in replacement for the original make_training_sets() in
# 9.2_add_contact_data.R. Knowns are split into train/val using a
# neighbour-joining tree built from their peptide sequences:
#
#   * "close"    pairs (below `close_q` quantile of cophenetic distance) MUST
#                stay together - same set.
#   * "moderate" and "distant" pairs are split ~50/50 across train and val.
#
# Controls (ctr_dat) are still drawn randomly as in the original function.
# The returned structure and naming match the original so downstream code
# (e.g. nn_input$C$train, nn_input$chem$val) keeps working.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ape); library(igraph)
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Build an amino-acid matrix (rows = peptides, cols = aligned AA positions)
# from a k_dat subset. Uses meta_data$AA for each row; peptides must be the
# same length, otherwise we right-pad with "-".
.build_aa_mat <- function(k_sub, id_col = "pep_id") {
  aas <- purrr::map(k_sub$meta_data, ~ as.character(.x$AA))
  aas <- purrr::map(aas, ~ ifelse(is.na(.x), "-", .x))
  L   <- max(lengths(aas))
  aas <- purrr::map(aas, ~ c(.x, rep("-", L - length(.x))))
  mat <- do.call(rbind, aas)
  rownames(mat) <- k_sub[[id_col]]
  mat
}

# Pairwise cophenetic distances from a seq matrix via dist.aa -> nj -> ladderize
.pair_distances <- function(aa_mat) {
  if (nrow(aa_mat) < 4) {
    # NJ needs >= 4 tips; fall back to raw pairwise AA distance.
    D <- as.matrix(ape::dist.aa(aa_mat))
    return(D)
  }
  D   <- ape::dist.aa(aa_mat)
  tre <- ape::ladderize(ape::nj(D))
  as.matrix(ape::cophenetic.phylo(tre))
}

# Given a distance matrix and quantile cutoffs, return a data.frame of pairs
# labelled close / moderate / distant.
.label_pairs <- function(D, close_q = 0.10, mod_q = 0.66) {
  labs <- rownames(D)
  idx  <- which(upper.tri(D), arr.ind = TRUE)
  vals <- D[upper.tri(D)]
  close_cut <- as.numeric(quantile(vals, close_q))
  mod_cut   <- as.numeric(quantile(vals, mod_q))
  data.frame(
    a        = labs[idx[, 1]],
    b        = labs[idx[, 2]],
    distance = vals,
    relation = cut(vals,
                   breaks = c(-Inf, close_cut, mod_cut, Inf),
                   labels = c("close", "moderate", "distant")),
    stringsAsFactors = FALSE
  )
}

# Assign peptides to "train"/"val" obeying close-pair groups and balancing
# cross-set moderate+distant edges near 50%. Returns a named vector:
#   names = peptide id, values = "train" | "val".
.assign_train_val <- function(pep_ids, pairs, train_frac = 0.5,
                              w_size = 0.6, w_cross = 0.4,
                              seed = 42) {
  set.seed(seed)
  close_edges <- pairs[pairs$relation == "close", c("a", "b")]
  g    <- igraph::graph_from_data_frame(close_edges, directed = FALSE,
                                        vertices = pep_ids)
  comp <- igraph::components(g)
  gmap <- setNames(comp$membership, names(comp$membership))
  gsz  <- as.integer(table(gmap))

  md   <- pairs[pairs$relation %in% c("moderate", "distant"), ]
  md$ga <- gmap[md$a]; md$gb <- gmap[md$b]
  md   <- md[md$ga != md$gb, ]

  n_total      <- length(pep_ids)
  target_train <- round(train_frac * n_total)
  total_md     <- nrow(md)

  order_idx <- order(-gsz)
  assign    <- setNames(rep(NA_character_, length(gsz)), seq_along(gsz))
  train_n   <- 0L; val_n <- 0L; cross_md <- 0L

  score <- function(side, gid, size) {
    t_n <- train_n + if (side == "train") size else 0L
    nbrs <- md[md$ga == gid | md$gb == gid, , drop = FALSE]
    other <- ifelse(nbrs$ga == gid, nbrs$gb, nbrs$ga)
    oside <- assign[as.character(other)]
    known <- !is.na(oside)
    adds  <- if (any(known)) sum(oside[known] != side) else 0L
    new_cross <- cross_md + adds
    size_pen  <- abs(t_n - target_train) / n_total
    cross_pen <- if (total_md > 0) abs(new_cross / total_md - 0.5) else 0
    w_size * size_pen + w_cross * cross_pen
  }

  for (i in order_idx) {
    gid <- as.character(i); size <- gsz[i]
    s_t <- score("train", gid, size); s_v <- score("val", gid, size)
    side <- if (isTRUE(all.equal(s_t, s_v))) sample(c("train","val"), 1)
            else if (s_t < s_v) "train" else "val"
    assign[gid] <- side
    if (side == "train") train_n <- train_n + size else val_n <- val_n + size
    nbrs <- md[md$ga == gid | md$gb == gid, , drop = FALSE]
    other <- ifelse(nbrs$ga == gid, nbrs$gb, nbrs$ga)
    oside <- assign[as.character(other)]
    known <- !is.na(oside) & as.character(other) != gid
    if (any(known)) cross_md <- cross_md + sum(oside[known] != side)
  }

  setNames(assign[as.character(gmap[pep_ids])], pep_ids)
}

# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
# Args:
#   k_dat      : tibble of known peptides (as produced earlier in the script)
#   ctr_dat    : tibble of control windows
#   banned     : vector of gene symbols to exclude from controls
#   aa_mats    : optional named list of aa matrices, e.g.
#                list(N = mat_N, C = mat_C, chem = mat_chem). If NULL, matrices
#                are built from k_dat$meta_data$AA for each subset.
#   ctrl_n     : number of controls per set (default 400, 100 for chem train)
#   close_q    : quantile cutoff for "close" pairs
#   mod_q      : quantile cutoff for "moderate" pairs
#   train_frac : fraction of knowns to put in train (default 0.5)
#   seed       : RNG seed

make_training_sets <- function(k_dat, ctr_dat,
                               banned     = character(),
                               aa_mats    = NULL,
                               ctrl_n     = 400,
                               close_q    = 0.10,
                               mod_q      = 0.66,
                               train_frac = 0.5,
                               seed       = 42) {

  targets_map <- list(N = c("N", "loop_N"), C = c("C", "loop_C"))

  # Helper that produces train/val/all for one (non-chem) target group
  split_group <- function(tag) {
    targ <- targets_map[[tag]]
    k_sub <- k_dat %>% filter(target %in% targ)
    if (nrow(k_sub) < 2) {
      return(list(train = k_sub[0, ], val = k_sub[0, ], all = k_sub))
    }

    aa_mat <- if (!is.null(aa_mats[[tag]])) aa_mats[[tag]] else .build_aa_mat(k_sub)
    # Keep aa_mat aligned to k_sub order
    aa_mat <- aa_mat[k_sub$pep_id, , drop = FALSE]

    D      <- .pair_distances(aa_mat)
    pairs  <- .label_pairs(D, close_q = close_q, mod_q = mod_q)
    assign <- .assign_train_val(k_sub$pep_id, pairs,
                                train_frac = train_frac, seed = seed)

    k_train <- k_sub %>% filter(pep_id %in% names(assign)[assign == "train"])
    k_val   <- k_sub %>% filter(pep_id %in% names(assign)[assign == "val"])

    c_pool  <- ctr_dat %>%
      filter(win_type == "db", !gene %in% banned, target %in% targ)

    out <- list(
      train = bind_rows(k_train,
                        c_pool %>% slice_sample(n = min(ctrl_n, nrow(c_pool)))),
      val   = bind_rows(k_val,
                        c_pool %>% slice_sample(n = min(ctrl_n, nrow(c_pool)))),
      all   = bind_rows(k_sub,
                        ctr_dat %>% filter(target %in% targ))
    )
    # Stash diagnostics as attributes so the list shape matches the original.
    attr(out, "split") <- tibble(pep_id = names(assign), set = unname(assign))
    attr(out, "pairs") <- pairs
    out
  }

  tmp <- list(N = split_group("N"), C = split_group("C"))

  # ---- chem block (chemokines targeting N) --------------------------------
  chem_dat <- k_dat %>%
    filter(grepl("^(CCL|CXCL|XCL|CX3CL)", gene), target %in% "N")

  if (nrow(chem_dat) >= 2) {
    aa_mat <- if (!is.null(aa_mats[["chem"]])) aa_mats[["chem"]] else .build_aa_mat(chem_dat)
    aa_mat <- aa_mat[chem_dat$pep_id, , drop = FALSE]
    D      <- .pair_distances(aa_mat)
    pairs  <- .label_pairs(D, close_q = close_q, mod_q = mod_q)
    assign <- .assign_train_val(chem_dat$pep_id, pairs,
                                train_frac = train_frac, seed = seed + 1)
    k_train <- chem_dat %>% filter(pep_id %in% names(assign)[assign == "train"])
    k_val   <- chem_dat %>% filter(pep_id %in% names(assign)[assign == "val"])
  } else {
    k_train <- chem_dat[0, ]; k_val <- chem_dat
    pairs <- NULL; assign <- setNames(character(), character())
  }

  c_pool_chem <- ctr_dat %>%
    filter(!gene %in% banned, target %in% "N")

  chem_out <- list(
    train = bind_rows(k_train,
                      c_pool_chem %>% slice_sample(n = min(ctrl_n, nrow(c_pool_chem)))),
    val   = bind_rows(k_val,
                      c_pool_chem %>% slice_sample(n = min(ctrl_n, nrow(c_pool_chem))))
  )
  attr(chem_out, "split") <- tibble(pep_id = names(assign), set = unname(assign))
  attr(chem_out, "pairs") <- pairs
  tmp[["chem"]] <- chem_out

  tmp
}


# ---------------------------------------------------------------------------
# Example usage (mirrors the original call):
# ---------------------------------------------------------------------------
#
# nn_input <- make_training_sets(
#   k_dat      = known_dat,
#   ctr_dat    = c_dat,
#   banned     = banned,
#   train_frac = 0.5
# )
#
# # Inspect how the known peptides were split for the C group:
# attr(nn_input$C, "split")
# # And the relation counts across train/val:
# sp <- attr(nn_input$C, "split")
# attr(nn_input$C, "pairs") %>%
#   mutate(set_a = sp$set[match(a, sp$pep_id)],
#          set_b = sp$set[match(b, sp$pep_id)],
#          location = ifelse(set_a == set_b, paste0("within_", set_a),
#                            "cross_train_val")) %>%
#   count(relation, location)
