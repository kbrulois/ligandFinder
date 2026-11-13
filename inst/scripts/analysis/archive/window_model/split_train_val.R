# split_train_val.R
# Split nodes into training / validation sets using pair classifications.
#
# Rules:
#   - "close" pairs MUST stay together (both train or both val).
#   - "moderate" and "distant" pairs should be split ~evenly between the two
#     sets (i.e. we want roughly half of those pairs to cross the train/val
#     boundary).
#
# Approach:
#   1. Build a graph whose edges are only "close" pairs.
#   2. Connected components of that graph are indivisible groups.
#   3. Assign components to train/val aiming for a 50/50 split by node count,
#      while minimising the imbalance of moderate+distant pairs that end up
#      within-set vs cross-set (i.e. keeping cross-set share near 0.5).
#
# Inputs : ~/AF2_analysis/tre_pairs.csv  (from group_nodes.R)
# Outputs: ~/AF2_analysis/split_assignments.csv   (node, set)
#          ~/AF2_analysis/split_pair_summary.csv  (relation x location counts)

set.seed(42)

pairs_path  <- "~/AF2_analysis/tre_pairs.csv"
assign_path <- "~/AF2_analysis/split_assignments.csv"
summary_path <- "~/AF2_analysis/split_pair_summary.csv"
target_train_frac <- 0.50   # change if you want a different ratio

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

pairs <- read.csv(pairs_path, stringsAsFactors = FALSE)
stopifnot(all(c("node_a", "node_b", "relation") %in% names(pairs)))

nodes <- sort(unique(c(pairs$node_a, pairs$node_b)))
cat("Nodes:", length(nodes), "  Pairs:", nrow(pairs), "\n")
cat("Relation counts:\n"); print(table(pairs$relation))

# ---- 1. Indivisible groups from "close" pairs -----------------------------
close_edges <- pairs[pairs$relation == "close", c("node_a", "node_b")]
g <- graph_from_data_frame(close_edges, directed = FALSE, vertices = nodes)
comp <- components(g)
group_of <- setNames(comp$membership, names(comp$membership))
group_sizes <- as.integer(table(group_of))
cat("Indivisible groups:", length(group_sizes),
    "  largest:", max(group_sizes), "\n")

# ---- 2. Moderate + distant edges expressed at the group level --------------
md <- pairs[pairs$relation %in% c("moderate", "distant"), ]
md$ga <- group_of[md$node_a]
md$gb <- group_of[md$node_b]
md <- md[md$ga != md$gb, ]   # within-group pairs are fixed, ignore here

# ---- 3. Greedy balanced assignment of groups to train/val -----------------
# Process groups largest-first. For each group, pick the side (train/val) that
# (a) keeps node counts near the target ratio, and (b) keeps the share of
# cross-set moderate+distant pairs near 0.5. We combine both into a single
# score and pick the side with lower score.

n_total <- length(nodes)
target_train_n <- round(target_train_frac * n_total)
target_val_n   <- n_total - target_train_n
target_cross_frac <- 0.50   # we want ~50% of md pairs to cross train<->val

gid_order <- order(-group_sizes)   # largest groups first
assign <- setNames(rep(NA_character_, length(group_sizes)),
                   seq_along(group_sizes))

# running tallies
train_n <- 0L; val_n <- 0L
cross_md <- 0L; total_md <- nrow(md)

score_side <- function(side, gid, size) {
  # Hypothetical counts if we put gid on `side`
  t_n <- train_n + if (side == "train") size else 0L
  v_n <- val_n   + if (side == "val")   size else 0L

  # Cross-set md pairs contributed by this group's already-assigned neighbours
  nbrs <- md[md$ga == gid | md$gb == gid, ]
  other_g <- ifelse(nbrs$ga == gid, nbrs$gb, nbrs$ga)
  other_side <- assign[as.character(other_g)]
  known <- !is.na(other_side)
  if (any(known)) {
    adds_cross <- sum(other_side[known] != side)
  } else {
    adds_cross <- 0L
  }
  new_cross <- cross_md + adds_cross

  # Normalised penalties
  size_penalty <- abs(t_n - target_train_n) / n_total
  # denominator for cross uses total md pairs; this is an approximation since
  # we haven't seen every pair yet, but works well greedily.
  cross_penalty <- if (total_md > 0) {
    abs(new_cross / total_md - target_cross_frac)
  } else 0
  # weight: node-balance gets a bit more weight than cross-edge balance
  0.6 * size_penalty + 0.4 * cross_penalty
}

for (i in gid_order) {
  gid <- as.character(i)
  size <- group_sizes[i]
  s_train <- score_side("train", gid, size)
  s_val   <- score_side("val",   gid, size)
  # Break ties randomly
  if (isTRUE(all.equal(s_train, s_val))) {
    side <- sample(c("train", "val"), 1)
  } else {
    side <- if (s_train < s_val) "train" else "val"
  }
  assign[gid] <- side
  if (side == "train") train_n <- train_n + size else val_n <- val_n + size

  # Update cross_md with edges to already-assigned neighbours
  nbrs <- md[md$ga == gid | md$gb == gid, ]
  other_g <- ifelse(nbrs$ga == gid, nbrs$gb, nbrs$ga)
  other_side <- assign[as.character(other_g)]
  known <- !is.na(other_side) & as.character(other_g) != gid
  if (any(known)) {
    cross_md <- cross_md + sum(other_side[known] != side)
  }
}

# ---- 4. Write node -> set table -------------------------------------------
node_set <- data.frame(
  node = nodes,
  set  = assign[as.character(group_of[nodes])],
  group = as.integer(group_of[nodes]),
  stringsAsFactors = FALSE
)
write.csv(node_set, assign_path, row.names = FALSE)

# ---- 5. Pair-level diagnostics --------------------------------------------
pairs$set_a <- assign[as.character(group_of[pairs$node_a])]
pairs$set_b <- assign[as.character(group_of[pairs$node_b])]
pairs$location <- ifelse(pairs$set_a == pairs$set_b,
                         paste0("within_", pairs$set_a),
                         "cross_train_val")
tab <- as.data.frame(table(relation = pairs$relation,
                           location = pairs$location))
write.csv(tab, summary_path, row.names = FALSE)

cat("\n=== Results ===\n")
cat(sprintf("train: %d nodes  (%.1f%%)\n",
            sum(node_set$set == "train"),
            100 * mean(node_set$set == "train")))
cat(sprintf("val  : %d nodes  (%.1f%%)\n",
            sum(node_set$set == "val"),
            100 * mean(node_set$set == "val")))
cat("\nPair distribution:\n")
print(xtabs(Freq ~ relation + location, data = tab))
cat("\nWrote:\n  ", assign_path, "\n  ", summary_path, "\n")
