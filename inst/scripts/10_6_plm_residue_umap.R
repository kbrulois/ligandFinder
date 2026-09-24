#!/usr/bin/env Rscript
## ---- per-residue UMAP of the PLM embedding channels ----------------------------
## Is the per-residue structure the window models read off the embeddings visible
## in the embedding space itself, before any training?
##
## Takes the PCA-reduced per-residue channels `python -m lf_plm reduce` wrote
## (plm_01..plm_k, one row per residue of the 5,182 secretome proteins), embeds a
## background sample plus EVERY labelled residue in one UMAP per backend, and
## colours the labelled ones by their per-residue class.
##
## The labels come from the known peptide windows of BOTH termini
## (`known_idx`: pep_pocket, pep_other, DB, NT/CT_cleavage_context, gap), mapped
## from window-local positions to (accession, residue index). 17 of ~1,355
## labelled residues are covered by two known windows with DIFFERENT labels --
## all in PDYN, where adjacent dynorphin peptides overlap, so a residue is
## `pep_other` to the upstream peptide and `NT_cleavage_context` to the
## downstream one. Both are right from their own window, so they get their own
## `ambiguous` class rather than an arbitrary winner.
##
## Nothing here is trained: the coordinates are the foundation model's view.
##
## READ THE AA PANEL FIRST. The dominant axis of a per-residue PLM embedding is
## the residue's own identity: a 20-cluster split of the ProtT5 UMAP is ~84%
## pure by amino acid (ESM C ~32%). So a class that is defined by its residues
## -- DB is K/R by construction -- looks "separated" for trivial reasons. The
## kNN table therefore reports both the raw purity and a null in which labels
## are shuffled WITHIN each amino acid, which holds residue identity fixed and
## asks whether the embedding carries anything about the peptide class beyond it.
##
##   /usr/local/bin/Rscript inst/scripts/10_6_plm_residue_umap.R
##   ... --backends prot_t5 --n-background 50000
##
## The hand-built per-residue channels (relASA, cons_rs_n, min/mean_afm, DSSP
## secondary structure, backbone angles, H-bond energies) are joined on from
## secretome_aa by (accession, index), so the same space can be read against the
## features the model already had. secretome_aa is SPARSE -- only residues with
## DSSP/conservation coverage -- so roughly half the background residues have no
## metrics; they stay in the embedding and are NA in the metric columns.
##
## Output under ~/AF2_analysis, stamped:
##   plm_residue_umap_<stamp>.{svg,png}          classes + amino acid, both backends
##   plm_residue_umap_metrics_<stamp>.{svg,png}  the same space by relASA / cons / AFM
##   plm_residue_umap_byterm_<stamp>.{svg,png}   classes split by which terminus inserts
##   plm_residue_umap_labelled_<stamp>.csv       coords + label + metrics, labelled residues
##   plm_residue_umap_all_<stamp>.csv.gz         the same for every embedded residue
##   plm_residue_knn_<stamp>.csv                 kNN label purity vs the within-AA null
##   plm_residue_metric_knn_<stamp>.csv          how much of each metric the embedding encodes
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
backends   <- strsplit(.opt("--backends", "prot_t5,esm_c"), ",")[[1]]
n_bg       <- as.integer(.opt("--n-background", "120000"))
seed       <- as.integer(.opt("--seed", "42"))
k_nn       <- as.integer(.opt("--knn", "25"))
n_perm     <- as.integer(.opt("--n-perm", "200"))   # within-AA label shuffles for the null
nn_cache   <- path.expand(.opt("--nn-input", "~/AF2_analysis/lf_dcnn_compare_nn_input.rds"))
seq_parq   <- path.expand(.opt("--sequences", "~/AF2_analysis/lf_plm/sequences.parquet"))
out_dir    <- path.expand(.opt("--out-dir", "~/AF2_analysis"))
cache_rds  <- path.expand(.opt("--cache", "~/AF2_analysis/lf_plm/residue_umap_cache.rds"))
aa_rds     <- path.expand(.opt("--residue-metrics", "~/peptide_alg/build_residue_db/processed/secretome_aa.rds"))
met_cache  <- path.expand(.opt("--metric-cache", "~/AF2_analysis/lf_plm/residue_metrics_subset.rds"))
refresh    <- "--refresh" %in% .args
plm_of     <- c(prot_t5 = "~/AF2_analysis/lf_plm/prot_t5_pca32.parquet",
                esm_c   = "~/AF2_analysis/lf_plm/esm_c_pca32.parquet")
nice_of    <- c(prot_t5 = "ProtT5-XL-U50", esm_c = "ESM C 600M")
stamp      <- format(Sys.time(), "%Y%m%d_%H%M%S")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "plm_features.R"))          # lf_read_parquet
source(file.path(ROOT, "R", "known_peptide_ends.R"))    # lf_known_peptide_ends
if (!requireNamespace("uwot", quietly = TRUE)) stop("install.packages('uwot')")

## ---- per-residue labels from the known windows ---------------------------------
.cc <- readRDS(nn_cache)
kn  <- bind_rows(lapply(names(.cc$nn_input), function(t)
  .cc$nn_input[[t]]$all %>% filter(known == 1) %>% mutate(term = t)))
message(sprintf("known windows: %d (%s)", nrow(kn),
                paste(sprintf("%s=%d", names(table(kn$term)), table(kn$term)), collapse = " ")))

lab <- bind_rows(lapply(seq_len(nrow(kn)), function(i) {
  md <- kn$meta_data[[i]]; cl <- as.character(kn$known_idx[[i]])
  tibble(peps = kn$peps[i], gene = kn$gene[i], term = kn$term[i],
         index = md$index, AA = as.character(md$AA), label = cl)
})) %>%
  filter(!is.na(index), label != "padding")

## a residue in two known windows with two different labels keeps both facts
lab <- lab %>%
  group_by(gene, index, AA) %>%
  summarise(label = if (n_distinct(label) == 1) first(label) else "ambiguous",
            peps  = paste(sort(unique(peps)), collapse = "; "),
            term  = paste(sort(unique(term)),  collapse = "/"), .groups = "drop")
message(sprintf("labelled residues: %d unique (%d ambiguous across two known windows)",
                nrow(lab), sum(lab$label == "ambiguous")))
print(table(lab$label))

## Which terminus of the peptide inserts into the receptor pocket, per known
## window (R/known_peptide_ends.R). Every known window IS an insertion event --
## 9.2 required a terminal insertion to call a peptide a training positive --
## so the window's terminus is the inserting one, and `insertion` says whether
## the pocket-contacting residue sits at the very end or a few residues in.
.klist <- file.path(ROOT, "inst", "extdata", "ligand_list.rds")
.idmap <- file.path(ROOT, "data", "id_mapping.rds")
ins_of <- NULL
if (file.exists(.klist) && file.exists(.idmap)) {
  .kn2 <- suppressMessages(lf_known_peptide_ends(kn, .klist, "~/AF2_analysis/knowns.rds", readRDS(.idmap)))
  ins_of <- stats::setNames(.kn2$insertion, .kn2$peps)
  message("known windows by terminus x insertion:")
  print(table(terminus = sub("^loop_", "", as.character(.kn2$target)), insertion = .kn2$insertion))
}

## gene -> accession (2 genes carry 2 accessions; the AA check below picks the right one)
seqs <- lf_read_parquet(seq_parq)

class_cols <- c(CT_cleavage_context = "#FED439FF", DB = "#370335FF", gap = "#8A9197FF",
                NT_cleavage_context = "#D2AF81FF", pep_other = "#D5E4A2FF",
                pep_pocket = "#197EC0FF", ambiguous = "#D62728")
lab$label <- factor(lab$label, levels = names(class_cols))

## ---- one UMAP per backend --------------------------------------------------------
## Cached: the embedding is minutes of work and the plot is seconds, so tweaking
## the figure should not recompute it. --refresh forces a recompute.
.fp <- list(backends = backends, n_bg = n_bg, seed = seed, k_nn = k_nn, v = 2)
if (!refresh && file.exists(cache_rds) && identical(readRDS(cache_rds)$fingerprint, .fp)) {
  .cache <- readRDS(cache_rds)
  emb <- .cache$emb; knn_tbl <- .cache$knn_tbl
  message("reusing ", cache_rds, " (--refresh to recompute)")
  print(as.data.frame(bind_rows(knn_tbl)), digits = 3, row.names = FALSE)
} else {
emb <- list(); knn_tbl <- list()
for (bk in backends) {
  f <- path.expand(plm_of[[bk]])
  if (!file.exists(f)) { message("skipping ", bk, ": no ", f); next }
  message("\n== ", bk, " ==")
  tbl  <- lf_read_parquet(f)
  cols <- grep("^plm_", names(tbl), value = TRUE)
  tbl$row <- seq_len(nrow(tbl))

  ## map each labelled residue to its row, checking the residue letter -- a
  ## mismatch means a different sequence version, not a different accession
  cand <- lab %>% inner_join(seqs %>% select(accession, gene), by = "gene", relationship = "many-to-many")
  keyed <- tbl %>% select(accession, index, AA, row) %>%
    inner_join(cand %>% select(accession, index, AA, label, peps, term),
               by = c("accession", "index", "AA"))
  lost <- nrow(lab) - n_distinct(keyed$accession, keyed$index)
  if (lost) message(sprintf("  %d labelled residues had no matching (accession, index, AA) row", lost))
  keyed <- keyed %>% distinct(accession, index, .keep_all = TRUE)

  set.seed(seed)
  bg_rows <- sample(setdiff(tbl$row, keyed$row), min(n_bg, nrow(tbl) - nrow(keyed)))
  rows <- c(keyed$row, bg_rows)
  X <- as.matrix(tbl[rows, cols])
  message(sprintf("  UMAP on %d residues (%d labelled + %d background) x %d dims ...",
                  nrow(X), nrow(keyed), length(bg_rows), length(cols)))

  set.seed(seed)
  u <- uwot::umap(X, n_neighbors = 15, min_dist = 0.1, metric = "euclidean",
                  n_threads = max(1, parallel::detectCores() - 1), verbose = FALSE)
  emb[[bk]] <- tibble(backend = bk, UMAP1 = u[, 1], UMAP2 = u[, 2],
                      is_bg = c(rep(FALSE, nrow(keyed)), rep(TRUE, length(bg_rows))),
                      label = c(as.character(keyed$label), rep(NA_character_, length(bg_rows))),
                      peps  = c(keyed$peps, rep(NA_character_, length(bg_rows))),
                      term  = c(keyed$term, rep(NA_character_, length(bg_rows))),
                      accession = tbl$accession[rows], index = tbl$index[rows], AA = tbl$AA[rows])

  ## The quantitative read, in the ORIGINAL 32-d space (a UMAP picture can
  ## flatter or hide structure): of each labelled residue's k nearest
  ## neighbours among all labelled residues, what share carry its own label,
  ## against the share expected if labels were scattered at random.
  L <- as.matrix(tbl[keyed$row, cols])
  lv <- as.character(keyed$label); aa <- keyed$AA
  nn_idx <- FNN::get.knn(L, k = min(k_nn, nrow(L) - 1))$nn.index
  purity_of <- function(v) vapply(seq_along(v), function(i) mean(v[nn_idx[i, ]] == v[i]), numeric(1))
  obs <- purity_of(lv)

  ## Null: shuffle labels WITHIN each amino acid, so every AA keeps its label
  ## composition and only the assignment to individual residues moves. Purity
  ## above this null is class information the embedding holds beyond residue
  ## identity -- which is what the 20 islands of the ProtT5 map are made of.
  set.seed(seed)
  perm <- replicate(n_perm, {
    v <- lv
    for (a in unique(aa)) { i <- which(aa == a); v[i] <- sample(lv[i]) }
    tapply(purity_of(v), lv, mean)
  })
  knn_tbl[[bk]] <- tibble(backend = bk, label = lv, purity = obs) %>%
    group_by(backend, label) %>%
    summarise(n = n(), knn_purity = mean(purity), .groups = "drop") %>%
    mutate(baseline    = n / sum(n),
           null_withinAA = as.numeric(rowMeans(perm)[label]),
           enrichment  = knn_purity / baseline,
           vs_null     = knn_purity / null_withinAA,
           p_perm      = as.numeric(vapply(seq_along(label), function(i)
             (sum(perm[label[i], ] >= knn_purity[i]) + 1) / (n_perm + 1), numeric(1))))
  print(as.data.frame(knn_tbl[[bk]]), digits = 3, row.names = FALSE)
}
dir.create(dirname(cache_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(fingerprint = .fp, emb = emb, knn_tbl = knn_tbl), cache_rds)
}
if (!length(emb)) stop("no backend produced an embedding")

d <- bind_rows(emb) %>%
  mutate(backend = factor(nice_of[backend], levels = unname(nice_of[backends])),
         label   = factor(label, levels = names(class_cols)),
         ## the terminus whose insertion this residue's window describes
         terminus = factor(ifelse(is.na(term), NA_character_,
                                  ifelse(grepl("/", term), "both", term)),
                           levels = c("N", "C", "both")),
         insertion = if (is.null(ins_of)) NA_character_ else
           unname(ins_of[sub(";.*$", "", peps)]))

## ---- the hand-built per-residue channels, joined by (accession, index) ----------
## Cached subset: secretome_aa is 3 GB and has 1,209 columns; only the residues
## in this embedding and the channels the model uses are kept.
met_cols <- c("cons_rs", "cons_rs_n", "min_afm", "mean_afm", "max_afm", "relASA", "SS",
              "Phi", "Psi", "Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin",
              "NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy")
keys <- d %>% distinct(accession, index)
if (!file.exists(met_cache) && file.exists(aa_rds)) {
  message("building ", met_cache, " from ", aa_rds, " (3 GB, one-off) ...")
  sa <- readRDS(aa_rds)
  saveRDS(sa %>% semi_join(keys, by = c("accession", "index")) %>%
            select(any_of(c("accession", "index", "AA", met_cols))), met_cache)
  rm(sa); invisible(gc())
}
if (file.exists(met_cache)) {
  ## Join on (accession, index, AA), not just (accession, index): secretome_aa
  ## and sequence_uni disagree on the residue letter for ~1% of positions (a
  ## sequence-version mismatch), and 9.2 takes the same precaution -- a
  ## mismatch must surface as a missing metric, never as a silently misaligned
  ## one attached to the wrong residue.
  met <- readRDS(met_cache) %>% select(any_of(c("accession", "index", "AA", met_cols))) %>%
    distinct(accession, index, .keep_all = TRUE)
  n_key <- d %>% distinct(accession, index) %>% inner_join(met, by = c("accession", "index")) %>% nrow()
  d <- d %>% left_join(met, by = c("accession", "index", "AA"))
  n_ok <- d %>% filter(!is.na(relASA)) %>% distinct(accession, index) %>% nrow()
  message(sprintf("per-residue metrics joined on (accession, index, AA): %.0f%% of embedded residues covered",
                  100 * mean(!is.na(d$relASA))))
  message(sprintf("  %d residues matched on (accession, index) but were dropped by the AA check (sequence-version mismatch)",
                  max(0L, n_key - n_ok)))
} else {
  message("no ", met_cache, " and no ", aa_rds, " -- skipping the per-residue metrics")
  for (m in met_cols) d[[m]] <- NA
}

## How much of each hand-built channel does the embedding already encode? For a
## held-out half of the residues, predict the channel as the mean over its k
## nearest neighbours among the other half IN THE 32-d EMBEDDING, and score it
## as R^2 against the channel's own variance. High R^2 = that channel is largely
## redundant with the embedding; low = it carries something the PLM does not.
num_met <- setdiff(met_cols, "SS")
met_knn <- bind_rows(lapply(names(emb), function(bk) {
  x <- emb[[bk]] %>% left_join(met, by = c("accession", "index", "AA"))
  f <- path.expand(plm_of[[bk]])
  tbl <- lf_read_parquet(f); cols <- grep("^plm_", names(tbl), value = TRUE)
  X <- as.matrix(tbl[match(paste(x$accession, x$index), paste(tbl$accession, tbl$index)), cols])
  set.seed(seed); half <- sample(c(TRUE, FALSE), nrow(X), replace = TRUE)
  idx <- FNN::get.knnx(X[half, , drop = FALSE], X[!half, , drop = FALSE], k = k_nn)$nn.index
  bind_rows(lapply(num_met, function(m) {
    tr <- x[[m]][half]; te <- x[[m]][!half]
    pred <- rowMeans(matrix(tr[idx], nrow = nrow(idx)), na.rm = TRUE)
    ok <- is.finite(pred) & is.finite(te)
    tibble(backend = bk, metric = m, n = sum(ok),
           r2 = 1 - mean((te[ok] - pred[ok])^2) / stats::var(te[ok]),
           spearman = suppressWarnings(cor(te[ok], pred[ok], method = "spearman")))
  }))
}))
message("\nhow much of each hand-built channel the PLM embedding encodes (kNN R^2, held-out half):")
print(met_knn %>% select(backend, metric, n, r2) %>%
        tidyr::pivot_wider(names_from = backend, values_from = c(n, r2)) %>%
        as.data.frame(), digits = 3, row.names = FALSE)

## ---- plot ------------------------------------------------------------------------
bg  <- d %>% filter(is.na(label))
fg  <- d %>% filter(!is.na(label))
## one panel per backend, composed with patchwork: the two UMAPs have
## unrelated coordinate systems (free scales), which facet_wrap cannot combine
## with coord_equal
panels <- lapply(levels(d$backend), function(nm) {
  b <- bg %>% filter(backend == nm); f <- fg %>% filter(backend == nm)
  g <- geom_point(data = b, colour = "grey86", size = 0.15, alpha = 0.5)
  ggplot(mapping = aes(UMAP1, UMAP2)) +
    (if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 200) else g) +
    geom_point(data = f, aes(colour = label), size = 0.9, alpha = 0.9) +
    scale_colour_manual(values = class_cols, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)) +
    coord_equal() +
    labs(title = nm, x = "UMAP1", y = "UMAP2") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
})
## second row: the SAME embeddings coloured by amino acid, because that is what
## the islands are -- without it the class panel invites over-reading
aa_lv   <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
aa_cols <- stats::setNames(scales::hue_pal()(length(aa_lv)), aa_lv)
aa_panels <- lapply(levels(d$backend), function(nm) {
  x <- d %>% filter(backend == nm, AA %in% aa_lv) %>% mutate(AA = factor(AA, levels = aa_lv))
  g <- geom_point(data = x, aes(colour = AA), size = 0.15, alpha = 0.5)
  ggplot(mapping = aes(UMAP1, UMAP2)) +
    (if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 200) else g) +
    scale_colour_manual(values = aa_cols, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) +
    coord_equal() +
    labs(title = paste0(nm, " \u2014 coloured by amino acid"), x = "UMAP1", y = "UMAP2") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
})
p <- patchwork::wrap_plots(
       patchwork::wrap_plots(panels,    nrow = 1, guides = "collect") & theme(legend.position = "bottom"),
       patchwork::wrap_plots(aa_panels, nrow = 1, guides = "collect") & theme(legend.position = "bottom"),
       ncol = 1) +
  patchwork::plot_annotation(
    title = "Per-residue PLM embedding space: the known peptides' residue classes, and what the space is actually organised by",
    subtitle = sprintf("32 PCA channels per residue, no training involved. Grey = %s background residues, coloured = the %d labelled residues of the %d known N/C windows.\nBottom row: the same embeddings coloured by residue identity -- the dominant axis, so read the class panel against the within-AA null in the kNN table.",
                       format(nrow(bg) / nlevels(d$backend), big.mark = ","),
                       nrow(fg) / nlevels(d$backend), nrow(kn)))

for (ext in c("svg", "png"))
  ggsave(file.path(out_dir, sprintf("plm_residue_umap_%s.%s", stamp, ext)), p,
         width = 6 * nlevels(d$backend) + 1, height = 15, dpi = 150)

## ---- the same space, coloured by the hand-built channels ------------------------
show_met <- c(relASA = "relative solvent accessibility", cons_rs_n = "conservation (normalised)",
              mean_afm = "AlphaMissense (mean)", min_afm = "AlphaMissense (min)")
show_met <- show_met[vapply(names(show_met), function(m) any(!is.na(d[[m]])), logical(1))]
if (length(show_met)) {
  mp <- lapply(names(show_met), function(m) {
    lapply(levels(d$backend), function(nm) {
      x <- d %>% filter(backend == nm, !is.na(.data[[m]]))
      g <- geom_point(data = x, aes(colour = .data[[m]]), size = 0.2, alpha = 0.6)
      ggplot(mapping = aes(UMAP1, UMAP2)) +
        (if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 200) else g) +
        scale_colour_viridis_c(name = NULL, option = "magma") +
        coord_equal() +
        labs(title = sprintf("%s \u2014 %s", nm, show_met[[m]]), x = NULL, y = NULL) +
        theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
    })
  })
  pm <- patchwork::wrap_plots(do.call(c, mp), ncol = nlevels(d$backend)) +
    patchwork::plot_annotation(
      title = "Per-residue PLM embedding space, coloured by the hand-built channels",
      subtitle = sprintf("%.0f%% of embedded residues have secretome_aa coverage; the rest are omitted from these panels",
                         100 * mean(!is.na(d$relASA))))
  for (ext in c("svg", "png"))
    ggsave(file.path(out_dir, sprintf("plm_residue_umap_metrics_%s.%s", stamp, ext)), pm,
           width = 5.5 * nlevels(d$backend), height = 5 * length(show_met), dpi = 150)
}

## ---- classes split by which terminus inserts -----------------------------------
## The two termini are not the same problem: the C-terminal windows are mostly
## clean END insertions (the pocket residue within 2 of the cut) while the
## N-terminal ones are mostly LOOP insertions, so a class can be well separated
## on one side and not the other. Panel counts say how many residues each holds.
term_lv <- c("N", "C")

## facet_grid, not patchwork: one shared legend without a composition layer to
## fight. The price is no coord_equal (free scales and a fixed aspect ratio are
## mutually exclusive in ggplot2), which costs nothing here -- UMAP axes are
## arbitrary units and the two backends are not on a common scale anyway.
fg_t <- d %>% filter(!is.na(label), terminus %in% term_lv) %>%
  mutate(terminus = factor(terminus, levels = term_lv))
## the background is the same cloud in every column, so repeat it per terminus
bg_t <- bind_rows(lapply(term_lv, function(tm)
  d %>% filter(is.na(label)) %>% mutate(terminus = factor(tm, levels = term_lv))))

## strip labels carry the counts: how many residues, and how those windows'
## peptides insert into the receptor
strip <- fg_t %>% filter(backend == levels(d$backend)[1]) %>%
  group_by(terminus) %>%
  summarise(txt = sprintf("%s-terminal insertion\n%d residues | %s", first(terminus), n(),
                          {ins <- table(insertion[!is.na(insertion)])
                           if (length(ins)) paste(sprintf("%s: %d", names(ins), ins), collapse = ", ")
                           else "insertion type unavailable"}), .groups = "drop")
lab_t <- stats::setNames(strip$txt, as.character(strip$terminus))

g_bg <- geom_point(data = bg_t, colour = "grey88", size = 0.15, alpha = 0.45)
pt <- ggplot(mapping = aes(UMAP1, UMAP2)) +
  (if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g_bg, dpi = 200) else g_bg) +
  geom_point(data = fg_t, aes(colour = label), size = 1.1, alpha = 0.9) +
  scale_colour_manual(values = class_cols, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)) +
  facet_grid(rows = vars(backend), cols = vars(terminus), scales = "free",
             labeller = labeller(terminus = lab_t)) +
  labs(title = "Per-residue PLM embedding space, stratified by which terminus inserts into the receptor",
       subtitle = "Same embeddings and background as the main figure; only the labelled residues are split, by the terminus of the known window they came from.",
       x = "UMAP1", y = "UMAP2") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9), strip.background = element_rect(fill = "grey95"))
for (ext in c("svg", "png"))
  ggsave(file.path(out_dir, sprintf("plm_residue_umap_byterm_%s.%s", stamp, ext)), pt,
         width = 6 * length(term_lv), height = 5.6 * nlevels(d$backend), dpi = 150)

## kNN label purity within each terminus: does the embedding separate the
## classes better where the peptide inserts by its end than by a loop?
term_knn <- bind_rows(lapply(names(emb), function(bk) {
  f <- path.expand(plm_of[[bk]]); tbl <- lf_read_parquet(f)
  cols <- grep("^plm_", names(tbl), value = TRUE)
  x <- emb[[bk]] %>% filter(!is.na(label)) %>%
    mutate(terminus = ifelse(grepl("/", term), "both", term))
  X <- as.matrix(tbl[match(paste(x$accession, x$index), paste(tbl$accession, tbl$index)), cols])
  bind_rows(lapply(term_lv, function(tm) {
    i <- which(x$terminus == tm); if (length(i) < k_nn + 1) return(NULL)
    nn <- FNN::get.knn(X[i, , drop = FALSE], k = k_nn)$nn.index
    lv <- as.character(x$label)[i]
    tibble(backend = bk, terminus = tm, label = lv,
           purity = vapply(seq_along(lv), function(j) mean(lv[nn[j, ]] == lv[j]), numeric(1)))
  }))
})) %>% group_by(backend, terminus, label) %>%
  summarise(n = n(), knn_purity = mean(purity), .groups = "drop") %>%
  mutate(baseline = n / ave(n, backend, terminus, FUN = sum), enrichment = knn_purity / baseline)
message("\nkNN label purity within each terminus (32-d space, among that terminus's residues only):")
print(term_knn %>% tidyr::pivot_wider(names_from = backend, values_from = c(knn_purity, enrichment)) %>%
        as.data.frame(), digits = 3, row.names = FALSE)
write.csv(term_knn, file.path(out_dir, sprintf("plm_residue_knn_byterm_%s.csv", stamp)), row.names = FALSE)

write.csv(fg, file.path(out_dir, sprintf("plm_residue_umap_labelled_%s.csv", stamp)), row.names = FALSE)
data.table::fwrite(d, file.path(out_dir, sprintf("plm_residue_umap_all_%s.csv.gz", stamp)))
write.csv(bind_rows(knn_tbl), file.path(out_dir, sprintf("plm_residue_knn_%s.csv", stamp)), row.names = FALSE)
write.csv(met_knn, file.path(out_dir, sprintf("plm_residue_metric_knn_%s.csv", stamp)), row.names = FALSE)
message("\nwrote ", file.path(out_dir, sprintf("plm_residue_umap_%s.{svg,png}", stamp)),
        "\n      plm_residue_umap_labelled_", stamp, ".csv  (coords per labelled residue)",
        "\n      plm_residue_knn_", stamp, ".csv  (kNN label purity in the 32-d space)")
