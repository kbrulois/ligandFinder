## Per-gene card notes for the CV gallery index.
## A page where every end was held out is not interchangeable with one where
## some were trained on, so the card says which it is. The index CSS greens any
## note containing "all", so the never-seen wording is chosen to carry it.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
u  <- fread(file.path(S2, "cv5_unseen_ends.csv"))
kn <- rbindlist(lapply(c("N","C"), function(t) {
  d <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", t)))
  ## ANXA1_w-5-30 and MT-RNR2_w-6-29: a window can start BEFORE residue 1,
  ## so the range may carry a minus sign. The sign-blind pattern left those
  ## two training genes labelled "prediction only" on their cards.
  d[, .(gene = sub("_w-?[0-9]+--?[0-9]+$", "", peps), term = t, anchor, fold)] }))

## the 14 labelled precursors the window-FREE model never trained on
L <- reticulate::import("numpy")$load(file.path(A2, "lf_resid_labels.npz"), allow_pickle = TRUE)
lab <- as.integer(py_to_r <- reticulate::py_to_r(L$`__getitem__`("labelled")))
pg  <- as.integer(reticulate::py_to_r(L$`__getitem__`("prec_group")))
pi_ <- as.integer(reticulate::py_to_r(L$`__getitem__`("prot_idx")))
mat <- as.integer(reticulate::py_to_r(L$`__getitem__`("mature")))
gn  <- as.character(reticulate::py_to_r(L$`__getitem__`("gene")))
set.seed(NULL)
rng  <- reticulate::import("numpy")$random$default_rng(42L)
labp <- which(lab == 1L) - 1L
ctrl <- as.integer(reticulate::py_to_r(rng$choice(reticulate::np_array(which(lab == 0L) - 1L),
                                                  size = 100L, replace = FALSE)))
rows <- which(pi_ %in% c(labp, ctrl) & mat == 1L)
rng2 <- reticulate::import("numpy")$random$default_rng(42L)
key  <- pg[pi_[rows] + 1L]
units <- sort(unique(key))
units <- as.integer(reticulate::py_to_r(rng2$permutation(reticulate::np_array(units))))
sz <- table(key); target <- 0.4 * length(rows); val_u <- c(); acc <- 0
for (uu in units) { if (acc >= target) break; val_u <- c(val_u, uu); acc <- acc + sz[[as.character(uu)]] }
val_prot <- unique(pi_[rows[key %in% val_u]])
resid_held <- sort(unique(gn[val_prot + 1L][lab[val_prot + 1L] == 1L]))
cat("per-residue held-out labelled precursors:", length(resid_held), "\n")

## ligand_list carries the compound symbol "DEFB4A; DEFB4B"; the secretome
## and the plot bundle both call it DEFB4A, so the note would never match
## its page. Map to the first symbol.
u[, gene := sub(";.*$", "", gene)]
ns <- u[set == "never seen (peptide, < 60 aa)", .(ends = .N), by = gene]
tr <- kn[, .(ends = .N, folds = paste(sort(unique(fold)), collapse = ", ")), by = gene]
genes <- readLines(file.path(S2, "pages_all.txt"))
note <- vapply(genes, function(g) {
  if (g %in% ns$gene) {
    n <- ns$ends[ns$gene == g]
    s <- sprintf("all %d end(s) never seen — no fold trained on this gene", n)
  } else if (g %in% tr$gene) {
    r <- tr[gene == g]
    s <- sprintf("training gene — %d known end(s), held out in fold %s", r$ends, r$folds)
  } else s <- "no annotated ends — prediction only"
  if (g %in% resid_held) s <- paste0(s, "; also held out of the per-residue model")
  s
}, "")
saveRDS(tibble(gene = genes, note = unname(note)), file.path(S2, "cv5_notes.rds"))
cat(sprintf("notes: %d never-seen, %d training, %d prediction-only\n",
            sum(genes %in% ns$gene), sum(genes %in% tr$gene),
            sum(!genes %in% c(ns$gene, tr$gene))))
