#!/usr/bin/env Rscript
## ---- per-RESIDUE binary labels, for the window-free residue models ------------
## The window model asks "does a peptide end at this anchor" and reads a 36-mer
## to answer. This exports the labels for a different question: taking each
## residue on its own, is it a POCKET residue, and is it C-terminal CLEAVAGE
## CONTEXT. No window, so the labels have to be precursor-global rather than
## window-relative, which `lf_pepend_labels()` is not -- it answers for one
## peptide inside one 36-residue frame.
##
## Both labels come from the SAME source the window labels do (`knowns.rds`, the
## docked models), through the same funnel 10_7 uses -- relevant site, rank-1
## model, per peptide the model whose inserting residue sits closest to a
## terminus, chemokines dropped -- so the two label sets cannot drift apart.
##
##   pep_pocket  residue is in the pocket of ANY docked peptide of this precursor
##               (9.2's clean_contacts: area > 1, dist < 6, then in_pocket)
##   pep_other   residue is INSIDE a docked peptide but not in its pocket, which
##               is exactly how `lf_pepend_labels()` splits the peptide's own
##               span. Pocket wins where two peptides overlap, as it does there.
##   ct_context  residue lies 1..CTX downstream of ANY peptide's C terminus.
##               CTX defaults to 8, which is exactly the span a C-anchored window
##               labels CT_cleavage_context: the anchor sits at window position
##               28 of 36, so positions 29..36 are anchor+1..anchor+8.
##   nt_context  residue lies 1..CTX_N UPSTREAM of ANY peptide's N terminus.
##               CTX_N is 7, NOT 8, and the asymmetry is real: an N-anchored
##               window puts its anchor at position 8 of 36, so the residues
##               before the peptide are positions 1..7 = anchor-7..anchor-1.
##               Using 8 here would label one residue the window never sees.
##
## Rows are emitted in the row order of lf_pepend_residues.npz -- row
## `offset[i] + r - 1` is residue r of precursor i -- so python can use the
## labels and the features by index with no join at all.
##
##   Rscript inst/scripts/10_9p_export_residue_labels.R
##   ... --ctx 8 --out ~/AF2_analysis/lf_resid_labels.npz
##
## Output: an npz of flat arrays, plus a csv of per-precursor counts.
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
docked_p <- path.expand(.opt("--docked", "~/AF2_analysis/knowns.rds"))
cache_p  <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
res_npz  <- path.expand(.opt("--residues-npz", "~/AF2_analysis/lf_pepend_residues.npz"))
out_p    <- path.expand(.opt("--out", "~/AF2_analysis/lf_resid_labels.npz"))
CTX      <- as.integer(.opt("--ctx", "8"))
CTX_N    <- as.integer(.opt("--ctx-n", "7"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "dcnn_bridge.R"))
np <- reticulate::import("numpy", convert = FALSE)

## ---- 1. the docked peptides, through 10_7's funnel ----------------------------
id_map <- readRDS(file.path(ROOT, "data", "id_mapping.rds"))
e2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])
is_chemokine <- function(gene) grepl("^(CCL|CXCL|XCL|CX3CL)", gene)

con <- readRDS(docked_p) %>%
  filter(location == "relevant") %>%
  filter(rank == 1, .by = afpd_dir_name) %>%
  mutate(gene      = unname(e2s[p2_name]),
         accession = p2_id,
         pep_id    = paste0(gene, "_", p2_range),
         pep_start = as.integer(stringr::str_extract(p2_range, "^\\d+")),
         pep_end   = as.integer(stringr::str_extract(p2_range, "\\d+$")),
         ins_ind   = as.integer(stringr::str_extract(lig1_end, "\\d+"))) %>%
  group_by(pep_id) %>% arrange(ins_ind, desc(iptm), .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
  filter(!is_chemokine(gene))

con$pocket <- lapply(con$contacts, function(ct) {
  if (is.null(ct) || !nrow(ct)) return(integer(0))
  ct <- ct %>% filter(area > 1 & dist < 6) %>%
    group_by(lig1_residue) %>% summarise(in_pocket = first(in_pocket), .groups = "drop")
  as.integer(stringr::str_extract(ct$lig1_residue[ct$in_pocket %in% TRUE], "\\d+$"))
})
con <- con %>% select(gene, accession, pep_id, pep_start, pep_end, pocket)
message(sprintf("docked peptides: %d over %d precursors",
                nrow(con), dplyr::n_distinct(con$accession)))

## ---- 2. project onto residues, in npz row order --------------------------------
fc <- readRDS(cache_p)
z  <- np$load(res_npz)
offset <- as.integer(reticulate::py_to_r(z[["offset"]]))
n_prot <- as.integer(reticulate::py_to_r(z[["n_prot"]]))
c_prot <- as.integer(reticulate::py_to_r(z[["c_prot"]]))
n_rows <- as.integer(reticulate::py_to_r(z[["feat"]]$shape[[0]]))
stopifnot(length(offset) == nrow(fc$prec))
## The npz carries no accession array, so its row blocks are identified only by
## position. Check the lengths line up with the cache before trusting that.
lens <- nchar(fc$prec$seq)
stopifnot(identical(as.integer(c(offset[-1], n_rows) - offset), as.integer(lens)))
message("npz row blocks match the cache's precursor lengths")

y_pocket <- integer(n_rows); y_ct <- integer(n_rows); y_nt <- integer(n_rows)
in_pep   <- integer(n_rows); prot_idx <- integer(n_rows); resno <- integer(n_rows)
for (i in seq_len(nrow(fc$prec))) {
  rows <- offset[i] + seq_len(lens[i])
  prot_idx[rows] <- i - 1L
  resno[rows]    <- seq_len(lens[i])
}
by_acc <- split(seq_len(nrow(con)), con$accession)
hit <- 0L
for (acc in names(by_acc)) {
  i <- match(acc, fc$prec$accession)
  if (is.na(i)) next
  hit <- hit + 1L
  base <- offset[i]; Lp <- lens[i]
  for (j in by_acc[[acc]]) {
    pk <- con$pocket[[j]]
    pk <- pk[pk >= 1L & pk <= Lp]
    if (length(pk)) y_pocket[base + pk] <- 1L
    ## 1..CTX downstream of this peptide's C terminus, clipped to the precursor
    ct <- (con$pep_end[j] + 1L):(con$pep_end[j] + CTX)
    ct <- ct[ct >= 1L & ct <= Lp]
    if (length(ct)) y_ct[base + ct] <- 1L
    ## and the mirror on the N side, CTX_N wide rather than CTX -- see the header
    nt <- (con$pep_start[j] - CTX_N):(con$pep_start[j] - 1L)
    nt <- nt[nt >= 1L & nt <= Lp]
    if (length(nt)) y_nt[base + nt] <- 1L
    ps <- con$pep_start[j]:con$pep_end[j]
    ps <- ps[ps >= 1L & ps <= Lp]
    if (length(ps)) in_pep[base + ps] <- 1L
  }
}
message(sprintf("precursors with docked peptides: %d of %d", hit, nrow(fc$prec)))

## The peptide's own span minus its pocket. Derived rather than accumulated in
## the loop because the pocket of ONE peptide must suppress pep_other for every
## overlapping peptide too, which a per-peptide write cannot see.
y_pep_other <- as.integer(in_pep == 1L & y_pocket == 0L)

## A residue outside the mature range is signal peptide or past the end; the
## window model calls it `padding` and never scores it, so mark it and let the
## python side drop it rather than learning from it.
mature <- integer(n_rows)
for (i in seq_len(nrow(fc$prec))) {
  rows <- offset[i] + seq_len(lens[i])
  mature[rows] <- as.integer(resno[rows] >= n_prot[i] & resno[rows] <= c_prot[i])
}
## which precursors have any label at all -- the rest are controls (all-negative)
labelled <- integer(nrow(fc$prec))
labelled[match(intersect(names(by_acc), fc$prec$accession), fc$prec$accession)] <- 1L

message(sprintf("\nresidues: %s total, %s mature", format(n_rows, big.mark = ","),
                format(sum(mature), big.mark = ",")))
for (nm in c("pep_pocket", "pep_other", "ct_context", "nt_context", "in_peptide")) {
  v <- switch(nm, pep_pocket = y_pocket, pep_other = y_pep_other, ct_context = y_ct,
              nt_context = y_nt, in_peptide = in_pep)
  message(sprintf("  %-11s positives %6s  (%.3f%% of mature, %.2f%% of labelled precursors' mature)",
                  nm, format(sum(v), big.mark = ","), 100 * sum(v) / sum(mature),
                  100 * sum(v) / sum(mature[labelled[prot_idx + 1L] == 1L])))
}

np$savez_compressed(out_p,
  y_pocket    = np$asarray(y_pocket,     dtype = "int8"),
  y_pep_other = np$asarray(y_pep_other,  dtype = "int8"),
  y_ct        = np$asarray(y_ct,         dtype = "int8"),
  y_nt        = np$asarray(y_nt,         dtype = "int8"),
  in_pep   = np$asarray(in_pep,   dtype = "int8"),
  mature   = np$asarray(mature,   dtype = "int8"),
  prot_idx = np$asarray(prot_idx, dtype = "int32"),
  resno    = np$asarray(resno,    dtype = "int32"),
  labelled = np$asarray(labelled, dtype = "int8"),
  ## channel names travel with the labels so xgboost's importance table
  ## names real features instead of c0..c25
  channel_names = np$asarray(fc$all_params3),
  accession = np$asarray(fc$prec$accession),
  gene      = np$asarray(fc$prec$gene),
  ctx       = np$asarray(CTX),
  ctx_n     = np$asarray(CTX_N))
message("\nlabels -> ", out_p)
