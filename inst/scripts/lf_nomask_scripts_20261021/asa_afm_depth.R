## relASA vs AlphaMissense for (1) 100,000 random residues, (2) every peptide
## residue, (3) the pocket residues -- now with each peptide residue's INSERTION
## DEPTH into the receptor from its docked model.
##
## Residue sets: the window-free models' labels (lf_resid_labels.npz), exactly as
## the 2026-10-19 asa_vs_afm.png drew them: in_pep / y_pocket residues of the 66
## labelled precursors, minus the residues excluded by design (cons_rs_n,
## min_afm, mean_afm and relASA all 0: intracellular / TM). -> 2,571 / 501.
##
## Depth: 9.2_add_contact_data.R's definition, as the current pipeline applies
## it. 9.2 itself is NOT run any more (10_7 / 10_9u read knowns.rds directly and
## repeat its contact filter); this does the same:
##   * the model per peptide = 10_7's rule (10_9u's funnel): location ==
##     "relevant", rank 1 per run, then the inserting residue closest to a
##     terminus, ties by iptm -- the model the pocket labels were built from
##   * contacts with area > 1 and dist < 6, one row per peptide residue,
##     IC_dist = first(IC_dist)  (9.2's clean_contacts)
##   * IC_dist (R/afpd_metric_processing.R) is the residue's Ca projected on the
##     receptor's EC->IC axis, as distance to the INTRACELLULAR reference over
##     receptor length: 1 = extracellular surface, 0 = intracellular end, > 1 =
##     above the receptor. Plotted as depth = 1 - IC_dist (deeper = larger,
##     <= 0 = outside). A residue with no qualifying contact has no depth (NA).
suppressMessages({library(dplyr); library(data.table); library(reticulate)})
A2 <- path.expand("~/AF2_analysis")
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
out <- file.path(A2, "lf_headcmp", "asa_afm_depth.csv")
np <- import("numpy", convert = TRUE)

## ---- every residue: features + labels ---------------------------------------
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
L  <- np$load(file.path(A2, "lf_resid_labels.npz"), allow_pickle = TRUE)
lab_acc <- as.character(L[["accession"]])
lab <- data.table(accession = lab_acc[L[["prot_idx"]] + 1L], resno = as.integer(L[["resno"]]),
                  in_pep = as.logical(L[["in_pep"]]), pocket = as.logical(L[["y_pocket"]]))
res <- rbindlist(lapply(names(fc$feats), function(a) {
  m <- fc$feats[[a]]
  data.table(accession = a, resno = seq_len(nrow(m)), mean_afm = m[, "mean_afm"], relASA = m[, "relASA"],
             excluded = m[, "cons_rs_n"] == 0 & m[, "min_afm"] == 0 & m[, "mean_afm"] == 0 & m[, "relASA"] == 0)
}))
res <- merge(res, lab, by = c("accession", "resno"), all.x = TRUE)
res[is.na(in_pep), `:=`(in_pep = FALSE, pocket = FALSE)]
ok <- res[excluded == FALSE]
message(sprintf("residues: %s total, %s not excluded", format(nrow(res), big.mark = ","),
                format(nrow(ok), big.mark = ",")))

## ---- insertion depth, through 10_7's model choice ---------------------------
id_map <- readRDS(lf_src("data/id_mapping.rds"))
e2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])
con <- readRDS(file.path(A2, "knowns.rds")) %>%
  filter(location == "relevant") %>%
  filter(rank == 1, .by = afpd_dir_name) %>%
  mutate(gene = unname(e2s[p2_name]), accession = p2_id, pep_id = paste0(gene, "_", p2_range),
         pep_start = as.integer(stringr::str_extract(p2_range, "^\\d+")),
         pep_end   = as.integer(stringr::str_extract(p2_range, "\\d+$")),
         ins_ind   = as.integer(stringr::str_extract(lig1_end, "\\d+"))) %>%
  group_by(pep_id) %>% arrange(ins_ind, desc(iptm), .by_group = TRUE) %>% slice(1) %>% ungroup()
dep <- rbindlist(lapply(seq_len(nrow(con)), function(i) {
  ct <- con$contacts[[i]]
  if (is.null(ct) || !nrow(ct)) return(NULL)
  ct <- ct %>% filter(area > 1 & dist < 6) %>% group_by(lig1_residue) %>%
    summarise(IC_dist = first(IC_dist), in_pocket = first(in_pocket), .groups = "drop")
  data.table(accession = con$accession[i], pep_id = con$pep_id[i],
             resno = as.integer(stringr::str_extract(ct$lig1_residue, "\\d+$")),
             IC_dist = ct$IC_dist, in_pocket = ct$in_pocket %in% TRUE)
}))
## a residue covered by two docked peptides (overlapping products) keeps its
## DEEPEST placement
n_multi <- dep[, .N, by = .(accession, resno)][N > 1, .N]
dep <- dep[order(IC_dist)][, .SD[1], by = .(accession, resno)]
message(sprintf("depth: %d docked peptides, %s contacting residues (%d covered by >1 peptide; deepest kept)",
                nrow(con), format(nrow(dep), big.mark = ","), n_multi))

## ---- the three sets -----------------------------------------------------------
set.seed(42)
rnd <- ok[sample.int(nrow(ok), 100000L)][, set := "100,000 random residues"]
pep <- ok[in_pep == TRUE][, set := "all peptide residues"]
pkt <- ok[pocket == TRUE][, set := "pocket residues only"]
d <- rbindlist(list(rnd, pep, pkt))
d <- merge(d, dep[, .(accession, resno, IC_dist, in_pocket_docked = in_pocket)],
           by = c("accession", "resno"), all.x = TRUE, sort = FALSE)
d[set == "100,000 random residues", IC_dist := NA_real_]      # not peptides: no docked model
d[, depth := 1 - IC_dist]

## the labels and the docked contacts must agree on the pocket, or depth is
## being read off a different model than the one the labels came from
chk <- d[set == "pocket residues only"]
message(sprintf("pocket residues: %d; with a depth %d; docked in_pocket TRUE for %d (%.1f%%)",
                nrow(chk), sum(!is.na(chk$depth)), sum(chk$in_pocket_docked %in% TRUE),
                100 * mean(chk$in_pocket_docked %in% TRUE)))
for (s in unique(d$set)) {
  x <- d[set == s]
  message(sprintf("  %-26s n %6s | with depth %5d | depth median %s range %s",
                  s, format(nrow(x), big.mark = ","), sum(!is.na(x$depth)),
                  if (all(is.na(x$depth))) "-" else sprintf("%.2f", median(x$depth, na.rm = TRUE)),
                  if (all(is.na(x$depth))) "-" else paste(sprintf("%.2f", range(x$depth, na.rm = TRUE)), collapse = " to ")))
}
fwrite(d[, .(set, accession, resno, mean_afm, relASA, IC_dist, depth)], out)
message("-> ", out)
