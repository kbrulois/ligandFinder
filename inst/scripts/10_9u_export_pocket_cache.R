#!/usr/bin/env Rscript
## ---- per-peptide pocket residues, as a plain rds the plot pages can read -------
## The protein page colours each known peptide's own span by whether the residue
## inserts into the receptor pocket -- `peptide (inserting)` vs
## `peptide (non-inserting)`, the same two classes the window models predict. That
## needs the docked contacts, which live in `knowns.rds` behind a nested
## `contacts` column: too heavy to open inside a page, and needing the same
## funnel 10_7/10_9p use or the colours would disagree with the labels the models
## were fitted against.
##
## So it is funnelled ONCE here and written as a small lookup. Chemokines are
## kept, unlike 10_9p: that script drops them because the window TRAINING set
## does, but a page for a chemokine should still show its pocket.
##
##   Rscript inst/scripts/10_9u_export_pocket_cache.R
##
## Output: ~/AF2_analysis/lf_pocket_by_peptide.rds
##   accession, pep_id, pep_start, pep_end, pocket (list of integer residues)
## ------------------------------------------------------------------------------
suppressMessages({ library(dplyr); library(purrr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
docked_p <- path.expand(.opt("--docked", "~/AF2_analysis/knowns.rds"))
out_p    <- path.expand(.opt("--out", "~/AF2_analysis/lf_pocket_by_peptide.rds"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
id_map <- readRDS(file.path(ROOT, "data", "id_mapping.rds"))
e2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])

con <- readRDS(docked_p) %>%
  filter(location == "relevant") %>%
  filter(rank == 1, .by = afpd_dir_name) %>%
  mutate(gene      = unname(e2s[p2_name]),
         accession = p2_id,
         pep_id    = paste0(gene, "_", p2_range),
         pep_start = as.integer(stringr::str_extract(p2_range, "^\\d+")),
         pep_end   = as.integer(stringr::str_extract(p2_range, "\\d+$")),
         ins_ind   = as.integer(stringr::str_extract(lig1_end, "\\d+"))) %>%
  ## one model per peptide: the one whose inserting residue sits closest to a
  ## terminus, ties broken by iptm -- 10_7's rule, so the pocket a page draws is
  ## the pocket the labels were built from
  group_by(pep_id) %>% arrange(ins_ind, desc(iptm), .by_group = TRUE) %>%
  slice(1) %>% ungroup()

con$pocket <- lapply(con$contacts, function(ct) {
  if (is.null(ct) || !nrow(ct)) return(integer(0))
  ct <- ct %>% filter(area > 1 & dist < 6) %>%
    group_by(lig1_residue) %>% summarise(in_pocket = first(in_pocket), .groups = "drop")
  as.integer(stringr::str_extract(ct$lig1_residue[ct$in_pocket %in% TRUE], "\\d+$"))
})

out <- con %>% select(accession, gene, pep_id, pep_start, pep_end, pocket)
saveRDS(out, out_p)
message(sprintf("%d peptide(s) over %d precursor(s); %d with a non-empty pocket\n  -> %s",
                nrow(out), n_distinct(out$accession),
                sum(lengths(out$pocket) > 0), out_p))
