## =============================================================================
## Where do known peptides get lost between knowns.rds and the red rings?
##
## IMPORTANT: most of this attrition is INTENTIONAL. The training set is
## deliberately selective. It wants peptides that
##   (a) insert an N- or C-TERMINUS into the receptor pocket -- peptides that
##       insert their middle are excluded on purpose, being harder to predict;
##   (b) already have a cognate receptor-ligand model, because the per-residue
##       pep_pocket / pep_other channels are read off that complex.
## A peptide with no docked complex cannot be a positive no matter how well
## UniProt annotates it. So `knowns.rds` is not an incomplete peptide list, it
## IS the selection criterion.
##
## The point of this audit is therefore NOT "how many are missing" but "which
## losses are the intended selectivity and which are accidental". Every peptide
## is tagged:
##   by_design  the selection criteria doing their job
##   technical  an implementation detail, worth knowing about
##   bug        a genuinely distinct peptide silently collapsed
##
## Reads only cached/derived artefacts and writes one CSV. Nothing feeds back
## into the pipeline; this is a diagnostic.
##
## Usage:  source("inst/scripts/audit_known_peptides.R")
## =============================================================================

suppressPackageStartupMessages({library(dplyr); library(tibble)})

knowns_rds  <- path.expand("~/AF2_analysis/knowns.rds")
cache_rds   <- path.expand("~/AF2_analysis/nn_dat_cache.rds")
id_map_rds  <- "/Users/kbrulois/R_projects/ligandFinder/data/id_mapping.rds"
ligand_rds  <- "/Users/kbrulois/R_projects/ligandFinder/inst/extdata/ligand_list.rds"
out_csv     <- path.expand("~/AF2_analysis/known_peptide_audit.csv")

## Reference peptides are length-filtered rather than typed: ligand_list's
## `ligand_type` is NA for 345 of its 550 rows -- secretin among them -- so
## filtering on it silently drops real peptides. (`gene_containing_known_ligand`
## is TRUE for every row, so it is not a filter either.) The list is a broad
## GPCR-ligand table including whole proteins (BMP2, COL1A1, CTNNB1), and a
## length cut is what turns it into a cleaved-peptide reference.
ref_max_len <- 100

stopifnot(file.exists(knowns_rds), file.exists(cache_rds), file.exists(ligand_rds))

id_map <- readRDS(id_map_rds)
e2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])
a2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry"]])

k <- readRDS(knowns_rds) %>%
  mutate(gene = e2s[p2_name], pep_id = paste0(gene, "_", p2_range))

peps  <- k %>% distinct(gene, p2_range, pep_id)
audit <- peps %>% mutate(status = NA_character_, stage = NA_character_,
                         category = NA_character_, detail = NA_character_)

mark <- function(audit, ids, stage, category, detail) {
  hit <- audit$pep_id %in% ids & is.na(audit$status)      # FIRST stage only
  audit$status[hit]   <- paste0("lost_", stage)
  audit$stage[hit]    <- stage
  audit$category[hit] <- category
  audit$detail[hit]   <- detail
  audit
}

## ---- replay the con_dat funnel ---------------------------------------------
s1 <- k %>% filter(location == "relevant")
audit <- mark(audit, setdiff(k$pep_id, s1$pep_id), "location", "by_design",
              "docked pose is not at the relevant site")

s2 <- s1 %>% filter(rank == 1, .by = afpd_dir_name)
audit <- mark(audit, setdiff(s1$pep_id, s2$pep_id), "rank", "technical",
              "no rank-1 model for its docking run")

s2 <- s2 %>% mutate(lig1_end_ind = as.integer(stringr::str_extract(lig1_end, "\\d+")),
                    lig1_end_t   = stringr::str_remove(lig1_end, "\\d+"),
                    p2_start = stringr::str_extract(p2_range, "^\\d+"),
                    p2_end   = stringr::str_extract(p2_range, "\\d+$"))
s3 <- s2 %>% group_by(pep_id) %>% arrange(lig1_end_ind, desc(iptm), .by_group = TRUE) %>%
  slice(1) %>% ungroup()
audit <- mark(audit, setdiff(s2$pep_id, s3$pep_id), "pep_id_dedup", "technical",
              "duplicate pep_id")

## The two historical dedup passes, replayed so the audit still names which
## terminus caused a collapse even now that the pipeline keys on the full range.
d_end <- s3 %>% group_by(lig1_end_t, lig1_end_ind, p2_name, p2_end) %>%
  arrange(desc(iptm)) %>% slice(1) %>% ungroup()
audit <- mark(audit, setdiff(s3$pep_id, d_end$pep_id), "dedup_shared_end", "bug",
              "shares a C-terminal end with a higher-iptm peptide (fixed: dedup now keys on the full range)")

d_start <- d_end %>% group_by(lig1_end_t, lig1_end_ind, p2_name, p2_start) %>%
  arrange(desc(iptm)) %>% slice(1) %>% ungroup()
audit <- mark(audit, setdiff(d_end$pep_id, d_start$pep_id), "dedup_shared_start", "bug",
              "shares an N-terminal start with a higher-iptm peptide (fixed: dedup now keys on the full range)")

## Downstream stages are evaluated on the FULL-range survivors, i.e. what the
## corrected dedup keeps, so nothing is double-counted.
s4 <- s3 %>% group_by(lig1_end_t, lig1_end_ind, p2_name, p2_start, p2_end) %>%
  arrange(desc(iptm)) %>% slice(1) %>% ungroup()

## The terminal-insertion criterion, not a cleanup: lig1_end_ind is how far the
## pocket-inserting residue sits from the peptide terminus, so >= 16 means the
## peptide inserts its MIDDLE. Those are excluded deliberately as harder to
## predict. (Below that, 9.2 still separates < 3 as a clean terminal insertion
## from 3-15, which it labels loop_N / loop_C.)
s5 <- s4 %>% filter(lig1_end_ind < 16)
audit <- mark(audit, setdiff(s4$pep_id, s5$pep_id), "end_ind", "by_design",
              "middle-inserting peptide: pocket contact >= 16 residues from either terminus")

## ---- what actually became a window -----------------------------------------
kd  <- readRDS(cache_rds)$known_dat
win <- kd %>% distinct(pep_id, .keep_all = TRUE) %>%
  select(pep_id, window = peps, win_type, target)
audit <- audit %>% left_join(win, by = "pep_id")

audit <- mark(audit, audit$pep_id[is.na(audit$window) & is.na(audit$status)],
              "no_window", "technical",
              "no 36-residue window built (missing sequence, or truncated at a terminus)")

audit$in_scored_set <- !is.na(audit$win_type) & audit$win_type == "db"
audit <- mark(audit, audit$pep_id[!is.na(audit$win_type) & audit$win_type != "db" &
                                  is.na(audit$status)], "db_filter", "by_design",
              "window is chym/pep_end; training is restricted to dibasic-anchored windows")

audit$status[is.na(audit$status)]   <- "kept"
audit$category[audit$status == "kept"] <- "kept"

## ---- reference peptides with no cognate receptor model ----------------------
## Matched on gene AND overlapping precursor coordinates, so a gene that is
## docked for one peptide still shows its OTHER peptides as unrepresented.
ll <- readRDS(ligand_rds) %>%
  mutate(gene = coalesce(a2s[accession], e2s[paste0(uniprot_name, "_HUMAN")]),
         start = as.integer(start), end = as.integer(end),
         len = end - start + 1L) %>%
  filter(!is.na(gene), !is.na(start), !is.na(end), len <= ref_max_len)

dk <- peps %>% mutate(ds = as.integer(stringr::str_extract(p2_range, "^\\d+")),
                      de = as.integer(stringr::str_extract(p2_range, "\\d+$")))

overlaps <- function(g, s, e) {
  d <- dk[dk$gene %in% g, ]
  nrow(d) > 0 && any(pmin(d$de, e) >= pmax(d$ds, s), na.rm = TRUE)
}
ll$represented <- mapply(overlaps, ll$gene, ll$start, ll$end)

never <- ll %>% filter(!represented) %>%
  transmute(gene, p2_range = paste0(start, "x", end),
            pep_id = NA_character_,
            status = "no_receptor_model", stage = "no_receptor_model",
            category = "by_design",
            detail = paste0("in ligand_list (", coalesce(final_name, "unnamed"),
                            ", ", len, " aa) but no cognate receptor-ligand model"),
            window = NA_character_, win_type = NA_character_, target = NA_character_,
            in_scored_set = FALSE)

audit <- bind_rows(audit, never) %>%
  arrange(factor(category, c("bug", "technical", "by_design", "kept")), gene, p2_range)

if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fwrite(audit, out_csv)
} else {
  readr::write_csv(audit, out_csv)
}

message("wrote ", out_csv, "  (", nrow(audit), " rows)")
print(audit %>% count(category, stage, sort = TRUE), n = 30)
message(sprintf("\n%d of %d docked peptides reach the scored set (%d genes)",
                sum(audit$in_scored_set), nrow(peps),
                dplyr::n_distinct(audit$gene[audit$in_scored_set])))
message(sprintf("reference peptides (<= %d aa) with no cognate receptor model: %d across %d genes",
                ref_max_len, nrow(never), dplyr::n_distinct(never$gene)))
