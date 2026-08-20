## =============================================================================
## Select "active state" receptor conformations for initial-guess docking.
##
## Across the ~1.5M models (300 receptors x ~1000 peptides x 5 models), find the
## models where a peptide actually got INTO the orthosteric pocket. Those models
## carry a receptor conformation with an open pocket -- N-terminus and ECL2 out
## of the way -- which is exactly what the apo AF-DB monomer lacks.
##
## Selection is on POCKET OCCUPANCY, not ipTM: a peptide stuck on the
## extracellular surface or wrapped around the N-terminus can score well on ipTM
## while leaving the pocket occluded, which is the failure mode we're escaping.
## `num_res_in_pocket` is the geometric test (EC_dist < 1 & IC_dist < 1 & h_dist
## < 1, see R/afpd_metric_processing.R) and `lig1_location == "E"` keeps the
## extracellular-side binders.
##
## Ranking within a receptor: known pairs first (a real cognate ligand opening
## the pocket is the most trustworthy evidence the conformation is real), then
## best ipTM.
##
## Writes:
##   <out_stem>_manifest.csv  one row per selected model, full provenance;
##                            this is the input to ap_extract_open_models.R
##   <out_stem>_coverage.csv  per-receptor counts, incl. receptors with NO open model
## =============================================================================

## tidyverse the meta-package is awkward to install on Sherlock, and nothing
## here needs it: the verbs come from dplyr, and readr/data.table calls are
## namespaced. Added to the front of the search path rather than replacing it,
## and only when present, so this still runs off a laptop checkout.
.lf_lib <- "/home/groups/ebutcher/programs/pipeline/R_libs4.1"
if (dir.exists(.lf_lib)) .libPaths(c(.lf_lib, .libPaths()))

library(dplyr)
library(tibble)
library(readr)

## ---- options ---------------------------------------------------------------
metrics_csv     <- path.expand("~/AF2_analysis/Aug7_all.csv")
out_stem        <- path.expand("~/AF2_analysis/open_receptors")
n_per_receptor  <- 5        # conformers to keep per receptor
min_depth   <- 0
max_depth <- -0.45
max_radius <- 0.8
max_per_peptide <- NULL     # NULL = no limit. Set to e.g. 2 to force chemotype
                            # diversity, so one peptide can't define the pocket
                            # shape for a receptor (see the bias note below).



metrics <- data.table::fread(metrics_csv) %>% as_tibble()
message(sprintf("read %s models x %d columns",
                format(nrow(metrics), big.mark = ","), ncol(metrics)))

## ---- coverage: does every receptor actually have an open conformer? ---------
## Worth knowing BEFORE selecting -- receptors with zero open models need the
## trim/swing fallback on the apo AF-DB model instead.
is_open <- metrics$depth < min_depth &
           metrics$depth > max_depth &
           metrics$radius < max_radius

coverage <- metrics %>%
  mutate(open = is_open) %>%
  group_by(p1_name) %>%
  summarise(n_models       = n(),
            n_open         = sum(open, na.rm = TRUE),
            n_open_known   = sum(open & known_pair == "known", na.rm = TRUE),
            best_occupancy = max(num_res_in_pocket, na.rm = TRUE),
            best_pLDDT_lig1_mean_all_open = suppressWarnings(max(pLDDT_lig1_mean_all[open], na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(best_pLDDT_lig1_mean_all_open = ifelse(is.finite(best_pLDDT_lig1_mean_all_open), best_pLDDT_lig1_mean_all_open, NA_real_)) %>%
  arrange(n_open)

message(sprintf("receptors: %d total | %d with >=1 open model | %d with >=%d",
                nrow(coverage), sum(coverage$n_open > 0),
                sum(coverage$n_open >= n_per_receptor), n_per_receptor))
if (any(coverage$n_open == 0))
  message("  NO open model for: ",
          paste(coverage$p1_name[coverage$n_open == 0], collapse = ", "))

## ---- select the conformers -------------------------------------------------
sel <- metrics %>%
  filter(is_open) %>%
  mutate(known_first = as.integer(known_pair == "known")) %>%
  arrange(p1_name, desc(known_first), desc(pLDDT_lig1_mean_all), desc(num_res_in_pocket))

## optional: cap how many conformers any single peptide contributes
if (!is.null(max_per_peptide)) {
  sel <- sel %>%
    group_by(p1_name, p2_name) %>%
    slice_head(n = max_per_peptide) %>%
    ungroup() %>%
    arrange(p1_name, desc(known_first), desc(pLDDT_lig1_mean_all), desc(num_res_in_pocket))
}

sel <- sel %>%
  group_by(p1_name) %>%
  mutate(conformer = sprintf("c%02d", row_number())) %>%
  slice_head(n = n_per_receptor) %>%
  ungroup()

manifest <- sel %>%
  transmute(receptor = p1_name, receptor_id = p1_id, conformer,
            afpd_dir_name, code, model_c, model_e, rank,
            pdb_file = pdb_files, pae_file = pae_files,
            run_name, afpd_dir,
            peptide = p2_name, peptide_id = p2_id, peptide_range = p2_range,
            known_pair, iptm, radius, depth, pLDDT_lig1_mean_all,
            algorithm, complex_type)

## ---- write -----------------------------------------------------------------
## ap_extract_open_models.R reads the manifest directly: it carries run_name and
## afpd_dir_name, which is all that's needed to compute each archive's path
## (<models_root>/<run_name>/<afpd_dir_name>.tar).
readr::write_csv(manifest, paste0(out_stem, "_manifest.csv"))
readr::write_csv(coverage,  paste0(out_stem, "_coverage.csv"))

message(sprintf("selected %d models across %d receptors (%d from known pairs)",
                nrow(manifest), dplyr::n_distinct(manifest$receptor),
                sum(manifest$known_pair == "known")))
message("  ", paste0(out_stem, "_manifest.csv"), "   <- input to ap_extract_open_models.R")
message("  ", paste0(out_stem, "_coverage.csv"))

## How many receptors got the full quota, and how many are leaning on a single
## peptide? The second number is the bias check: if a receptor's conformers all
## came from one peptide, its pocket shape is defined by that chemotype.
manifest %>%
  group_by(receptor) %>%
  summarise(n = n(), n_peptides = dplyr::n_distinct(peptide), .groups = "drop") %>%
  count(n, n_peptides) %>%
  print(n = 40)
