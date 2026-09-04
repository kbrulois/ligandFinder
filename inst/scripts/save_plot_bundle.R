## =============================================================================
## Run ONCE, in the session that already has a trained nn_input_comb.
## Writes everything the protein plots need -- and nothing else -- to one file.
##
## Deliberately DROPS the `data` (36 x n_channels feature matrices) and
## `known_idx` list-columns: the plots never read them, and they are most of the
## weight. Nothing here needs keras, secretome_aa, c_dat, known_dat or models.
## =============================================================================
bundle_path <- "~/AF2_analysis/plot_bundle.rds"

stopifnot(exists("the_input"), exists("pep_input"), exists("nn_input_comb"))

keep <- c("peps", "gene", "target", "win_type", "pred", "pred_raw", "rank_cat",
          "end_type", "nn_closest_peptide", "nn_closest_sim",
          "per_index", "meta_data", "known")

saveRDS(
  list(the_input   = the_input,
       pep_input   = pep_input,
       pred_to_plot = nn_input_comb %>% dplyr::select(dplyr::any_of(keep)),
       mets        = if (exists("mets")) mets else NULL,
       all_mets    = if (exists("all_mets")) all_mets else NULL,
       species_dat = if (exists("species_dat")) species_dat else NULL),
  bundle_path)

message("wrote ", bundle_path, "  (",
        round(file.size(path.expand(bundle_path)) / 2^20), " MB)")
