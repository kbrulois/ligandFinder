

parse_pdb <- function(pdb) {

  if(!bio3d::is.pdb(pdb)) pdb <- bio3d::read.pdb(pdb)

  pdb[["atom"]]    %>%
    mutate(AA = bio3d::aa321(resid)) %>%
    select(all_of(c("AA", "x", "y", "z", "resno", "b", "chain", "elety", "elesy"))) %>%
    nest_by(chain, resno, .key = "atom_level") %>%
    ungroup %>%
    mutate(AA = map_chr(atom_level, ~unique(.[["AA"]])), .before = everything()) %>%
    mutate(pLDDT = map_dbl(atom_level, ~mean(.[["b"]])), .after = "AA") %>%
    mutate(map_df(atom_level, \(x) {bind_cols(x %>%
                                                filter(elety == "CA") %>%
                                                select(all_of(c("x", "y", "z"))) %>%
                                                rename_with(.fn = ~paste0("CA_", .x)),
                                              x %>%
                                                filter(!elety %in% c("CA", "N", "O", "C")) %>%
                                                summarise(across(all_of(c("x", "y", "z")), mean))
    )})) %>%
    mutate(
      x = coalesce(na_if(x, NaN), CA_x),
      y = coalesce(na_if(y, NaN), CA_y),
      z = coalesce(na_if(z, NaN), CA_z)
    )

}

import_raw_metrics <- function(input_data) {

  message("importing raw metrics")

  input_data %>%

    mutate(model = map(og_file_name, ~list.files(paste0(dir_path, "/", .)) %>%
                         stringr::str_extract(., "model_\\d_.*_pred_\\d") %>%
                         unique(.) %>%
                         .[!is.na(.)])) %>%
    mutate(iptm = map(og_file_name, \(x) {
      tryCatch({
        tmp <- jsonlite::read_json(paste(dir_path, x, "ranking_debug.json", sep = "/"))
        tmp %>% as_tibble %>% unnest(everything()) %>% select(-order)
      }, error = function(e) {rep(NA, 5)})
    })) %>%

    unnest(cols = c("model", "iptm")) %>%

    mutate(pdb = map2(.x = og_file_name,
                      .y = model,
                      ~bio3d::read.pdb(paste(dir_path, .x, paste0("unrelaxed_", .y, ".pdb"), sep = "/")))) %>%

    mutate(pdb.xyz = map(pdb, parse_pdb)) %>%

    mutate(pae = map2(.x = og_file_name, .y = model, \(.x, .y) {
      tmp <- jsonlite::read_json(paste(dir_path, .x, paste0("pae_", .y, ".json"), sep = "/"))
      tmp[[1]][["predicted_aligned_error"]]
    }))

}

process_metrics <- function(input_data) {

  input_data %>%

    mutate(pdb.xyz = map(pdb.xyz, .f = \(x) {
      pdb_chains <- unique(x[["chain"]])
      lut <- setNames(c("rec", paste0("lig", 1:(length(pdb_chains) - 1))), pdb_chains)
      x[["chain"]] <- lut[x[["chain"]]]
      x <- x %>%
        mutate(residue_name = stringr::str_c(chain, "_", AA, resno), .before = everything())
      x
    })) %>%

    mutate(bw_feat = map(`bw: full_table`, ~factor_to_uniprotFeature(.[["protein_segment"]]))) %>%

    mutate(tm_inds = map(bw_feat, \(x) {

      x %>%

        filter(grepl("^TM\\d", type)) %>%

        mutate(mid = (start + end) %/% 2) %>%

        pivot_longer(cols = c("start", "end", "mid")) %>%

        mutate(TM_type = as.numeric(stringr::str_remove(type, "^TM")) %% 2 == 0) %>%

        mutate(TM_name = case_when(TM_type & name == "start" ~ paste0(type, "_IC"),
                                   TM_type & name == "end" ~ paste0(type, "_EC"),
                                   !TM_type & name == "start" ~ paste0(type, "_EC"),
                                   !TM_type & name == "end" ~ paste0(type, "_IC"),
                                   name == "mid" ~ paste0(type, "_mid"))) %>%
        mutate(TM_type2 = stringr::str_remove(TM_name, "^TM\\d_"))

    })) %>%

    mutate(pdb_rec_seq = map_chr(pdb.xyz, \(x) {
      stringr::str_c(x$AA[x$chain == "rec"], collapse = "")
    })) %>%

    mutate(bw_seq = map_chr(`bw: full_table`, \(x) {stringr::str_c(x[["AA"]], collapse = "")})) %>%

    mutate(seq_match = if_else(pdb_rec_seq == bw_seq, "match", "different")) %>%

    mutate(plddt = map(pdb.xyz, \(x) {

      ligands <- unique(x[["chain"]])[-1]

      ligands <- factor_to_uniprotFeature(x[["chain"]]) %>%
        filter(type %in% ligands) %>%
        pivot_longer(c("start", "end")) %>%
        mutate(name = setNames(c("NT", "CT", "mid"), c("start", "end", "mid"))[name]) %>%
        rowwise() %>%
        mutate(pLDDT = case_when(name == "NT" ~ mean(x[["pLDDT"]][value:(value + 4)]),
                                 name == "CT" ~ mean(x[["pLDDT"]][(value - 4):value]))) %>%
        mutate(chain = paste0(type, "_", name)) %>%
        select(chain, pLDDT)

      bind_rows(
        tibble(chain = "all_residues",
               pLDDT = mean(x[["pLDDT"]])),
        x %>%
          group_by(chain) %>%
          summarise(pLDDT = mean(pLDDT)),
        ligands
      ) %>%
        pivot_wider(names_from = chain, names_prefix = "pLDDT_", values_from = "pLDDT")

    })) %>%

    unnest(plddt)
}

compute_RLdists <- function(input_data) {

  input_data %>%
    mutate(RLdists = map2(.x = pdb.xyz, .y = tm_inds, \(.x, .y) {

      ligands <- unique(.x[["chain"]])[-1]
      dist_dat <- .x %>% select(CA_x, CA_y, CA_z)

      ligands <- factor_to_uniprotFeature(.x[["chain"]]) %>%
        filter(type %in% ligands) %>%
        mutate(mid = (start + end) %/% 2) %>%
        pivot_longer(c("start", "end", "mid")) %>%
        mutate(value_plus4 = value + 4,
               value_minus4 = value - 4) %>%
        mutate(value_minus4 = if_else(value_minus4 < 0, 1, value_minus4)) %>%
        rowwise %>%
        mutate(case_when(name == "start" ~ dist_dat %>%
                           slice(value:value_plus4) %>%
                           summarise(across(everything(), mean)),
                         name == "end" ~ dist_dat %>%
                           slice(value_minus4:value) %>%
                           summarise(across(everything(), mean)),
                         name == "mid" ~ dist_dat %>%
                           slice(value))) %>%
        mutate(name = setNames(c("NT", "CT", "mid"), c("start", "end", "mid"))[name]) %>%
        ungroup() %>%
        mutate(ligand_name = paste(type, name, sep = "_"))

      dists_to_comp

      receptors <- .y %>%
        mutate(dist_dat %>% slice(value))

      all_dists <- bind_rows(
        map(1:nrow(dists_to_comp), \(i)
            cross_join(receptors %>% filter(TM_type2 == dists_to_comp$receptor[i]),
                       ligands %>% filter(name == dists_to_comp$ligand[i]), suffix = c(".rec", ".lig"))
        )) %>%
        rowwise() %>%
        mutate(distance = distance(s = c_across(all_of(c("CA_x.rec", "CA_y.rec", "CA_z.rec"))),
                                   p = c_across(all_of(c("CA_x.lig", "CA_y.lig", "CA_z.lig"))))) %>%
        mutate(final_name = paste0(TM_name, "_", ligand_name))

      bind_rows(
        all_dists %>%
          group_by(TM_type2, name.lig, type.lig) %>%
          summarise(distance = mean(distance), .groups = "drop") %>%
          mutate(final_name = paste0(TM_type2, "_", type.lig, "_", name.lig)) %>%
          select(final_name, distance),
        all_dists %>%
          select(final_name, distance)
      ) %>%
        pivot_wider(names_from = final_name, values_from = distance)

    })) %>%

    unnest(RLdists)
}

get_contacts <- function(pw_dist, bw, pdb.xyz, pae) {

  chains <- unname(pdb.xyz[["chain"]])
  uni_chains <- unique(chains)
  rec_lig_pairs <- c(paste("rec", uni_chains[-1], sep = "_"),
                     paste(uni_chains[-1], "rec", sep = "_"))


  contacts <- cbind(expand.grid(residue1 = 1:nrow(pw_dist), residue2 = 1:nrow(pw_dist)), dist = c(pw_dist)) %>%
    as_tibble %>%
    filter(!is.na(dist)) %>%
    filter(dist > 2 & dist < 5) %>%
    mutate(type1 = chains[residue1]) %>%
    mutate(type2 = chains[residue2]) %>%
    mutate(type = paste(type1, type2, sep = "_")) %>%
    filter(type %in% rec_lig_pairs) %>%
    rowwise %>%
    mutate(paeR = pae[[residue1]][[residue2]]) %>%
    mutate(paeL = pae[[residue2]][[residue1]]) %>%
    ungroup %>%
    mutate(BW = case_when(type1 == "rec" ~ bw[["BW"]][residue1],
                          type2 == "rec" ~ bw[["BW"]][residue2])) %>%
    mutate(AA1 = pdb.xyz[["AA"]][residue1]) %>%
    mutate(AA2 = pdb.xyz[["AA"]][residue2]) %>%
    rowwise %>%
    mutate(favorability = case_when(any(map_lgl(residue_pairs[["favorable"]], \(x) sum(c(AA1, AA2) %in% x) == 2)) ~ 1,
                                    any(map_lgl(residue_pairs[["unfavorable"]], \(x) sum(c(AA1, AA2) %in% x) == 2)) ~ -1,
                                    TRUE ~ 0))


  contacts <- left_join(contacts, bw_align, by = "BW")
  contacts <- left_join(contacts, bw[!is.na(bw$BW), ], by = "BW")

  pae_summary <- contacts %>%
    ungroup %>%
    summarise(across(all_of(c("paeR", "paeL", "favorability")), mean, .names = "{.col}_mean"))

  grped_contacts <- contacts %>%
    group_by(protein_segment) %>%
    summarise(count = n())

  grped_contacts <- left_join(bw_segs, grped_contacts, by = "protein_segment") %>%
    mutate(count = replace_na(count, 0)) %>%
    pivot_wider(names_from = protein_segment, values_from = count, names_prefix = "ligContacts_")

  binary_contact <- bw_align %>%
    mutate(contact = if_else(BW %in% contacts[["BW"]], 1, 0)) %>%
    select(name, contact) %>%
    pivot_wider(names_from = name, values_from = contact)

  bind_cols(pae_summary, grped_contacts, binary_contact) %>%
    mutate(totalCP = sum(contacts[["CP"]] == "CP", na.rm = TRUE), .before = everything())



}
