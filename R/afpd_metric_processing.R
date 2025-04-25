

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

import_raw_metrics <- function(dir_name,
                               run_name = run_id,
                               algorithm = alg,
                               pdb = "_ark_[ru]_.*.pdb$",
                               pae = "_pae_.*.json$",
                               code_model_rank = "_[A-Z]{4}[a-z]{7}_s\\d+m\\d+p\\d+_r\\d+") {

  message("importing raw metrics for ", dir_name)

  files <- list.files(paste(input_path_models, dir_name, sep = "/"))

  models <- tibble(pdb_files = grep(pdb, files, value = TRUE),
                   pae_files = grep(pae, files, value = TRUE),
                   code_model_rank = stringr::str_extract(pdb_files, code_model_rank),
                   rlx = stringr::str_extract(pdb_files, "_ark_[ru]_") %>%
                         stringr::str_remove(., "_ark_") %>%
                         stringr::str_remove(., "_$"),
                   model_c = stringr::str_extract(code_model_rank, "s\\d+m\\d+p\\d+"),
                   rank = stringr::str_extract(code_model_rank, "_r\\d+$") %>% stringr::str_remove(., "^_r"),
                   model_num = stringr::str_extract(model_c, "m\\d+") %>% stringr::str_remove(., "^m"),
                   pred_num = stringr::str_extract(model_c, "p\\d+") %>% stringr::str_remove(., "^p"),
                   model_e = paste0("model_", model_num, "_multimer_v3_", "pred_", pred_num)) %>%
                   mutate(code = stringr::str_extract(code_model_rank, "[A-Z]{4}[a-z]{7}"), .before = "pdb_files") %>%
                   mutate(algorithm = algorithm, .before = "pdb_files") %>%
                   mutate(run_name = run_name, .before = "pdb_files") %>%
              select(-code_model_rank) %>%


  dat <- jsonlite::read_json(paste(input_path_models, dir_name, "ranking_debug.json", sep = "/"))

  dat <- tibble(model_e = names(dat[["iptm"]]),
                iptm = unlist(dat[["iptm"]]),
                `iptm+ptm` = unlist(dat[["iptm+ptm"]]))

  if(sum(models[["model_e"]] %in% dat[["model_e"]]) != nrow(models)) warning("file name models do not match ranking_debug.json")

  pairing_info <- parse_dirname(pairing_dir = dir_name,
                                delim_proteins = "_",
                                delim_ranges = "x",
                                delim_start_end = "x") %>%
                  mutate(parsed_pair = map(parsed_pair,
                                           ~pivot_wider(.,
                                                        names_from=c("protein", "annotation"),
                                                        values_from=value))) %>%
                  unnest(parsed_pair)

  dat <- bind_cols(dat, pairing_info)

  dat <- left_join(models, dat, by = "model_e")

  dat %>%
    mutate(pdb_files = setNames(models[["pdb_files"]], models[["model_e"]])[model_e]) %>%

    mutate(pdb = map(pdb_files, ~bio3d::read.pdb(paste(input_path_models, dir_name, ., sep = "/")))) %>%

    mutate(pdb.xyz = map(pdb, parse_pdb)) %>%

    mutate(pae = map(setNames(models[["pae_files"]], models[["model_e"]])[model_e],
                     \(x) {
      tmp <- jsonlite::read_json(paste(input_path_models, dir_name, x, sep = "/"))
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

      ligand_dat <- factor_to_uniprotFeature(.x[["chain"]]) %>%
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


      receptors <- .y %>%
        mutate(dist_dat %>% slice(value))

      all_dists <- bind_rows(
        map(1:nrow(dists_to_comp), \(i)
            cross_join(receptors %>% filter(TM_type2 == dists_to_comp$receptor[i]),
                       ligand_dat %>% filter(name == dists_to_comp$ligand[i]), suffix = c(".rec", ".lig"))
        )) %>%
        rowwise() %>%
        mutate(distance = distance(s = c_across(all_of(c("CA_x.rec", "CA_y.rec", "CA_z.rec"))),
                                   p = c_across(all_of(c("CA_x.lig", "CA_y.lig", "CA_z.lig"))))) %>%
        mutate(final_name = paste0(TM_name, "_", ligand_name))

      all_dists <- bind_rows(
        all_dists %>%
          group_by(TM_type2, name.lig, type.lig) %>%
          summarise(distance = mean(distance), .groups = "drop") %>%
          mutate(final_name = paste0(TM_type2, "_", type.lig, "_", name.lig)) %>%
          select(final_name, distance),
        all_dists %>%
          select(final_name, distance)
      ) %>%
        pivot_wider(names_from = final_name, values_from = distance) %>%
        mutate(lig1_location = case_when(EC_lig1_mid < IC_lig1_mid ~ "E",
                                         EC_lig1_mid > IC_lig1_mid ~ "I",
                                         TRUE ~ NA), .before = everything())


      receptor_mid <- receptors %>%
                        filter(TM_type2 == "mid") %>%
                        summarise(across(all_of(c("CA_x", "CA_y", "CA_z")), mean))

      if(length(ligands) > 1) {
        stop
      } else {
      lig_dat <- .x %>%
        filter(chain == ligands) %>%
        rowwise %>%
        mutate(mid_dist = distance(s = c_across(all_of(c("CA_x", "CA_y", "CA_z"))),
                                   p = receptor_mid)) %>%
        ungroup

      closest_res <- which.min(lig_dat[["mid_dist"]])

      lig_end <- c(N = closest_res - 1,
                   C = nrow(lig_dat) - closest_res)

      lig_end <- lig_end[which.min(lig_end)]

      lig_end <- tibble(lig1_end = paste0(lig_end, names(lig_end)))

      }

      bind_cols(lig_end, all_dists)

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



