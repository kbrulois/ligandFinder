
compute_coev <- function(a3m_file = "~/peptide_alg/boltz_results_CCR9_ccl25_short/msa/CCR9_ccl25_short_paired_tmp_pairgreedy-env/pair.a3m", 
                         pdb_path = paste0("~/peptide_alg/models/CCR9_CCL25xP1_keep_pkl/ranked_", 0:4, ".pdb"), 
                         pdb_names = stringr::str_extract(pdb_path, ".{8}(?=\\.pdb$)") %>% stringr::str_remove_all(., "_"),
                         out_path = dirname(a3m_file),
                         pq_path = "~/peptide_alg/residue_db",
                         chain_names = c("CCR9", "CCL25")) {
  
  names(pdb_path) <- pdb_names
  
  pdb_dat <- Map(\(x) {
    parse_pdb(pdb_path = x) %>% 
      mutate(residue_name = paste0(AA, resno)) %>% 
      split(., .[["chain"]])
  }, pdb_path)
  
  #get residue data
  
  res_db <- arrow::open_dataset(source = pq_path)
  
  repeat {
    
    residue_anno <- res_db %>% 
      filter(gene %in% chain_names) %>%
      group_by(gene_grp) %>%
      collect() %>%
      ungroup() %>%
      mutate(gene = factor(gene, levels = chain_names)) %>%
      arrange(gene)
    
    for(x in names(pdb_dat)) {
      residue_anno <- residue_anno %>%
        mutate(!!paste0("pdb_dat_", x) := furrr::future_pmap(.l = list(seq2 = sequence_uni, to_map = map(pdb_dat[[x]], ~select(., -atom_level))),
                                                             .f = map_table))
    }
    
    pdb_cols <- paste0("pdb_dat_", names(pdb_dat))
    
    residue_anno <- residue_anno %>%
      mutate(across(all_of(pdb_cols), 
                    .fns = list(tagtoremove = ~map(.x, .f = ~`[[`(., "ms")),
                                score = ~map_dbl(.x, .f = ~`[[`(., "score"))),
                    .unpack = TRUE)) %>%
      select(-all_of(pdb_cols)) %>%
      rename_with(.cols = ends_with("_tagtoremove"), .fn = ~sub("_tagtoremove", "", .))
    
    
    chain_order_good <- residue_anno %>%
      summarise(across(all_of(paste0(pdb_cols, "_score")), .f = \(x) !all(x < 0)))
    
    if(all(chain_order_good)) {
      message("chain order is correct")
      break
    } else if(all(!chain_order_good)) {
      message("chain order is incorrect. swapping order")
      chain_names <- chain_names[c(2,1)]
    } else {
      message("chain order does not match between pdb files. Stopping")
    }
    
  }
  
  gpcrdb_gene_names <- setNames(c(bellestros_mapping[chain_names[1]], "CXCL14"), chain_names) ##temp patch
  
  bw_res <- Map(\(x) getBallesterosFromGPCRDB(gene_name = x), gpcrdb_gene_names)
  
  
  residue_anno <- left_join(residue_anno, tibble(gene = names(bw_res),
                                                 bw = bw_res), by = "gene")
  
  residue_anno <- residue_anno %>%
    mutate(bw = pmap(.l = list(seq2 = sequence_uni, to_map = bw),
                     .f = map_table))
  
  
  
  residue_anno <- residue_anno %>%
    mutate(across(all_of("bw"), 
                  .fns = list(tagtoremove = ~map(.x, .f = ~`[[`(., "ms")),
                              score = ~map_dbl(.x, .f = ~`[[`(., "score"))),
                  .unpack = TRUE)) %>%
    select(-bw) %>%
    rename_with(.cols = ends_with("_tagtoremove"), .fn = ~sub("_tagtoremove", "", .))
  
  
  res_sub <- residue_anno %>%
    select(-af_xyz) %>%
    unnest(cols = all_of(c("cons", "dssp", "af_missense", pdb_cols, "bw")), names_sep = "_") %>%
    group_by(gene) %>%
    mutate(residue_name = paste0(gene, "_", af_missense_AA, 1:n())) %>% #index by pdb #index by uniprot
    ungroup() %>%
    mutate(has_na = rowSums(is.na(select(., all_of(paste0(pdb_cols, "_x"))))) == length(pdb_cols)) %>%
    filter(!is.na(has_na))
  #mutate(NearestDifferentNeighbor(data = tibble(pdb_dat_x,pdb_dat_y,pdb_dat_z), 
  #                                g = gene, 
  #                                k = 3))
  
  
  
  
  
  
  
  
  
  
  a3m <- readLines(a3m_file, skipNul = TRUE)
  
  a3m_len <- length(a3m)
  
  a3m <- tibble(id = a3m[seq(1, a3m_len, by = 2)] %>% sub("^>", "", .),
                sequence = a3m[seq(2, a3m_len, by = 2)] %>% gsub("[a-z]", "", .),
                gene = stringr::str_extract(id, "^\\d{3}")) %>%
    fill(gene, .direction = "down") %>%
    mutate(gene = setNames(chain_names, c("101", "102"))[gene]) %>%
    split(., .[["gene"]])
  
  a3m_anno_file <- paste0(dirname(a3m_file), "/", "a3m_anno.rds")
  
  if(file.exists(a3m_anno_file)) {
    a3m_anno <- readRDS(a3m_anno_file)
  } else {
    a3m_anno <- Map(\(x) do_species_mapping(x[["id"]]), a3m)
    saveRDS(a3m_anno, a3m_anno_file)
  }
  
  a3m <- Map(\(x) {
    
    residue_seq <- a3m[[x]][["sequence"]][1]
    row_names <- a3m[[x]][["id"]]
    data <- a3m[[x]] %>% 
      select(sequence) %>%
      rowwise() %>%
      mutate(sequence = strsplit(sequence, split = "")) %>%
      {do.call(rbind, .[["sequence"]])}
    
    matrix(data = data,
           nrow = length(row_names), 
           ncol = nchar(residue_seq), 
           dimnames = list(row_names, 
                           paste0(x, "_", paste0(strsplit(residue_seq, split = "")[[1]], 1:nchar(residue_seq)))))
    
  }, names(a3m)) %>%
    {do.call(cbind, .)}
  
  
  
  res_names <- list(grep(chain_names[1], colnames(a3m), value = TRUE),
                    grep(chain_names[2], colnames(a3m), value = TRUE))
  
  
  a3m <- cbind(a3m[,grepl(chain_names[1], colnames(a3m))], a3m[,grepl(chain_names[2], colnames(a3m))])
  
  for(ind in seq_along(pdb_cols)) {
    
    pdb_col <- paste0(pdb_cols[ind], "_residue_name")
    
    
    message("trying to match a3m colnames with ", pdb_col)
    
    res_sub <- res_sub %>%
      mutate(residue_name_a3m = paste0(gene, "_", !!sym(pdb_col))) %>%
      mutate(residue_name_a3m = if_else(residue_name_a3m %in% colnames(a3m), residue_name_a3m, NA))
    
    is_match <- identical(colnames(a3m), res_sub %>% filter(!is.na(residue_name_a3m)) %>% pull(residue_name_a3m))
    
    if(is_match) {
      message("Sucessfully matched a3m colnames with ", pdb_col)
      break
    } else if(ind < length(pdb_cols)) {
      message("Failed to match. Trying ", pdb_cols[ind + 1])
    } else {
      message("Failed to match a3m data with any of the pdb files")
      
      res_sub <- left_join(tibble(residue_name_a3m = colnames(a3m)), res_sub, by = "residue_name_a3m")
    }
    
  }
  
  
  
  future::plan(strategy = future::sequential())
  
  co_evol_mat <- expand.grid(chain1 = do.call(c, res_names),
                             chain2 = res_names[[2]], 
                             stringsAsFactors = FALSE) %>%
    as_tibble %>%
    mutate(MIog = furrr::future_map2_dbl(.x = chain1, 
                                         .y = chain2, 
                                         .f = \(x, y) {
                                           DescTools::MutInf(x = a3m[, x], 
                                                             y = a3m[, y])
                                         })) 
  
  unique_combos <- apply(expand.grid(res_names[[1]], res_names[[2]]), 1, paste, collapse = "_")
  unique_combos2 <- apply(combn(res_names[[2]], m = 2), 2, paste, collapse = "_")
  
  unique_combos <- c(unique_combos, unique_combos2)
  
  co_evol_mat <- co_evol_mat %>% 
    mutate(unique_pairs = paste0(chain1, "_", chain2) %in% unique_combos)
  
  to_exclude <- is.na(res_sub[["residue_name_a3m"]])
  
  for(pdb_col in pdb_cols) {
    coord_cols <- paste0(pdb_col, "_", c("x", "y", "z"))
    co_evol_mat <- co_evol_mat %>%
      mutate(!!sub("pdb_dat_", "Dist_", pdb_col) := furrr::future_map2_dbl(.x = chain1, 
                                                                           .y = chain2, 
                                                                           .f = \(x, y) {
                                                                             
                                                                             distance(s = res_sub[res_sub[["residue_name_a3m"]] == x & !to_exclude, coord_cols] %>% 
                                                                                        unlist(., use.names=FALSE), 
                                                                                      p = res_sub[res_sub[["residue_name_a3m"]] == y & !to_exclude, coord_cols] %>% 
                                                                                        unlist(., use.names=FALSE))
                                                                             
                                                                           }))
    
  }
  
  library(infotheo)
  
  start <- Sys.time()
  
  co_evol_mat <- co_evol_mat %>%
    mutate(MI = furrr::future_map2_dbl(.x = chain1, 
                                       .y = chain2, 
                                       .f = \(x, y) {
                                         mi_value <- mutinformation(a3m[, x],
                                                                    a3m[, y])
                                         hx <- entropy(a3m[, x])
                                         hy <- entropy(a3m[, y])
                                         mi_value / sqrt(hx * hy)
                                       })) 
  
  end <- Sys.time()
  
  end - start
  
  
  aa_map <- setNames(1:21, c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
                             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"))
  
  cov_mat <- apply(a3m, 2, function(col) as.numeric(aa_map[col]))
  
  cov_mat <- cov(cov_mat)
  to_remove <- apply(cov_mat, 1, \(x) all(is.na(x)))
  cov_inv <- MASS::ginv(cov_mat[!to_remove,!to_remove])
  
  dimnames(cov_inv) <- list(colnames(cov_mat)[!to_remove],
                            colnames(cov_mat)[!to_remove])
  
  cov_inv_cols <- colnames(cov_inv)
  
  
  co_evol_mat <- co_evol_mat %>%
    mutate(COVinv = furrr::future_map2_dbl(.x = chain1, 
                                           .y = chain2, 
                                           .f = \(x, y) {
                                             if(x %in% cov_inv_cols & y %in% cov_inv_cols) {
                                               return(abs(cov_inv[x, y]))
                                             } else {
                                               return(NA)
                                             }
                                           })) 
  
  
  
  
  
  
  
  
  
  
  
  co_evol_mat <- co_evol_mat %>%
    mutate(across(starts_with("dist_"), ~factor(case_when(
      .x < 2  ~ "0-2A",
      .x < 4  ~ "2-4A",
      .x < 6  ~ "4-6A",
      .x < 8  ~ "6-8A",
      is.na(.x) ~ NA,
      TRUE       ~ "> 8A"
    ), levels = c("0-2A", "2-4A", "4-6A", "6-8A", "> 8A")), .names = "disc{.col}", 
    .unpack = TRUE)) 
  
  
  
  co_evol_mat <- left_join(co_evol_mat, res_sub %>% 
                             rename(chain1 = residue_name_a3m), 
                           by = "chain1")
  
  co_evol_mat <- left_join(co_evol_mat, res_sub %>% 
                             filter(gene == chain_names[2]) %>% 
                             rename(chain2 = residue_name_a3m), 
                           by = "chain2")
  
  
  co_evol_mat <- co_evol_mat %>%
    mutate(PWcons = cons_frequency.x * cons_frequency.y) %>%
    mutate(PWcons = log(PWcons + 0.1))  %>%
    mutate(PWcons = PWcons - min(PWcons, na.rm = TRUE) + 0.1) %>%
    mutate(COVnorm = COVinv/PWcons) %>%
    mutate(dPWcons = factor(case_when(
      PWcons < 1  ~ "< 1",
      PWcons < 2  ~ "1-2",
      PWcons < 3  ~ "2-3",
      PWcons < 4  ~ "3-4",
      TRUE         ~ "> 4"
    ), levels = c("< 1", "1-2", "2-3", "3-4", "> 4")))
  
  co_evol_mat <- co_evol_mat %>%
    mutate(pair_type = paste(stringr::str_extract(chain1, "^[^_]+"),
                             stringr::str_extract(chain2, "^[^_]+"), sep = ":"))
  
  
  
  
  gate_on_dist <- function(var, dist_desc) {
    var_name <- deparse(substitute(var))
    tibble(`0-2A` = if_else(dist_desc == "0-2A", var, NA),
           `0-4A` = if_else(dist_desc %in% c("0-2A", "2-4A"), var, NA),
           `0-6A` = if_else(dist_desc %in% c("0-2A", "2-4A", "4-6A"), var, NA),
           `0-8A` = if_else(dist_desc %in% c("0-2A", "2-4A", "4-6A", "6-8A"), var, NA))
  }
  
  for(dist in grep("^discDist_", colnames(co_evol_mat), value = TRUE)) {
    suffix <- sub("^discDist_", "", dist)
    co_evol_mat <- co_evol_mat %>%
      mutate(across(all_of(c("MI", "COVinv", "COVnorm", "dPWcons")), 
                    .fns = ~gate_on_dist(var = ., !!sym(dist)), 
                    .names = "{.col}_{suffix}",
                    .unpack = TRUE))
  }
  
  co_evol_mat <- co_evol_mat %>% mutate(chain1 = factor(paste0(chain1),
                                                        levels = c(do.call(c, res_names))),
                                        chain2 = factor(chain2,
                                                        levels = res_names[[2]]))
  
  
  
  coev_res <- list(co_evol_mat = co_evol_mat,
                   a3m = a3m,
                   a3m_anno = a3m_anno,
                   res_sub = res_sub,
                   pdb_names = pdb_names,
                   out_path = out_path)
  
  saveRDS(coev_res, paste0(out_path, "/coev_res.rds"))
  
  return(coev_res)
  
}



