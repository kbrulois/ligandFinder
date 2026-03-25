


url_exists <- function(url) {
  con <- tryCatch(
    url(url, open = "rb"),
    error = function(e) NULL
  )
  if (is.null(con)) return(FALSE)
  close(con)
  TRUE
}



assembleCXC_anno_feats <- function(features, peps) {

  tmp <- features %>%
    as_tibble %>%
    mutate(type = if_else(type == "disulfide bond", "unidsb", type)) %>%
    filter(case_when(type %in% c("afdsb", "unidsb") ~ type %in% c("afdsb", "unidsb"),
                     source %in% c("phs", "gpcrdb_gtp", "sven", "top200NC", "disha", "uni_pep") ~ source %in% c("phs", "gpcrdb_gtp", "sven", "top200NC", "disha", "uni_pep"))) %>%
    mutate(feature = case_when(type %in% c("afdsb", "unidsb") ~ type,
                               source %in% c("phs", "gpcrdb_gtp", "sven", "top200NC", "disha", "uni_pep") ~ source), .before = everything()) %>%
    mutate(end = case_when(feature == "unidsb" & is.na(end) & stringr::str_detect(description, "Interchain") ~ stringr::str_extract(description, "C-\\d+") %>% stringr::str_remove(., "C-") %>% as.integer(.),
                           TRUE ~ end)) %>%
    filter(!is.na(end)) %>%
    select(feature, start, end)


  dsb <- tmp %>%
    filter(feature %in% c("afdsb", "unidsb")) %>%
    mutate(group = paste0(start, "_", end)) %>%
    group_by(group) %>%
    summarise(feature = if(n() == 2) {"both"} else {first(feature)},
              start = first(start),
              end = first(end))

  tmp <- bind_rows(tmp %>% filter(!feature %in% c("afdsb", "unidsb")),
                   dsb)


  if("start" %in% colnames(peps)) {
  tmp <- bind_rows(tmp, peps)
  }

   tmp <- tmp %>%
    mutate(label = paste0(feature, "_", start, "-", end))

   if("selection" %in% colnames(tmp)) {
     tmp <- tmp %>%
    mutate(label = if_else(feature == "new_pep", paste0(label, "_", selection), label))
     }

   tmp %>%
    mutate(color = anno_feat_colors[feature]) %>%
    select(feature, start, end, color, label) %>%
    tidyr::drop_na()

}


assembleCXC_res_mod <- function(features, sequence_uni) {

  AA_sequence <- sequence_uni %>% stringr::str_split(., "", simplify = TRUE) %>% `c`

  seq_dat <- tibble(AA = AA_sequence) %>%
              mutate(start = row_number())

  tmp <- features %>%
    filter(type %in% c("glycosylation site", "modified residue")) %>%
    mutate(label = case_when(grepl("^Phospho", description) ~ "phosphorylation",
                                   grepl("^Sulfo", description) ~ "sulfation",
                                   grepl(" amide$", description) ~ "amidation",
                                   grepl("^Deamidated ", description) ~ "deamidated",
                                   grepl("hydroxy", description, ignore.case = TRUE) ~ "hydroxylation",
                                   grepl("methyl", description, ignore.case = TRUE) ~ "methylation",
                                   grepl("N-acetyl", description) ~ "acetylation",
                                   type == "glycosylation site" ~ "glycosylation",
                                   TRUE ~ "other"), .before = everything()) %>%
    mutate(color = all_desc_colors[["modification"]][res_mod_lut[label]]) %>%
    mutate(end = start) %>%
    mutate(feature = "modification") %>%
    select(feature, start, end, color, label) %>%
    distinct(.keep_all = TRUE) %>%
    tidyr::drop_na()

  if(nrow(tmp) > 0) {
    tmp <- left_join(tmp, seq_dat, by = "start") %>%
      mutate(label = paste0(AA, start, "-", label))
  }

  tmp %>%
    select(feature, start, end, color, label)

}

assembleCXC_dom <- function(features) {

  tmp <- features %>%
    filter(type == "domain") %>%
    mutate(feature = "domain", label = stringr::str_remove(description, " \\d+")) %>%
    mutate(label = stringr::str_replace_all(label, " ", "_")) %>%
    select(feature, start, end, label)

  if(nrow(tmp) > 0) {
  all_doms <- unique(tmp[["label"]])

  desc_col_scale <- all_desc_colors[["domain"]]
  common_doms <- all_doms[all_doms %in% names(all_desc_colors[["domain"]])]
  uncommon_doms <- all_doms[!all_doms %in% names(all_desc_colors[["domain"]])]
  desc_col_scale <- all_desc_colors[["domain"]][names(all_desc_colors[["domain"]]) %in% common_doms]
  uncommon_doms <- color_scalify(uncommon_doms, colors = desc_colors[!desc_colors %in% desc_col_scale])
  desc_col_scale <- c(desc_col_scale, uncommon_doms)

  tmp <- tmp %>%
          mutate(color = desc_col_scale[label]) %>%
          select(feature, start, end, color, label)
  }
  tmp
}

assembleCXC_SV <- function(features) {

  features %>%
    filter(type == "sequence variant") %>%
    mutate(feature = "SV", label = description, color = "#FED439FF") %>%
    mutate(label = stringr::str_replace_all(label, " ", "_")) %>%
    mutate(label = stringr::str_replace_all(label, ";", "")) %>%
    mutate(label = stringr::str_replace_all(label, ".", "")) %>%
    mutate(label = if_else(label == "", "unknown-SV", label)) %>%
    mutate(end = start) %>%
    select(feature, start, end, color, label) %>%
    tidyr::drop_na()

}

assembleCXC_DBC <- function(features, sequence_uni) {

  tmp <- features %>%
    filter(source == "dbc") %>%
    group_by(type, start) %>%
    summarise(start = start[!is.na(start)][1]) %>%
    mutate(type = stringr::str_extract(type, "^NTC|^CTC")) %>%
    mutate(source = type) %>%
    mutate(end = start) %>%
    ungroup()


}



assembleCXC_topo <- function(topo) {

if(is.na(topo)) {
  return(NULL)
} else {
tmp <- stringr::str_split(topo, "", simplify = TRUE) %>% `c`

tmp <- ligandFinder::factor_to_uniprotFeature(tmp)

if(nrow(tmp) == 0) {
  return(NULL)
} else {
tmp <- tmp %>%
  mutate(color = all_desc_colors[["topo"]][type]) %>%
  mutate(label = topo_lut[type]) %>%
  mutate(label = stringr::str_replace_all(label, " ", "_")) %>%
  mutate(feature = "topo") %>%
  select(feature, start, end, color, label)
return(tmp)
}
}

}

assembleCXC_SS <- function(features) {

  features %>%
    filter(source == "alpha fold") %>%
    mutate(color = all_desc_colors[["SS"]][type]) %>%
    mutate(label = ss_features[type]) %>%
    filter(!label %in% ss_features[c(1,3)]) %>%
    mutate(label = stringr::str_replace_all(label, " ", "_")) %>%
    mutate(feature = "SS") %>%
    select(feature, start, end, color, label)

}

assembleCXC_Dibasic <- function(features) {

  features %>%
    filter(type == "Dibasic" & source == "sites") %>%
    mutate(color = "#8030CC") %>%
    mutate(label = "DB") %>%
    mutate(feature = "Dibasic") %>%
    select(feature, start, end, color, label)

}



generate_cm_script <- function(uni_id,
                               gene,
                               disc_mets,
                               cont_mets = do.call(c, mets) %>% unname,
                               cont_colors = viridis::turbo(n = 20) %>% paste(., collapse = ":"),
                               defattr_path = "~/defattr",
                               defattr_type = "local") {


  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uni_id, "-F1-model_v6.pdb")
  cm_sct <- list(paste0("open ", url))

  def_files <- tibble(mets = cont_mets,
                      files = paste0(defattr_path, "/", uni_id, "/", cont_mets, ".defattr")) %>%
    mutate(exists = if(defattr_type == "url") {map_lgl(files, url_exists)} else {map_lgl(files, file.exists)}) %>%
    filter(exists) %>%
    mutate(good = if(defattr_type == "url") {TRUE} else {map_lgl(files, \(x) {
      tryCatch({tmp <- data.table::fread(x)
      !anyNA(tmp$V3)}, error = function(e) FALSE)
    })}) %>%
    filter(good)

  cm_sct <- c(cm_sct,
              map(def_files[["mets"]], ~paste("alias", ., "color byattribute", ., "palette", cont_colors))
  )

  db_sel <- disc_mets %>%
    filter(feature == "Dibasic") %>%
    mutate(cmd = paste0("select :", start, "-", end, "; ", "color sel purple; ", "show sel atoms; ", "style sel sphere")) %>%
    add_row(cmd = "select clear") %>%
    summarise(cmd = paste(cmd, collapse = "; ")) %>%
    pull(cmd) %>%
    paste("alias", "db", .) %>%
    as.list(.)


  dsb_mets <- disc_mets %>%
                filter(feature %in% c("both", "afdsb", "unidsb"))

  dsb_sel <- dsb_mets %>%
    mutate(cmd = paste0("select :", start, ",", end, "; ", "color sel ",  color, "; ", "show sel atoms; ", "style sel stick;")) %>%
    mutate(cmd = if_else(grepl("^both_", label), cmd, paste0(cmd, "select :", start, ",", end, "; ", "label sel text ", feature, "; ", "label height 0.5"))) %>%
    add_row(cmd = "select clear") %>%
    summarise(cmd = paste(cmd, collapse = "; ")) %>%
    pull(cmd) %>%
    paste("alias dsb ", .) %>%
    as.list(.)



  dsb_dist <- dsb_mets %>%
      mutate(cmd = paste0("select :", start, ",", end, "; ", "distance :", start, "@SG :", end, "@SG")) %>%
      add_row(cmd = "select clear") %>%
      summarise(cmd = paste(cmd, collapse = "; ")) %>%
      pull(cmd) %>%
      paste("alias dsb_dist", .) %>%
      as.list(.)


  disc_mets <- disc_mets %>%
    filter(!feature %in% c("Dibasic", "afdsb", "unidsb", "both")) %>%
    split(., .[["feature"]])


  pep_sel <- purrr::map2(disc_mets, names(disc_mets), \(d_met, feat_name) {

    d_met %>%
      mutate(color = paste0("select :", start, "-", end, "; ", "color sel ", color),
             label = paste0("select :", (start + end) %/% 2, "; ", "label sel text ", label, "; ", "label height 4")) %>%
      add_row(color = "select all; color #D2AF81FF; db; dsb", label = "label delete", .before = 1) %>%
      add_row(color = "select clear", label = "select clear") %>%
      select(color, label) %>%
      summarise(across(everything(), ~paste(., collapse = "; "))) %>%
      mutate(cmd = paste0(c_across(c(color, label)), collapse = "; ")) %>%
      pull(cmd) %>%
      paste("alias", feat_name, .)

  })

  cm_sct <- c(cm_sct, db_sel, dsb_sel, dsb_dist, pep_sel)

  cm_sct <- c(do.call(`c`, cm_sct), "db", "dsb", "view")

  c(cm_sct, map_chr(def_files[["files"]], ~paste("open", .)), "dsb_dist")


}





###tests

#input <- the_input[[1]]

#pep_tp <- pep_input[[1]]

do_chimera_scripts <- function(input, pep_tp) {


  features <- input %>% pull(features) %>% `[[`(1)
  sequence_uni <- input %>% pull(sequence_uni)


  topo <- input %>% pull(topo) %>% `[[`(1)

  peps <- pep_tp %>% `[[`("data") %>% `[[`(1)

  gene_name <- input %>% pull(gene)


  disc_mets <- bind_rows(assembleCXC_Dibasic(features),
                         assembleCXC_anno_feats(features, peps),
                         assembleCXC_res_mod(features, sequence_uni),
                         assembleCXC_dom(features),
                         assembleCXC_SV(features),
                         assembleCXC_topo(topo),
                         assembleCXC_SS(features))

  disc_mets <- disc_mets %>%
                  tidyr::drop_na()


  cxc_file_path <- paste0("~/chimerax_scripts/", gene_name, ".cxc")

  cm_script <- generate_cm_script(uni_id = input %>% pull(accession),
                                  gene = gene_name,
                                  disc_mets = disc_mets)


  writeLines(cm_script,
             con = cxc_file_path)


}



start <- Sys.time()
#future::plan(strategy = future::multisession(workers = 6))
future::plan(strategy = future::sequential())


purrr::walk2(.x = the_input[c(1:5435, 5441:5614)],
             .y = pep_input[c(1:5435, 5441:5614)],
             do_chimera_scripts)

end <- Sys.time()

end - start






start <- Sys.time()
#future::plan(strategy = future::multisession(workers = 6))
future::plan(strategy = future::sequential())


purrr::walk2(.x = the_input[c(1:5435, 5441:5614)],
             .y = pep_input[c(1:5435, 5441:5614)],
             do_chimera_scripts)

end <- Sys.time()

end - start














