


.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)

#source("path/to/utils.R")


dir_path <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/pocForGrant"

getwd()
setwd("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts")

analysis_name <- "pocForGrant"
out_file_name <- paste0(analysis_name, "_qc_dists")


dir.create(analysis_name)
setwd(analysis_name)

dir.create("input")
setwd("input")

dists_to_comp <- tibble(receptor = c("EC", "IC", "mid", "mid"),
                        ligand = c("mid", "mid", "CT", "NT"))

#gpcr_list <- readRDS("~/oak/deorphan-AI-ze/uniprot/gpcrdb_receptor_list_KB_250304_final.rds")
gpcr_list <- readRDS("/oak/stanford/groups/ebutcher/deorphan-AI-ze/uniprot/gpcrdb_receptor_list_KB_250304_final.rds")

#id_map <- readRDS("~/oak/deorphan-AI-ze/uniprot/id_mapping.rds")
id_map <- readRDS("/oak/stanford/groups/ebutcher/deorphan-AI-ze/uniprot/id_mapping.rds")


comp_jobs <- tibble(file_name = list.files(dir_path))

comp_jobs <- comp_jobs %>%
  mutate(parse_proteins(file_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x"))


test <- comp_jobs %>%
  mutate(num_files = map_int(file_name, ~length(list.files(paste0(dir_path, "/", .))))) %>%
  filter(num_files == get_mode(num_files)) %>%
  filter(p1_id %in% (gpcr_list %>% filter(bw_avail == "available") %>% pull(uniprot_name)))

num_of_grps <- 16
test <- test %>%
  mutate(group = paste0("job", ntile(n = num_of_grps), ".rds"))



test %>%
  group_by(group) %>%
  group_walk(~ saveRDS(.x, file = .y$group))

jobs <- unique(test$group)

setwd('..')

future::plan(strategy = future::multicore(workers = 16))

start <- Sys.time()

furrr::future_map(jobs, \(job) {
  
  tryCatch({

test <- readRDS(file = paste0("input/", job))

test2 <- test %>%
  mutate(model = map(file_name, ~list.files(paste0(dir_path, "/", .)) %>% 
                                            stringr::str_subset(., "pred_0.pdb$"))) %>%
  mutate(iptm = map(file_name, \(x) {
    tryCatch({
    tmp <- jsonlite::read_json(paste(dir_path, x, "ranking_debug.json", sep = "/"))
    tmp %>% as_tibble %>% unnest(everything()) %>% select(-order)
    }, error = function(e) {rep(NA, 5)})
  })) %>%
  unnest(cols = c("model", "iptm")) %>%
  mutate(pdb = map2(.x = file_name, .y = model, ~bio3d::read.pdb(paste(dir_path, .x, .y, sep = "/")))) %>%
  mutate(pdb.xyz = map(pdb, parse_pdb)) %>%
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name) %>% select(p1_id, `bw: full_table`), by = "p1_id")} %>%
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

test3 <- test2 %>%
  filter(seq_match == "match") %>%
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

test3 <- bind_rows(test3, 
                   test2 %>% 
                     filter(seq_match == "different"))


saveRDS(test3, file = job)
message("completed ", job)

}, error = function(e) message("problem with ", job))

})

end <- Sys.time()

end - start

res <- bind_rows(
  lapply(jobs, readRDS)
)



umap_config <- umap::umap.defaults
umap_config$min_dist <- 0.5
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200

subsetter <- res %>%
  select(starts_with("iptm"), 
         starts_with("TM"), 
         starts_with("pLDDT"),
         starts_with("mid_"), 
         starts_with("EC_"), 
         starts_with("IC_")) %>%
  select(!where(is.list)) %>%
              rowwise() %>%
              mutate(subsetter = anyNA(c_across(everything()))) %>%
          pull(subsetter)

to_umap <- res %>%
  select(starts_with("iptm"), 
         starts_with("TM"), 
         starts_with("pLDDT"),
         starts_with("mid_"), 
         starts_with("EC_"), 
         starts_with("IC_")) %>%
  select(!where(is.list)) %>%
  filter(!subsetter) %>%
  as.matrix
  

umap_res <- umap::umap(d = to_umap, config = umap_config)

res[["UMAP1"]] <- NA
res[["UMAP2"]] <- NA

res[["UMAP1"]][!subsetter] <- umap_res[["layout"]][,1]
res[["UMAP2"]][!subsetter] <- umap_res[["layout"]][,2]

gpcr_cols <- c("p1_id", 
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily",
               "bw_avail")

res <- res %>% 
  relocate(starts_with("UMAP"), .after = "model") %>%
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)), 
             by = "p1_id")} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "model")


saveRDS(res, paste0(out_file_name, ".rds"))

data.table::fwrite(res %>% 
                     select(!where(is.list)), 
                   paste0(out_file_name, ".csv"))




