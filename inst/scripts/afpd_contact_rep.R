



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)

#source("path/to/utils.R")
#source("path/to/name_parsing.R")

#dir_path <- "~/oak/deorphan-AI-ze/models/pocForGrant"
dir_path <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/pocForGrantdT"

getwd()
#setwd("~/peptide_alg/contact_test")
setwd("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts")

analysis_name <- "pocForGrantdT3"
out_file_name <- paste0(analysis_name, "_contacts")


dir.create(analysis_name)
setwd(analysis_name)
dir.create("input")

dists_to_comp <- tibble(receptor = c("EC", "IC", "mid", "mid"),
                        ligand = c("mid", "mid", "CT", "NT"))

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

bw_align <- summarize_bw()

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))


comp_jobs <- tibble(file_name = list.files(dir_path))

if(any(grepl("_and_", comp_jobs[["file_name"]]))) { stop("files need to be renamed") }


if(FALSE) {

  all(grepl("_and_", comp_jobs[["file_name"]]))

  tmp <- comp_jobs %>%
    mutate(parse_proteins(file_name,
                          delim_proteins = "_and_",
                          delim_ranges = "_",
                          delim_start_end = "-",
                          p1_range_type = "dT")) %>%
    mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
    mutate(new_file_name = paste(paste0("h", p1_id, p1_range_type), paste(p2_id, p2_range, sep = "x"), sep = "_"))

  tmp %>% select(file_name, p1_id, p2_id, new_file_name)


  ###execute with caution
  tmp %>%
    mutate(file.rename(paste0(dir_path, "/", file_name),
                       paste0(dir_path, "/", new_file_name)))


  comp_jobs <- tibble(file_name = list.files(dir_path))


  rm(tmp)


  tmp <- comp_jobs %>%
    mutate(new_file_name = paste(paste0(p1_id), paste(p2_id, sub("-", "x", p2_range), sep = "x"), sep = "_"))


}

comp_jobs <- comp_jobs %>%
  mutate(parse_proteins(file_name,
                        delim_proteins = "_",
                        delim_ranges = "x",
                        delim_start_end = "x",
                        p1_range_type = "dT")) %>%
  mutate(p1_id = sub("^h", "", p1_id)) %>%
  mutate(p1_id = sub("dT$", "", p1_id))


comp_jobs_sub <- comp_jobs %>%
  mutate(num_files = map_int(file_name, ~length(list.files(paste0(dir_path, "/", .))))) %>%
  filter(num_files == get_mode(num_files)) %>%
  filter(p1_id %in% (gpcr_list %>% filter(bw_avail == "available") %>% pull(uniprot_name)))

num_of_grps <- 16
comp_jobs_sub <- comp_jobs_sub %>%
  mutate(group = paste0("job", ntile(n = num_of_grps), ".rds"))



comp_jobs_sub %>%
  group_by(group) %>%
  group_walk(~ saveRDS(.x, file = paste0("input/", .y$group)))

jobs <- unique(comp_jobs_sub$group)

future::plan(strategy = future::multicore(workers = 16))

start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(jobs, \(job) {

  tryCatch({

    to_do <- readRDS(file = paste0("input/", job))

    metrics <- import_raw_metrics(input_data = to_do)

    metrics <- left_join(metrics,
                             gpcr_list %>%
                               rename(p1_id = uniprot_name) %>%
                               select(p1_id, `bw: full_table`),
                             by = "p1_id")

    metrics <- process_metrics(input_data = metrics)

    metrics <- bind_rows(
                        compute_RLdists(input_data = metrics %>%
                                    filter(seq_match == "match")),
                        metrics %>%
                         filter(seq_match == "different"))

    metrics <- metrics %>%
                  mutate(pw_dist = map(pdb, bio3d::dm.pdb))



    metrics <- bind_rows(
                  metrics %>%
                    filter(seq_match == "match") %>%
                     mutate(contacts = pmap(list(pw_dist = pw_dist,
                                  bw = `bw: full_table`,
                                  pdb.xyz = pdb.xyz,
                                  pae = pae), get_contacts)) %>%
                                  unnest(contacts),
                    metrics %>%
                           filter(seq_match == "different"))


    saveRDS(metrics %>% select(!where(is.list)), file = job)
    message("completed ", job)

  }, error = function(e) message("problem with ", job))

})

end <- Sys.time()

end - start

yo()

nums <- list(1:4, 5:8, 9:12, 13:16)

for(x in seq_along(nums)) {

test <- lapply(paste0("job", x, ".rds"), readRDS)

test2 <- lapply(test, \(x) {x %>% select(!where(is.list))})

saveRDS(bind_rows(test2), paste0("consolidated", x, ".rds"))

}

res2 <- bind_rows(
  lapply(paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/pocForGrant7/consolidated", 0:3, ".rds"), readRDS)
)


res <- bind_rows(
  lapply(jobs, readRDS)
)

res <- bind_rows(res, res2)



known_pairs <- list(c("GPR25", "CXL17"),
                    c("CCR9", "CCL25"),
                    c("GPR15", "GP15L"),
                    c("CML1", "RARR2"),
                    c("CML2", "RARR2"),
                    c("CCRL2", "RARR2"))

res <- res %>%
  rowwise %>%
  mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_id, p2_id) %in% x) == 2)) ~ "known",
                                TRUE ~ "unknown"), .after = "model") %>%
  ungroup



res <- res %>%
  mutate(ligand_location = case_when(EC_lig1_mid < IC_lig1_mid ~ "EC",
                                     EC_lig1_mid > IC_lig1_mid ~ "IC",
                                     TRUE ~ "bw_not_available"), .after = "model") %>%
  mutate(totalCP = replace_na(totalCP, 0))







col_types <- list(all_bw = bw_align[["name"]],
                  cp_bw = grep("_CP$", bw_align[["name"]], value = TRUE),
                  orient = c("EC_lig\\d_mid", "IC_lig\\d_mid", "mid_lig\\d_CT", "mid_lig\\d_NT"),
                  af_qc = c("^pLDDT", "^pae", "iptm", "iptm+ptm"),
                  sum_contacts = c("^ligContacts", "totalCP"),
                  sumcon_cp_bw = c(grep("_CP$", bw_align[["name"]], value = TRUE), c("^ligContacts", "totalCP")))




pca_res <- list()
nmf_res <- list()
umap_res <- list()

umap_config <- umap::umap.defaults
umap_config$min_dist <- 0.5
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200



for(col_type in names(col_types)) {

subsetter <- res %>%
  select(ligand_location,
         matches(col_types[[col_type]])) %>%
  select(!where(is.list)) %>%
  rowwise() %>%
  mutate(subsetter = anyNA(c_across(-ligand_location))) %>%
  mutate(subsetter2 = ligand_location == "IC") %>%
  mutate(subsetter = subsetter | subsetter2) %>%
  pull(subsetter)

dim_red_input <- res %>%
  select(matches(col_types[[col_type]])) %>%
  select(!where(is.list)) %>%
  ungroup %>%
  dplyr::filter(!subsetter) %>%
  as.matrix


pca_res[[col_type]] <- prcomp(t(dim_red_input), rank. = 10)

for(i in 1:ncol(pca_res[[col_type]][["rotation"]])) {

  dim_red_name <- paste0("PC", i, "_", col_type)

  res[[dim_red_name]] <- NA

  res[[dim_red_name]][!subsetter] <- pca_res[[col_type]][["rotation"]][,i]

}



nmf_res[[col_type]] <- NMFN::nnmf(x = dim_red_input, k = 8)

for(i in 1:ncol(nmf_res[[col_type]][["W"]])) {

  dim_red_name <- paste0("NMF", i, "_", col_type)

  res[[dim_red_name]] <- NA

  res[[dim_red_name]][!subsetter] <- nmf_res[[col_type]][["W"]][,i]

}


umap_res[[col_type]] <- umap::umap(d = nmf_res[[col_type]][["W"]], config = umap_config)

for(i in 1:ncol(umap_res[[col_type]][["layout"]])) {

  dim_red_name <- paste0("nmfUMAP", i, "_", col_type)

  res[[dim_red_name]] <- NA

  res[[dim_red_name]][!subsetter] <- umap_res[[col_type]][["layout"]][,i]

}



}




gpcr_cols <- c("p1_id",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily",
               "bw_avail")

res <- res %>%
  relocate(starts_with("PC"), .after = "model") %>%
  relocate(starts_with("UMAP"), .after = "model") %>%
  relocate(starts_with("NMF"), .after = "model") %>%
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_id")} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "model")






saveRDS(pca_res, paste0(out_file_name, "_pca.rds"))
saveRDS(umap_res, paste0(out_file_name, "_umap.rds"))
saveRDS(nmf_res, paste0(out_file_name, "_nmf.rds"))


data.table::fwrite(res %>%
                     select(!where(is.list)),
                   paste0(out_file_name, ".csv"))


saveRDS(res, paste0(out_file_name, "_res.rds"))

message("scp kbrulois@dtn.sherlock.stanford.edu:", getwd(), "/", out_file_name, ".csv", " ~/Desktop/", out_file_name, ".csv")



