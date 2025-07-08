

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)
num_of_grps <- 16

future::plan(strategy = future::multicore(workers = num_of_grps))

run_dirs = c("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking",
             "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking_APACE")

run_dirs <- input_path_models <- run_dir <- "/scratch/groups/ebutcher/deorphan/models/benchmarking"

res <- bind_rows(map(run_dirs, ~get_metrics(run_dir = .,
                                            file = "metrics_v2.csv",
                                            reader = data.table::fread)))



gpcr_cols <- c("uniprot_name",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily")

res <- res %>%
  {left_join(., gpcr_list %>% select(all_of(gpcr_cols)),
             by = join_by(p1_name == uniprot_name))} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "run_name")


res <- res %>%
  mutate(contact_raw = furrr::future_map(afpd_dir_name, \(x) readRDS(paste(input_path_models, x, "metrics_v2c.rds", sep = "/"))))


known_pairs <- get_known_pairs()

res <- res %>%
  rowwise %>%
  mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_name, p2_name) %in% x) == 2)) ~ "known",
                                TRUE ~ "unknown"), .after = "iptm") %>%
  ungroup %>%
  mutate(known_pair2 = if_else(known_pair == "known", paste0(p1_name, ";", p2_name), NA))


res <- res %>%
  group_by(afpd_dir_name) %>%
  mutate(contact_raw = contact_raw[[1]][row_number()])

res <- res %>%
  rename(mean_pLDDT_lig1 = pLDDT_lig1,
         mean_pLDDT_rec = pLDDT_rec) %>%
  unnest(contact_raw)


score_cols <- c("paeL",
                "paeR",
                "pLDDT_rec",
                "pLDDT_lig1",
                "frequency_scaled_lig1",
                "frequency_scaled_rec",
                "mean_af_missense_rec",
                "mean_af_missense_lig1",
                "favorability",
                "CP",
                "sb",
                "ds",
                "area_scaled")

res <- res %>%
  filter(dist > 2) %>%
  #mutate(ligand_index = as.numeric(stringr::str_remove(ligand_index, "^L"))) %>%
  mutate(across(all_of(c("paeL", "paeR")), ~ (30 - .)/ 30)) %>%
  mutate(CP = if_else(CP == "CP", 1, 0)) %>%
  mutate(sb = if_else(tags == "sb", 1, 0)) %>%
  mutate(ds = if_else(tags == "ds", 1, 0)) %>%
  mutate(area_scaled = if_else(area > 30, 1, area/30))


if(FALSE) {
res_sub <- res %>%
              group_by(afpd_dir_name) %>%
              filter(known_pair == "known")


res_sub2 <- res %>%
  group_by(afpd_dir_name) %>%
  filter(known_pair == "unknown")


to_sample <- res_sub2 %>%
  distinct(afpd_dir_name) %>%
  pull(afpd_dir_name)

known_size <- res_sub %>%
                distinct(afpd_dir_name) %>%
                pull(afpd_dir_name) %>%
                length

unknown_samp <- sample(to_sample, known_size * 5)

res_sub2 <- res_sub2 %>%
  filter(afpd_dir_name %in% unknown_samp)

res <- bind_rows(res_sub, res_sub2)

}

data.table::fwrite(res, "/oak/stanford/groups/ebutcher/kevin/ds_cons.csv")

contact_rep <- res %>%
  filter(!is.na(BW)) %>%
  ungroup %>%
  mutate(all_con_score = CP * favorability * paeR * paeL * sb * area * cons_rec * frequency_scaled_lig1 * frequency_scaled_rec) %>%
  group_by(afpd_dir_name, rank, ligand_index, BW) %>%
  filter(all_con_score == max(all_con_score)) %>%
  distinct(afpd_dir_name, rank, ligand_index, BW, .keep_all = TRUE) %>%
  ungroup %>%
  select(all_of(c(colnames(res)[1:76], "BW", "ligand_index", score_cols))) %>%
  #group_by(afpd_dir_name, rank, ligand_index) %>%
  pivot_wider(names_from = "BW", values_from = score_cols) %>%
  mutate(across(!any_of(c(gpcr_cols, "code", "afpd_dir_name", "rank", "ligand_index")), ~replace_na(data = ., replace = 0)))

test %>%
  ungroup %>%
  dplyr::summarise(n = dplyr::n(), .by = c(afpd_dir_name, rank, ligand_index, BW)) %>%
  dplyr::filter(!is.na(BW)) %>%
  dplyr::filter(n > 1L)

data.table::fwrite(contact_rep, "/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_lig1-15.csv")

contact_rep2 <- res %>%
  filter(!is.na(BW)) %>%
  ungroup %>%
  mutate(across(all_of(score_cols), ~replace_na(data = ., replace = 0))) %>%
  mutate(all_con_score = CP * favorability * paeR * paeL * sb * area * cons_rec * frequency_scaled_lig1 * frequency_scaled_rec) %>%
  group_by(afpd_dir_name, rank, ligand_index) %>%
  filter(all_con_score == max(all_con_score)) %>%
  ungroup %>%
  distinct(afpd_dir_name, rank, ligand_index, .keep_all = TRUE)

contact_rep <- contact_rep %>%
                  mutate(ligand_res_id = paste0(afpd_dir_name, rank, ligand_index))

contact_rep2 <- contact_rep2 %>%
  mutate(ligand_res_id = paste0(afpd_dir_name, rank, ligand_index))



umap_res <- data.table::fread("/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_UMAP_eu_n30.csv")

umap_res <- bind_cols(ligand_res_id = contact_rep$ligand_res_id, umap_res)

to_save <- left_join(umap_res %>% as_tibble, contact_rep2, by = "ligand_res_id")


umap_res2 <- data.table::fread("/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_UMAP_eu.csv")
umap_res3 <- data.table::fread("/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_UMAP_eu_gated.csv")
umap_res4 <- data.table::fread("/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_UMAP_jac_gated.csv")

gated_umaps <- bind_cols(umap_res3, umap_res4)

subsetter <- contact_rep$ligand_index %in% paste0("L", 1:15)

contact_rep <- contact_rep %>%
  mutate(umap_res, .after = "run_name") %>%
  mutate(umap_res2, .after = "run_name")


for(i in colnames(gated_umaps)) {

  contact_rep[[i]] <- NA
  contact_rep[[i]][subsetter] <- gated_umaps[[i]]

}



umap_file <- "/oak/stanford/groups/ebutcher/kevin/con_rep_60kby3k_UMAP5.csv"

data.table::fwrite(to_save, umap_file)

message(paste('scp', paste0("kbrulois@dtn.sherlock.stanford.edu:", umap_file), "~/Desktop"))


dim_red_input <- contact_rep %>%
  ungroup %>%
  select(matches(score_cols)) %>%
  as.matrix




subsetter <- res %>%
  ungroup %>%
  select(all_of(score_cols)) %>%
  mutate(subsetter = if_any(everything(), is.na)) %>%
  pull(subsetter)

dim_red_input <- res %>%
  ungroup %>%
  select(all_of(score_cols)) %>%
  dplyr::filter(!subsetter) %>%
  as.matrix





umap_config <- umap::umap.defaults
umap_config$n_neighbors <- 30
umap_config$min_dist <- 0.5
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200



umap_res <- umap::umap(d = dim_red_input, config = umap_config)

for(i in 1:ncol(umap_res[["layout"]])) {

  dim_red_name <- paste0("UMAP", i)

  res[[dim_red_name]] <- NA

  res[[dim_red_name]][!subsetter] <- umap_res[["layout"]][,i]


}
