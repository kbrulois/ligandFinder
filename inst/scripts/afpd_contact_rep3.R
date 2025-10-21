

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

run_dirs <- input_path_models <- run_dir <- "/scratch/groups/ebutcher/deorphan/models/benchmarking_test"

run_dirs <- input_path_models <- run_dir <- "/scratch/groups/ebutcher/deorphan/models/top200NCnew"





res <- tibble(afpd_dir_name = fs::dir_ls(input_path_models) %>% stringr::str_remove(., ".tar$") %>% basename(),
               file_parts = map(afpd_dir_name, ~stringr::str_split(., "_", simplify = TRUE))) %>%
  mutate(file_part_len = map_int(file_parts, length)) %>%
  mutate(parse_proteins(afpd_dir_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  mutate(data_files = furrr::future_map(afpd_dir_name, ~fs::dir_ls(paste0(input_path_models, "/", .)) %>% basename())) %>%
  mutate(num_files = furrr::future_map_int(afpd_dir_name, ~length(list.files(paste0(input_path_models, "/", .))))) %>%
  mutate(complete = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "ranking_debug.json", sep = "/"))
  })) %>%
  mutate(complete2 = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "metrics_v1.csv", sep = "/"))
  }))

res <- res %>%
  filter(complete)


res <- res %>%
  mutate(num_E_models = map_int(data_files, \(x) {
    x_sub <- stringr::str_detect(x, ".pdb$")
    grep("_ark_.*_[A-Z]{4}[a-z]{7}_s\\d+m\\d+p\\d+_r\\d+", x[x_sub], value = TRUE) %>%
      stringr::str_detect(., "_E_") %>%
      sum})) %>%
  mutate(v2c_present = map_chr(data_files, \(x) {
    ifelse("metrics_v2c.rds" %in% x, "yes", "no")
  }))

res <- res %>%
  mutate(contact_raw = furrr::future_map(afpd_dir_name, \(x) {

    file <- paste(input_path_models, x, "metrics_v2c.rds", sep = "/")
    if(file.exists(file)) {
      return(readRDS(file))
    } else {
      return("none")
    }
  }))

res <- res %>%
  mutate(contact_good = furrr::future_map_int(contact_raw, \(x) {
    sum(!sapply(x, is.null))
  }))


res <- res %>%
  mutate(contact_exclude = map_lgl(contact_raw, \(x) {
    logi <- all(sapply(x, is.null)) | x[1] == "none"
    if(length(logi) == 0) { logi <- FALSE}
    return(logi)
  })) %>%
  filter(!contact_exclude)


res %>%
  filter(num_E_models == 1 & contact_good == 5) %>%
  pull(data_files)

res %>%
  slice(3) %>%
  pull(contact_raw)


res <- res %>%
  mutate(contact_raw2 = map2(data_files, contact_raw, \(x, y) {

    pdb = "_ark_.*.pdb$"
    code_model_rank = "_[A-Z]{4}[a-z]{7}_s\\d+m\\d+p\\d+_r\\d+"


    models <- tibble(pdb_files = grep(pdb, x, value = TRUE),
                     code_model_rank = stringr::str_extract(pdb_files, code_model_rank),
                     rlx = stringr::str_extract(pdb_files, "_[ru]_") %>%
                       stringr::str_remove(., "^_") %>%
                       stringr::str_remove(., "_$"),
                     model_c = stringr::str_extract(code_model_rank, "s\\d+m\\d+p\\d+"),
                     rank = stringr::str_extract(code_model_rank, "_r\\d+$") %>% stringr::str_remove(., "^_r"),
                     model_num = stringr::str_extract(model_c, "m\\d+") %>% stringr::str_remove(., "^m"),
                     pred_num = stringr::str_extract(model_c, "p\\d+") %>% stringr::str_remove(., "^p"))

    if(length(y) != nrow(models)) {y <- y[1:nrow(models)]
    message("contact_raw and num models don't match for ", models$pdb_files[1])
    message("length contact_raw: ", length(y))
    message("length models: ", nrow(models))}

    models %>%
      mutate(code = stringr::str_extract(code_model_rank, "[A-Z]{4}[a-z]{7}"), .before = "pdb_files") %>%
      mutate(algorithm = "AF2v3", .before = "pdb_files") %>%
      mutate(run_name = "bm", .before = "pdb_files") %>%
      select(-code_model_rank) %>%
      mutate(contact_raw = y) %>%
      filter(!is.null(contact_raw))


  }))


known_pairs <- get_known_pairs()

res <- res %>%
  rowwise %>%
  mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_id, p2_id) %in% x) == 2)) ~ "known",
                                TRUE ~ "unknown"), .after = "afpd_dir_name") %>%
  ungroup %>%
  mutate(known_pair2 = if_else(known_pair == "known", paste0(p1_id, ";", p2_id), NA))



gpcr_cols <- c("uniprot_name",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily")

res <- res %>%
  {left_join(., gpcr_list %>% select(all_of(gpcr_cols)),
             by = join_by(p1_id == uniprot_name))} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "afpd_dir_name")



res <- res %>%
        select(-contact_raw) %>%
  filter(!contact_exclude) %>%
        unnest(contact_raw2)

res <- res %>%
  filter(grepl("_E_", pdb_files)) %>%
  unnest(contact_raw)

yo()

















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

score_cols <- c("paeL",
                "paeR",
                "mean_af_missense_rec",
                "mean_af_missense_lig1",
                "favorability",
                "area_scaled")

lig_lut <- c(setNames(rep("L1-5", 5), paste0("L", 1:5)),
             setNames(rep("L6-10", 5), paste0("L", 6:10)),
             setNames(rep("L11-20", 10), paste0("L", 11:20)),
             setNames(rep("L21-50", 30), paste0("L", 21:50)),
             setNames(rep("L51-277", 250), paste0("L", 51:300)))

lig_lut2 <- c(setNames(paste0("L", 1:10), paste0("L", 1:10)),
              setNames(rep("L11x20", 10), paste0("L", 11:20)),
              setNames(rep("L21x50", 30), paste0("L", 21:50)),
              setNames(rep("L51x277", 250), paste0("L", 51:300)))

bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

res_pairs <- readRDS("/oak/stanford/groups/ebutcher/kevin/res_pairs.rds")

pair_lut <- setNames(res_pairs[["name"]], res_pairs[["pairs"]])

runs_c <- runs_c %>%
  mutate(AA_pair = paste0(AA_rec, AA_lig1)) %>%
  mutate(AA_pair_type = if_else(AA_pair %in% res_pairs[["pairs"]], pair_lut[AA_pair], "unclassified"))

test <- runs_c %>%
  #filter(mean_af_missense_rec > 0.8 & mean_af_missense_lig1 > 0.8) %>%
  mutate(quality = if_else(iptm > 0.8, "good", if_else(iptm < 0.5, "bad", "ok"))) %>%
  filter(!quality == "ok") %>%
  mutate(clash = dist < 2) %>%
  group_by(clash, AA_pair_type, quality, in_pocket, gpcr_family, known_pair) %>%
  summarise(total = n()/length(unique(code)),
            across(all_of(c("area", "dist", "paeL")), ~mean(., na.rm = TRUE), .names = "{.col}_mean"),
            )

test %>%
  filter(AA_pair_type == "hydrophobic_interactions") %>%
  print(n =500)

data.table::fwrite(test, "/oak/stanford/groups/ebutcher/kevin/int_tab.csv")

metrics <- c(13, 37, 54, 57, 58, 59, 100, 101, 118, 139, 182, 187, 188, 189, 190, 191, 325, 327)

runs_c <- left_join(runs_c, res3 %>% select(code, starts_with("data_split")), by = "code")

runs_c %>%
  #filter(data_split_1 %in% c("test", "train") | known_pair == "known") %>%
  select(!!metrics) %>%
  data.table::fwrite(., "/scratch/groups/ebutcher/deorphan/analysis/9M_bm_cons.csv")





test <- runs_c %>%
  filter(in_pocket & AA_pair_type == "hydrophobic_interactions" & protein_segment %in% paste0("TM", 1:7) & dist > 2) %>%
  mutate(critical = case_when(mean_af_missense_rec > 0.7 & mean_af_missense_lig1 > 0.7 ~ "good",
                              TRUE ~ "bad")) %>%
  group_by(code) %>%
  summarise(crital2 = (sum(critical == "good") + 1)/(sum(critical == "bad") + 1))

test2 <- left_join(runs_m, test, by = "code")

yardstick::roc_auc_vec(factor(test2$known_pair), test2$crital2)



res <- runs_c %>%
  filter(dist > 2) %>%
  #mutate(ligand_index = as.numeric(stringr::str_remove(ligand_index, "^L"))) %>%
  mutate(across(all_of(c("paeL", "paeR")), ~ (30 - .)/ 30)) %>%
  #mutate(CP = if_else(CP == "CP", 1, 0)) %>%
  #mutate(sb = if_else(tags == "sb", 1, 0)) %>%
  #mutate(ds = if_else(tags == "ds", 1, 0)) %>%
  mutate(area_scaled = if_else(area > 30, 1, area/30))







#file_path <- "/oak/stanford/groups/ebutcher/kevin/all_con.csv"
#file_path <- "/scratch/groups/ebutcher/deorphan/analysis/all_con.csv"

#data.table::fwrite(res, "/oak/stanford/groups/ebutcher/kevin/ds_cons.csv")

lig_inds <- unique(res[["ligand_index"]]) %>% .[gtools::mixedorder(.)]

data.table::fwrite(res, file_path)




res <- data.table::fread(file_path)


res <- as_tibble(res)

shap_bw <- readRDS("/oak/stanford/groups/ebutcher/kevin/shap_bw.rds")

res2 <- res %>%
  filter(!is.na(BW) & BW %in% bw_align[["BW"]]) %>%
  #filter(!is.na(BW) & BW %in% shap_bw) %>%
  mutate(across(all_of(score_cols), ~replace_na(data = ., replace = 0))) %>%
  mutate(all_con_score = (CP + favorability + paeR + paeL + sb + area + cons_rec + frequency_scaled_lig1 + frequency_scaled_rec) / 8) %>%
  mutate(ligand_grp = lig_lut2[ligand_index]) %>%
  group_by(afpd_dir_name, rank, ligand_grp, BW) %>%
  filter(all_con_score == max(all_con_score))

res2 <- res2 %>%
  filter(!ligand_index %in% paste0("L", 101:300))



res2 <- res2 %>%
  select(all_of(c("afpd_dir_name", "code", "complex_type", "p2_name", "gpcr_family", "known_pair", "p1_name",  "BW", "ligand_grp", score_cols))) %>%
  pivot_wider(names_from = all_of(c("BW", "ligand_grp")), values_from = all_of(score_cols), values_fill = 0)











###do dim reduction on contacts

sp_mat <- res2 %>%
  ungroup %>%
  select(37:ncol(.)) %>%
  as.matrix

col_variances <- apply(sp_mat, 2, var)

sp_mat <- sp_mat[, col_variances > 0.01]
nmf_res <- NMFN::nnmf(x = sp_mat, k = 8)

umap_config <- umap::umap.defaults
umap_config$min_dist <- 0.5
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200
umap_res <- umap::umap(d = nmf_res[["W"]], config = umap_config)

end <- Sys.time()
end


res2 <- res2 %>%
  tibble::add_column(UMAP1nmf = umap_res[["layout"]][,1],
                     UMAP2nmf = umap_res[["layout"]][,2], .after = "afpd_dir_name")


data.table::fwrite(res2, "/oak/stanford/groups/ebutcher/kevin/umap_shap_selected2.csv")


data.table::fwrite(res2 %>%
                     ungroup %>%
                     select(37:ncol(.)), "/scratch/groups/ebutcher/deorphan/analysis/con_rep.csv")










sp_mat <- res2 %>%
  ungroup %>%
  select(9:ncol(.))

sp_mat <- bind_rows(feat_names, sp_mat) %>%
  mutate(across(everything(), ~tidyr::replace_na(., replace = 0)))

sp_mat <- sp_mat %>%
  as.matrix %>%
  Matrix::Matrix(., sparse = TRUE)

anyNA(sp_mat)

Matrix::writeMM(sp_mat, "/scratch/groups/ebutcher/deorphan/analysis/xg_dat_new.mtx")

data.table::fwrite(res2[0, -c(1:8)], "/scratch/groups/ebutcher/deorphan/analysis/feat_names.csv")



generate_data_split <- function(group_data, unknown_fac = 1.8, iterations = 20) {

  for(iter in 1:iterations) {

    message("iteration ", iter)

    ds_name <- paste0("data_split_", iter)

    comps <- group_data %>%
      filter(known_pair == "known") %>%
      distinct(complex_type, .keep_all = TRUE)

    uni_comps <- comps[["complex_type"]]

    n <- length(uni_comps)

    if(n == 0) {
      group_data <- group_data %>%
        mutate(!!ds_name := "unused")
    } else {

      n_train <- case_when(
        n == 1 ~ 1,
        n == 2 ~ 1,
        n %in% 3:5 ~ 2,
        n %in% 6:7 ~ 3,
        TRUE ~ floor(0.6 * n)
      )

      train_knowns <- sample(uni_comps, size = n_train)

      group_data <- group_data %>%
        mutate(!!ds_name := case_when(known_pair == "known" & complex_type %in% train_knowns ~ "train",
                                      known_pair == "known" & !complex_type %in% train_knowns ~ "test",
                                      known_pair == "unknown" ~ "unused"))

      unknowns <- group_data %>%
        filter(known_pair == "known") %>%
        group_by(complex_type, p1_name, !!sym(ds_name)) %>%
        summarize(p1_name_c = length(unique(complex_type)))

      for(x in 1:nrow(unknowns)) {

        unknown <- unknowns[x, ]

        unknown_rec <- group_data %>%
          filter(known_pair == "unknown" & p1_name == unknown[["p1_name"]]) %>%
          distinct(complex_type, .keep_all = TRUE)

        uni_comps <- unknown_rec[["complex_type"]]

        unknown_rec <- sample(uni_comps, size = round(unknown[["p1_name_c"]] * unknown_fac, 0))

        group_data <- group_data %>%
          mutate(!!ds_name := case_when(known_pair == "unknown" & complex_type %in% unknown_rec ~ unknown[[ds_name]],
                                        TRUE ~ .data[[ds_name]]))

      }
    }
  }

  return(group_data)

}


plib <- readRDS("/oak/stanford/groups/ebutcher/kevin/plib.rds")

plib <- plib %>%
  rename(p1_name = p1_id) %>%
  nest_by(p1_name)

res2 <- left_join(res2 %>% select(-c(37:ncol(.))), plib, by = "p1_name")

res2 <- res2 %>%
  rowwise %>%
  mutate(known_pair = case_when(is.null(data) ~ "unknown",
                                any(stringr::str_c(p2_name, "_", p2_range) %in% stringr::str_c(data$p2_id, "_", data$p2_range)) ~ "known",
                                TRUE ~ "unknown"))



res3 <- runs_m %>%
  #select(1:8) %>%
  group_by(gpcr_family) %>%
  group_split() %>%
  map_df(generate_data_split)

res3 <- res3[match(res2[["code"]], res3[["code"]]), ]



data.table::fwrite(res3 %>% select(!where(is.list)), "/scratch/groups/ebutcher/deorphan/analysis/xg_dat_anno.csv")

data.table::fwrite(res2 %>% select(!where(is.list)), "/scratch/groups/ebutcher/deorphan/analysis/xg_dat_anno_new.csv")


res3 <- data.table::fread("/scratch/groups/ebutcher/deorphan/analysis/xg_dat_anno.csv")

feat_names <- data.table::fread("/scratch/groups/ebutcher/deorphan/analysis/feat_names.csv")




shap_names <- res2 %>%
  ungroup %>%
  select(-c(1:8)) %>%
  {tibble(feat_name = colnames(.))}

  shap_names <- tibble(feat_name = colnames(feat_names)) %>%
  mutate(ligand_index = stringr::str_extract(feat_name, "_L[^_]*$") %>% stringr::str_remove("^_")) %>%
  mutate(bw_index = stringr::str_remove(feat_name, "_L[^_]*$") %>%
           stringr::str_extract(., "[^_]*$")) %>%
  mutate(feature = stringr::str_extract(feat_name, paste0(score_cols, collapse = "|"))) %>%
  mutate(ligand_bw = paste0(bw_index, "_", ligand_index)) %>%
  mutate(ligand_bw = paste0(bw_index))

interaction_grps <- unique(shap_names[["ligand_bw"]])

interaction_cons <- map(interaction_grps, \(x) shap_names[["feat_name"]][shap_names[["ligand_bw"]] == x])

jsonlite::write_json(interaction_cons, "/oak/stanford/groups/ebutcher/kevin/interaction_cons.json", pretty = TRUE)




shap_data <- data.table::fread("/scratch/groups/ebutcher/deorphan/analysis/shap_tree_traintest.csv")

bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

shap_names <- left_join(shap_names, bw_align, by = join_by(bw_index == BW)) %>%
  mutate(CP = if_else(CP == "CP", "CP", "other"))

shap_names <- bind_cols(shap_names, shap_data)

shap_names <- shap_names %>%
  pivot_longer(starts_with("V"), names_to = "shap_name", values_to = "shap_value")

to_integrate <- res3 %>%
                  filter(data_split_6 %in% c("train", "test")) %>%
                  mutate(shap_name = paste0("V", 1:nrow(.))) %>%
                  select(shap_name, data_split_6, model_num, rank, code, known_pair, gpcr_family, num_E_models)

shap_names <- left_join(shap_names, to_integrate, by = "shap_name")

data.table::fwrite(shap_names, "/oak/stanford/groups/ebutcher/kevin/shap_final2.csv")










bm_update <- data.table::fread("~/Desktop/bm_update3.csv") %>% as_tibble()
bm_con <- data.table::fread("~/Desktop/bm_contact_score_latest.csv") %>% as_tibble()

bm_update <- bm_update %>%
  select(-c(196:541))

res <- left_join(bm_con, bm_update, by = "code")

res_clean <- res %>%
  select(-ends_with(".y")) %>%
  rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))







