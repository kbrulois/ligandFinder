

msa_names <- paste0("V", 1:22)

names(msa_names) <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
  'S', 'T', 'W', 'Y', 'V', 'X', '-')


library(tidyverse)
library(reticulate)
library(keras3)

np <- import("numpy")

af_dir <- "/Users/kbrulois/Desktop/BDKRB2_human_full_and_CXCL14-61-F"
pkl_dir <- paste0("~/Desktop/result_model_1_multimer_v3_pred_0_parsed_pkl")
pdb_path <- paste0(af_dir, "/unrelaxed_model_1_multimer_v3_pred_0.pdb")

seqs <- readLines(paste0(pkl_dir, "/seqs.txt"))

chain_names <- c("BDKRB2", "CXCL14") #get from file name once formalized

chain <- do.call(c, lapply(seq_along(seqs), \(x) rep(chain_names[x], nchar(seqs[x]))))

seq_cat <- stringr::str_split_1(stringr::str_flatten(seqs), "")


pdb_dat <- parse_pdb(pdb_path = pdb_path) %>%
  nest_by(chain, .key = "pdb_dat") %>%
  ungroup %>%
  mutate(chain = chain_names)


pkl_dat <- tibble(files = list.files(pkl_dir)) %>%
              filter(grepl(".npy$", files)) %>%
              mutate(data = map(files, ~ np$load(paste0(pkl_dir, "/", .))))

funcs <- list(mean = mean,
              sd = sd)

bind_cols(
lapply(names(funcs), \(x) {

apply(keras3::op_softmax(pkl_dat$data[pkl_dat$files == "masked_msa.logits.npy"][[1]]),
      MARGIN = c(2,3),
      funcs[[x]]) %>%
  as_tibble %>%
  rename(all_of(msa_names)) %>%
  mutate(AA = seq_cat,
         gene = chain, .before = everything()) %>%
  nest(.key = paste0("msa", "_", x), .by = gene)

}), .name_repair = "minimal"
) %>% select(which(!duplicated(names(.)))) -> msa_dat



pq_path <- paste0(ligandFinder::get_db_path(), "/residue_db")


res_db <- arrow::open_dataset(source = pq_path)

residue_anno <- res_db %>%
  filter(gene %in% chain_names) %>%
  group_by(gene_grp) %>%
  collect() %>%
  ungroup()


residue_anno <- left_join(residue_anno, msa_dat, by = "gene")
residue_anno <- left_join(residue_anno, pdb_dat %>% rename(gene = chain), by = "gene")



residue_anno <- residue_anno %>%
  mutate(msa_mean = furrr::future_pmap(.l = list(seq2 = sequence_uni, to_map = msa_mean),
                                   .f = map_table)) %>%
  mutate(msa_sd = furrr::future_pmap(.l = list(seq2 = sequence_uni, to_map = msa_sd),
                                  .f = map_table)) %>%
  mutate(pdb_dat = furrr::future_pmap(.l = list(seq2 = sequence_uni, to_map = map(pdb_dat, ~select(., -atom_level))),
                                     .f = map_table))

residue_anno <- residue_anno %>%
  mutate(across(all_of(c("msa_mean", "msa_sd", "pdb_dat")),
                .fns = list(tagtoremove = ~map(.x, .f = ~`[[`(., "ms")),
                            score = ~map_dbl(.x, .f = ~`[[`(., "score"))),
                .unpack = TRUE)) %>%
  select(-msa_mean, -msa_sd, -pdb_dat) %>%
  rename_with(.cols = ends_with("_tagtoremove"), .fn = ~sub("_tagtoremove", "", .))



residue_anno <- residue_anno %>%
  mutate(across(all_of(c("msa_mean", "msa_sd")),
                .fns = ~map(.x, \(x) {
                  x %>%
                    mutate(AA_mean = rowMeans(select(., all_of(names(msa_names)[1:20])), na.rm = TRUE))
                })))




res_sub <- residue_anno %>%
              select(-features, -af_xyz) %>%
              unnest(cols = c("cons", "dssp", "af_missense", "msa_mean", "msa_sd", "pdb_dat"), names_sep = "_") %>%
              #group_by(gene) %>%
              #mutate(residue_name = paste0(gene, "_", pdb_dat_AA, 1:n())) %>% #index by pdb #index by uniprot
              #ungroup() %>%
              filter(!is.na(pdb_dat_x)) %>%
              mutate(NearestDifferentNeighbor(data = tibble(pdb_dat_CA_x,pdb_dat_CA_y,pdb_dat_CA_z),
                         g = gene,
                         k = 3)) %>%
              group_by(gene) %>%
              mutate(residue_name = paste0(gene, "_", pdb_dat_AA, 1:n())) %>% #index by pdb
              ungroup()

stack_interacting_residues <- function(data, idx1, idx2, idx3) {
  interaction_name <- cur_data()[["residue_name"]]
  bind_rows(cur_data(), data[c(idx1, idx2, idx3),]) %>%
    mutate(interaction = interaction_name)
}

to_plot <- res_sub %>%
            filter(dist1 < 8 & gene == chain_names[2]) %>%
            rowwise() %>%
            reframe(stack_interacting_residues(data = res_sub, idx1 = idx1, idx2 = idx2, idx3 = idx3)) %>%
            dplyr::select(starts_with("msa_"),
                          starts_with("cons_frequency_scaled"),
                          starts_with("dssp_relASA"),
                          starts_with("af_missense_"),
                          starts_with("dist"),
                          starts_with("interaction"),
                          starts_with("residue_name")) %>%
            dplyr::select(-ends_with("_AA"), -ends_with("_score")) %>%
            reshape2::melt(id.vars = c("residue_name", "interaction"), variable.name = "params") %>%
            mutate(param_type = str_extract(params, "^[^_]+")) %>%
            mutate(params = str_remove(params, "^[^_]*_")) %>%
            mutate(value = as.numeric(value)) %>%
            group_by(param_type) %>%
            mutate(value = scales::rescale(value, to = c(0,1)))


interaction_plot <- ggplot2::ggplot(to_plot) +
  ggplot2::geom_tile(aes(x = residue_name, y = params, fill = value), show.legend = TRUE) +
  ggh4x::facet_nested(cols = vars(interaction), rows = vars(param_type), scales = "free", space = "free", switch = "y") +
  xlab("") +
  ylab("") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggplot2::scale_fill_viridis_c(option = "H")

ggsave(filename = "~/Desktop/interacting_residues.svg",
       plot = interaction_plot,
       device = svglite::svglite, width = 30, height = 14)










