
###extra stuff


five <- five %>% 
  as_tibble %>%
  mutate(new = if_else(EC_lig1_mid < IC_lig1_mid, "EC", "IC"))


order_by_dist <- res_sub %>%
  group_by(gene) %>%
  dplyr::arrange(dist1, .by_group = TRUE) %>%
  select(residue_name, dist1) %>%
  mutate(residue_name2 = str_extract(residue_name, "(?<=_).*"))

co_evol_mat <- co_evol_mat %>%
  mutate(chain1 = factor(chain1, levels = order_by_dist %>% 
                           filter(gene == chain_names[1]) %>% 
                           pull(residue_name2)),
         chain2 = factor(chain2, levels = order_by_dist %>% 
                           filter(gene == chain_names[2]) %>% 
                           pull(residue_name2)))

dist_mat <- res_sub %>%
  select(gene, pdb_dat_x,pdb_dat_y,pdb_dat_z) %>%
  {split(., .[["gene"]])} %>%
  {proxy::dist(x = .[[2]] %>% select(-gene), y = .[[1]] %>% select(-gene))} %>%
  {colnames(.) <- res_sub %>% filter(gene == chain_names[[2]]) %>% pull(residue_name); .}

dist_mat %>% as.table() %>% dim


ggplot(data = co_evol_mat, mapping = aes(x = chain1, y = chain2, fill = MI)) + 
  geom_tile() + 
  xlab(chain_names[1]) +
  ylab(chain_names[2]) +
  scale_fill_viridis_c(option = "H")








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


