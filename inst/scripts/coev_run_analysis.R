


###parse a3m AF3 from Irina

test <- jsonlite::read_json(path = "~/peptide_alg/models/hAGTR1dT_hAGTR2dT_hBKRB1dT_hBKRB2dT_lnkCXCL14x52x68w_af3_tscc_data/hagtr1dt_hagtr2dt_hbkrb1dt_hbkrb2dt_lnkcxcl14x52x68w_data.json")

msa_names <- c("AGTR1", "AGTR2", "BKRB1", "BKRB2", "CXCL14")

a3ms <- lapply(seq_along(msa_names), \(x) {
  
  message("Processing ", msa_names[x])
  a3m <- strsplit(test$sequences[[5]]$protein$pairedMsa, "\\n")[[1]]
  
  a3m_len <- length(a3m) 
  
  tibble(id = a3m[seq(1, a3m_len, by = 2)] %>% sub("^>", "", .),
         sequence = a3m[seq(2, a3m_len, by = 2)] %>% gsub("[a-z]", "", .),
         gene = msa_names[x])
  
})







### set input

bellestros_mapping = c(setNames("bkrb1", "BDKRB1"), #todo figure out what these gene symbols are
                       setNames("bkrb2", "BDKRB2"),
                       setNames("agtr1", "AGTR1"),
                       setNames("agtr2", "AGTR2"),
                       setNames("ccr9", "CCR9"))

bellestros_mapping_rev <- names(bellestros_mapping)
names(bellestros_mapping_rev) <- toupper(unname(bellestros_mapping))

pdb_location <- "~/peptide_alg/models/hAGTR1dT_hAGTR2dT_hBKRB1dT_hBKRB2dT_lnkCXCL14x52x68w_af3_tscc_data/"

to_iterate <- tibble(pdb_files = list.files(pdb_location)) %>%
              filter(grepl(".pdb$", pdb_files)) %>%
              mutate(gene = stringr::str_extract(pdb_files, "^[^_]+")) %>%
              mutate(gene = bellestros_mapping_rev[gene]) %>%
              mutate(pdb_files = paste0(pdb_location, "/", pdb_files)) %>%
              mutate(a3m_files = paste0("~/peptide_alg/models/", gene, "/CXCL14_68W_pairgreedy/pair.a3m")) %>%
              group_by(gene) %>%
              summarise(pdb_files = list(pdb_files),
                        a3m_files = unique(a3m_files))


if(x == "AGTR2") {
pdb_path <- paste0("~/peptide_alg/models/", x, "/CXCL14_68W_rlx_r1m4.pdb")
} else {
  pdb_path <- paste0("~/peptide_alg/models/", x, "/CXCL14_68W_rlx_r1m1.pdb")
}

pdb_path <- c(pdb_path, paste0("~/peptide_alg/models/", x, "/CXCL14_68W_BKRB1_lnkCXCL14x52x68W_rnk4_JCC.pdb"))

pdb_path <- paste0("~/peptide_alg/models/CCR9_CCL25xP1_keep_pkl/ranked_", 0:4, ".pdb")

pdb_path <- c("~/peptide_alg/models/BDKRB2/BDKRB2_CXCL14_68W_mult3_8seeds_Dropout_calciptm_Feb15_2025__5_83332_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_006.pdb #1 RELAXED separately_Feb15.pdb",
              "~/peptide_alg/models/BDKRB2/CXCL14_68W_rlx_r1m1.pdb",
              "~/peptide_alg/models/BDKRB2/CXCL14_68W_rlx_r2m4.pdb")

to_iterate %>%
  filter(gene == x) %>%
  mutate(pdb_files = paste0(pdb_location, "/", pdb_files)) %>%
  pull(pdb_files) -> pdb_path






all_gpcrs <- c("BDKRB2", "CCR9", "BDKRB1", "ANGTR1", "ANGTR2") #get list





coev_res <- compute_coev(a3m_file = "~/peptide_alg/boltz_results_CCR9_ccl25_short/msa/CCR9_ccl25_short_paired_tmp_pairgreedy-env/pair.a3m", 
                         pdb_path = paste0("~/peptide_alg/models/CCR9_CCL25xP1_keep_pkl/ranked_", 0:4, ".pdb"),
                         chain_names = c("CCR9", "CCL25"))

plot_coev_ov(coev_res = coev_res,
             pdb_names = coev_res[["pdb_names"]][1:2])

plot_coev_pp(coev_res = coev_res)

plot_coev_viol(coev_res = coev_res)




coev_res <- compute_coev(a3m_file = "~/peptide_alg/", 
                         pdb_path = paste0("~/peptide_alg/CXCL12_CXCR4/8u4o.pdb"),
                         chain_names = c("CXCR4", "CXCL12"))

plot_coev_ov(coev_res = coev_res,
             pdb_names = coev_res[["pdb_names"]][1:2])

plot_coev_pp(coev_res = coev_res)

plot_coev_viol(coev_res = coev_res)







for(i in 1:nrow(to_iterate)) {
  
  gene_name <- to_iterate$gene[i]
  
  chain_names = c(gene_name, "CXCL14")
  
  message("doing coev for ", gene_name)

coev_res <- compute_coev(a3m_file = to_iterate %>% filter(gene == gene_name) %>% pull(a3m_files), 
                         pdb_path = to_iterate %>% filter(gene == gene_name) %>% pull(pdb_files) %>% unlist,
                         chain_names = c(gene_name, "CXCL14"))

plot_coev_ov(coev_res = coev_res,
             pdb_names = coev_res[["pdb_names"]][1:2])

plot_coev_pp(coev_res = coev_res,
             pdb_names = coev_res[["pdb_names"]][1:2])

plot_coev_viol(coev_res = coev_res)

}




lapply(names(bellestros_mapping)[c(1,2,4)], do_coevol)




co_evol_mat <- lapply(names(bellestros_mapping), \(x) {
  tmp <- readRDS(paste0("~/peptide_alg/", x, ".rds"))[["co_evol_mat"]]
  tmp %>%
    mutate(pair = paste0(x, "_CXCL14"))
})

co_evol_mat <- bind_rows(co_evol_mat)




