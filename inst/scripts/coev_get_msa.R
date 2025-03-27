
msa_names <- paste0("V", 1:22)

names(msa_names) <- alpha_fold_AA_order



library(reticulate)


np <- import("numpy")

af_dir <- "/Users/kbrulois/Desktop/BDKRB2_human_full_and_CXCL14-61-F"
pkl_dir <- paste0("~/Desktop/result_model_1_multimer_v3_pred_0_parsed_pkl")
pdb_path <- paste0(af_dir, "/unrelaxed_model_1_multimer_v3_pred_0.pdb")

seqs <- readLines(paste0(pkl_dir, "/seqs.txt"))

chain_names <- c("BDKRB2", "CXCL14") #get from file name once formalized

chain <- do.call(c, lapply(seq_along(seqs), \(x) rep(chain_names[x], nchar(seqs[x]))))

seq_cat <- stringr::str_split_1(stringr::str_flatten(seqs), "")

