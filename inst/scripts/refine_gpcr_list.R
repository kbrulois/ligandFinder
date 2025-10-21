
ucsd_gpcrs <- c("ACKR1", "ACKR2", "ACKR3", "ACKR4", "ACTHR", "AGTR1", "AGTR2",
                "APJ", "BKRB1", "BKRB2", "BRS3", "C3AR", "C5AR1", "C5AR2", "CCKAR",
                "CCR10", "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7",
                "CCR8", "CCR9", "CCRL2", "CX3C1", "CXCR1", "CXCR2", "CXCR3",
                "CXCR4", "CXCR5", "CXCR6", "EDNRA", "EDNRB", "G32P1", "G37L1",
                "GALR1", "GALR2", "GALR3", "GASR", "GHSR", "GNRHR", "GNRR2",
                "GP101", "GP119", "GP132", "GP135", "GP139", "GP141", "GP142",
                "GP146", "GP148", "GP149", "GP150", "GP151", "GP152", "GP153",
                "GP160", "GP161", "GP162", "GP171", "GP173", "GP174", "GP176",
                "GP182", "GP183", "GPR12", "GPR15", "GPR17", "GPR18", "GPR19",
                "GPR20", "GPR21", "GPR22", "GPR25", "GPR26", "GPR27", "GPR31",
                "GPR32", "GPR33", "GPR34", "GPR35", "GPR37", "GPR39", "GPR3",
                "GPR42", "GPR45", "GPR4", "GPR52", "GPR55", "GPR61", "GPR62",
                "GPR63", "GPR6", "GPR75", "GPR78", "GPR82", "GPR83", "GPR84",
                "GPR85", "GPR87", "GPR88", "GRPR", "KISSR", "LGR4", "LGR5", "LGR6",
                "LSHR", "MAS1L", "MAS", "MC3R", "MC4R", "MC5R", "MCHR1", "MCHR2",
                "MRGRD", "MRGRE", "MRGRF", "MRGRG", "MRGX1", "MRGX2", "MRGX3",
                "MRGX4", "MSHR", "MTLR", "MTR1L", "NK1R", "NK2R", "NK3R", "NMBR",
                "NMUR1", "NMUR2", "NPBW1", "NPBW2", "NPFF1", "NPFF2", "NPSR1",
                "NPY1R", "NPY2R", "NPY42", "NPY4R", "NPY5R", "NPY6R", "NTR1",
                "NTR2", "OGR1", "OPRD", "OPRK", "OPRM", "OPRX", "OPSG2", "OPSG3",
                "OPSX", "OX1R", "OX2R", "OXYR", "P2RY8", "P2Y10", "PAR1", "PAR2",
                "PAR3", "PAR4", "PKR1", "PKR2", "PRLHR", "PSYR", "PTAFR", "QRFPR",
                "RGR", "RL3R1", "RL3R2", "RXFP1", "RXFP2", "SSR1", "SSR2", "SSR3",
                "SSR4", "SSR5", "TAAR2", "TAAR3", "TAAR5", "TAAR6", "TAAR8",
                "TAAR9", "TRFR", "TSHR", "UR2R", "V1AR", "V1BR", "V2R", "XCR1"
)


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


known_pairs <- get_known_pairs()

known_pairs_tbl <- do.call(rbind, known_pairs) %>% as_tibble

se_known <- data.table::fread("~/AF2_analysis/most_recent/2025_07_23_UCSD508_pairings.csv") %>% as_tibble()

se_receptors <- se_known$rec %>% stringr::str_remove(., "^h")

se_known <- se_known %>%
  mutate(parse_proteins(paste0(rec, "_", lig), delim_proteins = "_", delim_ranges = "x", delim_start_end = "x", num_proteins = 2))

se_known2 <- lapply(1:nrow(se_known), \(x) c(se_known[x, "p1_id"][[1]], se_known[x, "p2_id"][[1]]))

gpcr_sub <- gpcr_list %>%
  filter(grepl("^#", `ecb: Order of runs (priority)`) | uniprot_name %in% ucsd_gpcrs)


dat <- data.table::fread("~/AF2_analysis/bm_update3.csv") %>% as_tibble

dat <- dat %>%
      mutate(site = stringr::str_sub(code, 1, 2)) %>%
      filter(site != "AP")

dat <- dat %>%
 rowwise %>%
           mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_name, p2_name) %in% x) == 2)) ~ "known",
                                         any(map_lgl(se_known2, \(x) sum(c(p1_name, p2_name) %in% x) == 2)) ~ "known",
                                         TRUE ~ "unknown"))


test <- dat %>%
  group_by(p1_name, known_pair) %>%
  summarise(percent_E = 100 * sum(lig1_location == "E", na.rm = TRUE)/(sum(lig1_location %in% c("E", "I"), na.rm = TRUE)),
            mean_IPTM_Eonly = mean(iptm[lig1_location == "E"], na.rm = TRUE),
            max_IPTM_Eonly = max(iptm[lig1_location == "E"], na.rm = TRUE))

test <- test %>%
  pivot_wider(names_from = "known_pair", values_from = all_of(c("percent_E", "mean_IPTM_Eonly", "max_IPTM_Eonly")))

gpcr_sub <- left_join(gpcr_sub, test, by = join_by(uniprot_name == p1_name))


dat <- Matrix::readMM("~/AF2_analysis/most_recent/xg_dat.mtx")

dat <- as.matrix(dat)

feat_names <- data.table::fread("~/AF2_analysis/most_recent/feat_names.csv")

colnames(dat) <- colnames(feat_names)

dat2 <- data.table::fread("~/AF2_analysis/most_recent/xg_dat_anno.csv") %>% as_tibble

shap_bw <- readRDS("~/AF2_analysis/most_recent/shap_bw.rds")

shap_names <- tibble(feat_name = colnames(dat)) %>%
  mutate(ligand_index = stringr::str_extract(feat_name, "_L[^_]*$") %>% stringr::str_remove("^_")) %>%
  mutate(bw_index = stringr::str_remove(feat_name, "_L[^_]*$") %>%
           stringr::str_extract(., "[^_]*$")) %>%
  mutate(feature = stringr::str_extract(feat_name, paste0(score_cols, collapse = "|")))


shap_bw <- shap_bw[13:length(shap_bw)]

shap_names <- shap_names %>%
  filter(ligand_index %in% c("L1-5", "L6-10") & bw_index %in% shap_bw)

all_feats <- unique(shap_names$feature)

for(feat in all_feats) {

dat2[[feat]] <- rowMeans(dat[, shap_names %>% filter(feature == feat) %>% pull(feat_name)])
}


test2 <- dat2 %>%
  group_by(p1_name, known_pair) %>%
  summarise(across(all_feats, mean), n_contacts = nrow(pick(everything())))

test2 <- test2 %>%
  pivot_wider(names_from = "known_pair", values_from = all_of(c(all_feats, "n_contacts")))

gpcr_sub <- left_join(gpcr_sub, test2, by = join_by(uniprot_name == p1_name))




gpcr_sub <- gpcr_sub %>%
  mutate(bw_avail = map_lgl(`bw: full_table`, \(x) {if(is.null(nrow(x))) {FALSE} else {nrow(x) > 0}})) %>%
  mutate(ucsd_gpcrs = uniprot_name %in% ucsd_gpcrs)


gpcr_sub <- gpcr_sub %>%
            mutate(has_known_gpcr_lig = if_else(uniprot_name %in% c(known_pairs_tbl[["V2"]], se_known$p1_id), "yes", "no"), .before = "model_name")




umap_config <- umap::umap.defaults
umap_config$min_dist <- 0.3
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200

dim_red_dat <- gpcr_sub %>% select(c(84:87, 91:108)) %>%   mutate(across(everything(), ~replace_na(., 0))) %>% as.matrix

col_variances <- apply(dim_red_dat, 2, var)

dim_red_dat <- dim_red_dat[, col_variances > 0.1 | !is.nan(col_variances)]

umap_res <- umap::umap(d = dim_red_dat, config = umap_config)

gpcr_sub$umap1 <- umap_res$layout[,1]
gpcr_sub$umap2 <- umap_res$layout[,2]


to_plot <- c("uniprot_name",
             "gene_name_primary", "ecb: Order of runs (priority)",
             "ecb: Exclude due to N term >160AA", "ecb: Prioritization Notes",
             "ecb: Class or type", "percent_E", "mean_IPTM", "max_IPTM", "known",
             "bw_avail", "ucsd_gpcrs")

data.table::fwrite(gpcr_sub %>% select(!where(is.list)), "~/Desktop/gcpr_list_new4.csv")


ggplot2::theme_set(ggplot2::theme_bw())
p <- ggplot2::ggplot(gpcr_sub %>%
                       mutate()
                       mutate(`ecb: Order of runs (priority)` = "ecb: Order of runs (priority)")) +
        ggplot2::geom_tile(mapping = aes(fill = `ecb: Order of runs (priority)`, x = `ecb: Order of runs (priority)`, y = uniprot_name))



        ggnewscale::new_scale_fill() +
  ggplot2::geom_tile(mapping = aes(fill = bw_avail, x = constant, y = uniprot_name))


svglite::svglite(filename = "~/Desktop/gpcr_list.svg", width = 4, height = 40)
print(p)
dev.off()

















