



clean_coev_data <- function(x) {
    x %>%  
    filter(unique_pairs) %>%
    filter(!as.character(chain1) == as.character(chain2)) %>%
    mutate(pair_name = if_else(is.na(bw_BW.x),stringr::str_c(str_remove(chain1, "^[^_]*_"), ":", str_remove(chain2, "^[^_]*_")),
                               stringr::str_c("[", bw_BW.x, "] ", str_remove(chain1, "^[^_]*_"), ":", str_remove(chain2, "^[^_]*_")))) %>%
    mutate(pair_name_u = paste0(chain1, ":", chain2)) %>%
    pivot_longer(cols = starts_with("discDist"), 
                 names_to = "pdb_model", 
                 values_to = "discDist") %>%
    mutate(pdb_model = sub("^discDist_", "", pdb_model)) %>%
    mutate(interaction_status = case_when(
      discDist == "> 8A" ~ "NI",
      discDist == "0-2A" ~ "C",
      discDist %in% c("2-4A", "4-6A", "6-8A") ~ "I",
      TRUE ~ NA)) %>%
    mutate(PWcons = -1 * PWcons) %>%
    group_by(pair_type, pdb_model, interaction_status) %>%
    mutate(across(all_of(vars_tp2), list(rank = ~rank(-., ties.method = "first")), .unpack = TRUE)) %>%
    mutate(across(all_of(paste0(vars_tp2, "_rank")), list(top10 = \(x) {if_else(x <= 10, "top10", NA)}), .unpack = TRUE)) %>%
    mutate(PWcons = -1 * PWcons)

}

plot_coev_pp <- function(coev_res = coev_res,
                         params = c("MI", "COVinv", "COVnorm", "dPWcons"),
                         plot_dir = "~/Desktop/tmp_plots",
                         pdb_names = coev_res[["pdb_names"]],
                         extra_residues = "CXCL14_W68",
                         extra_residues_name = "interacting_68W",
                         output_filename = paste0(out_path, "/top10pairs/pp_coev_")) {
  
  list2env(coev_res[!names(coev_res) == "pdb_names"], envir = environment())
  
  pairs_input <- clean_coev_data(co_evol_mat)
  
  out_dir <- dirname(output_filename)
  dir.create(out_dir)

  top10pairs <- pairs_input %>%
  mutate(PWcons = -1 * PWcons) %>%
  pivot_longer(cols = all_of(vars_tp2),
               names_to = "method",
               values_to = "rank_input") %>%
  group_by(pair_type, pdb_model, interaction_status, method) %>%
  slice_max(order_by = rank_input, n = 10, with_ties = FALSE, na_rm = TRUE)



if(!is.null(extra_residues)) {
  
top10pairs <- bind_rows( 
  pairs_input %>%
    ungroup() %>%
    filter((chain1 == extra_residues | chain2 == extra_residues)) %>%
    group_by(pair_type, interaction_status, pair_name_u) %>%
    reframe(method = extra_residues_name,
            pdb_model2 = case_when(
              n() == length(pdb_cols) ~ "all",
              n() == 1 ~ pdb_model,
              TRUE ~ as.character(n())
            ),
            pdb_model = pdb_model,
            pair_name = unique(pair_name)) %>%
    distinct(pdb_model, pair_name_u, .keep_all = TRUE) %>%
    group_by(pair_type, pdb_model, interaction_status, method),
  top10pairs
)
}


rankedPairs <- pairs_input %>% 
  mutate(PWcons = -1 * PWcons) %>%
  pivot_longer(cols = all_of(vars_tp2),
               names_to = "method",
               values_to = "rank_input") %>%
  group_by(pair_type, pdb_model, interaction_status, method) %>%
  mutate(rank = rank(-rank_input))




D <- ape::dist.aa(a3m)
tre <- nj(D)
tre <- ladderize(tre)
tree_plot <- ggtree::ggtree(tre, layout = "roundrect")

msa_ord <- ggtree::get_taxa_name(tree_plot)

branch_names <- unique(a3m_anno[[1]][["class"]])

branch_colors <- extra_color[1:length(branch_names)]

names(branch_colors) <- branch_names

tree_plot <- tree_plot %<+% a3m_anno[[1]][match(msa_ord, a3m_anno[[1]][["query"]]),] + 
  ggtree::geom_tree(mapping = aes(color = class), 
                    layout = "roundrect") +
  ggplot2::scale_color_manual(values = branch_colors) +
  theme(legend.position = "top")

msa_ord <- rev(ggtree::get_taxa_name(tree_plot))

interact_mapper <- c(setNames("Interacting", "I"),
                     setNames("Non-interacting", "NI"))


to_iterate2 <- top10pairs %>%
  ungroup %>%
  select(pair_type, interaction_status, pdb_model) %>%
  distinct %>%
  filter(!is.na(interaction_status)) %>%
  filter(interaction_status == "I")

lapply(nrow(to_iterate2):1, \(z) {
  
  pt <- to_iterate2[["pair_type"]][z]
  IorNI <- to_iterate2[["interaction_status"]][z]
  pdbMod <- to_iterate2[["pdb_model"]][z]
  
  to_plot <- top10pairs %>%
    ungroup %>%
    filter(pair_type == pt,
           interaction_status == IorNI,
           pdb_model == pdbMod) %>%
    filter(!(interaction_status == "NI" & method == "interacting_68W"))
  
  
  plot_dir <- paste0("~/Desktop/tmp_plotz", z)
  
  dir.create(plot_dir)
  
  lapply(1:nrow(to_plot), \(y) {
    
    interacting <- to_plot[["interaction_status"]][y]
    nice_name <- to_plot[["pair_name"]][y]
    uni_pair_name <- to_plot[["pair_name_u"]][y]
    pair <- strsplit(uni_pair_name, ":")[[1]]
    
    message("plotting ", uni_pair_name)
    
    p_method <- ggplot2::ggplot(data = rankedPairs %>%
                                  ungroup %>%
                                  mutate(dRank = case_when(
                                    interaction_status != interacting ~ interact_mapper[interaction_status],
                                    is.na(interaction_status) ~ "not\nin\npdb model",
                                    rank <= 10 ~ "1-10",
                                    rank <= 20 ~ "11-20",
                                    rank <= 30 ~ "21-30",
                                    rank <= 50 ~ "31-50",
                                    rank <= 100 ~ "51-100",
                                    rank > 100 ~ "> 100",
                                  )) %>%
                                  filter(pair_name_u == uni_pair_name),
                                mapping = aes(x = method, y = pdb_model, label = dRank)) +
      ggplot2::geom_tile(fill = "#FFFFFF00", color = "black") +
      ggplot2::geom_text() +
      scale_x_discrete(position = "top") +
      ggplot2::ggtitle("Rank in Other Method/Model combinations") +
      theme(panel.grid = element_blank(),
            panel.background = element_blank())
    
    msa_toplot <- a3m %>%
      as.data.frame %>%
      as_tibble %>%
      select(all_of(pair)) %>%
      mutate(sequence_id = a3m_anno[[1]][["query"]], .before = everything()) %>%
      pivot_longer(-sequence_id) %>%
      mutate(name = factor(name, levels = colnames(a3m))) %>%
      mutate(sequence_id = factor(sequence_id, levels = msa_ord))
    
    p_msa <- ggplot(data = msa_toplot, 
                    mapping = aes(x = name, y = sequence_id, fill = value)) + 
      ggrastr::rasterise(geom_tile(),
                         dev = "ragg",
                         dpi = 100)  +
      ggrastr::rasterise(geom_text(mapping = aes(label = value), size = 2),
                         dev = "ragg",
                         dpi = 200) +
      scale_fill_manual(values = AA_colors) +
      xlab("") +
      ylab("") +
      scale_y_discrete(labels = a3m_anno[[1]][["species"]][match(msa_ord, a3m_anno[[1]][["query"]])]) +
      ggplot2::scale_color_manual(values = branch_colors) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(),
            axis.text.y = element_text(color = branch_colors[a3m_anno[[1]][["class"]]][match(msa_ord, a3m_anno[[1]][["query"]])]),
            panel.background = element_blank(),   
            panel.grid = element_blank(),         
            panel.border = element_blank())
    
    msa_qc <- lapply(a3m_anno, get_qc_data)
    
    msa_qc <- lapply(msa_qc, \(x) {
      x %>%
        select(query, qc, database, status)
    })
    
    msa_qc <- Map(\(x) {colnames(msa_qc[[x]]) <- paste0(x, "_", colnames(msa_qc[[x]])); msa_qc[[x]]}, names(msa_qc))
    
    msa_qc <- bind_cols(msa_qc) %>%
      rename(sequence_id = paste0(names(a3m_anno)[1], "_", "query")) %>%
      select(-ends_with("_query")) %>%
      pivot_longer(-sequence_id) %>%
      mutate(sequence_id = factor(sequence_id, levels = msa_ord)) %>%
      mutate(sequence_id = factor(sequence_id, levels = msa_ord))
    
    qc_names <- unique(msa_qc[["value"]])
    
    qc_colors <- extra_color[1:length(qc_names)]
    
    names(qc_colors) <- qc_names
    
    p_msa_qc <- ggplot(data = msa_qc, 
                       mapping = aes(x = name, y = sequence_id, fill = value)) + 
      geom_tile() + 
      #geom_tile(data = to_plot %>% filter(dist_desc != "> 8A"), 
      #          mapping = aes(x = chain1, y = chain2, color = dist_desc), 
      #          fill = NA, size = 1.5) +
      xlab("") +
      ylab("") +
      scale_fill_manual(values = qc_colors) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
            panel.background = element_blank(),   # Remove grey background
            panel.grid = element_blank(),         # Remove all grid lines
            panel.border = element_blank())
    
    
    to_iterate <- unique(a3m_anno[[1]][["class"]])
    
    if(sum(is.na(to_iterate)) < 1) {
      to_iterate <- c(NA, to_iterate)
    }
    
    p_cons2 <- Map(\(x) {
      message(x)
      subsetter <- a3m_anno[[1]][["class"]] == x
      subsetter[is.na(subsetter)] <- FALSE
      if(is.na(x)) { 
        subsetter <- TRUE
        x <- "all"}
      
      to_return <- tryCatch({
        ggplot() + 
          ggseqlogo::geom_logo(apply(a3m[subsetter,pair], 1, paste0, collapse = ""), 
                               col_scheme = ggseqlogo::make_col_scheme(chars = names(AA_colors),
                                                                       cols = unname(AA_colors))) + 
          ggtitle(x) +
          ggseqlogo::theme_logo() + 
          scale_x_discrete(expand = expansion(mult = buffer))}, error = function(e) {
            ggplot()
          })
      return(to_return)
    }, to_iterate)
    
    p_cons2 <- Reduce(`+`, p_cons2) +   
      patchwork::plot_layout(nrow = 1, axes = "collect_x")
    
    p_cons2 <- ggplotify::as.ggplot(p_cons2)
    
    
    
    
    p_cor <- ggplot(msa_toplot) +
      geom_line(mapping = aes(x = sequence_id, y = value, color = name, group = name)) +
      theme_bw() +
      ylab("") +
      theme(axis.text.x = element_blank(),
            panel.grid.major.x = element_blank(),  
            panel.grid.minor.x = element_blank(),
            legend.position = "left")
    
    
    
    design = "
1115555
666####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
234####
"
    
    
    dist_dat <- co_evol_mat %>% filter(chain1 == pair[1], chain2 == pair[2]) %>% select(starts_with("Dist_")) %>% unlist
    
    dist_dat <- paste(paste("", names(dist_dat), round(unname(dist_dat), 2), ""), collapse = "")
    
    p_al <- p_cons2 + tree_plot + p_msa_qc + p_msa + p_cor + p_method +
      patchwork::plot_layout(design = design, axes = "collect_x") +
      patchwork::plot_annotation(title = paste("Top", interact_mapper[interacting], "Residue Pairs for", pt, " from pdb file: ", pdbMod),
                                 subtitle = paste0("Current Pair Category: ", to_plot$method[y],
                                                   " Current pair: ", to_plot$pair_name[y],
                                                   " MI: ", co_evol_mat %>% filter(chain1 == pair[1], chain2 == pair[2]) %>% pull(MI) %>% round(., 2),
                                                   " PWcons: ", co_evol_mat %>% filter(chain1 == pair[1], chain2 == pair[2]) %>% pull(PWcons) %>% round(., 2),
                                                   " COVinv: ", co_evol_mat %>% filter(chain1 == pair[1], chain2 == pair[2]) %>% pull(COVinv) %>% round(., 2),
                                                   dist_dat),
                                 theme = theme(plot.title = element_text(face = 2,
                                                                         size = 22,
                                                                         hjust = 0),
                                               plot.subtitle = element_text(face = 1,
                                                                            size = 12,
                                                                            hjust = 0)))
    
    ggsave(paste0(plot_dir, "/", y, ".svg"), 
           plot = p_al, 
           device = svglite::svglite, width = 20, height = 55, limitsize = FALSE)
    gc()
    
  })
  
  
  
  html_slide_show(svg_directory = plot_dir,
                  output_file = paste0(output_filename, sub(":", "_", pt), "_", IorNI, "_", pdbMod, ".html"),
                  frames = 1:nrow(to_plot),
                  buttonNames = to_plot$pair_name,
                  categories = to_plot$method,
                  title = paste0("topPairs_", pt, "_", IorNI, "_", pdbMod),
                  columns = 2)
  
})

}



