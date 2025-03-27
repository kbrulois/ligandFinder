



plot_coev_ov <- function(coev_res = coev_res,
                         params = c("MI", "COVinv", "COVnorm", "dPWcons"),
                         plot_dir = "~/Desktop/tmp_plots",
                         pdb_names = coev_res[["pdb_names"]],
                         output_filename = paste0(out_path, "/overview_coev_", paste0(chain_names, collapse = "_"),".html")) {

  list2env(coev_res[!names(coev_res) == "pdb_names"], envir = environment())
  
  feats_toplot <- factor_to_uniprotFeature(res_sub 
                                         %>% pull(bw_protein_segment), 
                                         source = "GPCRdb")

to_plot_cons <- res_sub %>% 
  filter(!is.na(residue_name_a3m)) %>%
  select(cons_AA, cons_frequency, cons_doubleSmooth, bw_BW)

to_plot_cons$cons_index <- 1:nrow(to_plot_cons)

p_afm_dat <- res_sub %>% 
  filter(!is.na(residue_name_a3m)) %>%
  select(starts_with("af_missense_")) %>%
  dplyr::rename_with(.fn = ~sub("^af_missense_", "", .)) %>%
  mutate(index = 1:nrow(.), .before = everything()) %>%
  dplyr::rename(mean = mean_af_missense) %>%
  dplyr::select(-AA, -score) %>%
  reshape2::melt(id.vars = "index", variable.name = "afm")


p_size <- nrow(to_plot_cons)

buffer <- 1/(10 * p_size)


dists <- c(paste0("Dist_", pdb_names), paste0("discDist_", pdb_names))
params_gated <- expand.grid(paste0(pdb_names, "_0-8A"), params) %>%
  {paste(.[["Var2"]], .[["Var1"]], sep = "_")}

vars_tp <- c(params, dists, params_gated)

names(vars_tp) <- c(rep("parameters", length(params)), 
                    rep("distance", length(dists)),
                    rep("gated parameters", length(params_gated)))



 cons_cutoff <- 0.2
  is_conserved <- ifelse(to_plot_cons[["cons_frequency"]] < cons_cutoff, "conserved", "not conserved")
  is_conserved[is.na(is_conserved)] <- "not conserved"
  
  p_cons <- ggplot2::ggplot(to_plot_cons) +
    ggrastr::rasterise(ggplot2::geom_col(aes(x = cons_index, y = cons_frequency), width = 0.8),
                       dev = "ragg",
                       dpi = 100)  + 
    ggplot2::geom_line(aes(x = cons_index, y = cons_doubleSmooth), 
                       color = "#8A919799", 
                       linewidth = 2) + 
    geom_hline(yintercept = cons_cutoff, linetype = "dashed", color = "red") +
    geom_text(aes(x = cons_index, y = -0.3, label = cons_AA, color = is_conserved,
                  fontface = ifelse(is_conserved == "conserved", 2, 1)),
              size = 3) +
    ggplot2::scale_color_manual(name = "conservation",
                                values = c(setNames("red", "conserved"),
                                           setNames("black", "not conserved"))) +
    ggplot2::expand_limits(x = c(0, p_size)) +
    scale_x_continuous(expand = expansion(mult = buffer),
                       breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
                       minor_breaks = seq(0, p_size, by = 10)) +
    ggplot2::xlab("") +
    ggplot2::ylab("conservation (Aminode)") +
    ggplot2::theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),  
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "solid"),
      panel.grid.minor.x = element_line(linetype = "dashed")
    )
  
  
  p_afm <- ggplot2::ggplot(p_afm_dat) +
    ggplot2::geom_tile(aes(x = index, y = afm, fill = value), show.legend = FALSE) + 
    ggplot2::scale_fill_viridis_c(option = "H") +
    scale_x_continuous(expand = expansion(mult = buffer)) +
    ylab("per AA AF missense") +
    ggplot2::theme_bw()
  
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
  
  msa_toplot <- a3m %>%
    as.data.frame %>%
    as_tibble %>%
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
    theme(axis.text.x = element_blank(),
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
    ggrastr::rasterise(geom_tile(),
                       dev = "ragg",
                       dpi = 100)  +
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
  
  
  
  
  p_cons2 <- ggplot() + 
    ggseqlogo::geom_logo(apply(a3m, 1, paste0, collapse = ""), 
                         col_scheme = ggseqlogo::make_col_scheme(chars = names(AA_colors),
                                                                 cols = unname(AA_colors))) + 
    ggseqlogo::theme_logo() + 
    scale_x_discrete(expand = expansion(mult = buffer)) +
    theme(xis.text.x = element_blank())
  
  p_feats <- ggplot(data = feats_toplot) +
    ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_segment(aes(x = start, 
                                                                       xend = end, 
                                                                       y = source, 
                                                                       color = type), 
                                                                   linewidth = 4),
                                             colour = "black",
                                             sigma = 1, 
                                             expand = 5),
                       dev = "ragg",
                       dpi = 300)  +
    geom_text(data = to_plot_cons,  aes(x = cons_index, y = -0.3, label = bw_BW, color = is_conserved,
                                        fontface = ifelse(is_conserved == "conserved", 2, 1)),
              size = 3, angle = 90, hjust = 0) +
    scale_color_manual(values = extra_color) + 
    scale_x_discrete(expand = expansion(mult = buffer)) + 
    ggplot2::theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),  # Remove major horizontal grid lines
      panel.grid.minor.y = element_blank(), 
      panel.grid.major.x = element_line(linetype = "solid"),
      panel.grid.minor.x = element_line(linetype = "dashed"),
      legend.position = "right",
      legend.justification = "left",
      legend.title.position = "top"
    )
  
  
  
  design = "
##1
##2
##3
##3
##3
##3
##5
##8
674
674
674
674
674
674
674
674
674
674
674
674
674
674
674
674
674
674
674
"
  
  design2 = "
##1
##2
##3
##3
##3
##3
##5
##8
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
"

lapply(vars_tp, \(param) {
    
    message("plotting ", param)
  
  p <- ggplot(data = co_evol_mat, 
              mapping = aes(x = chain1, y = chain2, fill = !!sym(param))) + 
    ggrastr::rasterise(geom_tile(),
                       dev = "ragg",
                       dpi = 100)  + 
    #geom_tile(data = co_evol_mat %>% filter(dist_desc != "> 8A"), 
    #          mapping = aes(x = chain1, y = chain2, color = dist_desc), 
    #          fill = NA, size = 1.5) +
    xlab("") +
    ylab(chain_names[2]) +
    scale_x_discrete(expand = expansion(mult = buffer)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),   
          panel.grid = element_blank(),         
          panel.border = element_blank())
  
  if(grepl("^discDist_",param)) {
    p <- p + 
      scale_fill_manual(values = dist_colors)
  } else if(grepl("^dPWcons", param)) {
    p <- p + 
      scale_fill_manual(values = cons_colors)
  } else {
    p <- p +
      scale_fill_viridis_c(option = "H")
  }
  
  p_al <- p_cons + p_feats + p + p_msa + p_cons2 + tree_plot + p_msa_qc + p_afm +
    patchwork::plot_layout(design = design, ncol = 2, axes = "collect_x", widths = c(13, 2, 55)) +
    patchwork::plot_annotation(title = paste("Overview of Pairwise Coevolution and Distance Metrics"),
                               subtitle = paste(chain_names, collapse = " "),
                               theme = theme(plot.title = element_text(face = 2,
                                                                       size = 22,
                                                                       hjust = 0),
                                             plot.subtitle = element_text(face = 1,
                                                                          size = 12,
                                                                          hjust = 0)))
  
  ggsave(paste0(plot_dir, "/", param, ".svg"), 
         plot = p_al, 
         device = svglite::svglite, width = 60, height = 55, limitsize = FALSE)
})

html_slide_show(svg_directory = plot_dir,
                output_file = output_filename,
                frames = vars_tp,
                categories = names(vars_tp),
                title = paste0("overview_", paste0(chain_names, collapse = "_")),
                columns = 1)

}

