


file_name <- "~/AF2_analysis/db_analysis_N-term5.svg"

plot_title <- "N-terminal Cleavage Sites"

scg3 <- tibble(position = 83,
               terminus = "end",
               gene = "SCGB3A2",
               source = "other",
               has_db = TRUE,
               seq_dbcenter = "SEAVKKLLEA")


scg3 <- tibble(position = 70,
               terminus = "start",
               gene = "SCGB3A2",
               source = "other",
               has_db = TRUE,
               seq_dbcenter = "VEGLRKCVNE")

scg4 <- tibble(position = c(86, 83, 28),
               terminus = "end",
               gene = c("SCGB1D2", "SCGB2B2", "SCGB1A1"),
               source = "other",
               has_db = TRUE,
               seq_dbcenter = c("VKILKKCSV-", "SVVIKKILQS", "GAQLKKLVDT"))

scg3 <- bind_rows(scg3, scg4)




align_input <- test %>%
  filter(terminus %in% c("start")) %>%
  #filter(source != "uniprot") %>%
  filter(has_db) %>%
  distinct(seq_dbcenter, .keep_all = TRUE) %>%
  pull(seq_dbcenter) %>%
  unique()



a3m <- test %>%
  #{bind_rows(., scg3)} %>%
  filter(terminus %in% c("start")) %>%
  #filter(source != "uniprot") %>%
  filter(has_db) %>%
  distinct(seq_dbcenter, .keep_all = TRUE) %>%
  mutate(seq_name = paste(gene, terminus, position, sep = "_"))

a3m_mat <- do.call(rbind, strsplit(a3m[["seq_dbcenter"]], split = ""))
rownames(a3m_mat) <- a3m[["seq_name"]]
colnames(a3m_mat) <- c(paste0("P", 10:1), paste0("P", 1:10, "'"))





D <- ape::dist.aa(a3m_mat)
tre <- ape::nj(D)
tre <- ape::ladderize(tre)
tree_plot <- ggtree::ggtree(tre, layout = "roundrect")

msa_ord <- ggtree::get_taxa_name(tree_plot)

tree_plot <- tree_plot +
  ggtree::geom_tree(layout = "roundrect") +
  theme(legend.position = "top")

msa_ord <- rev(ggtree::get_taxa_name(tree_plot))


msa_toplot <- a3m_mat %>%
  as.data.frame %>%
  as_tibble %>%
  mutate(sequence_id = rownames(a3m_mat), .before = everything()) %>%
  pivot_longer(-sequence_id) %>%
  mutate(name = factor(name, levels = colnames(a3m_mat))) %>%
  mutate(sequence_id = factor(sequence_id, levels = msa_ord))

p_msa <- ggplot(data = msa_toplot,
                mapping = aes(x = name, y = sequence_id, fill = value)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(aes(label = value), size = 1.8, color = "black") +#%>%
  #ggfx::with_outer_glow(., colour = "white", sigma = 0.8, expand = 3, blend_type = "add") %>%
  #ggrastr::rasterise(., dev = "ragg", dpi = 300) +
  scale_fill_manual(values = AA_colors) +
  xlab("") +
  ylab("") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 90, size = 6, hjust = 0),
        plot.margin = margin(0, 0, 0, 0))


p_consensus <- ggplot() +
  ggseqlogo::geom_logo(align_input) +
  scale_y_continuous(limits = c(0, 1)) +
  ggseqlogo::theme_logo() +
  theme(legend.position = "right", axis.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0))


design = "
#2
13
13
13
13
13
13"

svglite::svglite(file_name, height = 22, width = 9)
tree_plot + p_consensus + p_msa +
  patchwork::plot_layout(design = design) +
  patchwork::plot_annotation(title = plot_title) &
  theme(plot.margin = margin(0,0,0,0))
dev.off()

