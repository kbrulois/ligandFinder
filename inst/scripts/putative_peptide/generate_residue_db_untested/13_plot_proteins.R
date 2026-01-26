



library(tidyverse)
library(ggplot2)
library(ggtext)
library(patchwork)
library(ligandFinder)



s_localDir <- "~/peptide_alg/build_residue_db"

species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))

sim_mats <- readRDS(system.file("data/sim_mats.rds", package = "ligandFinder"))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))



if(FALSE) {

all_topo <- secretome$topo %>%
                stringr::str_c(.) %>%
                stringr::str_split(., "", simplify = TRUE) %>%
                `c` %>%
                unique %>%
                `[`(!is.na(.))

all_SS <- secretome$dssp %>%
            map(., \(x) {
              if(is.data.frame(x)) return(stringr::str_c(x[["SS"]]))
              else return(NA)}) %>%
            do.call(`c`, .) %>%
            unique %>%
            `[`(!is.na(.))


all_feat_types <- secretome$features %>%
  map(., \(x) {
    if(is.data.frame(x)) {

      logi <- x[["type"]] == "modified residue"

      return(if_else(logi,
              stringr::str_c(x[["description"]]) %>% stringr::str_remove(., ";.*$"),
              stringr::str_c(x[["type"]])))
    } else {
      return(NA)
      }
    }) %>%
  do.call(`c`, .) %>%
  unique %>%
  `[`(!is.na(.)) %>%
  `[`(!grepl("^phs_", .))

all_feats <- bind_rows(map2(secretome$features, secretome$gene, \(x, y) x %>% mutate(gene = y))) %>% as_tibble

all_feats %>%
  filter(type %in% c("glycosylation site", "modified residue")) %>%
  mutate(description = case_when(grepl("^Phospho", description) ~ "phosphorylation",
                                 grepl("^Sulfo", description) ~ "sulfation",
                                 grepl(" amide$", description) ~ "amidation",
                                 grepl("^Deamidated ", description) ~ "deamidated",
                                 TRUE ~ description)) %>%
  group_by(gene) %>%
  summarise(n_unique = n_distinct(description), .groups = "drop") %>%
  arrange(desc(n_unique))

all_domains <- all_feats %>%
  filter(type == "domain") %>%
  mutate(description = stringr::str_remove(description, " \\d+")) %>%
  reframe(tibble_table(description)) %>%
  arrange(desc(freq))

all_residue_mods <- all_feats %>%
  filter(type %in% c("glycosylation site", "modified residue")) %>%
  mutate(description = case_when(grepl("^Phospho", description) ~ "phosphorylation",
                                 grepl("^Sulfo", description) ~ "sulfation",
                                 grepl(" amide$", description) ~ "amidation",
                                 grepl("^Deamidated ", description) ~ "deamidated",
                                 grepl("hydroxy", description, ignore.case = TRUE) ~ "hydroxylation",
                                 grepl("methyl", description, ignore.case = TRUE) ~ "methylation",
                                 grepl("N-acetyl", description) ~ "acetylation",
                                 type == "glycosylation site" ~ "glycosylation",
                                 TRUE ~ description)) %>%
  reframe(tibble_table(description)) %>%
  arrange(desc(freq))

}






#####rough input

peps_tp <- known_peps %>%
  mutate(gene = setNames(id_map$`Gene Names (primary)`, id_map$`Entry Name`)[uniprot_name]) %>%
  filter(!is.na(gene)) %>%
  {left_join(., secretome %>% select(ami_map_score, gene), by = "gene")} %>%
  filter(!is.na(ami_map_score)) %>%
  mutate(gene = stringr::str_remove(gene, ";.*$")) %>%
  mutate(feature = "new_pep") %>%
  nest(.by = "gene")
genes <- unname(peps_tp[["gene"]])



genes <- c("CCL22", "CXCL12", "BRINP1", "SCGB3A2", "TLR3", "IL1B", "BUB1", "CXCL3", "CEACAM4", "CCL18", "ACTC1", "RPS14", "SEC61G", "HOXC8", "DAD1", "CALM1")
names(genes) <- c(rep("goi", 4), rep("cons_mid_species_high", 3), rep("species_low", 3), rep("cons_high_species_high", 3), rep("cons_high_species_low", 3))

genes <- c("brinp2", "cort", "cxcl14", "cxcl17", "pomc") %>% toupper
names(genes) <- rep("goi", 5)

genes <- c("ANGPTL2", "ANO8", "COCH", "CRTAP", "DHH", "ENAM", "GRM4",
           "ITIH2", "KCNV2", "MCEMP1", "NUCB1", "PLVAP", "SHISA9", "SLC39A6",
           "SMOC2", "TMC7", "TRABD2A", "WNT1")[c(13,17)]


peps_tp <- peps %>%
  mutate(feature = "new_pep") %>%
  nest(.by = "gene")

genes <- peps_tp[["gene"]]

already_run <- c("~/AF2_analysis/new_peps_html", "~/AF2_analysis/known_peps_new") %>%
  lapply(., \(x) list.files(x) %>% stringr::str_remove(., ".html$")) %>%
  do.call(`c`, .)

genes <- genes[!genes %in% c("NPIPB15", already_run)]

#genes <- c("IL5", "KCNV2", "ANO8")








html_file_name <- "~/AF2_analysis/plot_test.html"

plot_dir <- "~/AF2_analysis/all_peps_new/"

plot_title <- "new_peps"

mets <- list(cons = c("blos_wt_all_n", "cons_rs", "blos_wt_mam", "blos_wt_all", "gran_wt_all"),
             af_missense = c("mean_afm", "min_afm"),
             dssp = c("relASA"),
             aa_scores = c("pep_xgb4c", "chem_xgb3c", "pep_nn4c", "chem_nn4c")
)

mets$aa_scores <- c(mets$aa_scores, paste0(mets$aa_scores, "_s6"))

all_mets <- do.call(`c`, mets) %>% unname

names(all_mets) <- rep(names(mets), sapply(mets, length))


secretome <- readRDS("~/AF2_analysis/secretome.rds")

peps_tp <- readRDS("~/AF2_analysis/peps_tp.rds")





dir.create(plot_dir)
unlink(plot_dir)


the_input <- secretome %>%
  filter(!accession %in% c("A0AAG2TCD0", "A0AAG2UXZ5")) %>%
  #filter(gene %in% !!genes) %>%
  mutate(aa_scores = map(aa_scores, \(x) x[, colnames(x) %in% mets[["aa_scores"]]])) %>%
  group_split(gene)

pep_input <- peps_tp %>%
  #filter(gene %in% genes) %>%
  group_split(gene)

input_gene <- map_chr(the_input, \(x) x$gene)

pep_input_gene <- map_chr(pep_input, \(x) x$gene)

identical(input_gene, pep_input_gene)


start <- Sys.time()
#future::plan(strategy = future::multisession(workers = 5))
future::plan(strategy = future::sequential())


furrr::future_walk2(.x = the_input[1:5],
                    .y = pep_input[1:5],
                    make_protein_plot,
                    .options = furrr::furrr_options(
                    packages = c("dplyr", "ggplot2", "patchwork", "ligandFinder")
                    ))

end <- Sys.time()

end - start






















if(FALSE) {
  options(future.globals.maxSize = Inf)

  furrr::future_walk(names(p), \(x) {
    p_width <- num_res[[x]]/7
    p_height <- (num_mets[[x]]  + p[[x]][["sub_title_len"]])/5
    message("saving ", x)
    ggsave(plot = p[[x]][["plot"]], filename = paste0(plot_dir, "/", x, ".svg"), device = "svg", limitsize = FALSE, width = , height = )
  })



ligandFinder::html_slide_show(svg_directory = plot_dir,
                              output_file = html_file_name,
                              frames = genes,
                              categories = names(genes),
                              title = plot_title,
                              columns = 1)







link_text <- c("UniProt",
               "Aminode",
               "chatGPT",
               "GWAS",
               "Literature",
               "Disease Associations",
               "AlphaFoldDB")


furrr::future_walk(names(p), \(x) {

  file_name <- paste0(plot_dir, x, ".html")

  xml <- xml2::read_html(file_name)

  xml2::xml_ns(xml)
  xml2::html_structure(xml)

  xml %>%
    xml2::xml_find_all("/html/body")

  xml %>%
    xml2::xml_find_all("//*[contains(text(), 'KCNV2')]")

  xml %>%
    xml2::xml_find_all(xpath="//p[contains(.,'Links: ')]") %>%
    keep(xml2::xml_text(.) %in% link_text) %>%
    xml2::xml_add_parent("a", "xlink:href" = p[[x]][["links"]][[xml2::xml_text(.)]], target = "_blank")

  xml2::write_xml(xml, file_name)

})




}


