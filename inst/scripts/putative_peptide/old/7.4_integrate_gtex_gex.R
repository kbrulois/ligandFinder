


dat <- data.table::fread("~/oak/gtexxx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")

dat <- as_tibble(dat)

dat[["Name"]] <- sub('\\.[0-9]*$', '', dat[["Name"]])

gene_lut <- ensembldb::select(EnsDb.Hsapiens.v79,
                              keys= dat[["Name"]],
                              keytype = "GENEID",
                              columns = c("GENEID", "SYMBOL", "GENEBIOTYPE"))

colnames(gene_lut)[1] <- "Name"

gene_lengths <- ensembldb::lengthOf(EnsDb.Hsapiens.v79)

gene_lut$gene_length <- unname(gene_lengths[gene_lut[["Name"]]])/1000

dat <- lazy_dt(dat)

dat <- dat %>%
  group_by(Name) %>%
  summarise(across(starts_with("GTEX-"), sum)) %>%
  as_tibble()

dat <- dplyr::left_join(dat, gene_lut, "Name")

tmp_dat <- dat %>%
  #mutate(across(starts_with("GTEX-"), ~ normTPM(x = .x, lengths = dat[["gene_length"]]))) %>%
  dplyr::filter(GENEBIOTYPE %in% c("protein_coding")) %>%
  dplyr::select(-c("Name", "gene_length", "GENEBIOTYPE"))

row_names <- tmp_dat[["SYMBOL"]]

tmp_dat <- tmp_dat %>%
  dplyr::select(starts_with("GTEX-")) %>%
  as.matrix

rownames(tmp_dat) <- row_names


meta_all <- meta_all[meta_all[["submitter_id"]] %in% colnames(tmp_dat), ]

meta_all <- meta_all[match(colnames(tmp_dat), meta_all[["submitter_id"]]),]


library(DESeq2)

meta_dat <- data.frame(tissue = meta_all$tissue_type,
                       donor = meta_all$subjects.id)

#run DESeq2

dds <- DESeqDataSetFromMatrix(countData = tmp_dat, colData = DataFrame(meta_dat), ~tissue + donor)
dds <- estimateSizeFactors(dds)
log2_dat <- DESeq2::normTransform(dds)
vst_dat <- vst(dds)



sce <- readRDS('~/GPR25_variants/gtex_log_Sept20_2024.rds')

DESeq2::DESeq(