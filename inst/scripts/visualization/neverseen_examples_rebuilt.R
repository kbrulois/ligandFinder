## Reconstructed from the inline heredoc that wrote lf_headcmp/esmonly36_neverseen_examples.png
## (session 3c9d38d3, 2026-10-24T23:11Z), parameterised so the first arm can be AAswap.
##   ARM1=""                -> the original: curated (26 hand-built) vs ESM-C   (byte-for-byte the old logic)
##   ARM1=_curated_AAswap   -> AAswap vs ESM-C
##   GENES=NMS,RLN2,MSMB,INS (default)   OUT=<png>
suppressMessages({library(data.table); library(ggplot2)})
A2 <- path.expand("~/AF2_analysis"); S <- path.expand(Sys.getenv("LF_RUN", file.path(A2, "lf_run")))
ARM1 <- Sys.getenv("ARM1", "")
LAB1 <- if (ARM1 == "") "curated (26 hand-built)" else if (ARM1 == "_curated_AAswap") "AAswap" else sub("^_", "", ARM1)
GEN <- strsplit(Sys.getenv("GENES", "NMS,RLN2,MSMB,INS"), ",")[[1]]
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
prec <- as.data.table(fc$prec)[gene %in% GEN, .(accession, gene)]; rm(fc)
un <- fread(file.path(S, "cv5_unseen_ends.csv"))[grepl("^never seen", set)]
D <- list(); E <- list()
for (tm in c("C","N")) {
  cur <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, ARM1))); setDT(cur)
  esm <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten_esmonly36.rds", tm))); setDT(esm)
  for (g in GEN) {
    a <- prec[gene == g, accession]
    c1 <- cur[accession == a, .(anchor, score = cnn__all50)][, `:=`(arm = LAB1, gene = g, term = tm)]
    e1 <- esm[accession == a, .(anchor, score = cnn__all50)][, `:=`(arm = "ESM-C only (36 PCA)",  gene = g, term = tm)]
    D[[paste(g,tm)]] <- rbind(c1, e1)
    ends <- un[accession == a & term == tm, .(anchor)]
    if (nrow(ends)) E[[paste(g,tm)]] <- ends[, `:=`(gene = g, term = tm)]
  }
  rm(cur, esm); invisible(gc())
}
D <- rbindlist(D); E <- rbindlist(E)
D[, panel := paste0(gene, " — ", term, " terminus")]
E[, panel := paste0(gene, " — ", term, " terminus")]
D <- D[panel %in% E$panel]
cols <- setNames(c("#4C72B0", "#DD8452"), c(LAB1, "ESM-C only (36 PCA)"))
p <- ggplot(D, aes(anchor, score, colour = arm)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, colour = "grey75", linetype = "22") +
  geom_vline(data = E, aes(xintercept = anchor), linewidth = 0.45, colour = "grey25") +
  geom_line(linewidth = 0.65) +
  facet_wrap(~ panel, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(title = sprintf("Never-seen peptide ends: %s vs ESM-C only", LAB1),
       subtitle = "CNN linear head, mean of all 50 members. Vertical line = the true annotated end; dashed = T 0.5",
       x = "residue", y = "end score") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
## an ARM1 run is named "_rebuilt": curated_AAswap_vs_esmonly36_neverseen_examples.png
## is neverseen_examples_pair.R's figure, the fair one, and this must not replace it
out <- Sys.getenv("OUT", file.path(A2, "lf_headcmp", sprintf("%s_vs_esmonly36_neverseen_examples%s.png",
  if (ARM1 == "") "curated" else sub("^_", "", ARM1), if (ARM1 == "") "" else "_rebuilt")))
ggsave(out, p, width = 11, height = 7, dpi = 170, bg = "white")
cat("wrote", out, "\n"); cat("panels:", length(unique(D$panel)), " ends marked:", nrow(E), "\n")
