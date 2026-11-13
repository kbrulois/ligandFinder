suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")
kg <- unique(x$knowns$accession)
f0 <- function(a) { m <- fc$feats[[a]]; mean(rowSums(abs(m[, ch4, drop=FALSE])) == 0) }
kn <- vapply(intersect(kg, prec$accession), f0, numeric(1))
cat(sprintf("KNOWN peptide precursors (n = %d): median no-data %.1f%%, mean %.1f%%, any gap %d, >25%% missing %d\n",
            length(kn), 100*median(kn), 100*mean(kn), sum(kn > 0), sum(kn > 0.25)))
set.seed(1); other <- setdiff(prec$accession, kg)
os <- vapply(sample(other, 1500), f0, numeric(1))
cat(sprintf("ALL OTHER precursors    (n = 1500 sample): median no-data %.1f%%, mean %.1f%%, >25%% missing %d (%.0f%%)\n",
            100*median(os), 100*mean(os), sum(os > 0.25), 100*mean(os > 0.25)))
cat("\nknown precursors with the most missing data:\n")
h <- sort(kn, decreasing = TRUE)[1:8]
print(data.frame(gene = prec$gene[match(names(h), prec$accession)],
                 accession = names(h), pct_no_data = round(100*h, 1),
                 length = nchar(prec$seq[match(names(h), prec$accession)])), row.names = FALSE)
