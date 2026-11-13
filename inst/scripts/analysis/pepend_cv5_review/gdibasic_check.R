suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus"); TOL <- 2L
source("R/pepend_windows.R")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")
seqs <- setNames(fc$prec$seq, fc$prec$accession)
sc <- as.data.table(readRDS(path.expand("~/AF2_analysis/lf_pepend_scan_C_near.rds"))$scan)[, .(accession, gene, anchor, score)]
sc[, motif := { s <- seqs[[.BY$accession]]
  if (is.null(s)) rep(NA_character_, .N) else lf_pepend_motif_vec(s, "C", anchor, 1L) }, by=accession]
cat("=== how many C anchors of each motif exist in the whole scan ===\n")
print(sc[, .(anchors = .N, pct = sprintf("%.1f%%", 100*.N/nrow(sc))), by=motif][order(-anchors)])

tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
g4 <- tst %>% mutate(m = mapply(function(a,e) lf_pepend_motif(list(seqs[[a]]),"C",list(e),list(1L)),
                                accession, end)) %>% filter(m == "G + dibasic")
ll <- fread(file.path(L,"ligand_list.csv")); up <- fread(file.path(L,"uniprot_peptides.csv")) %>%
  filter(type %in% c("peptide","propeptide"))
kn <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")$knowns %>% filter(term=="C")
allends <- bind_rows(kn %>% transmute(accession, anchor),
  tst %>% transmute(accession, anchor=end),
  ll %>% transmute(accession, anchor=as.integer(end)),
  up %>% transmute(accession, anchor=as.integer(end))) %>% filter(!is.na(anchor)) %>% distinct()
near <- CJ(d=-TOL:TOL, i=seq_len(nrow(allends)))[, .(k=paste(allends$accession[i], allends$anchor[i]+d))]
sc[, key := paste(accession, anchor)]
negs <- sc[!(key %in% near$k)]
cat(sprintf("\nshared negative pool: %s anchors, of which G+dibasic: %s\n",
            format(nrow(negs), big.mark=","), format(sum(negs$motif=="G + dibasic", na.rm=TRUE), big.mark=",")))
out <- merge(as.data.table(g4)[, .(accession, name, end)], sc[, .(accession, anchor, score)],
             by.x=c("accession","end"), by.y=c("accession","anchor"))
out[, rank_global  := sapply(score, function(v) sum(negs$score > v) + 1L)]
out[, rank_in_class := sapply(score, function(v) sum(negs$score > v & negs$motif == "G + dibasic") + 1L)]
out[, above_that_are_Gdib := sprintf("%.1f%%", 100*rank_in_class/rank_global)]
cat("\n=== the 4 never-seen G+dibasic peptides ===\n")
print(out[order(rank_global), .(name, end, score=round(score,4), rank_global, rank_in_class, above_that_are_Gdib)])
cat(sprintf("\nmedian rank: global %s, within the G+dibasic class %s\n",
            format(median(out$rank_global)), format(median(out$rank_in_class))))
cat("\n=== motif mix of the anchors outranking the MEDIAN of these 4 ===\n")
thr <- median(out$score); up2 <- negs[score > thr]
print(up2[, .(n=.N, pct=sprintf("%.1f%%", 100*.N/nrow(up2))), by=motif][order(-n)])
