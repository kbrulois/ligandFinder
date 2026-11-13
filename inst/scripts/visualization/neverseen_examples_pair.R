## Never-seen peptide ends, AAswap vs ESM-C, with examples chosen by a fixed rule in
## BOTH directions (the original esmonly36_neverseen_examples.png picked genes by
## ESM-C's wins only, and its script was never saved; neverseen_examples_rebuilt.R
## reproduces it).
##
## The rule, decided 2026-11-10 after an independent check of the first version:
##  * an arm CALLS an end when its best CNN linear-head score (all-50 mean) within
##    +/-TOL residues reaches the arm's own threshold (arms/<arm>.sh THRESH_DEFAULT:
##    AAswap 0.5, ESM-C 0.65). TOL = 4, the tolerance of segments_roc.R/anchors_roc.R;
##    the exact anchor alone turned UCN N83 into an "AAswap win" that ESM-C makes one
##    residue later.
##  * ends at a precursor TERMINUS (N end on n_prot = first residue after the signal
##    peptide, C end on c_prot = last residue) are counted but not drawn: their windows
##    are about half padding, so a call there needs no cleavage-site recognition.
##  * precursors with no ESM-C embedding are left out (ESM-C scores them on zeros);
##    the subtitle says how many and how many of those AAswap calls.
##  * the N_EACH genes with the largest interior wins each way are drawn; a win's
##    size = winner's margin over its own threshold + loser's distance below its own.
##
##   Rscript neverseen_examples_pair.R
##   N_EACH=3 TOL=4 OUT=<png>
suppressMessages({library(data.table); library(ggplot2)})
A2 <- path.expand("~/AF2_analysis"); S <- path.expand(Sys.getenv("LF_RUN", file.path(A2, "lf_run")))
ARMS <- list(list(key = "_curated_AAswap", lab = "AAswap (43 ch)", short = "AAswap", thr = 0.5),
             list(key = "_esmonly36", lab = "ESM-C only (37 ch)", short = "ESM-C", thr = 0.65))
N_EACH <- as.integer(Sys.getenv("N_EACH", "3")); TOL <- as.integer(Sys.getenv("TOL", "4"))
## the peptide set; HAMP C84 is listed twice in the file, so dedupe on the end
un <- unique(fread(file.path(S, "cv5_unseen_ends.csv"))[grepl("peptide", set)], by = c("accession", "term", "anchor"))
un[, terminus := (term == "N" & anchor == n_prot) | (term == "C" & anchor == c_prot)]
emb <- unique(arrow::read_parquet(file.path(A2, "lf_plm/esm_c_pca36.parquet"), col_select = "accession")$accession)
un[, embedded := accession %in% emb]
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
g2a <- unique(as.data.table(fc$prec)[, .(accession, gene)]); rm(fc)
TR <- list()
for (tm in c("C", "N")) for (a in ARMS) {
  x <- readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, a$key))); setDT(x)
  TR[[paste(tm, a$key)]] <- x[accession %in% un$accession, .(accession, anchor, score = cnn__all50)]
  rm(x); invisible(gc())
}
sc <- function(tm, k) TR[[paste(tm, k)]]
## best score within +/-TOL of each end, per arm
best <- function(e, tm, k) {
  w <- e[, .(accession, a0 = anchor, lo = anchor - TOL, hi = anchor + TOL)]
  r <- sc(tm, k)[w, on = .(accession, anchor >= lo, anchor <= hi), .(accession, a0 = i.a0, score), allow.cartesian = TRUE]
  r[, .(s = max(score, na.rm = TRUE)), by = .(accession, a0)]
}
for (i in 1:2) {
  b <- rbindlist(lapply(c("C", "N"), function(tm) best(un[term == tm], tm, ARMS[[i]]$key)[, term := tm]))
  un[b, on = .(accession, term, anchor = a0), (paste0("s", i)) := i.s]
}
un[, call1 := s1 >= ARMS[[1]]$thr][, call2 := s2 >= ARMS[[2]]$thr]
un <- g2a[un, on = "accession"]
left <- un[embedded == FALSE]; ends <- un[embedded == TRUE]
stopifnot(!anyNA(ends$s1), !anyNA(ends$s2))
k <- function(x) sum(x, na.rm = TRUE)
tab <- ends[, .(n = .N, both = k(call1 & call2), only1 = k(call1 & !call2), only2 = k(!call1 & call2),
                neither = k(!call1 & !call2), only1_term = k(call1 & !call2 & terminus),
                only2_term = k(!call1 & call2 & terminus))]
print(tab); cat("left out (no ESM-C embedding):", nrow(left), " of which AAswap calls", k(left$call1), "\n")
w1 <- ends[call1 & !call2 & !terminus][, gap := (s1 - ARMS[[1]]$thr) + (ARMS[[2]]$thr - s2)][order(-gap)]
w2 <- ends[call2 & !call1 & !terminus][, gap := (s2 - ARMS[[2]]$thr) + (ARMS[[1]]$thr - s1)][order(-gap)]
pick <- function(w) w[gene %in% head(unique(w$gene), N_EACH)]
sel <- rbind(pick(w1)[, winner := ARMS[[1]]$short], pick(w2)[, winner := ARMS[[2]]$short])
print(sel[, .(gene, term, anchor, s1 = round(s1, 2), s2 = round(s2, 2), gap = round(gap, 2), winner)])
P <- unique(sel[, .(accession, gene, term, winner)])
P[, panel := sprintf("%s  %s end  (%s wins)", gene, term, winner)]
D <- rbindlist(lapply(seq_len(nrow(P)), function(j) {
  r <- P[j]; rbindlist(lapply(ARMS, function(a)
    sc(r$term, a$key)[accession == r$accession, .(anchor, score)][, `:=`(arm = a$lab, panel = r$panel)]))
}))
E <- un[P, on = .(accession, term), .(anchor, panel = i.panel,
          won = paste(accession, term, anchor) %in% sel[, paste(accession, term, anchor)])]
TH <- rbindlist(lapply(ARMS, function(a) data.table(arm = a$lab, thr = a$thr)))
lv <- c(P[winner == ARMS[[1]]$short, panel], P[winner == ARMS[[2]]$short, panel])
D[, panel := factor(panel, lv)]; E[, panel := factor(panel, lv)]
cols <- setNames(c("#3c63a0", "#c96a35"), sapply(ARMS, `[[`, "lab"))
p <- ggplot(D, aes(anchor, score, colour = arm)) +
  geom_hline(data = TH, aes(yintercept = thr, colour = arm), linewidth = 0.35, linetype = "22", show.legend = FALSE) +
  geom_vline(data = E, aes(xintercept = anchor, linewidth = won), colour = "grey25", show.legend = FALSE) +
  scale_linewidth_manual(values = c(`TRUE` = 0.7, `FALSE` = 0.3)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ panel, ncol = 2, scales = "free_x", dir = "v") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(title = "Never-seen peptide ends: AAswap vs ESM-C only",
       subtitle = sprintf(paste0(
         "CNN linear head, mean of all 50 members. An arm calls an end when its best score within %d residues reaches its own threshold (dashed: AAswap 0.5, ESM-C 0.65).\n",
         "Of %d never-seen peptide ends: both call %d, ESM-C only %d, AAswap only %d, neither %d.\n",
         "%d of ESM-C's one-arm calls and %d of AAswap's sit at a precursor terminus (window half padding); away from the termini it is %d vs %d.\n",
         "Left out: %d ends on precursors with no ESM-C embedding, which ESM-C scores on zeros; AAswap calls %d of them.\n",
         "Shown: the %d genes with the largest wins each way among ends away from the precursor termini. Heavy line = the winning end."),
         TOL, tab$n, tab$both, tab$only2, tab$only1, tab$neither, tab$only2_term, tab$only1_term,
         tab$only2 - tab$only2_term, tab$only1 - tab$only1_term,
         nrow(left), k(left$call1), N_EACH),
       x = "residue", y = "end score") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
out <- Sys.getenv("OUT", file.path(A2, "lf_headcmp", "curated_AAswap_vs_esmonly36_neverseen_examples.png"))
ggsave(out, p, width = 12, height = 2.1 * ceiling(nrow(P) / 2) + 1.9, dpi = 170, bg = "white")
cat("wrote", out, " panels:", nrow(P), "\n")
