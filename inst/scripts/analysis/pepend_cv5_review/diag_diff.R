suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
old <- readRDS(file.path(S2, "sa_slim.rds"))    # Jan table (extracellular only)
new <- readRDS(file.path(S2, "sa1_slim.rds"))   # Aug table (all topologies)
cat("shared columns:", paste(intersect(names(old), names(new)), collapse=", "), "\n\n")
## compare on raw table values, matched by (accession, index) where old index valid
chk <- readRDS(file.path(S2, "sa1_chk.rds"))
good <- chk$acc[chk$idx_eq_pos %in% TRUE]
set.seed(1); samp <- sample(good, 300)
ro <- split(seq_len(nrow(old)), old$accession); rn <- split(seq_len(nrow(new)), new$accession)
cols <- intersect(c("cons_rs_n","min_afm","mean_afm","relASA","SS","Phi","cons_rs","species_limit"),
                  intersect(names(old), names(new)))
agg <- list()
for (a in samp) {
  if (is.null(ro[[a]]) || is.null(rn[[a]])) next
  o <- old[ro[[a]], ]; o <- o[!is.na(o$index), ]
  if (!nrow(o)) next
  n <- new[rn[[a]], ]
  if (max(o$index) > nrow(n)) next
  n <- n[o$index, ]
  for (cn in cols) {
    a1 <- o[[cn]]; b1 <- n[[cn]]
    d <- if (is.numeric(a1)) {
      both <- !is.na(a1) & !is.na(b1); c(n_both = sum(both), n_diff = sum(a1[both] != b1[both]),
        max_abs = if (any(both)) max(abs(a1[both]-b1[both])) else 0)
    } else {
      both <- !is.na(a1) & !is.na(b1); c(n_both = sum(both), n_diff = sum(a1[both] != b1[both]), max_abs = NA)
    }
    agg[[cn]] <- rbind(agg[[cn]], d)
  }
}
cat("per channel, over residues where BOTH tables have a value:\n")
for (cn in names(agg)) {
  m <- colSums(agg[[cn]][, c("n_both","n_diff"), drop=FALSE])
  cat(sprintf("  %-14s both %7d  differ %7d (%.3f%%)  max|diff| %s\n", cn, m[1], m[2],
      100*m[2]/max(m[1],1), if (all(is.na(agg[[cn]][,"max_abs"]))) "-" else
        format(max(agg[[cn]][,"max_abs"], na.rm=TRUE), digits=3)))
}
