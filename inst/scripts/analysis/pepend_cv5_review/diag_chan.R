suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")
ap <- fc$all_params3; prec <- fc$prec
scale_cols <- setdiff(ap, c("AA_hydro","AA_charge","AA_mw","AA_pI","padding"))
nd <- readRDS("~/AF2_analysis/nn_dat_cache.rds"); cr <- nd$chan_range; rm(nd); invisible(gc())
x <- readRDS(file.path(S2, "sa1_slim.rds"))
for (cn in scale_cols) x[[cn]] <- pmin(pmax((x[[cn]] - cr[1,cn])/(cr[2,cn]-cr[1,cn]), 0), 1)
rn <- split(seq_len(nrow(x)), x$accession)
set.seed(1); samp <- sample(intersect(names(fc$feats), names(rn)), 300)
tot <- setNames(numeric(length(scale_cols)), scale_cols)
dif <- tot; mx <- tot
for (a in samp) {
  o <- fc$feats[[a]]; n <- x[rn[[a]], ]
  if (nrow(o) != nrow(n)) next
  had <- rowSums(abs(o[, scale_cols, drop=FALSE])) > 0
  if (!any(had)) next
  for (cn in scale_cols) {
    ov <- o[had, cn]; nv <- n[[cn]][had]; nv[is.na(nv)] <- 0
    tot[cn] <- tot[cn] + length(ov); d <- abs(ov - nv)
    dif[cn] <- dif[cn] + sum(d > 1e-6); mx[cn] <- max(mx[cn], max(d))
  }
}
r <- data.frame(channel = scale_cols, n = tot, n_diff = dif,
                pct = round(100*dif/pmax(tot,1), 3), max_abs = round(mx, 4)) %>%
  arrange(desc(n_diff))
print(r, row.names = FALSE)
