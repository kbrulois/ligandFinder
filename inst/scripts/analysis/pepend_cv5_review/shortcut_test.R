suppressMessages({library(dplyr); library(purrr)})
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")
pr <- readRDS("~/AF2_analysis/lf_pepend_run_C_near_slack8.rds")$out$pred
all_t <- x$nn_input$C$all
stopifnot(nrow(all_t) == length(pr))
f0 <- map_dbl(all_t$data, function(d) {
  d <- as.data.frame(d); k <- d$padding == 0
  if (!any(k)) return(NA_real_); mean(rowSums(d[k, ch4] == 0) == 4)
})
kn <- all_t$known == 1
sp <- x$knowns %>% filter(term == "C") %>% distinct(peps, .keep_all = TRUE) %>% select(peps, split)
d <- tibble(peps = all_t$peps, known = kn, score = pr, f0 = f0,
            motif = all_t$motif) %>% left_join(sp, by = "peps")

cat(sprintf("OLD C arm, %d candidate windows\n", nrow(d)))
cat(sprintf("  Spearman(score, share of residues with no data) over the %d NEGATIVES: %.3f\n",
            sum(!kn), cor(d$score[!kn], d$f0[!kn], method = "spearman", use = "complete.obs")))
cat(sprintf("  negatives, median share missing: %.2f   knowns: %.2f\n",
            median(d$f0[!kn], na.rm=TRUE), median(d$f0[kn], na.rm=TRUE)))
q <- cut(d$f0[!kn], c(-.01,.001,.25,.5,.75,1.01), labels=c("0","<25%","25-50%","50-75%",">75%"))
cat("  mean old score of negatives by how much data they are missing:\n")
print(round(tapply(d$score[!kn], q, mean), 4))

## re-rank the held-out knowns against negatives that have FULL data only
held <- d %>% filter(known, split == "val")
rank_vs <- function(pool) vapply(held$score, function(v) sum(pool > v) + 1L, integer(1))
neg_all  <- d$score[!kn]
neg_full <- d$score[!kn & d$f0 == 0]
cat(sprintf("\n  negatives: %d all, %d with complete data (%.0f%%)\n",
            length(neg_all), length(neg_full), 100*length(neg_full)/length(neg_all)))
r_all <- rank_vs(neg_all); r_full <- rank_vs(neg_full)
## scale the complete-data ranks back to the full pool size, so the two are comparable
r_scaled <- r_full * length(neg_all) / length(neg_full)
cat(sprintf("  held-out median rank: %s against all negatives -> %s against complete-data negatives (rescaled to the same pool)\n",
            format(median(r_all)), format(round(median(r_scaled)))))
cat("  by motif (rescaled):\n")
tibble(motif = held$motif, all = r_all, full = round(r_scaled)) %>%
  group_by(motif) %>% summarise(n = n(), vs_all = median(all), vs_complete = median(full),
                                .groups="drop") %>% arrange(vs_all) %>% as.data.frame() %>% print(row.names = FALSE)
