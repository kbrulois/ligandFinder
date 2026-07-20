

nn_input_comb %>%
  filter(grepl("ANO8", peps)) %>%
  select(peps, win_type, target, pred_raw) %>%
  arrange(desc(pred_raw)) %>%
  print(n = 30)



## --- format helpers ---------------------------------------------------------
## per_index_cat now outputs K_real independent sigmoids (columns =
## real_class_names), not an 8-class softmax; the input has n_channels incl. the
## end_type/position channels. build_x replaces the removed per-model nn_in objects.
build_x <- function(dataset) {
  n <- length(dataset[["data"]])
  aperm(array(unlist(dataset[["data"]]), dim = c(seq_len, n_channels, n)), c(3, 1, 2))
}
## Optional diagnostic: per-index (multi-label) accuracy over val, per model.
## Needs the trained `models` in the session; skipped when you're plotting from a
## saved nn_input_comb only.
if (exists("models")) {
real_class_names <- names(classes)[real_cols]                                   # 6 real classes (for the breakdown)
cat_class_names  <- names(classes)[c(real_cols, which(names(classes) == "none"))]  # K_cat = 6 real + none (output cols)
for (term in names(nn_input)) {
  vp    <- predict(models[[term]], build_x(nn_input[[term]]$val), verbose = 0)[["per_index_cat"]]  # (n, seq, K_cat)
  truth <- do.call(rbind, nn_input[[term]]$val$known_idx)                                          # (n, seq) labels
  if (anyNA(vp))
    message(sprintf("[%s] %d/%d per-index predictions are NaN (per-index head may have diverged)",
                    term, sum(is.na(vp)), length(vp)))
  # NaN-safe argmax that always returns a matrix (apply()+which.max() returns a
  # LIST when a cell is all-NaN -> the 'non-numeric matrix extent' error).
  # softmax over K_cat incl none -> argmax is the predicted label directly.
  d    <- dim(vp)
  flat <- matrix(vp, nrow = d[1] * d[2], ncol = d[3])   # (n*seq) x K_cat
  flat[is.na(flat)] <- -Inf
  ti   <- max.col(flat, ties.method = "first")
  top_i <- matrix(ti, nrow = d[1])                      # (n, seq)
  pred_lab <- matrix(cat_class_names[top_i], nrow = d[1])
  keep <- !truth %in% c("none", "padding")
  cat(sprintf("[%s] masked per-index accuracy (exc none/padding): %.3f  (n=%d)\n",
              term, mean(pred_lab[keep] == truth[keep]), sum(keep)))
  imp_classes <- c("DB", "pep_pocket", "pep_other")   # the classes that matter
  cls_acc <- setNames(rep(NA_real_, length(real_class_names)), real_class_names)
  for (cls in real_class_names) {
    idx <- truth == cls
    if (sum(idx) == 0) next
    cls_acc[cls] <- mean(pred_lab[idx] == truth[idx])
    cat(sprintf("    %-22s acc=%.3f  n=%d\n", cls, cls_acc[cls], sum(idx)))
  }
  imp <- cls_acc[intersect(imp_classes, names(cls_acc))]
  cat(sprintf("    >> MACRO acc over important classes (%s): %.3f\n",
              paste(names(imp), collapse = ", "), mean(imp, na.rm = TRUE)))
}
}  # end if (exists("models"))


## Plot per-index predictions for the KNOWNS straight from nn_input_comb, which
## already carries each window's per_index predictions + meta_data. No re-predict
## and no `models` needed -- runs from a saved nn_input_comb alone.
stopifnot("per_index" %in% names(nn_input_comb))

## Extra (non-known) peptide windows to plot alongside the knowns. Give window
## ids exactly as they appear in nn_input_comb$peps (e.g. "CXCL14_w76-111").
## Leave as character(0) to plot knowns only.
extra_peps <- c("ANO8_w12-47", "ANO8_w1-32", "ANO8_w27-62", "ANO8_w883-918")

## facet-strip text colours by set (rendered via ggtext markdown in the title).
## Install ggtext for the colours: install.packages("ggtext"). Without it the
## script still runs, just with plain (uncoloured) strip titles.
set_cols <- c(trn = "#1b7837", tst = "#2166ac", other = "#b35806")
use_md   <- requireNamespace("ggtext", quietly = TRUE)
if (!use_md) message("per_index: install.packages('ggtext') to colour the strip titles by trn/tst/other")

for(term in c("C", "N")) {

# knowns for this terminus (incl. loop_*) + any user-supplied extra windows
sub <- nn_input_comb %>%
  filter((known == 1 | peps %in% extra_peps),
         target %in% c(term, paste0("loop_", term)))

if (nrow(sub) == 0) { message("per_index: nothing to plot for ", term); next }

trn_peps <- if (exists("nn_input")) nn_input[[term]]$train$peps else character(0)
sub <- sub %>%
  mutate(set = dplyr::case_when(
    peps %in% trn_peps ~ "trn",     # known, in training set
    known == 1         ~ "tst",     # known, held out (val)
    TRUE               ~ "other"    # user-supplied extra window
  ))

# facet title carries a ggtext colour span so the strip text is coloured by set
make_title <- function(pep, tgt, st) {
  if (use_md) sprintf("<span style='color:%s'>%s%s %s</span>", set_cols[[st]], pep, tgt, st)
  else        paste0(pep, tgt, " ", st)
}

# per-index prediction curves (one facet per window)
df2 <- purrr::pmap(
  list(sub$peps, sub$target, sub$set, sub$per_index),
  function(pep, tgt, st, pidx) {
    pidx %>%
      pivot_longer(cols = -index) %>%
      mutate(title = make_title(pep, tgt, st))
  }) %>% bind_rows

# sequence / known_idx track underneath (known_idx pulled from the top-level
# column so extra/candidate windows -- which lack it in meta_data -- still work)
df3 <- purrr::pmap(
  list(sub$peps, sub$target, sub$set, sub$meta_data, sub$known_idx),
  function(pep, tgt, st, md, ki) {
    md %>%
      mutate(index     = row_number(),
             known_idx = ki,
             title     = make_title(pep, tgt, st),
             name      = "seq",
             value     = -0.1)
  }) %>% bind_rows

df2 <- bind_rows(df3, df2)

p <- ggplot(data = df2 %>% filter(name != "seq"), aes(x = index, y = value, color = name)) +
  geom_point(pch = 21, stroke = 1) +
  geom_smooth(method = "loess", span = 0.6, se = FALSE) +
  scale_color_manual(values = class_cols) +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.text   = if (use_md) ggtext::element_markdown() else element_text()) +   # colour-coded title if ggtext present
  ggplot2::geom_tile(data = df2 %>% filter(name == "seq"), mapping = aes(x = index, y = value, fill = known_idx), height = 0.2, linewidth= 0.01) +
  geom_point(data = df2 %>% filter(name == "seq"), pch = 15, size = 2.5, color = "grey85") +
  ggplot2::geom_text(data = df2 %>% filter(name == "seq"), mapping = aes(x = index, y = value, label = AA), size = 1.8, fontface = "bold", color = "black") +
  ggplot2::facet_wrap(vars(title)) +
  scale_fill_manual(values = class_cols)





ggsave(paste0("~/AF2_analysis/per_index_preds_knowns_", term, ".svg"), plot = p, width = 18, height = 18)

}




#pep_names_tp <- stringr::str_subset(nn_input_comb$pep_name, "BRINP")


if (FALSE) {   # candidate (unknown) top-rank profiles -- OFF by default; you asked
               # for knowns only. Flip to TRUE to plot ranked candidates, but it
               # then needs its own per-window predict (like the knowns loop above),
               # not the knowns-only `pred`, to be fast and correctly indexed.
max_rank <- 300
step <- 50

for (start in seq(1, max_rank, by = step)) {

  end <- min(start + step - 1, max_rank)



pep_names_tp <- nn_input_comb %>%
  filter(target %in% c(term, paste0("loop_", term))) %>%
  filter(location %in% c("4l", "3l", "3t")) %>%
  filter(known == 0) %>%
  dplyr::slice(start:end) %>%
  pull(peps)

pep_names_tp <- pep_names_tp[pep_names_tp %in% nn_input[[term]]$val$peps]

df2 <- map(pep_names_tp[-c(15:17)],
          \(x) {
            ind <- which(nn_input[[term]]$val$peps == x)
            pep_title <- paste0(x, nn_input[[term]]$val$target[ind])
              df <- pred[ind,,] %>%
                    as_tibble
              colnames(df) <- cat_class_names    # per_index_cat = K_cat softmax cols (6 real + none)
              df <- df %>%
                #select(-none) %>%
                mutate(index = row_number()) %>%
                pivot_longer(cols = -index) %>%
                mutate(title = pep_title)
          }) %>% bind_rows





p <- ggplot(df2, aes(x = index, y = value, color = name)) +
  geom_point() +
  geom_smooth(method = "loess", span = 0.6, se = FALSE) +
  facet_wrap(vars(title)) +
  scale_color_manual(values = class_cols) +
  theme_bw() +
  theme(legend.title = element_blank())


ggsave(paste0("~/AF2_analysis/per_index_preds_",start, "-", end, ".svg"), plot = p, width = 18, height = 9)

}
}  # end if (FALSE) -- candidate loop













