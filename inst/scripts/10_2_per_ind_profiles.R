


val_preds <- model %>% predict(nn_in$val$x)

true_classes <- apply(nn_in$val$y_per_index_cat, c(1,2), which.max) - 1
pred_classes <- apply(val_preds$per_index_cat, c(1,2), which.max) - 1

# masked accuracy - exclude none positions
non_none <- true_classes != 8

cat("overall accuracy (inc none):", mean(true_classes == pred_classes), "\n")
cat("masked accuracy (exc none):", mean(true_classes[non_none] == pred_classes[non_none]), "\n")

# per class breakdown
for(cls in names(classes)[names(classes) != "none"]) {
  idx <- true_classes == classes[cls]
  if(sum(idx) == 0) next
  acc <- mean(pred_classes[idx] == true_classes[idx])
  cat(sprintf("  %-30s acc=%.3f  n=%d\n", cls, acc, sum(idx)))
}


for(term in c("C", "N")) {

x_val <- array(
  unlist(nn_input[[term]]$all$data),
  dim = c(seq_len, n_channels, length(nn_input[[term]]$all$data))
)

x_val <- aperm(x_val, c(3, 1, 2))


pred <- models[[term]] %>%
          predict(x_val) %>%
          `[[`("per_index_cat")

pep_names_tp <- known_dat %>% filter(target == term) %>% pull(peps)

nn_input[[term]]$all$train <- if_else(nn_input[[term]]$all$peps %in% nn_input[[term]]$train$peps, "trn", "tst")



pep_names_tp2 <- grep("CXCL14", nn_input[[term]][["all"]]$peps, value = T)

pep_names_tp <- c("CXCL14_w76-111", "CXCL14_w88-111", pep_names_tp)

pep_names_tp <- pep_names_tp[pep_names_tp %in% nn_input[[term]]$all$peps]

if(term == "C") {ind_tr <- 4} else {ind_tr <- 1}

ind_tr <- c(ind_tr, 8)

df2 <- map(pep_names_tp,
           \(x) {
             ind <- which(nn_input[[term]]$all$peps == x)[1]
             pep_title <- paste0(x, nn_input[[term]]$all$target[ind], " ", nn_input[[term]]$all$train[ind])
             df <- pred[ind,,] %>%
               as_tibble
             colnames(df) <- names(class_cols)
             df <- df %>%
               mutate(index = row_number()) %>%
               pivot_longer(cols = -index) %>%
               mutate(title = pep_title)
           }) %>% bind_rows


df3 <- nn_input[[term]][["all"]] %>%
  filter(peps %in% pep_names_tp) %>%
  select(peps, target, meta_data, train) %>%
  unnest(meta_data) %>%
  mutate(title = paste0(peps, target)) %>%
  group_by(title) %>%
  mutate(index = row_number()) %>%
  mutate(title = paste(title, train)) %>%
  mutate(name = "seq") %>%
  mutate(value = -0.1)

df2 <- bind_rows(df3, df2)

p <-ggplot(data = df2 %>% filter(name != "seq"), aes(x = index, y = value, color = name)) +
  geom_point(pch = 21, stroke = 1) +
  geom_smooth(method = "loess", span = 0.6, se = FALSE) +
  scale_color_manual(values = class_cols) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ggplot2::geom_tile(data = df2 %>% filter(name == "seq"), mapping = aes(x = index, y = value, fill = known_idx), height = 0.2, linewidth= 0.01) +
  geom_point(data = df2 %>% filter(name == "seq"), pch = 15, size = 2.5, color = "grey85") +
  ggplot2::geom_text(data = df2 %>% filter(name == "seq"), mapping = aes(x = index, y = value, label = AA), size = 1.8, fontface = "bold", color = "black") +
  ggplot2::facet_wrap(vars(title)) +
  scale_fill_manual(values = class_cols)





ggsave(paste0("~/AF2_analysis/per_index_preds15_", term, "__",start, "-", end, ".svg"), plot = p, width = 18, height = 18)

}




#pep_names_tp <- stringr::str_subset(nn_input_comb$pep_name, "BRINP")


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
              colnames(df) <- names(class_cols)[-8]
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


ggsave(paste0("~/AF2_analysis/per_index_preds_topC_22dbt_pepval__",start, "-", end, ".svg"), plot = p, width = 18, height = 9)

}













