


library(yardstick)


contact_scores <- data.table::fread("~/AF2_analysis/most_recent/bm_final.csv") %>% as_tibble

dat2 <- dat %>%
  {left_join(., contact_scores %>% select(starts_with("code"), starts_with("contact_score"), starts_with("new_metrics_")) %>% select("code"| where(is.numeric)), by = "code")}

metrics <- grep("^new_metrics|^contact_score|^favorability_|^area_|^tags_|^pLDDT_|^frequency_scaled|^dist_|^doubleSmooth_scaled|^paeL_|^paeR_|^score_|^mean_af_missense_",
                colnames(dat2), value = TRUE)

metrics <- c(metrics, c("depth", "radius", "iptm", "iptm+ptm"))

inverse_metrics <- c("depth", "radius", grep("RTCNN|FART|minPAE|avgPAE|RMSF|ligStrain|dist_", metrics, value = TRUE))


metric_rank <- dat2 %>%
                filter(run_name == "bm_sep28" & location == "relevant") %>%
                mutate(across(any_of(metrics), \(x) { if(is.numeric(x)) {return(x)} else {return(0)}})) %>%
                mutate(across(any_of(inverse_metrics), \(x) {x * -1})) %>%
                summarise(across(any_of(metrics), \(x) {
                  roc_auc_vec(.[["known_pair"]], x)
                }, .names = "{.col}_roc"),
                across(any_of(metrics), \(x) {
                  pr_auc_vec(.[["known_pair"]], x)
                }, .names = "{.col}_pr"))



metric_rank <- metric_rank %>%
  pivot_longer(cols = everything()) %>%
  mutate(metric = stringr::str_extract(name, "roc$|pr$")) %>%
  mutate(name = stringr::str_remove(name, "_roc$|_pr$")) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(desc(roc))



