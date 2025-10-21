

models <- list(good = c("hACKR3_hBRNP3x369x380", "hACKR3_hBRNP3x364x380", "hCCKAR_hBRNP3x331x380", "hAGRA1_hBRNP1x354x368"),
               incomplete = c("hCCR1_hBRNP2x347x397", "hCCKAR_hBRNP3x330x380", "hAGTR2_hBRNP1x342x368", "hACTHR_hBRNP1x356x368"),
               none = c("hCCR1_hBRNP1x356x368", "hCCKAR_hBRNP3x355x380", "hC5AR2_hBRNP3x364x380", "hAGTR1_hBRNP3x369x380"))

models <- out_dat %>%
  filter(model %in% (job_dat %>% filter(complete & wanted == "yes") %>% pull(model))) %>%
  pull(afpd_dir_name)


dest <- "/scratch/groups/ebutcher/deorphan/models/cxc17_gp15l"

source_dir <- "/scratch/groups/ebutcher/deorphan/models/add_bm"

dir.create(dest)


num_cores <- 16
future::plan(strategy = future::multicore(workers = num_cores))


furrr::future_map(models, ~fs::dir_copy(fs::path(source_dir, .),
                                        fs::path(dest, .)))


