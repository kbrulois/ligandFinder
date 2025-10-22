

models <- list(good = c("hACKR3_hBRNP3x369x380", "hACKR3_hBRNP3x364x380", "hCCKAR_hBRNP3x331x380", "hAGRA1_hBRNP1x354x368"),
               incomplete = c("hCCR1_hBRNP2x347x397", "hCCKAR_hBRNP3x330x380", "hAGTR2_hBRNP1x342x368", "hACTHR_hBRNP1x356x368"),
               none = c("hCCR1_hBRNP1x356x368", "hCCKAR_hBRNP3x355x380", "hC5AR2_hBRNP3x364x380", "hAGTR1_hBRNP3x369x380"))

models <- out_dat %>%
  filter(model %in% (job_dat %>% filter(complete & wanted == "yes") %>% pull(model))) %>%
  pull(afpd_dir_name)

models <- out_dat %>%
  slice(1:20) %>% pull(afpd_dir_name)


dest <- "/scratch/groups/ebutcher/deorphan/models/test_top200"

dest2 <- "/scratch/groups/ebutcher/deorphan/models/test_top200_test"


source_dir <- "/scratch/groups/ebutcher/deorphan/models/top200NC_ffinal"

dir.create(dest)


num_cores <- 16
future::plan(strategy = future::multicore(workers = num_cores))


furrr::future_map(unname(models), \(x) {fs::dir_copy(fs::path(source_dir, x),
                                        fs::path(dest, x))})


future::plan(strategy = future::sequential())

furrr::future_walk(models, \(x) {
  f <- fs::path(dest, x)
  f2 <- fs::path(dest2, x)
  tar_cmd <- glue::glue("tar --format=posix -cf {f2}.tar -C {f} .")

  tryCatch({
    system(tar_cmd)
  }, error = function(e) {
    message("❌ Failed to tar ", f, ": ", conditionMessage(e))
  })
})


