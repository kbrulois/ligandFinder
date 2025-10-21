

models <- list(good = c("hACKR3_hBRNP3x369x380", "hACKR3_hBRNP3x364x380", "hCCKAR_hBRNP3x331x380", "hAGRA1_hBRNP1x354x368"),
               incomplete = c("hCCR1_hBRNP2x347x397", "hCCKAR_hBRNP3x330x380", "hAGTR2_hBRNP1x342x368", "hACTHR_hBRNP1x356x368"),
               none = c("hCCR1_hBRNP1x356x368", "hCCKAR_hBRNP3x355x380", "hC5AR2_hBRNP3x364x380", "hAGTR1_hBRNP3x369x380"))

dest <- "/scratch/groups/ebutcher/deorphan/models/brinp_test/"

dir.create(dest)

fs::dir_copy(paste0("/scratch/groups/ebutcher/deorphan/models/brinp_final/",
                    do.call(c, models)),
             paste0(dest, do.call(c, models)))
