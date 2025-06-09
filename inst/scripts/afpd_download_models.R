

pos <- c("hAGTR2_hANGTx25x31.tar", "hCXCR4_hSDF1x22x89.tar", "hGP152_hCCL4x24x92.tar",
                 "hMC3R_hCOLIx77x87.tar", "hNK1R_hTKN1x58x68.tar")

neg <- c("hV1AR_hOREXx70x97.tar", "hTA2R3_hGLUCx98x127.tar", "hV1BR_hCXCL3x35x107.tar",
                 "hCXCR5_hCCL16x24x120.tar", "hV1BR_hCCL27x25x112.tar")

all_files <- c(pos, neg)

lapply(all_files, \(x) {
  fs::file_copy(paste0("~/oak/deorphan-AI-ze/models/benchmarking/", x), "~/peptide_alg/testing_set")
})
