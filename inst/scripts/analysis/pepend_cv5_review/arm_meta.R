## metadata for the arm comparison: which of the 41,870 `all` windows are known
## C ends, and which of those were held out
suppressMessages({library(dplyr); library(readr)})
S2 <- commandArgs(TRUE)[1]
x  <- readRDS(path.expand("~/AF2_analysis/lf_pepend_nn_input_near.rds"))
nn <- x$nn_input[["C"]]
all_t <- nn$all %>%
  left_join(x$knowns %>% filter(term == "C") %>% select(peps, split), by = "peps")
write_csv(all_t %>% transmute(row = row_number(), known = as.integer(known),
                              peps = as.character(peps),
                              split = ifelse(is.na(split), "", split),
                              motif = as.character(motif)),
          file.path(S2, "arm_all_meta.csv"))
message(sprintf("all=%d  known=%d  held out=%d", nrow(all_t), sum(all_t$known == 1),
                sum(all_t$known == 1 & all_t$split %in% "val", na.rm = TRUE)))
