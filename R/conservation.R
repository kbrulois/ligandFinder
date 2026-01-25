


do_cons_hu_ref <- function(cons,
                           sim_mat,
                           sim_cutoff = 0.85) {

  bind_cols(cons,
            map(names(sim_mat), \(x) {

              sim_mat[[x]] %>%
                mutate(index = row_number()) %>%
                pivot_longer(cols = -index, names_to = "aminode", values_to = "similarity") %>%
                {left_join(., species_dat, by = "aminode")} %>%
                filter(aminode != "Homo_sapiens") %>%
                mutate(wts = case_when(is.na(similarity) ~ NA,
                                       similarity < sim_cutoff ~ myo_dissim,
                                       similarity >= sim_cutoff ~ myo_sim)) %>%
                group_by(index) %>%
                summarise(!!paste0(x, "_wt_all") := weighted.mean(x = similarity, w = wts, na.rm = TRUE),
                          !!paste0(x, "_uw_all") := mean(similarity, na.rm = TRUE),
                          !!paste0(x, "_wt_mam") := weighted.mean(x = similarity[class == "Mammalia"], w = wts[class == "Mammalia"], na.rm = TRUE),
                          !!paste0(x, "_uw_mam") := mean(similarity[class == "Mammalia"], na.rm = TRUE)) %>%
                ungroup %>%
                select(-index)

            })
  )

}



do_cons <- function(cons, alignment_AA) {

  bind_cols(cons,
            map(names(sim_mats), \(x) {

              sm <- sim_mats[[x]]

              all_in <- alignment_AA %>%
                mutate(index = row_number()) %>%
                pivot_longer(cols = -index, names_to = "aminode")

              mam_in <- all_in %>%
                {left_join(., species_dat, by = "aminode")} %>%
                filter(class == "Mammalia") %>%
                select(index, aminode, value)

              inputs <- list(all = all_in,
                             mam = mam_in)

              bind_cols(
                map(names(inputs), \(y) {

                  tmp <- inputs[[y]] %>%
                    group_by(index) %>%
                    reframe(tibble_table(value)) %>%
                    pivot_wider(names_from = "index", values_from = freq, values_fill = 0) %>%
                    dplyr::slice(match(rownames(sm), AA)) %>%
                    {row_names <<- .[["AA"]]; .} %>%
                    select(-AA) %>%
                    as.matrix

                  rownames(tmp) <- row_names

                  sm <- sm[rownames(sm) %in% rownames(tmp), colnames(sm) %in% rownames(tmp)]

                  tmp <- tmp[!rownames(tmp) %in% c("-", "X"), ]

                  sm <- sm[!rownames(sm) %in% c("-", "X"), !colnames(sm) %in% c("-", "X")]

                  tmp <- apply(tmp, 2, \(y) y/sum(y))

                  tibble(!!paste0(x, "_nr_", y) := drop(apply(tmp, 2, \(z) crossprod(z, sm %*% z))))
                }))

            }))

}



