

library(ligandFinder)
library(tidyverse)

## ---- fast-iteration cache --------------------------------------------------
## known_dat + c_dat are the expensive part of this script: get_pep_data runs
## over ~40K windows, filtering the 2.8 GB secretome_aa per protein (~1 hr). But
## everything BELOW the cache (end_type channel + make_training_sets) is cheap.
## So cache known_dat/c_dat once; re-runs that only tweak the split / channel /
## val_frac reload in seconds and don't need secretome_aa in the session at all.
## Force a full rebuild by deleting the file or: options(lf.rebuild_nn_dat = TRUE)
nn_dat_cache <- "~/AF2_analysis/nn_dat_cache.rds"

if (file.exists(nn_dat_cache) && !isTRUE(getOption("lf.rebuild_nn_dat"))) {

  message("9.2: loading cached known_dat/c_dat from ", nn_dat_cache,
          "  (options(lf.rebuild_nn_dat = TRUE) to rebuild)")
  .cc         <- readRDS(nn_dat_cache)
  known_dat   <- .cc$known_dat
  c_dat       <- .cc$c_dat
  all_params3 <- .cc$all_params3

} else {

con_dat <- readRDS("~/AF2_analysis/knowns.rds")


con_dat <- con_dat %>%
  filter(location == "relevant") %>%
  filter(rank == 1, .by = afpd_dir_name) %>%
  mutate(p2_gene = setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])[p2_name]) %>%
  mutate(pep_id = paste0(p2_gene, "_", p2_range))

con_dat %>%
  group_by(pep_id, p1_name) %>%
  summarise(uniq_ins = paste0(unique(lig1_end), collapse = ";")) %>%
  print(n = 200)


con_dat <- con_dat %>%
  mutate(lig1_end_ind = stringr::str_extract(lig1_end, "\\d+") %>% as.integer, .before = everything()) %>%
  mutate(lig1_end_t = stringr::str_remove(lig1_end, "\\d+"), .before = everything()) %>%
  group_by(pep_id) %>%
  arrange(lig1_end_ind, desc(iptm), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

con_dat <- con_dat %>%
              mutate(tibble(p2_start = stringr::str_extract(p2_range, "^\\d+"),
                            p2_end = stringr::str_extract(p2_range, "\\d+$"))) %>%
              group_by(lig1_end_t, lig1_end_ind, p2_name, p2_end) %>%
              arrange(desc(iptm)) %>%
              slice(1) %>%
              ungroup() %>%
              group_by(lig1_end_t, lig1_end_ind, p2_name, p2_start) %>%
              arrange(desc(iptm)) %>%
              slice(1) %>%
              ungroup()



table(con_dat$lig1_end_t, con_dat$lig1_end_ind)


con_dat <- con_dat %>%
            filter(lig1_end_ind < 16) %>%
            mutate(target = if_else(lig1_end_ind < 3, lig1_end_t, paste0("loop_", lig1_end_t)))


secretome <- secretome %>%
  mutate(sp_ind = map_int(features, \(x) {
    n_trunc <- x %>%
      filter(type == "signal peptide") %>%
      pull(end)

    if(length(n_trunc) == 0) {
      n_trunc <- 1
    }
    if(is.na(n_trunc)) {
      n_trunc <- 1
    }
    if(length(n_trunc) > 1) {
      n_trunc <- max(n_trunc, na.rm = TRUE)
    }
    n_trunc

  }))

con_dat <- left_join(con_dat %>% mutate(accession = p2_id), secretome %>% select(accession, sequence_uni, sp_ind), by = "accession")

## Drop knowns whose precursor protein is absent from `secretome` (no
## sequence_uni / sp_ind): they cannot be windowed and would otherwise error in
## get_adj_db_sites at `if (sp_ind >= min_pep_ind)`.
{
  .missing <- con_dat %>% filter(is.na(sequence_uni) | is.na(sp_ind))
  if (nrow(.missing) > 0) {
    message(sprintf("Dropping %d known(s) not found in `secretome`: %s",
                    nrow(.missing),
                    paste(unique(.missing$p2_gene), collapse = ", ")))
  }
}
con_dat <- con_dat %>% filter(!is.na(sequence_uni), !is.na(sp_ind))





ind <- 1

get_last_loc <- function(x, pattern) {
  str_locate_all(x, pattern) %>%
    purrr::map(~ if (nrow(.x) > 0) .x[nrow(.x), ] else c(NA, NA)) %>%
    do.call(rbind, .)
}

ind <- 94
sequence_uni = con_dat$sequence_uni[ind]
start = con_dat$p2_start[ind] %>% as.numeric
end = con_dat$p2_end[ind] %>% as.numeric
sp_ind = con_dat$sp_ind[ind]
target = con_dat$target[ind]
dbn_ws = 20
dbc_ws = 20
pep_id = con_dat$pep_id[ind]

win_size <- list(N = list(start = -5,
                          end = 30),
                 C = list(start = -30,
                          end = 5)
)


get_adj_db_sites <- function(sequence_uni = con_dat$sequence_uni[ind],
                             start = con_dat$p2_start[ind] %>% as.numeric,
                             end = con_dat$p2_end[ind] %>% as.numeric,
                             sp_ind = con_dat$sp_ind[ind],
                             target = con_dat$target[ind],
                             dbn_ws = 20,
                             dbc_ws = 20,
                             term_tol = 2,
                             window_size = win_size,
                             pep_id = pep_id) {

start <- as.numeric(start)
end <- as.numeric(end)

pep_inds <- stringr::str_extract(pep_id, "\\d+x\\d+") %>% stringr::str_split(., "x", simplify = TRUE) %>% `c` %>% as.integer
min_pep_ind <- min(pep_inds)
if(sp_ind >= min_pep_ind) {sp_ind <- min_pep_ind - 1}
if(sp_ind == 1) {sp_ind <- 0}


n_dbw <- max(sp_ind, start - dbn_ws, na.rm = TRUE)
c_dbw <- min(nchar(sequence_uni), end + dbc_ws, na.rm = TRUE)


n_prot <- max(sp_ind + 1, 1, na.rm = TRUE)
c_prot <- nchar(sequence_uni)


if(target %in% c("C", "loop_C")) {
  db_seq <- stringr::str_sub(sequence_uni, start = end, end = c_dbw)
  db_ind <- stringr::str_locate(db_seq, "(?=(KK|KR|RK|RR))")[1, "start"]
  db_ind = end + db_ind
  db_spacer <- db_ind - end
  ws <- window_size[["C"]]
} else {
  db_seq <- stringr::str_sub(sequence_uni, start = n_dbw, end = start)
  db_ind <- get_last_loc(db_seq, "(?=(KK|KR|RK|RR))")[1, 1] + 1 ##start and end switched for this regex
  db_ind = n_dbw +  db_ind - 2
  db_spacer <- start - db_ind
  ws <- window_size[["N"]]
}

tmp <- tibble(db_seq = db_seq,
              db_ind = db_ind,
              db_spacer = db_spacer)

## Is the relevant peptide boundary at the precursor protein terminus?
## N boundary -> mature N-terminus (first residue after the signal peptide, n_prot);
## C boundary -> the protein's last residue (c_prot).
is_terminus <- if (target %in% c("C", "loop_C")) {
  isTRUE(end >= c_prot - term_tol)
} else {
  isTRUE(start <= n_prot + term_tol)
}

tmp <- tmp %>%
          ## end_type is the new split axis: "end" (boundary is the protein
          ## terminus) vs "cleavage" (db- or chym-anchored internal cut).
          mutate(end_type = if_else(is_terminus, "end", "cleavage")) %>%
          mutate(win_type = case_when(is_terminus                   ~ "pep_end",
                                      !is.na(db_ind) & db_spacer < 6 ~ "db",
                                      TRUE                           ~ "chym")) %>%
          mutate(window_origin = case_when(win_type == "db"             ~ db_ind,
                                           target %in% c("C", "loop_C") ~ end + 2,
                                           TRUE                         ~ start - 2))


tmp %>%
  rowwise() %>%
  mutate(wN = window_origin + ws[[1]],
         wC = window_origin + ws[[2]]) %>%
  mutate(across(all_of(c("wN")), ~n_prot - ., .names = "{.col}_pad")) %>%
  mutate(across(all_of(c("wC")), ~. - c_prot, .names = "{.col}_pad")) %>%
  mutate(across(all_of(c("wN")), ~max(., n_prot, na.rm = TRUE))) %>%
  mutate(across(all_of(c("wC")), ~min(., c_prot, na.rm = TRUE)))

}




con_dat <- con_dat %>%
  rowwise() %>%
  mutate(windows = pmap(list(sequence_uni = sequence_uni,
                             start = p2_start,
                             end = p2_end,
                             target = target,
                             sp_ind = sp_ind,
                             pep_id = pep_id),
                        .f = get_adj_db_sites))


res <- c("F", "W", "Y")

di_res <- as.vector(outer(res, res, paste0))


res1 <- c("L", "I", "V", "M")
res2 <- c("F", "Y", "W")

di_res2 <- as.vector(outer(res1, res2, paste0))
di_res3 <- as.vector(outer(res2, res1, paste0))
chymo_pat <- paste0("(", paste(c(di_res, di_res2, di_res3), collapse = "|"), ")")


clean_contacts <- function(contact_dat) {

  if(is.null(contact_dat)) return(contact_dat)

  contact_dat %>%
    filter(area > 1 & dist < 6) %>%
    mutate(p_seg = stringr::str_remove(protein_segment, "\\d+$")) %>%
    group_by(lig1_residue) %>%
    summarise(protein_segment = paste0(unique(p_seg) %>% .[gtools::mixedorder(.)], collapse = "_"),
              depth = first(IC_dist),
              in_pocket = first(in_pocket))

}

con_dat <- con_dat %>%
              ungroup() %>%
              mutate(clean_cons = map(contacts, clean_contacts))


if(FALSE) {
dput(c(all_params[c(1,15,30,31,32)],
       all_params[stringr::str_detect(all_params, "^SS_")],
       all_params[stringr::str_detect(all_params, "^AA_")]))
}

all_params3 <- c("cons_rs", "cons_rs_n", "min_afm", "mean_afm", "relASA", "SS_P",
                 "SS_S", "SS_E", "SS_-", "SS_T", "SS_G", "SS_B", "SS_H", "SS_I",
                 "AA_Q", "AA_P", "AA_V", "AA_L", "AA_T", "AA_S", "AA_A", "AA_G",
                 "AA_R", "AA_F", "AA_C", "AA_I", "AA_N", "AA_Y", "AA_W", "AA_K",
                 "AA_D", "AA_E", "AA_M", "AA_H", "AA_U", "padding")


get_pep_data <- function(window, p_id, gene, nn_params = all_params3, meta_dat = NULL) {

  dat <- secretome_aa %>%
    filter(accession == p_id) %>%
    mutate(padding = 0) %>%
    select(all_of(c(nn_params, "AA", "index", "topo2", "known"))) %>%
    mutate(across(everything(), ~ replace_na(.x, 0))) %>%
    mutate(index = row_number())

  pad <- dat %>% mutate(across(everything(), ~ .x[NA_integer_]))

  slice_in_data <- function(N, C, N_pad, C_pad) {

    if(is.na(N_pad)) N_pad <- 0
    if(is.na(C_pad)) C_pad <- 0

    tmp <- dat[dat$index %in% N:C, ]

    if(N_pad > 0) {
      tmp <- bind_rows(pad %>%
                         dplyr::slice(1:N_pad) %>%
                         mutate(padding = 1),
                       tmp
      )
    }

    if(C_pad > 0) {
      tmp <- bind_rows(tmp,
                       pad %>%
                         dplyr::slice(1:C_pad) %>%
                         mutate(padding = 1)
      )
    }

    if(!is.null(meta_dat)) {
    tmp <- tmp %>%
              mutate(lig1_residue = paste0(AA, index)) %>%
              {left_join(., meta_dat, by = "lig1_residue")}
    }

    return(tmp)

  }

  if(is.null(window)){
    tmp <- tibble(peps = character(), dat = list(), win_type = character(), end_type = character(), target = character())
  } else {
  if(nrow(window) == 0) {
    tmp <- tibble(peps = character(), dat = list(), win_type = character(), end_type = character(), target = character())
  } else {
    tmp <- window %>%
      ungroup() %>%
      mutate(peps = paste0(gene, "_w", wN, "-", wC)) %>%
      mutate(dat = pmap(list(wN, wC, wN_pad, wC_pad), slice_in_data)) %>%
      select(any_of(c("peps", "dat", "win_type", "end_type", "target")))


  }
  }

  return(tmp)

}



make_candidate_windows <- function(sequence_uni, sp_ind, window_size = win_size) {

  n_prot <- max(sp_ind, 1, na.rm = TRUE)
  c_prot <- nchar(sequence_uni)

  ws <- window_size

  tmp <- stringr::str_locate_all(sequence_uni, "KK|KR|RK|RR")[[1]] %>%
    as_tibble

  if(nrow(tmp) > 0) {
    tmp <- tmp %>%
      group_by(start, end) %>%
      reframe(target = c("N", "C"), win_type = "db")
  }

  tmp2 <- stringr::str_locate_all(sequence_uni, paste0(chymo_pat, "(?!.{0,12}(KK|KR|RK|RR))"))[[1]] %>%
    as_tibble %>%
    mutate(win_type = "chym") %>%
    mutate(target = "C")

  tmp3 <- stringr::str_locate_all(sequence_uni, paste0("(?<!(KK|KR|RK|RR).{0,12})", chymo_pat))[[1]] %>%
    as_tibble %>%
    mutate(win_type = "chym") %>%
    mutate(target = "N")

  ## Protein-terminus windows, so the "end" models have terminus-anchored
  ## negatives: mature N-terminus (first residue after the signal peptide)
  ## and the protein's C-terminus (last residue).
  term_tmp <- tibble(start    = c(sp_ind + 1L, c_prot),
                     end      = c(sp_ind + 1L, c_prot),
                     win_type = "pep_end",
                     target   = c("N", "C"))

  tmp <- bind_rows(tmp, tmp2, tmp3, term_tmp)

  tmp <- tmp %>%
    filter(!start < n_prot) %>%
    filter(!end < n_prot) %>%
    mutate(end_type = if_else(win_type == "pep_end", "end", "cleavage"))

  if(nrow(tmp) > 0) {
    return(
      tmp %>%
        rowwise() %>%
        mutate(if_else(target == "C",
                       list(tibble(wN = end + ws[["C"]][[1]],
                                   wC = end + ws[["C"]][[2]])),
                       list(tibble(wN = end + ws[["N"]][[1]],
                                   wC = end + ws[["N"]][[2]]))) %>% bind_rows) %>%
        mutate(across(all_of(c("wN")), ~n_prot - ., .names = "{.col}_pad")) %>%
        mutate(across(all_of(c("wC")), ~. - c_prot, .names = "{.col}_pad")) %>%
        mutate(across(all_of(c("wN")), ~max(., n_prot, na.rm = TRUE))) %>%
        mutate(across(all_of(c("wC")), ~min(., c_prot, na.rm = TRUE)))
    )

  } else {return(NULL)}
}

secretome <- secretome %>%
  mutate(windows = map2(sequence_uni, sp_ind, make_candidate_windows))


ctrl_dbw <- secretome

yo()



ctrl_dbw <- ctrl_dbw %>%
    rowwise() %>%
    reframe(get_pep_data(window = windows, p_id = accession, gene = gene),
            gene = gene)



known_dat <- con_dat %>%
                rowwise() %>%
                reframe(get_pep_data(window = windows, p_id = accession, gene = p2_gene, meta_dat = clean_cons),
                        gene = p2_gene,
                        target = target,
                        pep_id = pep_id)


known_dat <- known_dat %>%
              filter(map_lgl(dat, ~nrow(.) == 36))

known_dat <- known_dat %>%
              mutate(known = 1)

sep_meta_dat <- function(dat, nn_params = all_params3) {
  list(
    meta_data = select(dat, -any_of(nn_params)),
    data = select(dat, any_of(nn_params))
  )
}



known_dat <- known_dat %>%
          mutate(split = map(dat, sep_meta_dat)) %>%
          select(-dat) %>%
          unnest_wider(split)

make_known_idx <- function(meta_data, target, pep_id, end_type = "cleavage") {

  pep_backbone <- rep("gap", 36)

  pep_inds <- stringr::str_extract(pep_id, "\\d+x\\d+") %>% stringr::str_split(., "x", simplify = TRUE) %>% `c` %>% as.integer
  w_inds <- which(meta_data[["index"]] %in% pep_inds)

  ## For "end" peptides the relevant boundary is the protein terminus, so the
  ## window region on that side lies outside the protein: label it "padding"
  ## rather than the db-/cleavage-context motifs used for internal cuts.
  is_end <- end_type == "end"

  if(target %in% c("C", "loop_C")) {
    if(length(w_inds) == 1) {w_inds <- c(1, w_inds)}
    pep_backbone[30:36] <- if (is_end) "padding" else c(rep("DB", 2), rep("CT_cleavage_context", 5))
    pep_backbone[w_inds[1]:w_inds[2]] <- "pep_other"
    if(w_inds[1] != 1) {
    pep_backbone[1:(w_inds[1] - 1)] <- "NT_cleavage_context"
    }
  } else {
    if(length(w_inds) == 1) {w_inds <- c(w_inds, 36)}
    pep_backbone[1:7] <- if (is_end) "padding" else c(rep("NT_cleavage_context", 5), rep("DB", 2))
    pep_backbone[w_inds[1]:w_inds[2]] <- "pep_other"
    if(w_inds[2] != 36) {
      pep_backbone[(w_inds[2] + 1):36] <- "CT_cleavage_context"
    }
  }

  pad_inds <- which(is.na(meta_data$AA))
  pep_backbone[pad_inds] <- "padding"

  if(!"in_pocket" %in% colnames(meta_data)) {
    meta_data <- meta_data %>% mutate(in_pocket = FALSE)
  }
  meta_data %>%
    mutate(known_idx = pep_backbone) %>%
    mutate(pocket = if_else(in_pocket, "pep_pocket", "pep_other")) %>%
    mutate(known_idx = if_else(is.na(in_pocket), known_idx, pocket))

}





known_dat <- known_dat %>%
                mutate(meta_data = pmap(list(meta_data, target, pep_id, end_type), make_known_idx)) %>%
                mutate(known_idx = map(meta_data, \(x) x[["known_idx"]])) %>%
                mutate(known_idx_detailed = known_idx) %>%
                mutate(known_idx2 = map(known_idx_detailed, ~case_when(. == "pep_pocket" ~ "pep_pocket",
                                             TRUE ~ "none")))


known_dat <- known_dat %>%
  mutate(data = map(data, \(x) { x %>% mutate(across(everything(), ~replace_na(., 0)))}))



filter_window_dat <- function(data) {

  map_lgl(data[["topo2"]], ~. %in% c("i", "t")) %>% any %>% `!` & map_lgl(data[["known"]], ~. == 1) %>% any %>% `!`

}

ctrl_dbw$dat[ctrl_dbw$gene == "CXCL14"] <- ctrl_dbw$dat[ctrl_dbw$gene == "CXCL14"] %>% map(., \(x) { x$known <- 0; return(x)})


c_dat <- ctrl_dbw

c_dat <- c_dat %>%
  filter(!is.na(peps)) %>%
  #filter(map_lgl(dat, filter_window_dat)) %>%
  filter(map_lgl(dat, ~nrow(.) == 36))


sep_meta_dat <- function(data, nn_params = all_params3) {
  list(
    meta_data = select(data, -any_of(nn_params)),
    data = select(data, any_of(nn_params))
  )
}

c_dat <- c_dat %>%
  mutate(split = map(dat, sep_meta_dat)) %>%
  select(-dat) %>%
  unnest_wider(split) %>%
  mutate(known = 0)

c_dat <- c_dat %>%
          mutate(known_idx = list(c(rep("none", 36))))

c_dat <- c_dat %>%
          mutate(data = map(data, \(x) { x %>% mutate(across(everything(), ~replace_na(., 0)))}))

saveRDS(list(known_dat = known_dat, c_dat = c_dat, all_params3 = all_params3),
        nn_dat_cache)
message("9.2: cached known_dat/c_dat to ", nn_dat_cache)

}  # end cache-rebuild branch (everything below is cheap)

yo()


gpcr_ligs <- data.table::fread(system.file("extdata/GPCRdb_known_pairings_human_plus2more_unique.csv",
                                             package = "ligandFinder")) %>% as_tibble %>%
  mutate(gene = stringr::str_remove(lig, "^h") %>% stringr::str_remove(., "x\\d+x\\d+"))


uni_pep <- data.table::fread("~/R_projects/ligandFinder/inst/extdata/Candidate_peps_KB_from_known_pep_unknownGPCR_run#1.txt") %>% as_tibble

uni_pep <- uni_pep %>%
  filter(!ALso_run %in% c("no", "not_now", "??", "did we run?"))

uni_pep <- uni_pep %>%
  mutate(uniprot_name = setNames(id_map$`Entry Name`, id_map$Entry)[Accession]) %>%
  mutate(model = paste0(uniprot_name, ",", Start, "-", End), .before = everything())

banned <- c(gpcr_ligs$gene, uni_pep$Gene) %>% unique





## Grouped, diversity-spread train/val split of the known windows.
## Groups near-identical windows (families, overlapping/same-gene peptides) by
## aligned-sequence similarity so no cluster straddles train/val (prevents the
## leakage that inflates val AUC), then sends every ~1/val_frac-th cluster -- in
## dendrogram-leaf order -- to val so each set spans the diversity. `sim_h` is the
## cut height in ape::dist.aa's scaled units (fraction of differing positions):
## 0.30 groups windows that are >~70% identical. Falls back to a random split if
## `ape` is unavailable or clustering fails.
split_knowns <- function(k_sub, val_frac = 0.34, sim_h = 0.30) {

  n <- nrow(k_sub)
  if (n <= 1) return(list(train = seq_len(n), val = integer(0)))

  random_split <- function() {
    n_val <- max(1L, min(n - 1L, round(n * val_frac)))
    val   <- sample(seq_len(n), n_val)
    list(train = setdiff(seq_len(n), val), val = val)
  }

  if (!requireNamespace("ape", quietly = TRUE)) {
    message("split_knowns: 'ape' not installed -> random split")
    return(random_split())
  }

  tryCatch({
    ## aligned AA matrix (n x window_len); padding / missing -> "X"
    aa <- do.call(rbind, lapply(k_sub$meta_data, function(m) {
      a <- as.character(m[["AA"]]); a[is.na(a)] <- "X"; a
    }))

    d   <- ape::dist.aa(aa, scaled = TRUE)          # fraction of differing sites [0,1]
    hc  <- hclust(as.dist(d), method = "average")
    cl  <- cutree(hc, h = sim_h)                    # windows within sim_h divergence -> same cluster

    ord  <- unique(cl[hc$order])                    # cluster ids in dendrogram-leaf order
    step <- max(2L, round(1 / val_frac))
    val_clusters <- ord[seq(2L, length(ord), by = step)]  # spread val across the ordering
    val   <- which(cl %in% val_clusters)
    train <- setdiff(seq_len(n), val)

    if (length(train) == 0L || length(val) == 0L) random_split()
    else list(train = train, val = val)
  }, error = function(e) {
    message("split_knowns: clustering failed (", conditionMessage(e), ") -> random split")
    random_split()
  })
}


## Two models: N and C, each POOLING end + cleavage windows. The 4-way split
## (N/C x end/cleavage) starved C_end (too few known C-terminus peptides), so we
## pool per terminus for statistical strength and instead hand the model the
## regime as an input channel (end_type_ch, injected below) so it can still
## specialize. end_type stays a per-window column for the outputs/plots.
## Sample `n_ctrl` control windows whose pep_end-vs-cleavage (end_type) mix matches
## the positive set `pos`, so positives and negatives share the same end_type
## proportions and the model can't separate them on end_type composition alone.
## Falls back to a plain random sample when `pos` is empty.
balanced_ctrl_sample <- function(pos, ctrl, n_ctrl) {
  n_take <- min(n_ctrl, nrow(ctrl))
  if (nrow(pos) == 0 || n_take == 0) return(dplyr::slice_sample(ctrl, n = n_take))

  ets  <- c("end", "cleavage")
  prop <- prop.table(table(factor(pos$end_type, levels = ets)))
  n_by <- round(n_take * as.numeric(prop)); names(n_by) <- ets
  d    <- n_take - sum(n_by)                       # fix rounding so counts sum to n_take
  if (d != 0) n_by[which.max(n_by)] <- n_by[which.max(n_by)] + d

  purrr::imap(n_by, function(n_i, et) {
    pool <- dplyr::filter(ctrl, end_type == et)
    if (n_i <= 0 || nrow(pool) == 0) return(pool[0, ])
    dplyr::slice_sample(pool, n = n_i, replace = nrow(pool) < n_i)  # replace only if pool too small
  }) %>% dplyr::bind_rows()
}

make_training_sets <- function(k_dat, ctr_dat,
                               n_ctrl   = 400,
                               val_frac = 0.34,
                               sim_h    = 0.30) {

  targets <- list(N = c("N", "loop_N"), C = c("C", "loop_C"))

  Map(function(targ) {

    k_sub <- k_dat %>%
      filter(target %in% targ)

    c_sub <- ctr_dat %>%
      filter(target %in% targ, !gene %in% banned)

    ## grouped, diversity-spread split (see split_knowns) instead of a random draw
    sp         <- split_knowns(k_sub, val_frac = val_frac, sim_h = sim_h)
    samp_train <- sp$train
    samp_val   <- sp$val

    k_train <- dplyr::slice(k_sub, samp_train)
    k_val   <- dplyr::slice(k_sub, samp_val)

    ## controls sampled to MATCH each positive set's pep_end/cleavage mix
    train_dat <- bind_rows(k_train, balanced_ctrl_sample(k_train, c_sub, n_ctrl))
    val_dat   <- bind_rows(k_val,   balanced_ctrl_sample(k_val,   c_sub, n_ctrl))

    all_dat <- bind_rows(
      k_sub,
      ctr_dat %>% filter(target %in% targ)
    )

    list(train = train_dat,
         val   = val_dat,
         all   = all_dat)
  }, targets)
}

## Add end_type as a per-window input channel (1 = "end"/protein terminus,
## 0 = "cleavage"), constant across the 36 positions. Appended as the LAST data
## column so the padding channel keeps its index (which(all_params3=="padding"));
## the model script's n_channels then auto-updates 36 -> 37 with no edits. This
## lets the pooled N/C models condition on the end-vs-cleavage regime.
inject_end_type_ch <- function(d, et) {
  d %>% dplyr::mutate(end_type_ch = if_else(et == "end", 1, 0))
}
known_dat <- known_dat %>% mutate(data = map2(data, end_type, inject_end_type_ch))
c_dat     <- c_dat     %>% mutate(data = map2(data, end_type, inject_end_type_ch))

nn_input <- make_training_sets(k_dat = known_dat,
                               ctr_dat = c_dat)

## nn_input now holds two datasets: N and C, each list(train, val, all) pooling
## end + cleavage windows. Each window carries its regime in the `end_type_ch`
## input channel (and the `end_type` metadata column for outputs/plots).













