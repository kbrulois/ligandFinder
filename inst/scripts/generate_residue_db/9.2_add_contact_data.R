

library(ligandFinder)
library(tidyverse)

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
#options(lf.rebuild_nn_dat = TRUE)
## ---- fast-iteration cache --------------------------------------------------
## known_dat + c_dat are the expensive part of this script: get_pep_data runs
## over ~40K windows, filtering the 2.8 GB secretome_aa per protein (~1 hr). But
## everything BELOW the cache (win_type filter + make_training_sets) is cheap.
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
  chan_range  <- .cc$chan_range

} else {

## Only the rebuild path needs these (~4.8 GB combined); a cached run doesn't.
## Skipped if they're already in the session from an earlier pipeline step.
s_localDir <- "~/peptide_alg/build_residue_db"

if (!exists("secretome")) {
  message("9.2: loading secretome ...")
  secretome <- readRDS("~/AF2_analysis/secretome_latest.rds")
}
if (!exists("secretome_aa")) {
  message("9.2: loading secretome_aa (~3 GB) ...")
  secretome_aa <- readRDS(file.path(s_localDir, "processed/secretome_aa.rds"))
}

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
                            p2_end = stringr::str_extract(p2_range, "\\d+$")))

## Collapse duplicate docking entries for the SAME peptide, keeping the best
## iptm. The key is the full range.
##
## It used to be two passes, one keyed on p2_end alone and one on p2_start
## alone, which also collapsed peptides that merely SHARE an end -- and for this
## class of ligand that is the normal case, not an edge case: a long and a short
## form sharing an amidated C-terminus (apelin-13/-17/-36, GRP/neuromedin C,
## kisspeptin-10/-54, NPY, PYY, TAC1) or sharing an N-terminus (PACAP-27/-38,
## POMC products, PTHLH, VIP). 27 real peptides were being dropped that way,
## with no message -- ADCYAP1 kept PACAP-27 and silently lost PACAP-38.
##
## Keying on start AND end keeps genuinely distinct peptides and still removes
## true duplicates. Set dedup_shared_termini = TRUE to restore the old
## behaviour for comparison.
dedup_shared_termini <- FALSE

if (isTRUE(dedup_shared_termini)) {
  con_dat <- con_dat %>%
              group_by(lig1_end_t, lig1_end_ind, p2_name, p2_end) %>%
              arrange(desc(iptm)) %>% slice(1) %>% ungroup() %>%
              group_by(lig1_end_t, lig1_end_ind, p2_name, p2_start) %>%
              arrange(desc(iptm)) %>% slice(1) %>% ungroup()
} else {
  .before <- nrow(con_dat)
  con_dat <- con_dat %>%
              group_by(lig1_end_t, lig1_end_ind, p2_name, p2_start, p2_end) %>%
              arrange(desc(iptm)) %>% slice(1) %>% ungroup()
  message(sprintf("9.2: range dedup %d -> %d peptides (old shared-terminus rule kept 118)",
                  .before, nrow(con_dat)))
}



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

## ---- amino-acid encoding: physicochemical properties, not one-hot ----------
## The 21 one-hot AA_* channels are replaced by 4 continuous property channels,
## looked up from the residue letter (secretome_aa$AA). This drops the input
## from 36 to 19 channels and lets the model generalise across chemically
## similar residues (D~E, K~R, I~L~V) instead of treating all 20 as unrelated.
##
## hydro  Kyte-Doolittle hydropathy index (J Mol Biol 157:105, 1982)
## charge net charge of the side chain at pH 7 (His = +0.1, partial protonation)
## mw     residue (in-chain) molecular weight, Da = free AA minus water
## pI     isoelectric point of the free amino acid
## Sec (U) is included; any other letter (X/B/Z) falls through to NA -> 0.
aa_props_raw <- tibble::tribble(
  ~AA, ~hydro, ~charge,   ~mw,   ~pI,
  "A",    1.8,     0.0,  71.08,  6.00,
  "R",   -4.5,     1.0, 156.19, 10.76,
  "N",   -3.5,     0.0, 114.10,  5.41,
  "D",   -3.5,    -1.0, 115.09,  2.77,
  "C",    2.5,     0.0, 103.14,  5.07,
  "Q",   -3.5,     0.0, 128.13,  5.65,
  "E",   -3.5,    -1.0, 129.12,  3.22,
  "G",   -0.4,     0.0,  57.05,  5.97,
  "H",   -3.2,     0.1, 137.14,  7.59,
  "I",    4.5,     0.0, 113.16,  6.02,
  "L",    3.8,     0.0, 113.16,  5.98,
  "K",   -3.9,     1.0, 128.17,  9.74,
  "M",    1.9,     0.0, 131.19,  5.74,
  "F",    2.8,     0.0, 147.18,  5.48,
  "P",   -1.6,     0.0,  97.12,  6.30,
  "S",   -0.8,     0.0,  87.08,  5.68,
  "T",   -0.7,     0.0, 101.10,  5.60,
  "W",   -0.9,     0.0, 186.21,  5.89,
  "Y",   -1.3,     0.0, 163.18,  5.66,
  "V",    4.2,     0.0,  99.13,  5.96,
  "U",    2.5,    -1.0, 150.04,  5.47)

## Min-max scale each property to [0,1] across the residues above, so the four
## channels share the dynamic range of relASA and the SS one-hots rather than
## letting mw (57-186) dominate. Padding positions carry 0 in all four, exactly
## as the old all-zero one-hot vector did, and the `padding` channel flags them.
aa_props <- aa_props_raw %>%
  mutate(across(-AA, ~ (.x - min(.x)) / (max(.x) - min(.x)))) %>%
  rename(AA_hydro = hydro, AA_charge = charge, AA_mw = mw, AA_pI = pI)

aa_prop_cols <- c("AA_hydro", "AA_charge", "AA_mw", "AA_pI")

## Backbone geometry + DSSP H-bond energetics, straight from secretome_aa
## (same names 9_expand_by_residue.R builds them under, where they are the
## `angles`/`energy` blocks of the nn2 feature set).
## Phi/Psi are supplied sin/cos-encoded so the model sees the angle as a point
## on the unit circle rather than a value with a wrap-around discontinuity.
angles <- c("Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin")
energy <- c("NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy")

## NOTE: "cons_rs" is deliberately excluded from the input channels; the
## normalised variant "cons_rs_n" is retained.
all_params3 <- c("cons_rs_n", "min_afm", "mean_afm", "relASA", "SS_P",
                 "SS_S", "SS_E", "SS_-", "SS_T", "SS_G", "SS_B", "SS_H", "SS_I",
                 angles, energy, aa_prop_cols, "padding")


## ---- global [0,1] scaling ---------------------------------------------------
## Every channel is min-max scaled onto [0,1] so none dominates by magnitude
## (raw ranges differ wildly: relASA ~[0,1], Phi/Psi sin-cos [-1,1], DSSP H-bond
## energies ~[-3,0]). Ranges are computed ONCE over all of secretome_aa -- they
## must be global, not per-protein, or the same residue would scale differently
## in different proteins -- and are cached so later data gets the same transform.
##
## The AA_* property channels are already min-max scaled in aa_props, and
## `padding` is a 0/1 flag, so both are excluded here.
##
## NOTE this changes what a padded position means for the signed channels. Raw 0
## used to be both "no H-bond" and "padding"; after scaling a real no-H-bond
## residue sits at the top of the energy range while padding stays 0, so the two
## are no longer conflated. Padding is still flagged by the `padding` channel.
scale_cols <- setdiff(all_params3, c(aa_prop_cols, "padding"))

if (!exists("chan_range")) {
  message("9.2: computing global [0,1] scaling ranges over secretome_aa ...")
  chan_range <- vapply(scale_cols, function(cn) {
    v <- secretome_aa[[cn]]
    r <- range(v, na.rm = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[2] <= r[1]) c(0, 1) else r  # constant/empty -> identity
  }, numeric(2))
  colnames(chan_range) <- scale_cols
  print(t(chan_range))
}

scale_01 <- function(d, cols = scale_cols, rng = chan_range) {
  for (cn in cols) {
    lo <- rng[1, cn]; hi <- rng[2, cn]
    d[[cn]] <- pmin(pmax((d[[cn]] - lo) / (hi - lo), 0), 1)   # clamp out-of-range to [0,1]
  }
  d
}

get_pep_data <- function(window, p_id, gene, nn_params = all_params3, meta_dat = NULL,
                         seq_uni = NULL) {

  ## fail loudly (and once) if a requested channel isn't in secretome_aa,
  ## instead of surfacing as an opaque select() error deep in the rowwise loop
  .want <- setdiff(nn_params, c(aa_prop_cols, "padding"))
  .miss <- setdiff(.want, colnames(secretome_aa))
  if (length(.miss))
    stop("secretome_aa is missing requested channel(s): ", paste(.miss, collapse = ", "))

  ## ---- residue scaffold ------------------------------------------------
  ## secretome_aa is SPARSE for most proteins: 2914/4907 accessions either start
  ## past residue 1 or have interior gaps (ANO8 = 719 rows spanning 2..1232, 512
  ## missing). The old `mutate(index = row_number())` overwrote secretome_aa's
  ## authoritative `index` with a positional count, which is only correct for the
  ## 40% of proteins whose table is contiguous from 1. For the rest every residue
  ## was renumbered, and since slice_in_data() selects on `dat$index %in% N:C`
  ## with N/C in true sequence_uni coordinates, each window then picked up the
  ## WRONG residues -- wrong feature channels, not just wrong letters on a plot.
  ##
  ## Instead build one row per precursor residue and left-join the sparse table
  ## onto it (the same scaffold pattern 10_1dcnn_new6.R uses for aa_scores).
  ## Joining on index AND AA means a sequence-version mismatch surfaces as NA
  ## rather than as silently misaligned data. Residues with no secretome_aa row
  ## keep their true coordinate and fall through to 0 like any other missing
  ## value; AA and the AA_* property channels are always correct because they
  ## come from the scaffold.
  sa <- secretome_aa %>% filter(accession == p_id)

  scaffold <- if (!is.null(seq_uni) && !is.na(seq_uni) && nchar(seq_uni) > 0) {
    tibble(index = seq_len(nchar(seq_uni)),
           AA    = stringr::str_split(seq_uni, "", simplify = TRUE) %>% c())
  } else {
    ## no sequence supplied: fall back to the table's own span (still keeps the
    ## authoritative index, just cannot fill residues absent from both).
    sa %>% select(index, AA) %>% filter(!is.na(index)) %>% arrange(index)
  }

  dat <- scaffold %>%
    left_join(sa, by = intersect(c("index", "AA"), colnames(sa))) %>%
    mutate(padding = 0) %>%
    ## property channels from the residue letter (unmatched letters -> NA -> 0)
    left_join(aa_props, by = "AA") %>%
    select(all_of(c(nn_params, "AA", "index", "topo2", "known"))) %>%
    scale_01() %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

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
    tmp <- tibble(peps = character(), dat = list(), win_type = character(), target = character())
  } else {
  if(nrow(window) == 0) {
    tmp <- tibble(peps = character(), dat = list(), win_type = character(), target = character())
  } else {
    tmp <- window %>%
      ungroup() %>%
      mutate(peps = paste0(gene, "_w", wN, "-", wC)) %>%
      mutate(dat = pmap(list(wN, wC, wN_pad, wC_pad), slice_in_data)) %>%
      select(any_of(c("peps", "dat", "win_type", "target")))


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
    filter(!end < n_prot)

  if(nrow(tmp) > 0) {
    return(
      tmp %>%
        rowwise() %>%
        ## C windows anchor on the site's last residue, N windows on its first.
        ## This matches the db_ind convention in get_adj_db_sites(), which puts
        ## the two-residue site at window positions 30-31 (C) and 6-7 (N) --
        ## the same slots make_known_idx() hard-codes as "DB".  Anchoring N on
        ## `end` shifts every candidate N window one residue relative to the
        ## knowns.  Single-residue pep_end rows have start == end, so they are
        ## unaffected either way.
        mutate(if_else(target == "C",
                       list(tibble(wN = end + ws[["C"]][[1]],
                                   wC = end + ws[["C"]][[2]])),
                       list(tibble(wN = start + ws[["N"]][[1]],
                                   wC = start + ws[["N"]][[2]]))) %>% bind_rows) %>%
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
    reframe(get_pep_data(window = windows, p_id = accession, gene = gene,
                         seq_uni = sequence_uni),
            gene = gene)



known_dat <- con_dat %>%
                rowwise() %>%
                reframe(get_pep_data(window = windows, p_id = accession, gene = p2_gene, meta_dat = clean_cons,
                                     seq_uni = sequence_uni),
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

make_known_idx <- function(meta_data, target, pep_id) {

  pep_backbone <- rep("gap", 36)

  pep_inds <- stringr::str_extract(pep_id, "\\d+x\\d+") %>% stringr::str_split(., "x", simplify = TRUE) %>% `c` %>% as.integer
  w_inds <- which(meta_data[["index"]] %in% pep_inds)

  if(target %in% c("C", "loop_C")) {
    if(length(w_inds) == 1) {w_inds <- c(1, w_inds)}
    pep_backbone[30:36] <- c(rep("DB", 2), rep("CT_cleavage_context", 5))
    pep_backbone[w_inds[1]:w_inds[2]] <- "pep_other"
    if(w_inds[1] != 1) {
    pep_backbone[1:(w_inds[1] - 1)] <- "NT_cleavage_context"
    }
  } else {
    if(length(w_inds) == 1) {w_inds <- c(w_inds, 36)}
    pep_backbone[1:7] <- c(rep("NT_cleavage_context", 5), rep("DB", 2))
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
                mutate(meta_data = pmap(list(meta_data, target, pep_id), make_known_idx)) %>%
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

saveRDS(list(known_dat = known_dat, c_dat = c_dat, all_params3 = all_params3,
             chan_range = chan_range),
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

## `banned` removes known-ligand genes from the CONTROL pool in
## make_training_sets(). It is matched against c_dat$gene / known_dat$gene, which
## are HGNC symbols (p2_gene comes from id_map "Gene Names (primary)").
##
## But gpcr_ligs$gene is stripped out of the CP-style `lig` id ("hGALAx33x62"),
## so it is a UniProt ENTRY NAME, not a symbol -- 538/548 of its values resolve
## as Entry Names, only 331 as symbols. Where the two differ the gene was never
## filtered: CCK/CCKN, CRH/CRF, GAL/GALA, GHRH/SLIB, TAC3/TKNK, PPY/PAHO,
## PTH2/TIP39, QRFP/OX26, PRLH/PRRP, APELA/ELA, MT-RNR2/HUNIN all leaked into the
## negatives while the same window was also a labelled positive.
##
## Resolve BOTH directions and union, so an id given either way is caught.
.entry2sym <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])
.sym2entry <- setNames(id_map[["Entry Name"]], id_map[["Gene Names (primary)"]])

banned <- c(gpcr_ligs$gene, uni_pep$Gene)
banned <- c(banned,
            unname(.entry2sym[banned]),    # entry name -> symbol
            unname(.sym2entry[banned])) %>% # symbol -> entry name
  na.omit() %>% unique()

message(sprintf("9.2: banned = %d genes (both namespaces resolved)", length(banned)))





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


## Two models: N and C. Training is restricted to win_type == "db" windows
## (see the filter below), so every window is a dibasic-anchored internal cut
## and the old end/cleavage (end_type) split axis no longer applies: it was
## constant across the training set and has been removed throughout.
## Controls are a plain random sample of the db windows for that terminus.
balanced_ctrl_sample <- function(pos, ctrl, n_ctrl) {
  n_take <- min(n_ctrl, nrow(ctrl))
  dplyr::slice_sample(ctrl, n = n_take)
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

    ## controls: random sample of same-terminus db windows
    train_dat <- bind_rows(k_train, balanced_ctrl_sample(k_train, c_sub, n_ctrl))
    val_dat   <- bind_rows(k_val,   balanced_ctrl_sample(k_val,   c_sub, n_ctrl))

    ## `all` = every window for this terminus, scored downstream as nn_input_comb.
    ## De-duplicated on `peps`: the knowns path (get_adj_db_sites) and the
    ## candidate scanner (make_candidate_windows) independently emit the SAME
    ## coordinates for a known peptide's cut site, so e.g. CCK_w76-111 arrived
    ## twice -- once known = 1, once known = 0. Keep the known copy: k_sub is
    ## bound first and distinct() keeps the first occurrence. The candidate
    ## scanner also re-emits a window under more than one win_type/target, which
    ## this collapses too.
    all_dat <- bind_rows(
      k_sub,
      ctr_dat %>% filter(target %in% targ)
    ) %>%
      dplyr::distinct(peps, .keep_all = TRUE)

    list(train = train_dat,
         val   = val_dat,
         all   = all_dat)
  }, targets)
}

## ---- db-only training set --------------------------------------------------
## Restrict positives AND controls to dibasic-anchored windows. Applied here,
## below the cache, so changing it is a seconds-long re-run rather than a full
## ~1 hr rebuild. Data channels are unchanged (no end_type_ch), so the model
## script's n_channels stays at length(all_params3).
message(sprintf("9.2: win_type filter -> db. knowns %d -> %d, controls %d -> %d",
                nrow(known_dat), sum(known_dat$win_type == "db"),
                nrow(c_dat),     sum(c_dat$win_type == "db")))

known_dat <- known_dat %>% filter(win_type == "db")
c_dat     <- c_dat     %>% filter(win_type == "db")

nn_input <- make_training_sets(k_dat = known_dat,
                               ctr_dat = c_dat)

## nn_input holds two datasets: N and C, each list(train, val, all), built only
## from db (dibasic-anchored) windows.













