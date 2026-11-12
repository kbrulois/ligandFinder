## ---- MSA conservation along the precursor -------------------------------------
## The `curated_AAswap_msacons*` arms replace the Aminode conservation channel
## (`cons_rs_n`) with a per-residue score from each precursor's ColabFold MSA
## (inst/scripts/msa_cons/msa_conservation.py), min-max scaled into the channel
## `plm_cons_msa` (build_msacons_parquet.py). These functions draw that channel
## along a precursor beside the old one, with the known peptides and their ends,
## and summarise how both channels sit at known ends against every other
## residue. Driver: inst/scripts/msa_cons/plot_msa_conservation.R.
##
## Four things every figure here has to keep straight:
##   * residues EXCLUDED BY DESIGN (cons_rs_n, min_afm, mean_afm and relASA all
##     exactly 0: not extracellular in the topology table, or past AF fragment
##     1) are 0 in BOTH channels. That zero encodes topology, not conservation,
##     so it is drawn as a band, never as a value.
##   * cons_rs_n is also exactly 0 on residues that are NOT excluded: no Aminode
##     value, filled with 0 by 9.2 (122,170 residues on 2026-11-12, all of
##     MT-RNR2 among them). The model saw that 0, so it stays in every statistic,
##     but a plot marks it at the baseline rather than drawing "unconserved".
##     (Checked 2026-11-12 against processed/aminode.rds: POMC's table skips 197,
##     MT-RNR2 is not in it, and only 1 of 10.7M Aminode rates lies past the
##     rate at which 9.2's scaling clips to 0 -- so a 0 is a missing value.)
##   * cons_rs_n's LEVEL does not compare across precursors. Within one it is an
##     exact linear map of the Aminode rate (3,835 of 4,357 precursors in
##     aminode.rds fit to 1e-6; most of the rest number residues differently),
##     but there are 18 different maps, all reaching 0 at rate 0.75-0.78 and
##     differing in their top value: an invariant residue (rate 0) reads 1.00 in
##     INS, NPY, CALCA; 0.97 in POMC, ADCYAP1, PENK; 0.74 in GCG; 0.013-1.00
##     overall (a quarter of precursors below 0.38). The top value is the
##     precursor's species_limit: cons_rs_n = cons_rs x species_limit
##     (putative_peptide/generate_residue_db_untested/9_expand_by_residue.R),
##     where cons_rs is the rate reversed and rescaled to 0-1 across genes and
##     species_limit is the largest myo_sim among the species in the gene's
##     Aminode alignment (2.1_parse_aminode_new.R) -- how far back in time the
##     alignment reaches, 1 = the deepest species in inst/extdata/species_dat.rds.
##     (Reviewer, 2026-11-12: the per-precursor tops sit at species_dat's myo_sim
##     levels, 0.967 / 0.82 / 0.743 / 0.373 / 0.23 / ...) So the old level says how
##     deep the alignment is, not how conserved the residue is; the flat
##     stretches of the old line are invariant residues at that ceiling. That
##     goes with the old channel's 87% per-protein variance share, and is why the
##     summary also reports each AUC against the same precursors' own residues.
##   * a C end followed by G + a dibasic pair sits ON THE G in the _padg windows
##     (10_7 --amid-g-anchor), one past the mature end. Of the end tables only
##     lf_pepend_cv5_knowns_<T>_padg.csv is already in that convention; the
##     never-seen ends are moved by arm_conv.R's arm_anchor() rule.
## ------------------------------------------------------------------------------

#' Where the MSA-conservation inputs live
#' @keywords internal
LF_MSACONS_PATHS <- list(
  raw         = "~/AF2_analysis/lf_msa/msa_cons_raw.parquet",
  per_protein = "~/AF2_analysis/lf_msa/msa_cons_per_protein.tsv",
  msa_dir     = "~/AF2_analysis/lf_msa",
  channel     = "~/AF2_analysis/lf_plm/aa_onehot_msacons.parquet",
  cache       = "~/AF2_analysis/lf_pepend_residue_cache.rds",
  knowns      = "~/AF2_analysis/lf_pepend_cv5_knowns_%s_padg.csv",
  unseen      = "~/AF2_analysis/lf_nomask_scripts_20261021/cv5_unseen_ends_padg.csv",
  python      = "~/.virtualenvs/lf-plm/bin/python",
  ## the two _padg arms' scan metrics (old channel vs new), for the summary's caption
  arm_metrics = "~/AF2_analysis/lf_headcmp/unseen_padg/final_metrics.csv"
)

#' The end sets, in the order every figure and table here uses
#' @keywords internal
LF_MSACONS_SETS <- c("trained", "never-seen peptide", "never-seen large product")

#' Colours of the conservation figures, all taken from the shared tables
#'
#' New score = categorical slot 1; the old one in secondary ink, because it is
#' the reference the new channel is read against, and because the terminus
#' colours (green N, orange C, as on the protein pages) are already on the
#' panel as end marks -- a second orange would say "C end".
#' @keywords internal
lf_msacons_cols <- function() {
  tc <- getOption("lf.term_cols", LF_TERM_COLS)
  list(new = LF_VIZ$s1, old = LF_VIZ$ink2, N = tc[["N"]], C = tc[["C"]],
       peptide = nn_class_cols[["pep_other"]], excluded = LF_VIZ$grid,
       other_end = grDevices::adjustcolor(LF_VIZ$ink2, alpha.f = 0.55),
       ink = LF_VIZ$ink, ink2 = LF_VIZ$ink2)
}

#' Every cached residue with both conservation channels
#'
#' One row per residue of every precursor in the residue cache: the OLD channel
#' (`cons_rs_n`, exactly as 10_7 slices it into the windows), the NEW one
#' (`plm_cons_msa`, exactly as the msacons arms read it), the raw MSA score the
#' new one was scaled from, and whether the residue is excluded by design.
#'
#' Refuses to return a table the figures would mislabel: the residue letter
#' must agree between the channel parquet, the raw score and the cached
#' sequence, and the channel must follow from the raw score by the min-max
#' recorded in the parquet's `.scale.json` (0 on excluded residues), because
#' the y axis is labelled with that transform.
#' @param paths `LF_MSACONS_PATHS`, or a copy with some entries changed.
#' @return list: `res` (data.table: accession, index, AA, old, new, excluded,
#'   jsd, gap_frac, n_cov), `prec` (accession, gene, seq, n_prot, c_prot),
#'   `per` (per-protein MSA depth, neff, source), `scale` (lo, hi, n_covered).
#' @export
lf_msacons_residues <- function(paths = LF_MSACONS_PATHS) {
  if (!requireNamespace("arrow", quietly = TRUE)) stop("lf_msacons_residues needs the arrow package")
  p <- lapply(paths, path.expand)
  fc <- readRDS(p$cache)
  excl_ch <- c("cons_rs_n", "min_afm", "mean_afm", "relASA")
  res <- data.table::rbindlist(lapply(names(fc$feats), function(a) {
    f <- fc$feats[[a]]
    data.table::data.table(accession = a, index = seq_len(nrow(f)),
                           old = unname(f[, "cons_rs_n"]),
                           excluded = rowSums(f[, excl_ch, drop = FALSE] != 0) == 0)
  }))
  prec <- data.table::as.data.table(fc$prec)[, c("accession", "gene", "seq", "n_prot", "c_prot")]
  rm(fc)

  ch <- data.table::as.data.table(arrow::read_parquet(
    p$channel, col_select = c("accession", "index", "AA", "plm_cons_msa")))
  raw <- data.table::as.data.table(arrow::read_parquet(
    p$raw, col_select = c("accession", "index", "AA", "jsd", "gap_frac", "n_cov")))
  data.table::setnames(ch, c("AA", "plm_cons_msa"), c("AA_ch", "new"))
  res <- merge(res, ch, by = c("accession", "index"), all.x = TRUE)
  res <- merge(res, raw, by = c("accession", "index"), all.x = TRUE)
  if (anyNA(res$new) || anyNA(res$jsd))
    stop(sum(is.na(res$new)), " cached residues have no plm_cons_msa and ", sum(is.na(res$jsd)),
         " no raw MSA score; the parquets do not cover the residue cache")
  ## the letter, three ways: the channel's, the raw score's, the cache sequence's
  seq_aa <- prec[, .(index = seq_len(nchar(seq)), AA_seq = strsplit(seq, "")[[1]]), by = accession]
  res <- merge(res, seq_aa, by = c("accession", "index"), all.x = TRUE)
  bad <- res[AA != AA_ch | AA != AA_seq | is.na(AA_seq)]
  if (nrow(bad))
    stop(nrow(bad), " residues disagree on the residue letter, e.g. ",
         bad$accession[1], " ", bad$index[1], ": raw ", bad$AA[1], ", channel ", bad$AA_ch[1],
         ", cache ", bad$AA_seq[1])
  res[, c("AA_ch", "AA_seq") := NULL]

  sc <- jsonlite::fromJSON(paste0(p$channel, ".scale.json"))
  expect <- ifelse(res$excluded, 0, pmin(pmax((res$jsd - sc$lo) / (sc$hi - sc$lo), 0), 1))
  dev <- max(abs(res$new - expect))
  ## float32 in the parquet: agreement to ~1e-7 is exact
  if (dev > 1e-6) stop(sprintf("plm_cons_msa is not the recorded min-max of the raw score (max |dev| %.3g)", dev))
  if (sum(!res$excluded) != sc$n_covered)
    stop(sprintf("%d residues are not excluded by design, but the scale was fitted over %d",
                 sum(!res$excluded), sc$n_covered))

  per <- data.table::fread(p$per_protein)
  data.table::setkey(res, accession, index)
  list(res = res, prec = prec, per = per,
       scale = list(lo = sc$lo, hi = sc$hi, n_covered = sc$n_covered, max_dev = dev))
}

#' A mature C end moved onto the G of a G + dibasic site
#'
#' The _padg convention (10_7 --amid-g-anchor), with arm_conv.R's
#' `arm_anchor()` rule: a C end followed by G + KK/KR/RK/RR moves one residue
#' right, onto the G. N ends never move. Vectorised.
#' @param seq_of named character vector, sequence by accession.
#' @return integer anchors in the windows' convention.
#' @keywords internal
lf_msacons_g_anchor <- function(seq_of, accession, term, mature) {
  mature <- as.integer(mature)
  s <- unname(seq_of[as.character(accession)])
  nxt <- substr(s, mature + 1L, mature + 3L)
  mature + as.integer(term == "C" & !is.na(s) & grepl("^G(KK|KR|RK|RR)", nxt))
}

#' Known peptide ends, at the anchors the window models score
#'
#' \describe{
#'   \item{`trained`}{every known end of the five _padg folds
#'     (`lf_pepend_cv5_knowns_<T>_padg.csv`: each end is `val` in exactly one
#'     fold, so together they are all of them), already G-anchored.}
#'   \item{`never-seen peptide`, `never-seen large product`}{the two parts of
#'     `cv5_unseen_ends_padg.csv` -- annotated ends of genes no fold trained on
#'     -- moved onto the G where the _padg convention puts them.}
#' }
#' `motif` is [lf_pepend_motif_vec()] at the model anchor under the G
#' convention, so "G + dibasic" marks exactly the C ends that moved.
#' @param prec the `prec` table of [lf_msacons_residues()].
#' @return data.table: accession, gene, term, anchor, mature, set, motif; one
#'   row per distinct end and set.
#' @export
lf_msacons_ends <- function(prec, paths = LF_MSACONS_PATHS) {
  if (!"amid_g" %in% names(formals(lf_pepend_motif_vec)))
    stop("lf_pepend_motif_vec() predates the G-anchor convention; source this checkout's R/pepend_windows.R")
  seq_of <- stats::setNames(prec$seq, prec$accession)
  kn <- data.table::rbindlist(lapply(c("N", "C"), function(tm) {
    d <- data.table::fread(sprintf(path.expand(paths$knowns), tm))
    ## an N row's pep_end_mature is the peptide's OTHER end; N anchors never move
    data.table::data.table(accession = d$accession, term = tm, anchor = as.integer(d$anchor),
                           mature = if (tm == "C") as.integer(d$pep_end_mature) else as.integer(d$anchor),
                           set = "trained")
  }))
  ## the rule below must BE the training convention, or the never-seen ends
  ## would be marked on a different residue from the trained ones
  moved <- lf_msacons_g_anchor(seq_of, kn$accession, kn$term, kn$mature)
  if (any(moved != kn$anchor))
    stop(sum(moved != kn$anchor), " trained ends disagree with the G-anchor rule, e.g. ",
         kn$accession[which(moved != kn$anchor)[1]])
  un <- data.table::fread(path.expand(paths$unseen))
  un <- data.table::data.table(
    accession = un$accession, term = un$term, mature = as.integer(un$anchor),
    set = ifelse(grepl("^never seen", un$set), "never-seen peptide", "never-seen large product"))
  un$anchor <- lf_msacons_g_anchor(seq_of, un$accession, un$term, un$mature)
  ## the never-seen table has a row per PEPTIDE, so an end two peptides share is
  ## listed twice (hepcidin-20/-25 at C 84, CXCL12's four isoforms at N 22): 88
  ## rows are 87 peptide ends, 585 are 578. Count and select on distinct ends.
  un <- unique(un, by = c("accession", "term", "anchor", "set"))
  e <- rbind(kn, un[, names(kn), with = FALSE])
  e <- merge(e, prec[, .(accession, gene, seq, n_prot)], by = "accession")
  e[, motif := lf_pepend_motif_vec(seq[[1]], term[[1]], anchor, n_prot[[1]], amid_g = TRUE),
    by = .(accession, term)]
  e[, set := factor(set, levels = LF_MSACONS_SETS)]
  e[, c("seq", "n_prot") := NULL][order(accession, term, anchor)]
}

#' Known peptide spans on some precursors
#'
#' The curated knowns (`ligand_list.rds`, names as listed there) plus, on a gene
#' ligand_list does not cover, the never-seen set's own peptides (an N row
#' spans `anchor .. anchor + len - 1`, a C row `anchor - len + 1 .. anchor`).
#' One row per distinct (accession, start, end): ligand_list's name wins.
#' Spans are MATURE peptides; a G-anchored C end sits one residue past its span.
#' @param ll `inst/extdata/ligand_list.rds`.
#' @return data.table: accession, start, end, name.
#' @export
lf_msacons_spans <- function(accessions, ll, paths = LF_MSACONS_PATHS) {
  ll <- data.table::as.data.table(ll)
  a <- ll[accession %in% accessions & !is.na(start) & !is.na(end),
          .(accession, start = as.integer(start), end = as.integer(end), name = as.character(final_name))]
  un <- data.table::fread(path.expand(paths$unseen))[accession %in% accessions]
  b <- un[, .(accession,
              start = as.integer(ifelse(term == "N", anchor, anchor - len + 1L)),
              end   = as.integer(ifelse(term == "N", anchor + len - 1L, anchor)),
              name  = as.character(name))]
  b <- b[!accession %in% a$accession]
  s <- rbind(a, b)
  s[!duplicated(s[, .(accession, start, end)])][order(accession, start, -end)]
}

#' Listed peptides with an end that neither end table holds
#'
#' The peptides whose ends [lf_msacons_protein_plot()] draws in grey: the N end
#' is the span's start, the C end its end moved onto the G by
#' [lf_msacons_g_anchor()], and an end is listed when `ends` has that terminus
#' at that anchor. On a trained gene this is a peptide never docked -- the
#' trained ends are 10_7's docked peptides (knowns.rds, relevant site) --
#' e.g. POMC's beta-endorphin and GCG's GLP-2.
#' @param spans [lf_msacons_spans()]; `ends` [lf_msacons_ends()]; `prec` the
#'   `prec` table of [lf_msacons_residues()].
#' @return the rows of `spans` with an unlisted end.
#' @keywords internal
lf_msacons_unlisted <- function(spans, ends, prec) {
  if (!nrow(spans)) return(spans)
  seq_of <- stats::setNames(prec$seq, prec$accession)
  key <- paste(ends$accession, ends$term, ends$anchor)
  c_at <- lf_msacons_g_anchor(seq_of, spans$accession, "C", spans$end)
  spans[!paste(accession, "N", start) %in% key | !paste(accession, "C", c_at) %in% key]
}

#' Pick the example precursors, by a stated rule
#'
#' \enumerate{
#'   \item `n_classic` classics, from `classics`: the most G + dibasic trained C
#'     ends (where the _padg convention moved the end onto the G), then the
#'     most trained ends, then the deeper MSA (neff).
#'   \item one never-seen gene: from the never-seen peptide set (< 60 aa), the
#'     gene with the most ends there, ties to the deeper MSA.
#'   \item one more never-seen peptide gene, chosen so the excluded-by-design
#'     band has an example: the most ends among genes whose precursor carries
#'     excluded residues, ties to the shorter precursor.
#'   \item one shallow MSA: the lowest neff of any gene carrying a trained or
#'     never-seen end -- to show what the evidence prior does.
#' }
#' @param ends [lf_msacons_ends()]; `res`, `prec`, `per` from [lf_msacons_residues()].
#' @return data.table: accession, gene, group, why (one line, for the panel).
#' @export
lf_msacons_select <- function(ends, res, prec, per,
                              classics = c("ADCYAP1", "GCG", "POMC", "NPY", "CALCA", "PENK"),
                              n_classic = 3L) {
  g <- ends[, .(n_trained = sum(set == "trained"),
                n_gdb = sum(set == "trained" & term == "C" & motif == "G + dibasic"),
                n_unseen_pep = sum(set == "never-seen peptide")), by = .(accession, gene)]
  g <- merge(g, per[, .(accession, neff, depth)], by = "accession")
  g <- merge(g, res[, .(n_excl = sum(excluded), L = .N), by = accession], by = "accession")

  cl <- g[gene %in% classics & n_trained > 0][order(-n_gdb, -n_trained, -neff)][seq_len(n_classic)]
  cl[, `:=`(group = "classic",
            why = sprintf("classic precursor: %d G + dibasic C end%s among %d trained ends",
                          n_gdb, ifelse(n_gdb == 1, "", "s"), n_trained))]
  up <- g[n_unseen_pep > 0 & n_trained == 0]
  u1 <- up[order(-n_unseen_pep, -neff)][1]
  u1[, `:=`(group = "never-seen",
            why = sprintf("never-seen gene with the most never-seen peptide ends (%d; ties to the deeper MSA)",
                          n_unseen_pep))]
  u2 <- up[n_excl > 0 & accession != u1$accession][order(-n_unseen_pep, L)][1]
  u2[, `:=`(group = "never-seen",
            why = sprintf(paste0("never-seen peptide gene whose precursor has excluded-by-design residues (%d of %d): ",
                                 "the most such ends (%d), ties to the shorter precursor"), n_excl, L, n_unseen_pep))]
  taken <- c(cl$accession, u1$accession, u2$accession)
  sh <- g[!accession %in% taken][order(neff)][1]      # every gene in g carries an end
  sh[, `:=`(group = "shallow MSA",
            why = sprintf("shallowest MSA of any gene with a known end (%s)",
                          if (sh$n_trained > 0) sprintf("%d trained ends", sh$n_trained) else "never-seen"))]
  out <- rbind(cl, u1, u2, sh)
  if (anyNA(out$accession)) stop("an example rule matched no gene")
  out[, gene := prec$gene[match(accession, prec$accession)]]
  out[, .(accession, gene, group, why, n_trained, n_gdb, n_unseen_pep, n_excl, L, depth, neff)]
}

#' The raw MSA score of one precursor recomputed, e.g. without the prior
#'
#' Runs msa_conservation.py's own `score_one()` -- so the score cannot drift
#' from the one the channel was built from -- on the precursor's MSA with
#' `prior` background pseudo-sequences. `prior = 10` must reproduce
#' msa_cons_raw.parquet; the driver checks that before it draws `prior = 0`.
#' Python runs isolated (`-I`: neither the working directory nor PYTHON* env
#' vars reach its path; numpy comes from the user site, added back explicitly)
#' and writes no bytecode (`-B`) into the scripts directory.
#' @param source the per-protein table's `source` ("hf_colabfold",
#'   "hf_colabfold (via P04001)", ...): the directory, and the accession whose
#'   file holds an identical sequence's MSA.
#' @param script_dir the checkout's inst/scripts/msa_cons.
#' @return data.table: index, jsd, gap_frac, n_cov.
#' @export
lf_msacons_rescore <- function(accession, seq, source, prior, script_dir,
                               paths = LF_MSACONS_PATHS, max_depth = 50000L) {
  py <- path.expand(paths$python)
  if (!file.exists(py)) stop("no python at ", py)
  src <- sub(" .*$", "", source)
  via <- regmatches(source, regexpr("(?<=via )[A-Za-z0-9_-]+", source, perl = TRUE))
  fa <- if (length(via)) via else accession
  cand <- file.path(path.expand(paths$msa_dir), src, paste0(fa, c(".a3m.gz", ".a3m")))
  msa <- cand[file.exists(cand)][1]
  if (is.na(msa)) stop("no MSA for ", accession, " under ", src)
  code <- paste(
    "import site, sys",
    "sys.path.append(site.getusersitepackages())",
    "acc, seq, path, prior, depth, sdir = sys.argv[1:7]",
    "sys.path.insert(0, sdir)",
    "from msa_conservation import score_one",
    "_, jsd, gap, ncov, _, _ = score_one(acc, seq, path, int(depth), float(prior))",
    "print('\\n'.join('%d\\t%.9g\\t%.9g\\t%d' % (i + 1, j, g, n) for i, (j, g, n) in enumerate(zip(jsd, gap, ncov))))",
    sep = "\n")
  out <- system2(py, c("-I", "-B", "-c", shQuote(code), shQuote(accession), shQuote(seq), shQuote(msa),
                       prior, max_depth, shQuote(script_dir)), stdout = TRUE)
  if (!is.null(attr(out, "status"))) stop("msa_conservation.score_one failed for ", accession)
  d <- data.table::fread(text = out, header = FALSE, col.names = c("index", "jsd", "gap_frac", "n_cov"))
  if (nrow(d) != nchar(seq)) stop("rescore of ", accession, " returned ", nrow(d), " residues, not ", nchar(seq))
  d
}

#' Pack spans into lanes so neither bars nor their labels collide
#'
#' Greedy by start, longer first at a tie: a peptide takes the lower lane and
#' the shorter ones nested in it stack above. A label is drawn from the bar's
#' start, so a span occupies `max(end, start + label width)` -- the label width
#' estimated from the character count at `size` mm over a panel `panel_in`
#' inches wide. A label that would run past the precursor's last residue is
#' drawn ending at the bar's end instead (`right`), and that span occupies
#' `min(start, end - label width) .. end`. (Reviewer, 2026-11-12: PMCH's
#' "melanin-concentrating-hormone", 147-165 of 165, ran off the panel.)
#' @return data.frame: `lane`, `right` (label right-aligned on the bar's end).
#' @keywords internal
lf_msacons_lanes <- function(s, L, panel_in, size = 2.5) {
  if (!nrow(s)) return(data.frame(lane = integer(0), right = logical(0)))
  res_per_in <- L / panel_in
  chr_in <- size * 0.6 / 25.4                    # mean glyph width ~0.6 em
  lab_w <- (nchar(s$name) + 1.5) * chr_in * res_per_in
  right <- s$start - 0.5 + lab_w > L + 0.5
  left  <- ifelse(right, pmin(s$start - 0.5, s$end + 0.5 - lab_w), s$start - 0.5)
  reach <- ifelse(right, s$end + 0.5, pmax(s$end + 0.5, s$start - 0.5 + lab_w))
  o <- order(left, -s$end)
  lane_end <- numeric(0); lane <- integer(nrow(s))
  for (i in o) {
    free <- which(lane_end < left[i])
    k <- if (length(free)) free[1] else length(lane_end) + 1L
    lane[i] <- k; lane_end[k] <- reach[i]
  }
  data.frame(lane = lane, right = right)
}

#' Runs of TRUE in a logical vector, as start/end indices
#' @keywords internal
lf_msacons_runs <- function(x) {
  r <- rle(as.logical(x))
  e <- cumsum(r$lengths); s <- e - r$lengths + 1L
  data.frame(start = s[r$values], end = e[r$values])
}

#' The key entries of the per-protein figure, shared by its panels and its legend
#' @keywords internal
lf_msacons_keys <- function() list(
  series = c(new = "NEW: MSA conservation, plm_cons_msa (what the msacons arms read)",
             old = "OLD: Aminode cons_rs_n (the channel it replaced)",
             noprior = "NEW without the evidence prior (shallow example only)"),
  mark   = c(N = "N end", C = "C end",
             G = "C end on the G of a G + dibasic site (_padg anchor = mature end + 1)",
             u = "grey: end of a listed peptide in neither end table",
             z = "cons_rs_n = 0 on a scored residue: no Aminode value"),
  region = c(pep = "known peptide (mature span)",
             excl = "excluded by design (not extracellular): 0 in both channels"))

#' The scales every panel and the legend share, so a key cannot drift from its marks
#' @keywords internal
lf_msacons_scales <- function(cols = lf_msacons_cols()) {
  k <- lf_msacons_keys()
  list(
    ggplot2::scale_colour_manual(values = c(new = cols$new, old = cols$old, noprior = cols$new),
                                 labels = k$series, limits = names(k$series), name = NULL),
    ggplot2::scale_linetype_manual(values = c(new = "solid", old = "solid", noprior = "42"),
                                   labels = k$series, limits = names(k$series), name = NULL),
    ggplot2::scale_linewidth_manual(values = c(new = 0.65, old = 0.45, noprior = 0.5),
                                    labels = k$series, limits = names(k$series), name = NULL),
    ## all characters: one numeric shape would turn the whole vector into text
    ggplot2::scale_shape_manual(values = c(N = "[", C = "]", G = "G", u = "]", z = "x"),
                                labels = k$mark, limits = names(k$mark), name = NULL),
    ggplot2::scale_fill_manual(values = c(pep = cols$peptide, excl = cols$excluded),
                               labels = k$region, limits = names(k$region), name = NULL))
}

#' The per-protein figure's one legend, as a grob to stack above its panels
#'
#' Built from keys alone, so it lists every mark whether or not a given set of
#' proteins draws it -- patchwork's own guide collection is not used, since it
#' misplaces the guides under ggplot2 4.x here.
#' @export
lf_msacons_legend <- function(cols = lf_msacons_cols(), text_size = 7.5) {
  k <- lf_msacons_keys()
  kl <- data.frame(x = rep(0:1, 3), y = 0.5,
                   series = factor(rep(names(k$series), each = 2), levels = names(k$series)))
  km <- data.frame(x = 0, y = 0.5, mark = factor(names(k$mark), levels = names(k$mark)))
  kr <- data.frame(xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                   region = factor(names(k$region), levels = names(k$region)))
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = kl, ggplot2::aes(x, y, colour = series, linetype = series, linewidth = series)) +
    ggplot2::geom_point(data = km, ggplot2::aes(x, y, shape = mark)) +
    ggplot2::geom_rect(data = kr, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = region)) +
    lf_msacons_scales(cols) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, ncol = 1, override.aes = list(linewidth = c(0.9, 0.7, 0.7))),
      linetype = ggplot2::guide_legend(order = 1, ncol = 1),
      linewidth = ggplot2::guide_legend(order = 1, ncol = 1),
      shape = ggplot2::guide_legend(order = 2, ncol = 1, override.aes = list(
        colour = c(cols$N, cols$C, cols$C, cols$other_end, cols$old), size = c(3.4, 3.4, 3.4, 3.4, 2.2))),
      fill = ggplot2::guide_legend(order = 3, ncol = 1)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "top", legend.justification = "left",
                   legend.text = ggplot2::element_text(size = text_size, colour = cols$ink2),
                   legend.key.height = ggplot2::unit(10, "pt"), legend.key.width = ggplot2::unit(22, "pt"),
                   legend.box = "horizontal", legend.box.just = "top",
                   legend.spacing.x = ggplot2::unit(4, "pt"),
                   legend.box.spacing = ggplot2::unit(0, "pt"),
                   legend.background = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_blank())
  g <- ggplot2::ggplotGrob(p)
  g$grobs[[which(g$layout$name == "guide-box-top")]]
}

#' One precursor: both conservation channels on one 0-1 axis
#'
#' Above the 0-1 score area: the end marks (`[` N end, `]` C end, `G` a C end
#' on the G of a G + dibasic site -- the _padg anchor, one past the mature end;
#' in grey where the end is a listed peptide's but in neither end table) and
#' the known peptides in lanes, each labelled. In the score area: the
#' peptides shaded, dotted guides down from every end, residues excluded by
#' design as a band with no value drawn, and residues the old channel has no
#' value for (cons_rs_n exactly 0 on a scored residue) as crosses on the
#' baseline. Nothing is smoothed: both channels are per residue, as the model
#' reads them. The key is [lf_msacons_legend()].
#' @param r one precursor's rows of `res` ([lf_msacons_residues()]).
#' @param spans its [lf_msacons_spans()] rows; `ends` its [lf_msacons_ends()] rows.
#' @param meta list(gene, accession, depth, neff, source, why).
#' @param noprior optional data.frame(index, value): the score without the
#'   prior, already on the channel's scale, drawn dashed.
#' @param panel_in the panel's width in inches, for packing the labels.
#' @param x_title x axis title, or NULL (a stacked figure titles its last panel).
#' @return a ggplot; `attr(, "y_top")` is its y extent, for panel heights.
#' @export
lf_msacons_protein_plot <- function(r, spans, ends, meta, noprior = NULL, panel_in = 11,
                                    x_title = NULL, cols = lf_msacons_cols()) {
  r <- r[order(r$index), ]
  L <- nrow(r)
  lv <- names(lf_msacons_keys()$series)
  ln <- rbind(
    data.frame(index = r$index, value = ifelse(r$excluded, NA, r$new), series = "new"),
    data.frame(index = r$index, value = ifelse(r$excluded | r$old == 0, NA, r$old), series = "old"),
    if (!is.null(noprior))
      data.frame(index = noprior$index, value = ifelse(r$excluded, NA, noprior$value), series = "noprior"))
  ln$series <- factor(ln$series, levels = lv)
  zero <- r[!r$excluded & r$old == 0, ]
  excl <- lf_msacons_runs(r$excluded)

  ## peptides: their union shaded in the score area, each one a bar in a lane
  sp <- as.data.frame(spans)
  ln_sp <- lf_msacons_lanes(sp, L, panel_in)
  sp$lane <- ln_sp$lane
  ## a label starts just inside its bar's start, or -- where it would run off
  ## the panel -- ends just inside its bar's end
  pad <- 0.05 * L / panel_in
  sp$lab_x <- ifelse(ln_sp$right, sp$end + 0.5 - pad, sp$start - 0.5 + pad)
  sp$lab_hjust <- ifelse(ln_sp$right, 1, 0)
  n_lanes <- if (nrow(sp)) max(sp$lane) else 0L
  mark_y <- 1.08; lane0 <- 1.21; lane_h <- 0.135
  y_top <- lane0 + max(n_lanes - 1L, 0L) * lane_h + 0.07
  sh <- if (nrow(sp)) lf_msacons_runs(vapply(seq_len(L), function(i) any(sp$start <= i & sp$end >= i), logical(1)))
        else data.frame(start = integer(0), end = integer(0))

  ## end marks, one per residue and terminus (ACTH and alpha-MSH share N 138):
  ## the end tables' ends in the terminus colours, then -- in grey, so every
  ## shaded span shows its ends -- a listed peptide's ends that neither table
  ## holds (POMC's beta-endorphin, GCG's GLP-2), moved onto the G by the same
  ## rule
  em <- as.data.frame(ends)[, c("term", "anchor", "motif")]
  em$col <- ifelse(em$term == "N", cols$N, cols$C)
  if (nrow(sp)) {
    sq <- paste(r$AA, collapse = "")
    c_at <- lf_msacons_g_anchor(stats::setNames(sq, "p"), rep("p", nrow(sp)), "C", sp$end)
    oe <- data.frame(term = rep(c("N", "C"), each = nrow(sp)), anchor = c(sp$start, c_at))
    oe$motif <- ifelse(oe$term == "C", lf_pepend_motif_vec(sq, "C", oe$anchor, 1L, amid_g = TRUE), "")
    oe$col <- cols$other_end
    em <- rbind(em, oe)
  }
  em <- em[!duplicated(em[, c("term", "anchor")]), ]
  em$mark <- ifelse(em$term == "N", "N", ifelse(em$motif == "G + dibasic", "G", "C"))

  brk <- if (L <= 40) c(1, seq(5, L, 5)) else if (L <= 200) c(1, seq(25, L, 25)) else c(1, seq(50, L, 50))
  ttl <- sprintf("%s   %s, %d aa   |   MSA depth %s, neff %s, %s", meta$gene, meta$accession, L,
                 format(meta$depth, big.mark = ","), format(round(meta$neff), big.mark = ","),
                 lf_msacons_source_label(meta$source))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = sh, ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5, ymin = 0, ymax = 1),
                       fill = cols$peptide, alpha = 0.45) +
    ggplot2::geom_rect(data = excl, ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5, ymin = 0, ymax = 1),
                       fill = cols$excluded) +
    (if (nrow(excl)) ggplot2::geom_text(
      data = excl[excl$end - excl$start + 1 >= 0.08 * L, , drop = FALSE],
      ggplot2::aes(x = (start + end) / 2, y = 0.5, label = "excluded\nby design"),
      size = 2.6, colour = cols$ink2, lineheight = 0.9)) +
    (if (nrow(em)) ggplot2::geom_segment(
      data = em, ggplot2::aes(x = anchor, xend = anchor, y = 0, yend = mark_y - 0.04),
      colour = em$col, linewidth = 0.3, linetype = "22", alpha = 0.8)) +
    ggplot2::geom_line(data = ln, ggplot2::aes(index, value, colour = series, linetype = series,
                                               linewidth = series), na.rm = TRUE) +
    ggplot2::geom_point(data = zero, ggplot2::aes(index, 0), shape = "x", size = 2.2, colour = cols$old) +
    (if (nrow(em)) ggplot2::geom_point(data = em, ggplot2::aes(anchor, mark_y, shape = mark),
                                       colour = em$col, size = 3.4)) +
    ggplot2::geom_rect(data = sp, ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                                               ymin = lane0 + (lane - 1) * lane_h - 0.05,
                                               ymax = lane0 + (lane - 1) * lane_h + 0.05),
                       fill = cols$peptide, colour = NA) +
    ggplot2::geom_text(data = sp, ggplot2::aes(x = lab_x, y = lane0 + (lane - 1) * lane_h, label = name,
                                               hjust = lab_hjust),
                       size = 2.5, colour = cols$ink) +
    lf_msacons_scales(cols) +
    ggplot2::scale_x_continuous(breaks = brk, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", ".25", ".5", ".75", "1"),
                                expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0.5, L + 0.5), ylim = c(-0.03, y_top), clip = "off") +
    ggplot2::labs(title = ttl, subtitle = meta$why, x = x_title, y = "channel value") +
    lf_viz_theme(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9.5, face = "bold", colour = cols$ink),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = cols$ink2,
                                                         margin = ggplot2::margin(b = 2)),
                   axis.title.y = ggplot2::element_text(size = 7.5, hjust = 0.3),
                   axis.title.x = ggplot2::element_text(size = 8),
                   plot.margin = ggplot2::margin(12, 8, 2, 4),
                   legend.position = "none")
  attr(p, "y_top") <- y_top + 0.03
  p
}

#' How a per-protein table's MSA source reads on a figure
#' @keywords internal
lf_msacons_source_label <- function(source) {
  s <- sub(" .*$", "", source)
  lab <- c(hf_colabfold = "ColabFold MSA (LiteFold HF)", api_colabfold = "ColabFold MSA (API, env)")[s]
  via <- regmatches(source, regexpr("via [A-Za-z0-9_-]+", source))
  paste0(ifelse(is.na(lab), source, lab), if (length(via)) paste0(", ", via, "'s file") else "")
}

#' Both channels at known ends against every other residue
#'
#' For each end set: the residues within `k` of any of its ends (at the model's
#' anchors), against the background of every other scored residue -- scored =
#' not excluded by design, and "other" = not within `k` of an end of ANY set.
#' The SAME residues go into both channels, including the residues the old
#' channel has no value for (it is 0 there, and the model saw that 0).
#'
#' `auc` is the probability that a residue near an end outranks a background
#' residue (ties count half): 0.5 = no better than a random residue. Its
#' interval resamples GENES, since the residues of one precursor are not
#' independent. `auc_same` uses only the background residues of the same
#' precursors (does the channel pick the ends out within their own protein?);
#' `auc_cov` drops the no-Aminode residues from both channels.
#' @param res,ends from [lf_msacons_residues()] / [lf_msacons_ends()].
#' @return list(stats = one row per set x channel, long = the near-end and
#'   background values for drawing, k).
#' @export
lf_msacons_end_stats <- function(res, ends, k = 3L, n_boot = 2000L, seed = 1L) {
  sc <- res[excluded == FALSE, .(accession, index, old, new)]
  near <- unique(ends[, .(index = anchor + (-k:k)), by = .(accession, set, anchor)][, .(accession, set, index)])
  near <- merge(near, sc, by = c("accession", "index"))
  any_near <- unique(near[, .(accession, index)])
  bg <- sc[!any_near, on = c("accession", "index")]
  set.seed(seed)
  u_of <- function(x, b) {
    ## P(b < x) + P(b == x) / 2 for every x, against sorted b
    lt <- findInterval(x, b, left.open = TRUE); le <- findInterval(x, b)
    (lt + (le - lt) / 2) / length(b)
  }
  st <- data.table::rbindlist(lapply(levels(ends$set), function(s) {
    nn <- near[set == s]
    data.table::rbindlist(lapply(c("old", "new"), function(ch) {
      b <- sort(bg[[ch]])
      u <- u_of(nn[[ch]], b)
      g <- data.table::data.table(acc = nn$accession, u = u)[, .(S = sum(u), N = .N), by = acc]
      boot <- replicate(n_boot, { i <- sample.int(nrow(g), replace = TRUE); sum(g$S[i]) / sum(g$N[i]) })
      bs <- sort(bg[accession %in% nn$accession][[ch]])
      cov_n <- nn$old > 0; cov_b <- bg$old > 0
      data.table::data.table(
        set = s, channel = ch, n_res = nrow(nn), n_ends = nrow(ends[set == s]),
        n_genes = nrow(g), auc = mean(u),
        lo = unname(stats::quantile(boot, 0.025)), hi = unname(stats::quantile(boot, 0.975)),
        auc_same = mean(u_of(nn[[ch]], bs)),
        auc_cov = mean(u_of(nn[[ch]][cov_n], sort(bg[[ch]][cov_b]))),
        median_near = stats::median(nn[[ch]]), median_bg = stats::median(bg[[ch]]))
    }))
  }))
  long <- rbind(
    data.table::melt(near[, .(group = set, old, new)], id.vars = "group",
                     variable.name = "channel", value.name = "value"),
    data.table::melt(bg[, .(group = "other", old, new)], id.vars = "group",
                     variable.name = "channel", value.name = "value"))
  ## how much of each channel's variance is a per-PROTEIN level rather than
  ## residue-to-residue (share of the variance carried by the precursor means)
  between <- vapply(c(old = "old", new = "new"), function(ch) {
    v <- sc[[ch]]; stats::var(stats::ave(v, sc$accession)) / stats::var(v)
  }, numeric(1))
  list(stats = st, long = long, k = k, n_bg = nrow(bg), n_scored = nrow(sc),
       n_old0 = sum(sc$old == 0), between = between,
       spearman = stats::cor(sc$old, sc$new, method = "spearman"))
}

#' The summary panel: each channel's values near known ends vs elsewhere
#'
#' One panel per channel -- they sit on different scales (the old channel is
#' mostly a per-protein level, see `between` in [lf_msacons_end_stats()]), so
#' only the contrast WITHIN a panel means anything: violins of the channel at
#' residues within `k` of each end set and at every other scored residue, with
#' the AUC and its gene-bootstrap 95% interval over each end set.
#' @param es [lf_msacons_end_stats()].
#' @export
lf_msacons_summary_plot <- function(es, cols = lf_msacons_cols()) {
  st <- data.table::copy(es$stats)
  ch_lab <- c(old = "OLD: Aminode cons_rs_n", new = "NEW: MSA conservation, plm_cons_msa")
  grp_lv <- c(LF_MSACONS_SETS, "other")
  n_of <- st[channel == "new", stats::setNames(sprintf("%s ends\n%s residues", format(n_ends, big.mark = ","),
                                                         format(n_res, big.mark = ",")), set)]
  grp_lab <- c(`trained` = "trained known\nends", `never-seen peptide` = "never-seen\npeptide ends",
               `never-seen large product` = "never-seen large\ncleavage-product ends",
               other = "every other\nscored residue")
  xl <- paste0(grp_lab, "\n", c(n_of[LF_MSACONS_SETS], sprintf("%s residues", format(es$n_bg, big.mark = ","))))
  names(xl) <- grp_lv
  d <- data.table::copy(es$long)
  d[, `:=`(group = factor(group, levels = grp_lv), channel = factor(channel, levels = names(ch_lab)))]
  st[, `:=`(group = factor(set, levels = grp_lv), channel = factor(channel, levels = names(ch_lab)),
            lab = sprintf("AUC %.2f\n[%.2f, %.2f]", auc, lo, hi))]
  ggplot2::ggplot(d, ggplot2::aes(group, value)) +
    ggplot2::geom_violin(ggplot2::aes(fill = channel), colour = NA, scale = "width", width = 0.8,
                         alpha = 0.35, bw = 0.02) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = channel), width = 0.16, outlier.shape = NA,
                          linewidth = 0.4, fill = LF_VIZ$surface) +
    ggplot2::geom_text(data = st, ggplot2::aes(group, 1.1, label = lab), size = 2.6, colour = cols$ink,
                       lineheight = 0.9) +
    ggplot2::facet_wrap(~ channel, nrow = 1, labeller = ggplot2::as_labeller(ch_lab)) +
    ggplot2::scale_fill_manual(values = c(old = cols$old, new = cols$new), guide = "none") +
    ggplot2::scale_colour_manual(values = c(old = cols$old, new = cols$new), guide = "none") +
    ggplot2::scale_x_discrete(labels = xl) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", ".25", ".5", ".75", "1"),
                                limits = c(0, 1.18), expand = c(0.01, 0)) +
    ggplot2::labs(x = NULL, y = "channel value") +
    lf_viz_theme(base_size = 9) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", colour = cols$ink, size = 9),
                   axis.text.x = ggplot2::element_text(size = 7, lineheight = 0.95),
                   panel.spacing = ggplot2::unit(14, "pt"))
}
