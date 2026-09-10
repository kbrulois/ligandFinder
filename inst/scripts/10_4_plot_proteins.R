## =============================================================================
## Per-gene interactive protein plots (the <GENE>.html files).
##
## Prepares the plot inputs and renders them. Split out of 10_1dcnn_new6.R so the
## plots can be regenerated without retraining anything.
##
##   # in the session that has nn_input_comb (i.e. after 10_1dcnn_new6.R):
##   source("inst/scripts/10_4_plot_proteins.R")
##
##   # only some genes:
##   plot_genes <- c("NPY", "ANO8", "TAC1"); source("inst/scripts/10_4_plot_proteins.R")
##
## Needs in the session: nn_input_comb, secretome, secretome_aa, peps_tp.
## Writes the plot bundle at the end, so afterwards a FRESH session can re-render
## with no keras, no secretome_aa and no training:
##
##   Rscript inst/scripts/plot_only.R            # every gene in the bundle
##   Rscript inst/scripts/plot_only.R ANO8 NPY
##
## Knobs (set before sourcing, otherwise these defaults apply):
##   plot_dir     where the .html files go
##   plot_genes   character vector of gene symbols; default is every gene scored
##   write_bundle FALSE to skip writing plot_bundle.rds
## =============================================================================
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})

if (!exists("plot_dir"))     plot_dir     <- "~/AF2_analysis/new_meth_plot3"
if (!exists("write_bundle")) write_bundle <- TRUE

stopifnot(exists("nn_input_comb"), exists("secretome"))

## ---- score tracks to draw ---------------------------------------------------
mets <- list(cons        = c("blos_wt_mam", "blos_wt_all"),
             af_missense = c("mean_afm", "min_afm"),
             dssp        = c("relASA"),
             ## Smoothed peptide scores only: the chem_* tracks and the raw
             ## (unsmoothed) variants are deliberately excluded. These are already
             ## the _s6 columns, so nothing is appended -- the old
             ## `mets$aa_scores <- c(mets$aa_scores, paste0(mets$aa_scores, "_s6"))`
             ## would have asked for pep_nn4c_s6_s6.
             aa_scores   = c("pep_nn4c_s6", "pep_xgb4c_s6"))
all_mets <- do.call(`c`, mets) %>% unname
names(all_mets) <- rep(names(mets), sapply(mets, length))

## ---- attach per-residue scores to the secretome -----------------------------
## Idempotent AND atomic. Re-running this block used to left_join a second
## `aa_scores` onto a secretome that already had one, so dplyr disambiguated to
## aa_scores.x / aa_scores.y and secretome[["aa_scores"]] became NULL --
## surfacing one line down as "Can't recycle `.x` (size 0)". So drop any existing
## copy (plus the .x/.y wreckage from a previous failed run) first.
##
## Drop + join must be ONE pipeline: as two assignments, a failure in the join
## (e.g. secretome_aa not in the session) leaves secretome with the column
## already removed, and the next block dies with "object 'aa_scores' not found".
## Written this way, secretome is only reassigned once the join has succeeded.
stopifnot(exists("secretome_aa"))

secretome <- secretome %>%
  dplyr::select(-dplyr::any_of(c("aa_scores", "aa_scores.x", "aa_scores.y"))) %>%
  dplyr::left_join(
    secretome_aa %>% dplyr::select(!matches("_lead|_lag")) %>%
      dplyr::group_by(accession) %>% tidyr::nest(.key = "aa_scores"),
    by = "accession"
  )
stopifnot("aa_scores" %in% names(secretome))

secretome[["aa_scores"]] <- purrr::map2(
  secretome[["aa_scores"]], secretome[["sequence_uni"]],
  \(x, y) {
    if (is.null(x)) return(x)
    tmp <- tibble(index = 1:nchar(y), AA = stringr::str_split(y, "", simplify = TRUE) %>% c)
    dplyr::left_join(tmp, x, by = c("index", "AA"))
  }
)

## ---- which genes ------------------------------------------------------------
if (!exists("plot_genes")) plot_genes <- unique(nn_input_comb$gene)
genes <- plot_genes

## per-gene peptide tracks for the window plots (gene + nested `data`).
## Not produced anywhere in this script -- it comes from 13_plot_proteins.R /
## 12_extract_peptides.R -- so load it here unless it's already in the session.
if (!exists("peps_tp")) peps_tp <- readRDS("~/AF2_analysis/peps_tp.rds")

dir.create(path.expand(plot_dir), showWarnings = FALSE, recursive = TRUE)

## ---- one protein per gene ---------------------------------------------------
## group_split(gene) yields a MULTI-ROW group wherever a symbol maps to several
## accessions, and make_protein_plot_win assumes a single protein: it would
## expand several sequences onto one residue axis, emit several AlphaFold URLs,
## and write them all to the same <gene>.html. Longest sequence as a canonical
## proxy; with_ties = FALSE so the choice is deterministic.
the_input <- secretome %>%
  filter(!accession %in% c("A0AAG2TCD0", "A0AAG2UXZ5")) %>%
  filter(gene %in% !!genes) %>%
  slice_max(nchar(sequence_uni), n = 1, by = gene, with_ties = FALSE) %>%
  mutate(aa_scores = map(aa_scores, \(x) x[, colnames(x) %in% mets[["aa_scores"]]])) %>%
  group_split(gene)

if (!length(the_input)) stop("no genes matched: ", paste(utils::head(genes, 10), collapse = ", "))

## x[, colnames(x) %in% ...] is a set intersection with no fallback: when the
## score names drift -- 10_score_AA_xgboost.R derives the model suffix from
## nn[["neural_net"]][8], so it moves with the model list -- this silently yields
## a zero-column tibble and the score tracks disappear from the plots with no
## error at all. Fail loudly instead.
local({
  scored <- function(v) grep("^(pep|chem)_(nn|xgb)", v, value = TRUE)
  in_aa   <- if (exists("secretome_aa")) scored(names(secretome_aa)) else "<no secretome_aa>"
  nested  <- Filter(Negate(is.null), secretome[["aa_scores"]])
  in_sec  <- if (length(nested)) scored(names(nested[[1]])) else character(0)
  have    <- names(the_input[[1]][["aa_scores"]][[1]])
  want    <- mets[["aa_scores"]]

  message("score cols in secretome_aa       : ", paste(in_aa,  collapse = ", "))
  message("score cols in secretome$aa_scores: ", paste(in_sec, collapse = ", "))
  message("kept in the_input$aa_scores      : ", paste(have,   collapse = ", "))

  if (!length(have)) {
    hint <- if (!length(in_aa) || identical(in_aa, "<no secretome_aa>"))
      "secretome_aa has no score columns -- reload it: secretome_aa <- readRDS('~/peptide_alg/build_residue_db/processed/secretome_aa_1.rds')"
    else if (!length(in_sec))
      "secretome_aa has scores but secretome$aa_scores does not -- re-run the left_join block above"
    else
      paste0("names differ. secretome has: ", paste(in_sec, collapse = ", "),
             " -- update mets$aa_scores to match")
    stop("aa_scores came out empty. ", hint)
  }
  gone <- setdiff(want, have)
  if (length(gone))
    warning("score cols requested but not present: ", paste(gone, collapse = ", "))
})

## ---- pair the peptide tracks ------------------------------------------------
## Drop genes that have no peptide track. the two lists are paired POSITIONALLY,
## so a gene present in one but not the other shifts every later pair -- which
## would silently plot one gene's windows onto another gene's residues, for the
## whole rest of the run.
.tin_genes <- map_chr(the_input, \(x) as.character(x$gene)[1])
.pep_genes <- peps_tp %>% filter(gene %in% .tin_genes) %>% pull(gene) %>% unique()
if (any(!.tin_genes %in% .pep_genes))
  message("dropping ", sum(!.tin_genes %in% .pep_genes),
          " gene(s) with no peptide track: ",
          paste(utils::head(.tin_genes[!.tin_genes %in% .pep_genes], 10), collapse = ", "))
the_input <- the_input[.tin_genes %in% .pep_genes]

pep_input <- peps_tp %>%
  filter(gene %in% map_chr(the_input, \(x) as.character(x$gene)[1])) %>%
  group_split(gene)

## hard stop, not a printed TRUE/FALSE -- a misalignment here is invisible in the
## output and would corrupt every plot after the first mismatch
stopifnot(identical(map_chr(the_input,  \(x) as.character(x$gene)[1]),
                    map_chr(pep_input, \(x) as.character(x$gene)[1])))

pred_to_plot <- nn_input_comb

## make_protein_plot_win reads these from the calling environment rather than
## taking them as arguments -- they used to arrive incidentally from 9.2 and the
## training script. It catches its own errors per gene, so a missing one shows
## up as "no html written" rather than a stop(); load them explicitly.
if (!exists("species_dat"))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
if (!exists("id_map"))
  id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("classes"))
  classes <- setNames(seq_along(CLASS_NAMES) - 1L, CLASS_NAMES <- c(
    "CT_cleavage_context","DB","gap","NT_cleavage_context",
    "pep_other","pep_pocket","padding","none"))
for (.o in c("species_dat", "id_map", "all_mets", "classes"))
  if (!exists(.o)) stop("make_protein_plot_win needs `", .o, "` in the session", call. = FALSE)
rm(.o)

## ALWAYS source the repo copy of make_protein_plot_win -- never fall back to
## whatever is already attached. library(ligandFinder) exports an INSTALLED copy
## of this function, so an `if (!exists(...))` guard silently renders with stale
## code: it is the same name, it runs without error, and the only symptom is
## missing features in 5,000 html files. Sourcing unconditionally means the
## plots always match the source tree, installed package or not.
local({
  .cand <- c("R/plot_proteins_win_new.R",
             file.path(getwd(), "R", "plot_proteins_win_new.R"),
             "~/R_projects/ligandFinder/R/plot_proteins_win_new.R")
  .hit <- .cand[file.exists(path.expand(.cand))]
  if (!length(.hit))
    stop("R/plot_proteins_win_new.R not found -- run from the package root, or ",
         "devtools::load_all() the package before sourcing this script.",
         call. = FALSE)
  source(path.expand(.hit[[1]]), local = FALSE)
  message("make_protein_plot_win sourced from ", .hit[[1]])
})
stopifnot(is.function(make_protein_plot_win))

## Guard against the failure above ever recurring silently: if the ensemble
## columns are present, the plot function must be one that can draw them.
if ("per_index_sd" %in% names(nn_input_comb)) {
  ## the ribbon lives in make_detail_panel; the +/- sd in make_protein_plot_win
  .bad <- c(
    if (!grepl("geom_ribbon", paste(deparse(body(make_detail_panel)), collapse = "\n")))
      "make_detail_panel (no geom_ribbon)",
    if (!grepl("has_sd", paste(deparse(body(make_protein_plot_win)), collapse = "\n")))
      "make_protein_plot_win (no +/- sd)")
  if (length(.bad))
    stop("nn_input_comb carries per_index_sd but a stale copy is loaded -- ",
         paste(.bad, collapse = "; "), ". Probably the installed package ",
         "shadowing R/plot_proteins_win_new.R. The plots would silently omit ",
         "the ensemble spread.", call. = FALSE)
  rm(.bad)
}

## ---- render -----------------------------------------------------------------
message("plotting ", length(the_input), " genes -> ", plot_dir)
.t0 <- Sys.time()
for (i in seq_along(the_input)) {
  if (i %% 50 == 0 || i == 1)
    message(sprintf("  %d/%d  (%.1f min elapsed)", i, length(the_input),
                    as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
  make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]])
}
message(sprintf("done: %d genes in %.1f min", length(the_input),
                as.numeric(difftime(Sys.time(), .t0, units = "mins"))))

## ---- sanity-check one output ------------------------------------------------
local({
  f <- file.path(path.expand(plot_dir), paste0(the_input[[1]]$gene, ".html"))
  if (!file.exists(f)) { warning("no html written for ", the_input[[1]]$gene); return(invisible()) }
  html <- readLines(f, warn = FALSE)
  n <- function(p) sum(stringr::str_count(html, p))
  message(sprintf("%s.html: %d svg, %d nnwin_ ids, %d nn_show_panel hooks, %d status badges",
                  the_input[[1]]$gene, n("<svg"), n("data-id='nnwin_"),
                  n("nn_show_panel"), n("nn-status-badge")))
  if (n("data-id='nnwin_") == 0)
    warning("no nnwin_ data-ids in the svg -- the window click targets are missing")
})

## ---- bundle for re-rendering without a training session ---------------------
if (isTRUE(write_bundle)) source("inst/scripts/save_plot_bundle.R")
