## The same page, drawn from ONE ensemble member. C = seed 14, N = seed 7:
## the members that best retrieve the held-out ends. Both are oracle picks --
## chosen by scoring against the very ends being judged -- so these pages are an
## upper bound on single-seed behaviour, not a reproducible configuration. On C
## seed 14 is also what train-end selection picks; on N that rule picks seed 2.
suppressMessages({library(tidyverse); library(ggiraph); library(patchwork)})
suppressWarnings(suppressMessages(library(ligandFinder)))
## package source tree: $LF_SRC, or the checkout this script sits in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
A2 <- path.expand("~/AF2_analysis"); S2 <- Sys.getenv("SCRATCH2")
out_dir <- file.path(A2, "seed_pages"); dir.create(out_dir, showWarnings = FALSE)
.b <- readRDS(file.path(A2, "plot_bundle.rds"))
list2env(Filter(Negate(is.null), .b), .GlobalEnv)
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))
if (!exists("species_dat") || is.null(species_dat))
  species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))
source(lf_src("R/plot_proteins_win_new.R"))
source(lf_src("R/pepend_tracks.R"))
options(
  lf.track_scan_C   = file.path(A2, "lf_pepend_scan_C_seed14.rds"),
  lf.track_scan_N   = file.path(A2, "lf_pepend_scan_N_seed7.rds"),
  lf.track_winxgb_C = file.path(A2, "lf_winxgb_scan_C_seed14.npz"),
  lf.track_winxgb_N = file.path(A2, "lf_winxgb_scan_N_seed7.npz"),
  lf.track_member_C = 14L, lf.track_member_N = 7L,
  lf.class_tracks_by_term = TRUE)
gs <- purrr::map_chr(the_input, \(x) x$gene)
want <- readLines(file.path(S2, "val_genes.txt"))
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))$prec
want <- want[order(fc$c_prot[match(want, fc$gene)] - fc$n_prot[match(want, fc$gene)])]
done <- 0
for (g in want) {
  f <- file.path(out_dir, paste0(gsub("/", "_", g), ".html"))
  if (file.exists(f)) next
  i <- which(gs == g)[1]
  if (is.na(i)) { message("!!! ", g, " not in the bundle"); next }
  tm <- tryCatch(system.time(
    make_protein_plot_win(the_input[[i]], pred_to_plot, out_dir, pep_input[[i]]))[["elapsed"]],
    error = function(e) { message("!!! ", g, ": ", conditionMessage(e)); NA_real_ })
  done <- done + 1
  message(sprintf("=== %-10s %5.1f s  [%d]", g, tm, done))
}
message("FINISHED ", done)
