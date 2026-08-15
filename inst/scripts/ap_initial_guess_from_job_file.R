#!/usr/bin/env Rscript
## =============================================================================
## Build AlphaPulldown initial-guess PDBs for every pair in a screening job file.
##
## Batch version of ap_initial_guess_generate_test_input.R.  Same method: take a
## receptor conformation that some peptide already opened (selected by
## ap_initial_guess_active_state_receptors.R, extracted by
## ap_extract_open_models.R), strip the resident peptide, and superpose the test
## peptide's pocket-binding terminus onto the pose the resident one left behind.
##
## Input is the same .txt AP consumes: one input string per line, e.g.
##
##     NK2R;TKN1,72-107
##     GPR15;GP15L,24-81
##     ACKR3,1-362;SDF1,22-93
##
## Proteins are separated by --protein_delimiter (default ";"), a protein from
## its residue range by --range_delimiter (default ","). Both must match what
## you pass to AlphaPulldown. Names are UniProt entry names; accessions are
## resolved via ligandFinder's `id_mapping`.
##
## OUTPUT NAMING IS THE WHOLE POINT.  AlphaPulldown looks a job's description up
## in the initial-guess map with a plain dict .get() -- a miss returns None and
## the job runs UNGUIDED with no warning.  So each file is named exactly as AP
## will describe the job: one component per protein, "<accession>" or
## "<accession>_<start>-<end>", joined with "_and_".  Those separators are AP's
## own and do not follow the input delimiters:
##   NK2R;TKN1,72-107  ->  P21452_and_P20366_72-107.pdb
## Confirm on the first run against the "Initial guess PDB for <desc>:" line in
## AP's log -- that string is the key AP actually looks up.
##
## Usage:
##   Rscript inst/scripts/ap_initial_guess_from_job_file.R \
##     --jobs      /oak/.../scripts/mxrun/job1.txt \
##     --out       /scratch/groups/ebutcher/kevin/ig_seeds \
##     --open_dir  /oak/stanford/groups/ebutcher/deorphan-AI-ze/open_models
##
## Then point AP at --initial_guess_dir <out>.  Pairs with no open receptor
## model are skipped and listed in the summary; they need the geometric
## fallback (position_ligand_initial_guess) or a receptor added upstream.
## =============================================================================

.libPaths("/home/groups/ebutcher/programs/pipeline/R_libs4.1")
suppressPackageStartupMessages({
  library(optparse)
  library(ligandFinder)
  library(dplyr)
})

option_list <- list(
  make_option(c("-j", "--jobs"), type = "character", default = NULL,
              help = "Job .txt file: one AP input string per line. REQUIRED."),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Output directory for seed PDBs. REQUIRED."),
  make_option("--af_db", type = "character",
              default = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/af_db",
              help = "AlphaFold-DB cache root [default: %default]"),
  make_option("--open_dir", type = "character",
              default = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/open_models",
              help = "Directory of extracted open receptor models [default: %default]"),
  make_option("--manifest", type = "character", default = NULL,
              help = "Open-model manifest CSV [default: <open_dir>/open_receptors_manifest.csv]"),
  make_option("--terminus", type = "character", default = "C",
              help = "Default ligand terminus that enters the pocket [default: %default]"),
  make_option("--terminus_map", type = "character", default = NULL,
              help = "Optional TSV: ligand entry name <TAB> terminus, overriding --terminus per ligand."),
  make_option("--conformer", type = "character", default = "best",
              help = "Which open model to use per receptor: 'best' (highest ipTM) or 'random' [default: %default]"),
  make_option("--n_anchor", type = "integer", default = 8,
              help = "Terminal CA count used for superposition [default: %default]"),
  make_option("--clash_cut", type = "double", default = 2.5,
              help = "Heavy-atom clash threshold in Angstroms, reported not fixed [default: %default]"),
  make_option("--protein_delimiter", type = "character", default = ";",
              help = paste("Delimiter between proteins in an input string. Must match the",
                           "--protein_delimiter you pass to AlphaPulldown [default: %default]")),
  make_option("--range_delimiter", type = "character", default = ",",
              help = "Delimiter between a protein name and its residue range [default: %default]"),
  make_option("--overwrite", action = "store_true", default = FALSE,
              help = "Rebuild seeds that already exist in --out."),
  make_option("--seed", type = "integer", default = 42,
              help = "RNG seed, used only when --conformer random [default: %default]"),
  make_option("--id_mapping", type = "character", default = NULL,
              help = paste("Path to id_mapping.rds. Defaults to the copy installed with",
                           "ligandFinder; set this when running from a source checkout."))
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$jobs) || !nzchar(opt$jobs)) stop("--jobs is required")
if (is.null(opt$out)  || !nzchar(opt$out))  stop("--out is required")
if (!file.exists(opt$jobs)) stop("Job file not found: ", opt$jobs)
if (!opt$conformer %in% c("best", "random")) stop("--conformer must be 'best' or 'random'")

manifest_csv <- opt$manifest
if (is.null(manifest_csv)) manifest_csv <- file.path(opt$open_dir, "open_receptors_manifest.csv")
if (!file.exists(manifest_csv)) stop("Manifest not found: ", manifest_csv)

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

## ---- helpers (kept inline so this runs without reinstalling the package) ----

kabsch <- function(P, Q) {
  cP <- colMeans(P); cQ <- colMeans(Q)
  s  <- svd(t(sweep(P, 2, cP)) %*% sweep(Q, 2, cQ))
  d  <- sign(det(s$v %*% t(s$u)))          # guard against a reflection
  list(R = s$v %*% diag(c(1, 1, d)) %*% t(s$u), cP = cP, cQ = cQ)
}
apply_fit <- function(X, f) sweep(sweep(X, 2, f$cP) %*% t(f$R), 2, f$cQ, "+")

ca_xyz <- function(pdb, ch) {
  a <- pdb$atom[pdb$atom$chain == ch & pdb$atom$elety == "CA", ]
  a <- a[order(a$resno), ]
  list(resno = a$resno, xyz = as.matrix(a[, c("x", "y", "z")]))
}

min_cross <- function(A, B) {
  d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B)
  sqrt(max(0, min(d2)))
}

AA3 <- c(ALA="A",ARG="R",ASN="N",ASP="D",CYS="C",GLN="Q",GLU="E",GLY="G",HIS="H",
         ILE="I",LEU="L",LYS="K",MET="M",PHE="F",PRO="P",SER="S",THR="T",TRP="W",
         TYR="Y",VAL="V",SEC="U",PYL="O")

## Cleave a residue range straight out of the PDB TEXT.  Reads columns 13-16
## rather than selecting on `elety`, so a file whose atom-name column holds
## element symbols fails with a message naming the real problem instead of
## silently yielding zero CA atoms.  Numbering stays as UniProt positions.
cleave_pdb_text <- function(src, start, end, out) {
  raw <- grep("^(ATOM|HETATM)", readLines(src, warn = FALSE), value = TRUE)
  if (!length(raw)) stop("no atom records in ", src)

  aname <- trimws(substr(raw, 13, 16))
  rname <- trimws(substr(raw, 18, 20))
  rno   <- suppressWarnings(as.integer(substr(raw, 23, 26)))

  if (!any(aname == "CA"))
    stop("no CA atoms in ", src, " -- atom-name column holds element symbols; ",
         "delete the file so af_structure_to_pdb() rebuilds it")

  if (is.na(start)) start <- min(rno, na.rm = TRUE)
  if (is.na(end))   end   <- max(rno, na.rm = TRUE)

  keep <- !is.na(rno) & rno >= start & rno <= end
  if (!any(keep))
    stop("no residues in ", start, "-", end, "; file covers ",
         paste(range(rno, na.rm = TRUE), collapse = "-"))

  ca  <- keep & aname == "CA"
  got <- paste0(ifelse(is.na(AA3[rname[ca]]), "X", AA3[rname[ca]])[order(rno[ca])],
                collapse = "")
  writeLines(c(raw[keep], "END"), out)
  list(path = out, seq = got, n_res = sum(ca))
}

## Split one AP input string into its per-protein parts.  Everything is fixed
## string matching, so delimiters that are regex metacharacters ("+", ".") are
## safe to configure.
## "NK2R;TKN1,72-107" -> list(list(name="NK2R", start=NA, end=NA),
##                            list(name="TKN1", start=72,  end=107))
parse_ap_input <- function(line, delim, rdelim) {
  parts <- strsplit(trimws(line), delim, fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  lapply(parts, function(p) {
    seg <- strsplit(p, rdelim, fixed = TRUE)[[1]]
    nm  <- trimws(seg[1])
    if (length(seg) >= 2L && nzchar(trimws(seg[2]))) {
      se <- as.integer(strsplit(trimws(seg[2]), "-", fixed = TRUE)[[1]])
      if (length(se) != 2L || anyNA(se))
        stop("could not read a start-end range from '", p, "'")
      list(name = nm, start = se[1], end = se[2])
    } else {
      list(name = nm, start = NA_integer_, end = NA_integer_)
    }
  })
}

## One component of AP's job description.
desc_component <- function(acc, start, end) {
  if (is.na(start)) acc else sprintf("%s_%d-%d", acc, start, end)
}

## ---- lookups ---------------------------------------------------------------

## id_mapping ships as data/id_mapping.rds.  data() cannot read .rds (it only
## handles .R, .rda/.RData and the text formats) and LazyData does not expose
## .rds either, so it has to be read from its installed path.
load_id_mapping <- function(explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) stop("--id_mapping not found: ", explicit)
    return(readRDS(explicit))
  }
  for (sub in c("data", "extdata")) {
    p <- system.file(sub, "id_mapping.rds", package = "ligandFinder")
    if (nzchar(p) && file.exists(p)) return(readRDS(p))
  }
  # Running from a source checkout rather than an installed package.
  for (p in c("data/id_mapping.rds", "inst/extdata/id_mapping.rds")) {
    if (file.exists(p)) return(readRDS(p))
  }
  stop("could not locate id_mapping.rds; pass --id_mapping /path/to/id_mapping.rds")
}

id_mapping <- load_id_mapping(opt$id_mapping)
stopifnot(all(c("Entry", "Entry Name") %in% names(id_mapping)))
entry2acc <- setNames(id_mapping[["Entry"]], id_mapping[["Entry Name"]])

resolve_acc <- function(nm) {
  acc <- unname(entry2acc[nm])
  if (is.na(acc)) stop("no UniProt accession for entry name '", nm, "'")
  acc
}

man_all <- readr::read_csv(manifest_csv, show_col_types = FALSE) %>%
  mutate(path = file.path(opt$open_dir, pdb_file)) %>%
  filter(file.exists(path))

term_override <- character()
if (!is.null(opt$terminus_map) && nzchar(opt$terminus_map)) {
  tm <- utils::read.delim(opt$terminus_map, sep = "\t", header = FALSE,
                          stringsAsFactors = FALSE)
  term_override <- setNames(toupper(trimws(tm[[2]])), trimws(tm[[1]]))
}

## ---- build one seed --------------------------------------------------------

build_one <- function(line) {
  parts <- parse_ap_input(line, opt$protein_delimiter, opt$range_delimiter)
  if (length(parts) != 2L)
    stop("expected exactly 2 proteins per line, got ", length(parts))

  rec <- parts[[1]]; lig <- parts[[2]]
  rec_acc <- resolve_acc(rec$name)
  lig_acc <- resolve_acc(lig$name)

  description <- paste(
    desc_component(rec_acc, rec$start, rec$end),
    desc_component(lig_acc, lig$start, lig$end),
    sep = "_and_"
  )
  out_pdb <- file.path(opt$out, paste0(description, ".pdb"))

  if (file.exists(out_pdb) && !opt$overwrite) {
    return(tibble::tibble(input = line, description = description,
                          output = basename(out_pdb), status = "exists",
                          conformer = NA_character_, opened_by = NA_character_,
                          anchor_rmsd = NA_real_, closest = NA_real_,
                          n_clash = NA_integer_, note = "already built"))
  }

  man <- man_all %>% filter(receptor == rec$name | receptor_id == rec_acc)
  if (nrow(man) == 0)
    stop("no open model for ", rec$name, " (", rec_acc,
         ") -- needs the geometric fallback or an upstream receptor run")

  pick <- if (opt$conformer == "best") man[which.max(man$iptm), ]
          else man[sample(nrow(man), 1), ]

  cx <- bio3d::read.pdb(pick$path, verbose = FALSE)
  n_res <- tapply(cx$atom$resno, cx$atom$chain, function(x) length(unique(x)))
  ord   <- names(sort(n_res, decreasing = TRUE))
  if (length(ord) < 2) stop("expected >=2 chains in ", basename(pick$path))
  rec_ch <- ord[1]; res_ch <- ord[2]

  ## test peptide, cleaved fresh from AF-DB
  lig_src <- af_structure_to_pdb(lig_acc, cache_dir = opt$af_db)
  cl <- cleave_pdb_text(lig_src, lig$start, lig$end,
                        out = file.path(tempdir(), paste0(lig_acc, "_frag.pdb")))
  ligp <- bio3d::read.pdb(cl$path, verbose = FALSE)

  ## superpose onto the resident pose, anchored on the pocket-binding terminus
  terminus <- if (lig$name %in% names(term_override)) term_override[[lig$name]]
              else toupper(opt$terminus)
  res_ca <- ca_xyz(cx,   res_ch)
  new_ca <- ca_xyz(ligp, unique(ligp$atom$chain)[1])
  if (nrow(res_ca$xyz) < 3 || nrow(new_ca$xyz) < 3)
    stop("too few CA atoms to superpose (resident ", nrow(res_ca$xyz),
         ", test ", nrow(new_ca$xyz), ")")

  k   <- min(opt$n_anchor, length(res_ca$resno), length(new_ca$resno))
  idx <- function(n, k) if (terminus == "C") seq.int(n - k + 1, n) else seq_len(k)
  fit <- kabsch(new_ca$xyz[idx(nrow(new_ca$xyz), k), , drop = FALSE],
                res_ca$xyz[idx(nrow(res_ca$xyz), k), , drop = FALSE])

  lig_xyz_new <- apply_fit(as.matrix(ligp$atom[, c("x", "y", "z")]), fit)
  anchor_rmsd <- sqrt(mean(rowSums((
    apply_fit(new_ca$xyz[idx(nrow(new_ca$xyz), k), , drop = FALSE], fit) -
    res_ca$xyz[idx(nrow(res_ca$xyz), k), , drop = FALSE])^2)))

  ligp$atom[, c("x", "y", "z")] <- lig_xyz_new
  ligp$xyz <- bio3d::as.xyz(as.vector(t(lig_xyz_new)))

  recp    <- bio3d::trim.pdb(cx, bio3d::atom.select(cx, chain = rec_ch, verbose = FALSE))
  rec_xyz <- as.matrix(recp$atom[, c("x", "y", "z")])

  closest <- min_cross(lig_xyz_new, rec_xyz)
  n_clash <- sum(apply(lig_xyz_new, 1, function(p)
                   min(sqrt(colSums((t(rec_xyz) - p)^2))) < opt$clash_cut))

  bio3d::write.pdb(bio3d::cat.pdb(recp, ligp, rechain = TRUE), file = out_pdb)

  tibble::tibble(
    input = line, description = description, output = basename(out_pdb),
    status = if (n_clash > 0) "ok_with_clashes" else "ok",
    conformer = as.character(pick$conformer), opened_by = as.character(pick$peptide),
    anchor_rmsd = anchor_rmsd, closest = closest, n_clash = as.integer(n_clash),
    note = sprintf("%s terminus, %d-CA anchor, source %s", terminus, k, basename(pick$path))
  )
}

## ---- run over the job file -------------------------------------------------

lines <- readLines(opt$jobs, warn = FALSE)
lines <- trimws(lines)
lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
message(sprintf("Job file: %s  (%d pair(s))", opt$jobs, length(lines)))
message("Output:   ", opt$out)
message("Conformer selection: ", opt$conformer)

rows <- vector("list", length(lines))
for (i in seq_along(lines)) {
  rows[[i]] <- tryCatch(
    build_one(lines[i]),
    error = function(e) tibble::tibble(
      input = lines[i], description = NA_character_, output = NA_character_,
      status = "failed", conformer = NA_character_, opened_by = NA_character_,
      anchor_rmsd = NA_real_, closest = NA_real_, n_clash = NA_integer_,
      note = conditionMessage(e))
  )
  r <- rows[[i]]
  message(sprintf("  [%3d/%3d] %-10s %s%s", i, length(lines), r$status, lines[i],
                  if (identical(r$status, "failed")) paste0("  -- ", r$note) else ""))
}

summary_tbl <- dplyr::bind_rows(rows)
summary_csv <- file.path(opt$out, paste0(tools::file_path_sans_ext(basename(opt$jobs)),
                                         "_initial_guess_summary.csv"))
readr::write_csv(summary_tbl, summary_csv)

message("\n", paste(rep("-", 60), collapse = ""))
print(table(summary_tbl$status))
n_ok <- sum(summary_tbl$status %in% c("ok", "ok_with_clashes", "exists"))
message(sprintf("%d of %d seed(s) available in %s", n_ok, nrow(summary_tbl), opt$out))
message("Summary: ", summary_csv)
if (any(summary_tbl$status == "failed"))
  message("Failed pairs run UNGUIDED if you point AP at this directory -- ",
          "check the note column before submitting.")
