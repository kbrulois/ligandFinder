## Stage 2a of the candidate pipeline: the two FASTA files the MMseqs2 search
## needs, plus a target table carrying everything we want to report about a hit
## (Kevin, 2026-10-23). Run through candidate_pipeline.sh.
##
## TARGETS ("all known peptides"), deduplicated by sequence:
##   ligand_list.rds   the curated ligands                  type "ligand_list"
##   UniProt           peptide / propeptide features,       type "peptide" /
##                     sliced out of the precursor               "propeptide"
## A propeptide is a real cleavage product but NOT a signalling ligand, and it
## is 40% of this set, so the type travels with every hit.
##
## FASTA ids are surrogates (t000001...), not names: names contain characters
## FASTA and the blast-tab format would mangle, and two annotations can share a
## name. known.csv maps the surrogate to the gene, the coordinates and the type.
##
## known_id is <gene>x<start>x<end> -- the same shape as a candidate peptide_id,
## so a hit can be compared with a candidate directly.
suppressMessages({library(data.table)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
HEAD <- Sys.getenv("HEAD", "xgb_pm")
thr  <- as.numeric(Sys.getenv("THRESH", "0.5"))
TOL  <- as.integer(Sys.getenv("TOL", "4"))
O    <- file.path(A2, "mmseq_cand"); dir.create(O, showWarnings = FALSE, recursive = TRUE)

fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
prec <- as.data.table(fc$prec)[, .(accession, gene, pseq = seq)]; rm(fc)

ll <- as.data.table(readRDS(lf_src("inst/extdata/ligand_list.rds")))
ll <- ll[!is.na(sequence) & nzchar(sequence)]
ll[, nm := fifelse(!is.na(final_name) & nzchar(final_name), final_name, uniprot_name)]
k1 <- unique(ll[, .(accession, start, end, name = nm, seq = sequence, type = "ligand_list")])

up <- fread(file.path(A2, "pepend_locus", "uniprot_peptides.csv"))[type %in% c("peptide", "propeptide")]
up[prec, on = "accession", pseq := i.pseq]
up <- up[!is.na(pseq) & start >= 1 & end <= nchar(pseq)]
k2 <- up[, .(accession, start, end,
             name = fifelse(!is.na(description) & nzchar(description), description,
                            paste0(type, " ", start, "-", end)),
             seq = substr(pseq, start, end), type)]

kn <- rbind(k1, k2)
kn[prec, on = "accession", gene := i.gene]
kn[is.na(gene) | !nzchar(gene), gene := accession]      # fall back to the accession
kn <- kn[nchar(seq) >= 5]
## one row per sequence; a curated ligand wins over a UniProt feature for the
## same sequence, and a peptide over a propeptide
kn[, ord := match(type, c("ligand_list", "peptide", "propeptide"))]
setorder(kn, seq, ord, accession, start)
kn <- unique(kn, by = "seq")
kn[, ord := NULL]
kn[, known_id := fifelse(is.na(start) | is.na(end), gene, sprintf("%sx%dx%d", gene, start, end))]
kn[, tid := sprintf("t%06d", .I)]

writeLines(paste0(">", kn$tid, "\n", kn$seq), file.path(O, "known.fasta"))
fwrite(kn[, .(tid, known_id, name, gene, accession, start, end, type, seq)], file.path(O, "known.csv"))

cd <- fread(file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_T%s_tol%d.csv", HEAD, thr, TOL)))
u <- unique(cd[, .(sequence)]); u[, qid := sprintf("q%06d", .I)]
writeLines(paste0(">", u$qid, "\n", u$sequence), file.path(O, "cand.fasta"))
fwrite(u, file.path(O, "cand_ids.csv"))

cat(sprintf("targets: %s unique sequences (%s)\n", format(nrow(kn), big.mark = ","),
            paste(sprintf("%s %d", kn[, .N, by = type]$type, kn[, .N, by = type]$N), collapse = ", ")))
cat(sprintf("queries: %s candidates -> %s unique sequences\n",
            format(nrow(cd), big.mark = ","), format(nrow(u), big.mark = ",")))
