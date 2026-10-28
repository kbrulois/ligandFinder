## Pick candidates by a filter on the scored table, and the precursors to draw
## pages for (Kevin, 2026-10-28). candidate_pages.sh renders what this writes.
##
##   FILTER='curated_relation == "none" & mp_best_score > 0.6 & xgb_nomp_score > 0.45' \
##   OUT=~/AF2_analysis/pepend_pages_cand_nocurated_mp0.6_nomp0.45 \
##   HEAD=xgb_pm THRESH=0.5 TOL=4 ARM=_arms3 Rscript candidate_select.R
##
## FILTER is an R expression over the columns of candidate_peptides_<key>_scored.csv
## (stage 6 writes mp_best_class / mp_best_score there; rebuilt the same way if
## absent). FOCUS picks each gene's card anchor -- its best PASSING candidate by
## that column, ties to mp_best_score, then the id. Default: the first of
## xgb_nomp_score, xgb_score, mp_best_score, end_score_max, mmseqs_sim that
## FILTER names (in that order of preference, whatever order FILTER uses), else
## xgb_score.
##
## NOTE curated_relation compares a candidate with the CURATED knowns of 8-60 aa
## (ligand_list) only. "none" is not "unannotated": a UniProt peptide (label 1,
## e.g. gonadoliberin, orexin-A) can pass it. Add `& !label %in% 1L` to FILTER to
## drop those; the card note names any it keeps.
##
## Writes into OUT:
##   selected_candidates.csv  every passing candidate, all the scored table's columns,
##                            by gene (best FOCUS first) then FOCUS
##   selected_genes.csv       one row per precursor: passing candidates and groups,
##                            the anchor peptide and its scores, length, renderable?
##   genes.txt                the renderable genes, shortest precursor first
##   scan_genes.txt           the same as GENE@ACCESSION, for arm_gene_scans.sh (the
##                            accession each page draws: COL18A1 is P39060)
##   focus.csv                gene, peptide_id, rank, note, and the anchor's scores --
##                            the index's card anchors (candidate_gallery_index.R
##                            --focus), best FOCUS first
##   selection.txt            the filter, key and counts
##
## RENDERABLE = in plot_bundle.rds (final_pages.R stops on a gene it cannot
## find there) or in the extra bundle LF_EXTRA_BUNDLE names (built by
## inst/scripts/plot_bundle_extra.R for genes the old window model never
## scored), once, and with no "/" in the name (it becomes a filename).
suppressMessages({library(data.table)})
A2   <- path.expand("~/AF2_analysis")
HEAD <- Sys.getenv("HEAD", "xgb_pm"); ARM <- Sys.getenv("ARM", "_arms3")
thr  <- as.numeric(Sys.getenv("THRESH", "0.5")); TOL <- as.integer(Sys.getenv("TOL", "4"))
key  <- sprintf("%s_T%s_tol%d%s", HEAD, thr, TOL, ARM)
FILTER <- Sys.getenv("FILTER", "")
OUT    <- path.expand(Sys.getenv("OUT", ""))
if (!nzchar(FILTER) || !nzchar(OUT)) stop("set FILTER (an R expression) and OUT (a directory)")
f_sc <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_scored.csv", key))
if (!file.exists(f_sc)) stop("no scored table ", f_sc, " -- run candidate_pipeline.sh (stage 6)")

sc <- fread(f_sc)
## stage 6 writes mp_best_* into the scored table; rebuild them if absent
mp <- setdiff(grep("^mp_", names(sc), value = TRUE), c("mp_signalling_max", "mp_best_class", "mp_best_score"))
if (!"mp_best_score" %in% names(sc) && length(mp)) {
  M <- as.matrix(sc[, mp, with = FALSE]); k <- max.col(replace(M, is.na(M), -Inf), ties.method = "first")
  sc[, `:=`(mp_best_class = sub("^mp_", "", mp[k]), mp_best_score = M[cbind(seq_len(nrow(M)), k)])]
}
ex <- str2lang(FILTER)
miss <- setdiff(all.vars(ex), names(sc))
if (length(miss)) stop("FILTER names columns the scored table lacks: ", paste(miss, collapse = ", "))
keep <- sc[, eval(ex)]
if (!is.logical(keep) || length(keep) != nrow(sc)) stop("FILTER must give one TRUE/FALSE per candidate")
s <- sc[keep %in% TRUE]
if (!nrow(s)) stop("no candidate passes ", FILTER)

FOCUS <- Sys.getenv("FOCUS", "")
if (!nzchar(FOCUS)) {
  ## in this order of preference, whatever order FILTER names them in
  fv <- intersect(c("xgb_nomp_score", "xgb_score", "mp_best_score", "end_score_max", "mmseqs_sim"), all.vars(ex))
  FOCUS <- if (length(fv)) fv[1] else "xgb_score"
}
if (!FOCUS %in% names(s)) stop("FOCUS column not in the table: ", FOCUS)
## ---- which genes have a page, and the accession that page draws -----------------------
## plot_bundle.rds holds one entry per drawable gene (LF_EXTRA_BUNDLE adds more),
## and the page draws THAT entry's accession: under COL18A1 the residue cache
## also holds A0AAG2UXZ5, which no page draws
.bget <- function(l, f) vapply(l, function(x) as.character(x[[f]])[1], "")
b <- readRDS(file.path(A2, "plot_bundle.rds"))$the_input
bundle <- data.table(gene = .bget(b, "gene"), page_accession = .bget(b, "accession")); rm(b); invisible(gc())
xb <- Sys.getenv("LF_EXTRA_BUNDLE", "")
if (nzchar(xb) && !file.exists(path.expand(xb))) stop("LF_EXTRA_BUNDLE: no file ", xb)
extra_genes <- character(0)
if (nzchar(xb)) {
  x <- readRDS(path.expand(xb))$the_input
  extra_genes <- .bget(x, "gene")
  bundle <- rbind(bundle, data.table(gene = extra_genes, page_accession = .bget(x, "accession")))
}
dup <- unique(bundle$gene[duplicated(bundle$gene)])
s[bundle[!gene %in% dup], on = "gene", page_accession := i.page_accession]
## the anchor comes from the page's own accession whenever the gene passes there
s[, on_page := is.na(page_accession) | accession == page_accession]
setorderv(s, c("on_page", FOCUS, "mp_best_score", "peptide_id"), c(-1L, -1L, -1L, 1L), na.last = TRUE)

## ---- per gene: the anchor (best passing by FOCUS), counts, renderable ----------------
g <- s[, .(accession = accession[1], page_accession = page_accession[1], n_pass = .N,
           n_off_page = sum(!on_page), n_groups = uniqueN(group_id),
           focus_peptide = peptide_id[1], focus_length = length[1],
           focus_score = get(FOCUS)[1], focus_xgb = xgb_score[1], focus_mp = mp_best_score[1],
           focus_mp_class = mp_best_class[1], focus_mmseqs_hit = mmseqs_best_hit_name[1],
           focus_mmseqs_sim = mmseqs_sim[1],
           n_label1 = sum(label %in% 1L)), by = gene]
setorder(g, -focus_score, gene)
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
g[, precursor_length := nchar(fc$prec$seq)[match(fifelse(is.na(page_accession), accession, page_accession),
                                                 fc$prec$accession)]]
rm(fc)
g[, renderable := gene %in% bundle$gene & !gene %in% dup & !grepl("/", gene, fixed = TRUE)]
g[, why_not := fifelse(renderable, "",
                fifelse(!gene %in% bundle$gene, "not in plot_bundle.rds; see plot_bundle_extra.R",
                fifelse(grepl("/", gene, fixed = TRUE), "gene name contains /", "gene twice in plot_bundle.rds")))]
if (g[renderable == TRUE & n_off_page > 0, .N])
  message("passing candidates on an accession no page draws (kept in the csv, never an anchor): ",
          paste(g[renderable == TRUE & n_off_page > 0, sprintf("%s (%d)", gene, n_off_page)], collapse = ", "))
s[, c("page_accession", "on_page") := NULL]

## ---- write ---------------------------------------------------------------------------
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
s[, gene_rank := match(gene, g$gene)]
setorderv(s, c("gene_rank", FOCUS, "mp_best_score", "peptide_id"), c(1L, -1L, -1L, 1L))
s[, gene_rank := NULL]
fwrite(s, file.path(OUT, "selected_candidates.csv"))
fwrite(g, file.path(OUT, "selected_genes.csv"))
writeLines(g[renderable == TRUE][order(precursor_length, gene)]$gene, file.path(OUT, "genes.txt"))
writeLines(g[renderable == TRUE][order(precursor_length, gene), sprintf("%s@%s", gene, page_accession)],
           file.path(OUT, "scan_genes.txt"))
lab <- function(x) gsub("_", " ", sub("cellcellsignaling", "cell-cell signalling",
                                      sub("peptidehormone", "peptide hormone", x)))
fo <- g[, .(gene, peptide_id = focus_peptide, rank = seq_len(.N),
            note = sprintf("%d passing in %d group%s; best passing by %s: %s, %d aa, %.3f (MultiPep %.2f %s)",
                           n_pass, n_groups, fifelse(n_groups == 1L, "", "s"), FOCUS,
                           focus_peptide, focus_length, focus_score, focus_mp, lab(focus_mp_class)))]
## the anchor's own scores, so the index can list it when stage 7 did not pick it
## (stage 7 picks over ALL of a group's cuts, passing or not)
fcols <- intersect(c("group_id", "length", "xgb_score", "xgb_nomp_score", "end_score_max", "mp_best_score",
                     "mp_best_class", "mmseqs_sim", "mmseqs_best_hit_name", "label", "label_match",
                     "label_type", "curated_relation", "curated_match"), names(s))
fo <- cbind(fo, s[match(fo$peptide_id, s$peptide_id), fcols, with = FALSE])
## a known peptide at the anchor's range (label 1 = UniProt peptide or curated
## ligand; NA = propeptide): say so, since curated_relation "none" allows it
kn <- !is.na(fo$label_type) & nzchar(fo$label_type)
kind <- c(peptide = "UniProt peptide", propeptide = "UniProt propeptide", ligand_list = "curated ligand")
fo[kn, note := sprintf("%s -- matches %s%s", note, kind[label_type],
                       fifelse(label_type == "propeptide" | !nzchar(label_match), "", paste0(" ", label_match)))]
## and every known peptide the filter keeps for the gene, with how many cuts
km <- s[!is.na(label_type) & nzchar(label_type),
        .(k = fifelse(label_type == "propeptide" | !nzchar(label_match), unname(kind[label_type]), label_match)),
        by = gene][, .N, by = .(gene, k)][order(gene, -N)][
        , .(kept = paste(sprintf("%s (%d)", k, N), collapse = ", ")), by = gene]
fo[km, on = "gene", note := sprintf("%s. Known peptides among the passing cuts: %s", note, i.kept)]
fwrite(fo, file.path(OUT, "focus.csv"))
txt <- c(sprintf("key       %s", key),
         sprintf("table     %s", f_sc),
         sprintf("filter    %s", FILTER),
         sprintf("focus     %s (each gene's card opens at its best passing candidate by this)", FOCUS),
         sprintf("passing   %s candidates on %d precursors", format(nrow(s), big.mark = ","), nrow(g)),
         sprintf("pages     %d renderable", g[renderable == TRUE, .N]),
         if (length(extra_genes))
           sprintf("extra     %s (LF_EXTRA_BUNDLE): %s", xb, paste(intersect(extra_genes, g$gene), collapse = ", ")),
         if (g[renderable == FALSE, .N])
           sprintf("no page   %s", paste(g[renderable == FALSE, sprintf("%s (%s)", gene, why_not)], collapse = ", ")),
         sprintf("written   %s", format(Sys.time(), "%F %T")))
writeLines(txt, file.path(OUT, "selection.txt"))
cat(txt, sep = "\n")
