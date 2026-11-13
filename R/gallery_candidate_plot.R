## An index of a gallery's candidate peptides as one interactive scatter
## (Kevin, 2026-11-06, from his JMP view of the stage-6 table): end_score_max
## against mp_best_score, a panel per curated_relation, every point opening its
## candidate on its gene's page. Written as a single self-contained html from
## inst/templates/candidate_plot.html, which holds the drawing.

#' The best end of each precursor that has no candidate peptide
#'
#' One row per accession, in the columns of the scored table, so it can sit under
#' it: `start` and `end` are the best N and best C anchors, which need not be in
#' order, and `end_score_max` the higher of the two. A candidate's end score is
#' the highest of the arms at that anchor, so these are too. There is no peptide,
#' so no MultiPep score: `mp_best_score` holds where the better end lies along the
#' precursor instead, 0 to 1, and the plot labels its panel as such.
#' @param accs accessions.
#' @param prec the residue cache's precursor table.
#' @param key the candidate set, for its end segments.
#' @keywords internal
lf_best_ends <- function(accs, prec, key, a2) {
  arms <- c(curated = "", esmonly36 = "_esmonly36", curated_AAswap = "_curated_AAswap")
  ends <- lapply(c(N = "N", C = "C"), function(tm) {
    w <- Reduce(function(a, b) merge(a, b, by = c("accession", "anchor"), all = TRUE),
                lapply(names(arms), function(a) {
      b <- data.table::as.data.table(readRDS(file.path(a2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, arms[[a]]))))
      b <- b[b$accession %in% accs, c("accession", "anchor", "xgb_pm__all50"), with = FALSE]
      data.table::setnames(b, "xgb_pm__all50", a)
      b
    }))
    w$best <- do.call(pmax, c(w[, names(arms), with = FALSE], na.rm = TRUE))
    data.table::setorderv(w, c("accession", "best", "anchor"), c(1L, -1L, 1L))
    w[!duplicated(w$accession), ]
  })
  N <- ends$N; C <- ends$C[match(N$accession, ends$C$accession), ]
  i <- match(N$accession, prec$accession)
  seg <- data.table::fread(file.path(a2, "lf_headcmp", sprintf("end_segments_%s.csv", sub("_tol[0-9]+", "", key))),
                           select = c("accession", "term"))
  motif <- function(tm, an) vapply(seq_along(i), function(k)
    lf_pepend_motif_vec(prec$seq[i[k]], tm, an[k], prec$n_prot[i[k]]), "")
  L <- nchar(prec$seq[i])
  data.table::data.table(gene = prec$gene[i], accession = N$accession, start = N$anchor, end = C$anchor,
    curated_relation = "no candidate peptide", n_score_end = N$best, c_score_end = C$best,
    end_score_max = pmax(N$best, C$best), n_motif_end = motif("N", N$anchor), c_motif_end = motif("C", C$anchor),
    n_score_curated = N$curated, n_score_esmonly36 = N$esmonly36, n_score_curated_AAswap = N$curated_AAswap,
    c_score_curated = C$curated, c_score_esmonly36 = C$esmonly36, c_score_curated_AAswap = C$curated_AAswap,
    mp_best_score = ifelse(N$best >= C$best, N$anchor, C$anchor) / L,
    seg_n = tabulate(match(seg$accession[seg$term == "N"], N$accession), nrow(N)),
    seg_c = tabulate(match(seg$accession[seg$term == "C"], N$accession), nrow(N)))
}

#' Build the candidate scatter of a gallery
#'
#' @param dir the gallery directory of `<GENE>.html` pages.
#' @param template path to candidate_plot.html.
#' @param key the candidate set.
#' @param thr the end-score cut of the "end score type" colouring.
#' @param out the file to write: inside `dir`, unless an absolute path.
#' @param a2 the analysis directory.
#' @param base put before each link, for a page that is not beside the gene pages.
#' @param known an index of the known peptides to carry inside this file, opened
#'   at `#known`. Inside, it sits behind whatever gate this file does; a file of
#'   its own beside the public pages would not.
#' @return invisibly, the rows drawn.
#' @keywords internal
lf_gallery_candidate_plot <- function(dir, template, key = "xgb_pm_T0.5_tol4_arms3", thr = 0.6,
                                      out = "index_candidates.html", a2 = "~/AF2_analysis",
                                      base = "", known = NULL) {
  dir <- path.expand(dir); a2 <- path.expand(a2)
  cols <- c("gene", "accession", "start", "end", "sequence", "curated_relation", "curated_match",
            "uniprot_match", "n_score_end", "c_score_end", "end_score_max", "n_motif_end", "c_motif_end",
            "n_score_curated", "n_score_esmonly36", "n_score_curated_AAswap",
            "c_score_curated", "c_score_esmonly36", "c_score_curated_AAswap",
            "mp_best_score", "mp_best_class", "xgb_score", "xgb_nomp_score",
            "mmseqs_best_hit_name", "mmseqs_sim")
  d <- data.table::fread(file.path(a2, "lf_headcmp", sprintf("candidate_peptides_%s_scored.csv", key)),
                         select = cols)
  n_all <- nrow(d)
  pages <- sub("\\.html$", "", list.files(dir, "\\.html$"))
  ## a page draws one accession; where a name has two, scan_genes.txt names it
  pin <- grep("@", readLines(file.path(dir, "scan_genes.txt")), value = TRUE, fixed = TRUE)
  pin <- stats::setNames(sub(".*@", "", pin), sub("@.*", "", pin))
  keep <- d$gene %in% pages & (!d$gene %in% names(pin) | d$accession == pin[d$gene])
  d <- d[keep, ]

  prec <- readRDS(file.path(a2, "lf_pepend_residue_cache.rds"))$prec
  ## every page with no candidate joins as its best end, a row of its own
  page_g <- setdiff(grep("^index(_|$)", pages, value = TRUE, invert = TRUE), d$gene)
  page_a <- prec$accession[match(page_g, prec$gene)]
  hit <- page_g %in% names(pin); page_a[hit] <- pin[page_g[hit]]
  ## a page's gene name can be two genes' in the table ("DEFB4A; DEFB4B")
  miss <- is.na(page_a)
  if (any(miss)) page_a[miss] <- prec$accession[match(page_g[miss], sub(";.*$", "", prec$gene))]
  stopifnot(!anyNA(page_a))
  ncand <- nrow(d); ngc <- length(unique(d$gene))
  x <- lf_best_ends(page_a, prec, key, a2)
  x$gene <- page_g[match(x$accession, page_a)]
  stopifnot(nrow(x) == length(page_a))
  cand_seq <- d[, c("gene", "accession", "start", "end", "sequence")]
  d <- data.table::rbindlist(list(d, x), fill = TRUE)
  genes <- sort(unique(d$gene))
  acc <- d$accession[match(genes, d$gene)]
  ## the sequence is only for spelling out candidates
  seqs <- prec$seq[match(acc, prec$accession)]
  seqs[!genes %in% cand_seq$gene] <- ""
  g <- match(d$gene, genes)
  stopifnot(identical(substr(seqs[match(cand_seq$gene, genes)], cand_seq$start, cand_seq$end), cand_seq$sequence))

  ## fread reads an empty field as "", not NA: both mean "none" here
  dict <- function(x) { x[is.na(x)] <- ""; u <- sort(unique(c("", x))); list(u = u, i = match(x, u) - 1L) }
  ints <- function(x, k = 1000) { v <- as.integer(round(x * k)); v[is.na(v)] <- -1L
                                  paste0("[", paste(v, collapse = ","), "]") }
  js <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = FALSE))
  rel <- c("exact", "none", "overlap", "uniprot", "no candidate peptide")
  stopifnot(all(d$curated_relation %in% rel))
  motif <- c("dibasic", "G + dibasic", "monobasic", "terminus", "other")
  stopifnot(all(d$c_motif_end %in% motif), all(d$n_motif_end %in% motif))
  mpc <- dict(gsub("peptidehormone", "peptide hormone", d$mp_best_class, fixed = TRUE))
  cur <- dict(d$curated_match); uni <- dict(d$uniprot_match); hits <- dict(d$mmseqs_best_hit_name)
  ## the end score type is cut here, on the table's own values: the page's are
  ## rounded for size, which would move 0.6004 onto the wrong side of 0.6
  et <- ifelse(d$n_score_end > thr, ifelse(d$c_score_end > thr, 0L, 1L), ifelse(d$c_score_end > thr, 2L, 3L))
  f <- list(x = ints(d$mp_best_score), y = ints(d$end_score_max), n = ints(d$n_score_end), c = ints(d$c_score_end),
            et = ints(et, 1),
            nc = ints(d$n_score_curated), ne = ints(d$n_score_esmonly36), nw = ints(d$n_score_curated_AAswap),
            cc = ints(d$c_score_curated), ce = ints(d$c_score_esmonly36), cw = ints(d$c_score_curated_AAswap),
            xs = ints(d$xgb_score), xn = ints(d$xgb_nomp_score), si = ints(d$mmseqs_sim, 100),
            s = ints(d$start, 1), e = ints(d$end, 1), g = ints(g - 1L, 1),
            r = ints(match(d$curated_relation, rel) - 1L, 1),
            nm = ints(match(d$n_motif_end, motif) - 1L, 1), cm = ints(match(d$c_motif_end, motif) - 1L, 1),
            mc = ints(mpc$i, 1), cu = ints(cur$i, 1), un = ints(uni$i, 1), hi = ints(hits$i, 1),
            genes = js(genes), seqs = js(seqs), rel = js(rel), motif = js(motif),
            mpc = js(mpc$u), cur = js(cur$u), uni = js(uni$u), hits = js(hits$u),
            sn = ints(x$seg_n, 1), sc = ints(x$seg_c, 1), ncand = ncand, ngc = ngc,
            key = js(key), thr = thr, dropped = n_all - ncand,
            base = as.character(jsonlite::toJSON(base, auto_unbox = TRUE)))
  data <- paste0("{", paste(sprintf('"%s":%s', names(f), unlist(lapply(f, as.character))), collapse = ","), "}")
  data <- gsub("</", "<\\/", data, fixed = TRUE)
  tpl <- paste(readLines(template, warn = FALSE), collapse = "\n")
  kn <- if (is.null(known)) "null" else gsub("</", "<\\/", as.character(jsonlite::toJSON(
    paste(readLines(path.expand(known), warn = FALSE, encoding = "UTF-8"), collapse = "\n"), auto_unbox = TRUE)), fixed = TRUE)
  tpl <- sub("__KNOWN__", kn, tpl, fixed = TRUE)
  tpl <- strsplit(tpl, "__DATA__", fixed = TRUE)[[1]]
  stopifnot(length(tpl) == 2L)
  out <- path.expand(out)
  p <- if (startsWith(out, "/")) out else file.path(dir, out)
  con <- file(p, "w", encoding = "UTF-8")
  writeLines(c(tpl[1], data, tpl[2]), con, sep = "")
  close(con)
  message(sprintf("%s candidates on %d genes, %s genes as their best end -> %s (%.1f MB); %s left out, no page",
                  format(ncand, big.mark = ","), ngc, format(nrow(x), big.mark = ","), p, file.size(p) / 1e6,
                  format(n_all - ncand, big.mark = ",")))
  invisible(d)
}
