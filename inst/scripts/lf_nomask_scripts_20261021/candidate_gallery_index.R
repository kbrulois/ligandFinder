## The front page for a gallery of candidate pages (Kevin, 2026-10-27): every
## <GENE>.html in --dir gets a card that opens its page AT the gene's
## top-scoring group, with that group's picks (stage 7) listed on the card.
##
##   Rscript candidate_gallery_index.R --dir ~/AF2_analysis/pepend_pages_candidates_arms3 \
##     --key xgb_pm_T0.5_tol4_arms3 [--title "..."]
##
## The top group is the gene's group holding its highest xgb_score; the card
## links to <GENE>.html#<that group's xgb pick>, which the page opens with the N
## window at the pick's first residue, the C window at its last, centred
## between them. The card itself is one link, so the picks are listed as text:
## each distinct cut, the criteria it won (xgb / nomp / end / multipep / mmseqs) and the
## values that won them. The cards are built here and handed to
## make_gallery_index.R (--card-html, --card-href), which draws the page.
##
## --focus <csv> (gene, peptide_id[, rank, note, scores]; candidate_select.R
## writes it) gives a listed gene its card from THAT peptide instead: the link
## opens the page at it, the card leads with its note, then lists the picks of
## the peptide's own group -- with the peptide itself first, as "selected", when
## stage 7 picked other cuts on every criterion (its picks range over all of a
## group's cuts, passing the selection or not) -- and cards go in the csv's
## rank order. Only the csv's genes get a card. --selection "<text>" says in
## the blurb how the genes were chosen.
suppressMessages({library(data.table)})
.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) { i <- match(flag, .args); if (is.na(i)) default else .args[i + 1] }
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
A2  <- path.expand("~/AF2_analysis")
dir <- path.expand(.opt("--dir", file.path(A2, "pepend_pages_candidates_arms3")))
key <- .opt("--key", "xgb_pm_T0.5_tol4_arms3")
ttl <- .opt("--title", "LigandFinder — candidate peptides")
foc_p <- .opt("--focus", NA_character_)
sel_t <- .opt("--selection", "")
f_bg  <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_best_per_group.csv", key))
f_cd  <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s.csv", key))
f_seg <- file.path(A2, "lf_headcmp", sprintf("end_segments_%s.csv", sub("_tol[0-9]+", "", key)))
for (f in c(f_bg, f_cd)) if (!file.exists(f)) stop("missing ", f, " -- run candidate_pipeline.sh (stages 6-7) first")

genes <- sub("\\.html$", "", setdiff(list.files(dir, "\\.html$"), "index.html"))
if (!length(genes)) stop("no pages in ", dir)
bg  <- fread(f_bg)
cdg <- fread(f_cd, select = c("gene", "accession", "group_id", "peptide_id"))
## per accession: the residue cache holds two under COL18A1, and a page's
## candidate track shows its own accession's only
cnt <- cdg[, .(n_cand = .N, n_grp = uniqueN(group_id)), by = .(gene, accession)]
foc <- if (is.na(foc_p)) data.table(gene = character(), peptide_id = character()) else fread(path.expand(foc_p))
if (nrow(foc)) {
  if (!"note" %in% names(foc)) foc[, note := ""]
  if (!"rank" %in% names(foc)) foc[, rank := seq_len(.N)]
  ## overwritten from THIS set's table (focus.csv carries its own group_id), so
  ## a focus file from another candidate set is refused rather than drawn
  i <- match(foc$peptide_id, cdg$peptide_id)
  foc[, `:=`(group_id = cdg$group_id[i], accession = cdg$accession[i])]
  if (anyNA(i)) stop("--focus peptides not in ", basename(f_cd), ": ",
                     paste(foc$peptide_id[is.na(i)], collapse = ", "))
}
## focus mode: the csv's genes only -- a page left from an earlier selection in
## the same directory is not one of this one's, and says so
nopage <- nopage_why <- character(0)
if (nrow(foc)) {
  stale <- setdiff(genes, foc$gene)
  if (length(stale)) message(length(stale), " page(s) in ", dir, " are not in --focus, left off the index: ",
                             paste(stale, collapse = ", "))
  nopage <- setdiff(foc$gene, genes)
  ## why: not drawable at all (candidate_select.R says so), or the render failed
  sg <- file.path(dir, "selected_genes.csv")
  why <- if (file.exists(sg)) fread(sg, select = c("gene", "renderable", "why_not")) else NULL
  nopage_why <- if (is.null(why)) rep("no page", length(nopage)) else
    ifelse(nopage %in% why[renderable == FALSE, gene],
           paste0("no page: ", why$why_not[match(nopage, why$gene)]), "render failed, see .render/logs")
  genes <- intersect(genes, foc$gene)
  if (!length(genes)) stop("no page in ", dir, " for any --focus gene")
}
seg <- if (file.exists(f_seg)) fread(f_seg, select = "gene")[, .(n_seg = .N), by = gene] else data.table(gene = character(), n_seg = integer())

esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE)
                     x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
mp_lab <- function(x) gsub("_", " ", sub("cellcellsignaling", "cell-cell signalling",
                                         sub("peptidehormone", "peptide hormone", x)))
## what a pick won, and the value that won it
## (a "selected" row -- the focus peptide when stage 7 did not pick it -- shows
## every score it has)
won <- function(r) {
  txt <- c(xgb = sprintf("xgb %.3f", r$xgb_score),
           nomp = sprintf("xgb no-MultiPep %.3f", r$xgb_nomp_score),
           end = sprintf("end %.2f", r$end_score_max),
           multipep = sprintf("MultiPep %.2f %s", r$mp_best_score, mp_lab(r$mp_best_class)),
           mmseqs = sprintf("MMseqs %.2f %s", r$mmseqs_sim, r$mmseqs_best_hit_name))
  has <- !is.na(c(r$xgb_score, r$xgb_nomp_score, r$end_score_max, r$mp_best_score, r$mmseqs_sim))
  cr <- if (identical(r$picked_by, "selected")) names(txt)[has] else strsplit(r$picked_by, "+", fixed = TRUE)[[1]]
  paste(esc(txt[cr]), collapse = " &middot; ")
}

cards <- rbindlist(lapply(genes, function(g) {
  r <- bg[gene == g]
  if (!nrow(r)) {
    ns <- seg[gene == g, n_seg]
    return(data.table(gene = g, accession = NA_character_, href = paste0(g, ".html"), html = sprintf(
      '<div class="tg none-c">no candidates in this set%s</div>',
      if (length(ns) && ns > 0) sprintf(" &mdash; %d end segment%s, no N&ndash;C pair 8&ndash;60 aa",
                                        ns, if (ns == 1) "" else "s") else "")))
  }
  f <- foc[gene == g]
  if (nrow(f)) {                             # the gene's selected peptide and its group
    top <- f$group_id[1]; pick <- f$peptide_id[1]
    tr  <- r[group_id == top]
    if (!pick %in% tr$peptide_id && "xgb_score" %in% names(f))
      tr <- rbind(f[1, intersect(names(f), names(tr)), with = FALSE][, picked_by := "selected"], tr, fill = TRUE)
  } else {
    top <- r$group_id[1]                     # stage 7 orders groups by best xgb_score
    tr  <- r[group_id == top]
    pick <- tr[grepl("xgb", picked_by, fixed = TRUE)]$peptide_id[1]   # not "nomp": criterion names never nest
  }
  acc <- if (nrow(f)) f$accession[1] else tr$accession[1]
  cn  <- cnt[gene == g & accession == acc]
  known <- tr[label %in% 1L | curated_relation == "exact"]
  rows <- paste0('<tr', ifelse(nrow(f) & tr$peptide_id == pick, ' class="sel"', ''), '>',
                 '<td class="c">', esc(gsub("+", " + ", tr$picked_by, fixed = TRUE)), '</td>',
                 '<td class="p">', esc(tr$peptide_id), ' <span class="ln">', tr$length, ' aa</span>',
                 '<div class="w">', vapply(seq_len(nrow(tr)), function(i) won(tr[i]), ""), '</div></td></tr>',
                 collapse = "")
  data.table(gene = g, accession = acc, href = sprintf("%s.html#%s", g, pick), html = paste0(
    if (nrow(f) && nzchar(f$note[1])) sprintf('<div class="fs">%s</div>', esc(f$note[1])) else "",
    '<div class="tg">', if (nrow(f)) "its group" else "top group", ' <b>', esc(top), '</b> &middot; ',
    cn$n_grp, ' group', if (cn$n_grp == 1) "" else "s",
    ', ', format(cn$n_cand, big.mark = ","), ' candidates</div>',
    '<table class="pk">', rows, '</table>',
    {
      ## fread reads an empty field as "", not NA
      kn <- unique(c(known$label_match, known$curated_match))
      kn <- kn[!is.na(kn) & nzchar(kn)]
      if (length(kn)) sprintf('<div class="kn">known: %s</div>', esc(paste(kn, collapse = ", "))) else ""
    }))
}))

tmp <- tempfile("candidx_"); dir.create(tmp)
saveRDS(cards[, .(gene, html)], file.path(tmp, "cards.rds"))
saveRDS(cards[, .(gene, href)], file.path(tmp, "href.rds"))
if (nrow(foc)) saveRDS(foc[, .(gene, rank)], file.path(tmp, "order.rds"))
saveRDS(cards[!is.na(accession), .(gene, accession)], file.path(tmp, "acc.rds"))
writeLines(c(
  '<style>',
  ' .grid{grid-template-columns:repeat(auto-fill,minmax(340px,1fr)) !important}',
  ' .tg{font-size:11.5px;margin-top:6px;color:#0b0b0b} .tg.none-c{color:#52514e;font-style:italic}',
  ' table.pk{width:100%;border-collapse:collapse;margin-top:5px;font-size:11px}',
  ' table.pk td{padding:3px 4px;border-top:1px solid #efeee9;vertical-align:top}',
  ' table.pk td.c{color:#52514e;white-space:nowrap;width:1%} table.pk td.p{font-family:ui-monospace,Menlo,monospace}',
  ' .ln{color:#898781;font-family:ui-sans-serif,system-ui,sans-serif} .w{color:#52514e;font-family:ui-sans-serif,system-ui,sans-serif;font-size:10.5px}',
  ' .kn{font-size:11px;margin-top:5px;color:#1B7837;font-weight:600}',
  ' .fs{font-size:11.5px;margin-top:6px;color:#0b0b0b;font-weight:600} table.pk tr.sel td{background:#f6f5f1}',
  '</style>',
  'Each page carries the class-probability tracks, the window-free residue models, the',
  '<b>CNN linear head</b> and <b>xgb window head</b> for three arms (<b>curated</b>, <b>ESM-C only</b>,',
  '<b>AAswap</b>, curated with the AA properties &rarr; one-hot; thin lines per fold, heavy line all 50), the <b>end segments</b>',
  '(each arm at its own threshold; the digit is how many arms support a segment) and the <b>candidate',
  'peptides</b>: left half coloured by the best MultiPep score, right half by MMseqs similarity to the',
  'closest known 8&ndash;60 aa peptide, the symbol is the MultiPep class, end caps white / grey / black for',
  '1 / 2 / 3 supporting arms, and the tooltip carries the xgboost global scores.',
  if (nrow(foc)) c(
  if (nzchar(sel_t)) sprintf('<br><b>Genes</b>: every precursor with a candidate passing <code>%s</code>.', esc(sel_t)),
  if (length(nopage)) sprintf('%d passing precursor%s %s no card: %s.', length(nopage),
                              if (length(nopage) == 1) "" else "s", if (length(nopage) == 1) "has" else "have",
                              esc(paste(sprintf("%s (%s)", nopage, nopage_why), collapse = "; "))),
  '<br><b>A card opens its page at that gene&rsquo;s selected candidate</b> (its best passing cut, named on',
  'the card and shaded in its list): the N window at its first residue, the C window at its last. Listed',
  'with it: its group&rsquo;s best cut by xgb score, xgb score without MultiPep, end score, MultiPep and',
  'MMseqs, one line per distinct cut, chosen over all the group&rsquo;s cuts whether they pass or not',
  '&mdash; so when those are all other cuts, the selected one is listed first as <i>selected</i>.',
  'Cards are in order of the selected candidate&rsquo;s score.')
  else c(
  '<br><b>A card opens its page at the gene&rsquo;s top-scoring group</b> (the group holding its highest',
  'xgb score): the N window at the xgb pick&rsquo;s first residue, the C window at its last. Listed on the',
  'card: that group&rsquo;s best cut by xgb score, xgb score without MultiPep, end score, MultiPep and MMseqs (one line per distinct cut).')),
  file.path(tmp, "blurb.html"))
commit <- tryCatch(system2("git", c("-C", lf_src_root(), "rev-parse", "--short", "HEAD"), stdout = TRUE),
                   error = function(e) "?")
## pages drawn from a checkout with uncommitted page code are not that commit
dirty <- tryCatch(length(system2("git", c("-C", lf_src_root(), "status", "--porcelain", "--untracked-files=no",
                                          "--", "R", "inst/scripts"),
                                 stdout = TRUE)) > 0, error = function(e) FALSE)
if (dirty) commit <- paste(commit, "plus uncommitted changes")
writeLines(sprintf('Candidate set <code>%s</code> (candidate_pipeline.sh stages 1&ndash;7); pages drawn by <code>final_pages.R</code> at commit <code>%s</code>.',
                   esc(key), esc(commit)), file.path(tmp, "footer.html"))

rc <- system2(file.path(R.home("bin"), "Rscript"),
              c(lf_src("inst/scripts/make_gallery_index.R"), "--dir", shQuote(dir), "--title", shQuote(ttl),
                "--card-html", file.path(tmp, "cards.rds"), "--card-href", file.path(tmp, "href.rds"),
                "--blurb", file.path(tmp, "blurb.html"), "--footer", file.path(tmp, "footer.html"),
                if (nrow(foc)) c("--card-order", file.path(tmp, "order.rds"), "--only-ranked"),
                "--card-accession", file.path(tmp, "acc.rds"),
                "--no-stock-desc"))
if (rc != 0) stop("make_gallery_index.R failed")
print(cards[, .(gene, href)], row.names = FALSE)
