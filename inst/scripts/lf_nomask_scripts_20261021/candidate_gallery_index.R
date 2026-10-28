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
## each distinct cut, the criteria it won (xgb / end / multipep / mmseqs) and the
## values that won them. The cards are built here and handed to
## make_gallery_index.R (--card-html, --card-href), which draws the page.
suppressMessages({library(data.table)})
.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) { i <- match(flag, .args); if (is.na(i)) default else .args[i + 1] }
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
A2  <- path.expand("~/AF2_analysis")
dir <- path.expand(.opt("--dir", file.path(A2, "pepend_pages_candidates_arms3")))
key <- .opt("--key", "xgb_pm_T0.5_tol4_arms3")
ttl <- .opt("--title", "LigandFinder — candidate peptides")
f_bg  <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s_best_per_group.csv", key))
f_cd  <- file.path(A2, "lf_headcmp", sprintf("candidate_peptides_%s.csv", key))
f_seg <- file.path(A2, "lf_headcmp", sprintf("end_segments_%s.csv", sub("_tol[0-9]+", "", key)))
for (f in c(f_bg, f_cd)) if (!file.exists(f)) stop("missing ", f, " -- run candidate_pipeline.sh (stages 6-7) first")

genes <- sub("\\.html$", "", setdiff(list.files(dir, "\\.html$"), "index.html"))
if (!length(genes)) stop("no pages in ", dir)
bg  <- fread(f_bg)
cnt <- fread(f_cd, select = c("gene", "group_id"))[, .(n_cand = .N, n_grp = uniqueN(group_id)), by = gene]
seg <- if (file.exists(f_seg)) fread(f_seg, select = "gene")[, .(n_seg = .N), by = gene] else data.table(gene = character(), n_seg = integer())

esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE)
                     x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
mp_lab <- function(x) gsub("_", " ", sub("cellcellsignaling", "cell-cell signalling",
                                         sub("peptidehormone", "peptide hormone", x)))
## what a pick won, and the value that won it
won <- function(r) {
  cr <- strsplit(r$picked_by, "+", fixed = TRUE)[[1]]
  txt <- c(xgb = sprintf("xgb %.3f", r$xgb_score),
           end = sprintf("end %.2f", r$end_score_max),
           multipep = sprintf("MultiPep %.2f %s", r$mp_best_score, mp_lab(r$mp_best_class)),
           mmseqs = sprintf("MMseqs %.2f %s", r$mmseqs_sim, r$mmseqs_best_hit_name))
  paste(esc(txt[cr]), collapse = " &middot; ")
}

cards <- rbindlist(lapply(genes, function(g) {
  r <- bg[gene == g]
  if (!nrow(r)) {
    ns <- seg[gene == g, n_seg]
    return(data.table(gene = g, href = paste0(g, ".html"), html = sprintf(
      '<div class="tg none-c">no candidates in this set%s</div>',
      if (length(ns) && ns > 0) sprintf(" &mdash; %d end segment%s, no N&ndash;C pair 8&ndash;60 aa",
                                        ns, if (ns == 1) "" else "s") else "")))
  }
  top <- r$group_id[1]                       # stage 7 orders groups by best xgb_score
  tr  <- r[group_id == top]
  pick <- tr[grepl("xgb", picked_by, fixed = TRUE)]$peptide_id[1]
  cn  <- cnt[gene == g]
  known <- tr[label %in% 1L | curated_relation == "exact"]
  rows <- paste0('<tr><td class="c">', esc(gsub("+", " + ", tr$picked_by, fixed = TRUE)), '</td>',
                 '<td class="p">', esc(tr$peptide_id), ' <span class="ln">', tr$length, ' aa</span>',
                 '<div class="w">', vapply(seq_len(nrow(tr)), function(i) won(tr[i]), ""), '</div></td></tr>',
                 collapse = "")
  data.table(gene = g, href = sprintf("%s.html#%s", g, pick), html = paste0(
    '<div class="tg">top group <b>', esc(top), '</b> &middot; ', cn$n_grp, ' group', if (cn$n_grp == 1) "" else "s",
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
writeLines(c(
  '<style>',
  ' .grid{grid-template-columns:repeat(auto-fill,minmax(340px,1fr)) !important}',
  ' .tg{font-size:11.5px;margin-top:6px;color:#0b0b0b} .tg.none-c{color:#52514e;font-style:italic}',
  ' table.pk{width:100%;border-collapse:collapse;margin-top:5px;font-size:11px}',
  ' table.pk td{padding:3px 4px;border-top:1px solid #efeee9;vertical-align:top}',
  ' table.pk td.c{color:#52514e;white-space:nowrap;width:1%} table.pk td.p{font-family:ui-monospace,Menlo,monospace}',
  ' .ln{color:#898781;font-family:ui-sans-serif,system-ui,sans-serif} .w{color:#52514e;font-family:ui-sans-serif,system-ui,sans-serif;font-size:10.5px}',
  ' .kn{font-size:11px;margin-top:5px;color:#1B7837;font-weight:600}',
  '</style>',
  'Each page carries the class-probability tracks, the window-free residue models, the',
  '<b>CNN linear head</b> and <b>xgb window head</b> for three arms (<b>curated</b>, <b>ESM-C only</b>,',
  '<b>AAswap</b>, curated with the AA properties &rarr; one-hot; thin lines per fold, heavy line all 50), the <b>end segments</b>',
  '(each arm at its own threshold; the digit is how many arms support a segment) and the <b>candidate',
  'peptides</b>: left half coloured by the best MultiPep score, right half by MMseqs similarity to the',
  'closest known 8&ndash;60 aa peptide, the symbol is the MultiPep class, end caps white / grey / black for',
  '1 / 2 / 3 supporting arms, and the tooltip carries the xgboost global score.',
  '<br><b>A card opens its page at the gene&rsquo;s top-scoring group</b> (the group holding its highest',
  'xgb score): the N window at the xgb pick&rsquo;s first residue, the C window at its last. Listed on the',
  'card: that group&rsquo;s best cut by xgb score, end score, MultiPep and MMseqs (one line per distinct cut).'),
  file.path(tmp, "blurb.html"))
commit <- tryCatch(system2("git", c("-C", lf_src_root(), "rev-parse", "--short", "HEAD"), stdout = TRUE),
                   error = function(e) "?")
writeLines(sprintf('Candidate set <code>%s</code> (candidate_pipeline.sh stages 1&ndash;7); pages drawn by <code>final_pages.R</code> at commit <code>%s</code>.',
                   esc(key), esc(commit)), file.path(tmp, "footer.html"))

rc <- system2(file.path(R.home("bin"), "Rscript"),
              c(lf_src("inst/scripts/make_gallery_index.R"), "--dir", shQuote(dir), "--title", shQuote(ttl),
                "--card-html", file.path(tmp, "cards.rds"), "--card-href", file.path(tmp, "href.rds"),
                "--blurb", file.path(tmp, "blurb.html"), "--footer", file.path(tmp, "footer.html"),
                "--no-stock-desc"))
if (rc != 0) stop("make_gallery_index.R failed")
print(cards[, .(gene, href)], row.names = FALSE)
