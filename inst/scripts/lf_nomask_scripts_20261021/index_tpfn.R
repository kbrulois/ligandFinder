## The final gallery's index, with every never-seen end called TP or FN at one
## threshold (Kevin, 2026-10-22: "annotate TP and FN unseens at a threshold of
## 0.2"). Calls come from end_calls.R (THRESH must match): the segment rule on
## the all-50 profile -- the heavy squares on the page -- for both heads.
## Writes the card table (index_cards.rds), training-end notes without the
## never-seen spans (index_notes.rds; the card table lists those now) and a
## blurb with the counts, then runs the repo's make_gallery_index.R.
##   THRESH=0.2 Rscript index_tpfn.R
suppressMessages({library(data.table); library(tibble)})
A2 <- path.expand("~/AF2_analysis"); S <- file.path(A2, "lf_nomask_scripts_20261021")
W  <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7"
D  <- file.path(A2, "pepend_pages_flatten_cv5")
thr <- as.numeric(Sys.getenv("THRESH", "0.2"))
genes <- readLines(file.path(S, "page_genes.txt")); genes <- genes[nzchar(genes)]
esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
                     x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE) }
HN <- c(cnn = "linear head", xgb_pm = "xgb head")

ec <- fread(file.path(A2, "lf_headcmp", sprintf("end_calls_T%s.csv", thr)))[grepl("^never", set)]
ec[, `:=`(start = ifelse(term == "N", anchor, anchor - len + 1L),
          end   = ifelse(term == "N", anchor + len - 1L, anchor))]

tip <- function(r) {
  pk <- if (!is.na(r$peak_centre)) sprintf("nearest peak %.2f at %d", r$peak_height, r$peak_centre) else "no peak >= T on this precursor"
  sprintf("%s, %s %d: %s", HN[[r$head]], r$term, r$anchor,
          switch(r$reason,
                 "peak on the end"  = sprintf("TP -- peak %.2f centred at %d", r$peak_height, r$peak_centre),
                 "below T"          = sprintf("FN, below T -- end score %.2f; %s", r$score, pk),
                 "off-centre peak"  = sprintf("FN, inside a peak centred %d residues away (%.2f at %d); end score %.2f",
                                              abs(r$peak_centre - r$anchor), r$peak_height, r$peak_centre, r$score),
                 sprintf("FN, %s", r$reason)))
}
chip <- function(r) sprintf('<span class="c %s" title="%s">%s</span>',
                            tolower(r$call), esc(tip(r)), r$call)

cards <- rbindlist(lapply(genes, function(g) {
  e <- ec[gene == g]
  if (!nrow(e)) return(NULL)
  pep <- unique(e[, .(set, name, start, end)])[order(start, -end)]
  body <- vapply(seq_len(nrow(pep)), function(i) {
    p <- pep[i]
    hd <- sprintf('<div class="nsp">%s %d&ndash;%d%s</div>', esc(p$name), p$start, p$end,
                  if (grepl("large", p$set)) ' <span class="lp">large product</span>' else "")
    ends <- e[set == p$set & name == p$name & start == p$start & end == p$end]
    rows <- vapply(c("N", "C"), function(tm) {
      x <- ends[term == tm]; if (!nrow(x)) return("")
      sprintf('<div class="nsr"><span>%s %d</span>%s%s</div>', tm, x$anchor[1],
              chip(x[head == "cnn"][1]), chip(x[head == "xgb_pm"][1]))
    }, character(1))
    paste0(hd, paste(rows, collapse = ""))
  }, character(1))
  data.table(gene = g, html = paste0('<div class="ns"><div class="nsh"><span>never seen, T ', thr,
                                     '</span><span>linear</span><span>xgb</span></div>',
                                     paste(body, collapse = ""), '</div>'))
}))
saveRDS(as_tibble(cards), file.path(S, "index_cards.rds"))

## training-end notes, as final_notes.R writes them, minus the never-seen spans
kn <- rbindlist(lapply(c("N", "C"), function(t) {
  d <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", t)))
  d[, .(gene = sub("_w-?[0-9]+--?[0-9]+$", "", peps), term = t, anchor, fold)]
}))
note <- vapply(genes, function(g) {
  k <- kn[gene == g][order(term, anchor)]
  if (!nrow(k)) "no known end in any fold: all tracks out-of-sample" else
    paste0("training ends (4 of 5 folds saw each; held out by): ",
           paste(sprintf("%s %d (fold %d)", k$term, k$anchor, k$fold), collapse = ", "))
}, character(1))
saveRDS(tibble(gene = genes, note = unname(note)), file.path(S, "index_notes.rds"))

## counts for the blurb: unique ends (an end shared by two names counts once)
u <- unique(ec[, .(head, set, accession, term, anchor, call)])
cnt <- function(d, hd) { x <- d[head == hd]; sprintf("%d of %d", sum(x$call == "TP"), nrow(x)) }
pg <- u[unique(ec[gene %in% genes, .(accession)]), on = "accession", nomatch = 0]
pp <- function(d, st) d[set == st]
blurb <- paste(readLines(file.path(S, "final_blurb.html"), warn = FALSE), collapse = "\n")
old <- grep("magenta", strsplit(blurb, "\n")[[1]], value = TRUE)
stopifnot(length(old) == 1)
new <- sprintf(paste0(
  '<p>Peptides whose gene no fold trained on are drawn above the sequence as a <b>magenta &ldquo;never seen&rdquo; line</b> ',
  '(hover for name and span). Each card lists them with every end <b>called at T = %s</b> on the all-50 profile (the heavy ',
  'squares): <span class="c tp">TP</span> when a run of anchors &ge; T has its maximum within 2 residues of the end, ',
  '<span class="c fn">FN</span> otherwise &mdash; hover a call for why (below T, or inside a peak centred elsewhere) and the ',
  'scores. Left call = CNN linear head, right = xgb head.</p>\n',
  '<p>On these %d pages the linear head calls <b>%s</b> never-seen peptide ends TP and %s large-product ends ',
  '(xgb head: %s and %s). Over the whole never-seen set: peptides %s (xgb %s), large products %s (xgb %s).</p>'),
  thr, length(genes),
  cnt(pp(pg, "never seen: peptides"), "cnn"), cnt(pp(pg, "never seen: large products"), "cnn"),
  cnt(pp(pg, "never seen: peptides"), "xgb_pm"), cnt(pp(pg, "never seen: large products"), "xgb_pm"),
  cnt(pp(u, "never seen: peptides"), "cnn"), cnt(pp(u, "never seen: peptides"), "xgb_pm"),
  cnt(pp(u, "never seen: large products"), "cnn"), cnt(pp(u, "never seen: large products"), "xgb_pm"))
style <- paste0('<style>',
  '.ns{margin-top:7px;padding-top:6px;border-top:1px solid var(--grid);font-size:11px}',
  '.nsh,.nsr{display:grid;grid-template-columns:1fr 38px 38px;align-items:center;column-gap:4px}',
  '.nsh{color:var(--ink2);font-size:10.5px;margin-bottom:2px}.nsh span+span{text-align:center}',
  '.nsp{color:var(--ink);font-weight:600;margin-top:3px;overflow-wrap:anywhere}',
  '.nsr{color:var(--ink2);margin-top:1px}.lp{color:var(--ink2);font-weight:400}',
  '.c{display:inline-block;text-align:center;border-radius:3px;padding:0 4px;font-size:10.5px;font-weight:650;line-height:16px}',
  '.c.tp{background:#e3f1e5;color:#1B5E20}.c.fn{background:#f8e1e1;color:#8E1B1B}',
  '</style>')
blurb <- sub(old, paste0(style, new), blurb, fixed = TRUE)
writeLines(blurb, file.path(S, "final_blurb_tpfn.html"))

## the cards say which ends are FN now; let a "FN" search find them
args <- c(file.path(W, "inst/scripts/make_gallery_index.R"), "--dir", D,
          "--title", "LigandFinder — peptide-end pages, linear window head (5-fold CV)",
          "--notes", file.path(S, "index_notes.rds"), "--card-html", file.path(S, "index_cards.rds"),
          "--blurb", file.path(S, "final_blurb_tpfn.html"), "--footer", file.path(S, "final_footer.html"),
          "--no-stock-desc")
rc <- system2("/usr/local/bin/Rscript", shQuote(args))
cat(sprintf("%d cards with never-seen ends; index rc=%d\n", nrow(cards), rc))
cat(new, "\n")
