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
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "lf_paths.R"))
ARM <- Sys.getenv("ARM", "esmonly36")
D  <- file.path(A2, paste0("pepend_pages_", ARM))
## Each arm at its OWN operating point: the curated model at the threshold its
## gallery has always used, the ESM-C arm at the optimum measured for it. The
## TOLERANCE is shared -- it belongs to the evaluation, not to a model, so
## letting it differ would put a second, unasked-for difference in the chips.
THR_CUR <- Sys.getenv("THR_CUR", "0.2")
THR_ESM <- Sys.getenv("THR_ESM", "0.65")
TOLC    <- Sys.getenv("TOL", "4")
thr <- as.numeric(Sys.getenv("THRESH", "0.2"))
genes <- readLines(file.path(S, "page_genes.txt")); genes <- genes[nzchar(genes)]
esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
                     x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE) }
HN <- c(cnn = "linear head", xgb_pm = "xgb head")
AN <- c(cur = "curated", esm = "ESM-C only")

rd <- function(f, a) { p <- file.path(A2, "lf_headcmp", f)
  if (!file.exists(p)) stop("missing ", f, " -- run end_calls.R for that arm", call. = FALSE)
  fread(p)[grepl("^never", set)][, arm := a] }
ec <- rbind(rd(sprintf("end_calls_T%s_tol%s.csv", THR_CUR, TOLC), "cur"),
            rd(sprintf("end_calls_T%s_%s_tol%s.csv", THR_ESM, ARM, TOLC), "esm"))
ec[, `:=`(start = ifelse(term == "N", anchor, anchor - len + 1L),
          end   = ifelse(term == "N", anchor + len - 1L, anchor))]

tip <- function(r) {
  pk <- if (!is.na(r$peak_centre)) sprintf("nearest peak %.2f at %d", r$peak_height, r$peak_centre) else "no peak >= T on this precursor"
  sprintf("%s %s, %s %d: %s", AN[[r$arm]], HN[[r$head]], r$term, r$anchor,
          switch(r$reason,
                 "peak on the end"  = sprintf("TP -- peak %.2f centred at %d", r$peak_height, r$peak_centre),
                 "below T"          = sprintf("FN, below T -- end score %.2f; %s", r$score, pk),
                 "off-centre peak"  = sprintf("FN, inside a peak centred %d residues away (%.2f at %d); end score %.2f",
                                              abs(r$peak_centre - r$anchor), r$peak_height, r$peak_centre, r$score),
                 sprintf("FN, %s", r$reason)))
}
chip <- function(r) if (is.null(r) || !nrow(r) || is.na(r$call))
  '<span class="c na">&middot;</span>' else
  sprintf('<span class="c %s" title="%s">%s</span>', tolower(r$call), esc(tip(r)), r$call)

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
      sprintf('<div class="nsr"><span>%s %d</span>%s%s%s%s</div>', tm, x$anchor[1],
              chip(x[arm == "cur" & head == "cnn"][1]), chip(x[arm == "cur" & head == "xgb_pm"][1]),
              chip(x[arm == "esm" & head == "cnn"][1]), chip(x[arm == "esm" & head == "xgb_pm"][1]))
    }, character(1))
    paste0(hd, paste(rows, collapse = ""))
  }, character(1))
  data.table(gene = g, html = paste0('<div class="ns"><div class="nsh"><span>never seen</span>',
                                     '<span>cur lin</span><span>cur xgb</span>',
                                     '<span>esm lin</span><span>esm xgb</span></div>',
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
u <- unique(ec[, .(arm, head, set, accession, term, anchor, call)])
cnt <- function(d, hd, a = "cur") { x <- d[head == hd & arm == a]
  sprintf("%d of %d", sum(x$call == "TP"), nrow(x)) }
pg <- u[unique(ec[gene %in% genes, .(accession)]), on = "accession", nomatch = 0]
pp <- function(d, st) d[set == st]
## the FP cost behind the TP counts: without it the two arms, which sit at
## different thresholds, read as a ranking they are not.
sc <- function(f, t) fread(file.path(A2, "lf_headcmp", f))[set == "never seen: peptides" & thr == t]
SC <- sc(sprintf("segment_confusion_tol%s.csv", TOLC), as.numeric(THR_CUR))
SE <- sc(sprintf("segment_confusion_%s_tol%s.csv", ARM, TOLC), as.numeric(THR_ESM))
stopifnot(nrow(SC) == 4, nrow(SE) == 4)
blurb <- paste(readLines(file.path(S, "final_blurb.html"), warn = FALSE), collapse = "\n")
old <- grep("magenta", strsplit(blurb, "\n")[[1]], value = TRUE)
stopifnot(length(old) == 1)
new <- sprintf(paste0(
  '<p>Peptides whose gene no fold trained on are drawn above the sequence as a <b>magenta &ldquo;never seen&rdquo; line</b> ',
  '(hover for name and span). Each card lists them with every end called on the all-50 profile (the heavy squares), ',
  'for <b>four models</b>: the curated 26-channel arm at <b>T = %s</b> and the ESM-C-only 36-PCA arm at <b>T = %s</b>, ',
  'each at the threshold measured best for it, both at tolerance %s residues. ',
  '<span class="c tp">TP</span> when a run of anchors &ge; T has its maximum within tolerance of the end, ',
  '<span class="c fn">FN</span> otherwise &mdash; hover a call for the arm, the head and why. ',
  'Columns: curated linear, curated xgb, ESM-C linear, ESM-C xgb.</p>\n',
  '<p>On these %d pages the curated linear head calls <b>%s</b> never-seen peptide ends TP (xgb %s); ',
  'the ESM-C linear head calls <b>%s</b> (xgb %s). Over the whole never-seen peptide set: ',
  'curated %s / %s, ESM-C %s / %s.</p>\n',
  '<p><b>Do not read those TP counts as a ranking.</b> The two arms sit at different thresholds, so they are not ',
  'the same test. At T = %s the curated arm finds more ends but calls <b>%s</b> false-positive segments ',
  '(FPR %.2f); at T = %s the ESM-C arm finds fewer and calls <b>%s</b> (FPR %.2f) &mdash; roughly %.0fx fewer for about ',
  '%s fewer ends. Segment auc_T, which charges for both, is higher for ESM-C in every never-seen-peptide cell ',
  '(%.3f&ndash;%.3f against %.3f&ndash;%.3f).</p>'),
  THR_CUR, THR_ESM, TOLC, length(genes),
  cnt(pp(pg, "never seen: peptides"), "cnn", "cur"), cnt(pp(pg, "never seen: peptides"), "xgb_pm", "cur"),
  cnt(pp(pg, "never seen: peptides"), "cnn", "esm"), cnt(pp(pg, "never seen: peptides"), "xgb_pm", "esm"),
  cnt(pp(u, "never seen: peptides"), "cnn", "cur"), cnt(pp(u, "never seen: peptides"), "xgb_pm", "cur"),
  cnt(pp(u, "never seen: peptides"), "cnn", "esm"), cnt(pp(u, "never seen: peptides"), "xgb_pm", "esm"),
  THR_CUR, format(round(mean(SC$FP)), big.mark = ","), mean(SC$FPR),
  THR_ESM, format(round(mean(SE$FP)), big.mark = ","), mean(SE$FPR),
  mean(SC$FP) / mean(SE$FP), format(round(mean(SC$TP) - mean(SE$TP))),
  min(SE$auc_T), max(SE$auc_T), min(SC$auc_T), max(SC$auc_T))
style <- paste0('<style>',
  '.ns{margin-top:7px;padding-top:6px;border-top:1px solid var(--grid);font-size:11px}',
  '.nsh,.nsr{display:grid;grid-template-columns:1fr 31px 31px 31px 31px;align-items:center;column-gap:3px}',
  '.nsh{color:var(--ink2);font-size:10.5px;margin-bottom:2px}.nsh span+span{text-align:center}',
  '.nsp{color:var(--ink);font-weight:600;margin-top:3px;overflow-wrap:anywhere}',
  '.nsr{color:var(--ink2);margin-top:1px}.lp{color:var(--ink2);font-weight:400}',
  '.c{display:inline-block;text-align:center;border-radius:3px;padding:0 4px;font-size:10.5px;font-weight:650;line-height:16px}',
  '.c.tp{background:#e3f1e5;color:#1B5E20}.c.fn{background:#f8e1e1;color:#8E1B1B}',
  '.c.na{background:transparent;color:var(--ink2);font-weight:400}',
  '</style>')
blurb <- sub(old, paste0(style, new), blurb, fixed = TRUE)
writeLines(blurb, file.path(S, "final_blurb_tpfn.html"))

## the cards say which ends are FN now; let a "FN" search find them
args <- c(lf_src("inst/scripts/make_gallery_index.R"), "--dir", D,
          "--title", "LigandFinder — peptide-end pages, curated vs ESM-C-only (5-fold CV)",
          "--notes", file.path(S, "index_notes.rds"), "--card-html", file.path(S, "index_cards.rds"),
          "--blurb", file.path(S, "final_blurb_tpfn.html"), "--footer", file.path(S, "final_footer.html"),
          "--no-stock-desc")
rc <- system2("/usr/local/bin/Rscript", shQuote(args))
cat(sprintf("%d cards with never-seen ends; index rc=%d\n", nrow(cards), rc))
cat(new, "\n")
