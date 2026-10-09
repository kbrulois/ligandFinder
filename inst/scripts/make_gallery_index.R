#!/usr/bin/env Rscript
## ---- an index page for a directory of per-gene protein pages -------------------
## The per-gene pages are self-contained and know nothing about each other, so a
## set of them needs a front door: this writes an index.html listing every
## <GENE>.html beside it, with a filter box and each gene's accession and size.
##
## Regenerate it whenever pages are added -- it reads the directory, so nothing
## is hardcoded and a gene cannot go missing from the list.
##
##   Rscript inst/scripts/make_gallery_index.R
##   ... --dir ~/AF2_analysis/v2_test --title "..." --out index.html
## ------------------------------------------------------------------------------
suppressMessages({ library(dplyr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(f, d) { i <- match(f, .args); if (is.na(i) || i == length(.args)) d else .args[[i + 1L]] }
dir_p  <- path.expand(.opt("--dir", "~/AF2_analysis/v2_test"))
out_p  <- .opt("--out", file.path(dir_p, "index.html"))
ttl    <- .opt("--title", "LigandFinder — peptide-end predictions")
cache_p <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
## optional per-gene note, as an rds with a `gene` column and a `note` column.
## Used by the cross-validated gallery to say, on each card, how many of that
## gene known ends the page actually holds out -- a page where some ends were
## trained on is not interchangeable with one where none were, and that belongs
## on the card rather than in a readme.
note_p  <- .opt("--notes", NA_character_)

files <- sort(list.files(dir_p, pattern = "\\.html$", full.names = TRUE))
files <- files[basename(files) != basename(out_p)]
if (!length(files)) stop("no pages in ", dir_p, call. = FALSE)
genes <- sub("\\.html$", "", basename(files))

## accession and length come from the residue cache, which is what the pages
## themselves were built from -- absent cache just leaves the columns blank
acc <- len <- rep(NA_character_, length(genes))
if (file.exists(cache_p)) {
  fc <- readRDS(cache_p)
  i <- match(genes, fc$prec$gene)
  acc <- fc$prec$accession[i]
  len <- ifelse(is.na(i), NA_character_, as.character(nchar(fc$prec$seq[i])))
}

mb <- sprintf("%.1f MB", file.size(files) / 1048576)
note <- rep("", length(genes))
if (!is.na(note_p) && file.exists(path.expand(note_p))) {
  nd <- readRDS(path.expand(note_p))
  note <- ifelse(is.na(nd$note[match(genes, nd$gene)]), "",
                 nd$note[match(genes, nd$gene)])
}
esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE)
                     x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }

rows <- paste0(
  '<a class="card" href="', esc(basename(files)), '" data-note="', esc(note), '" data-gene="',
  esc(toupper(paste(genes, ifelse(is.na(acc), "", acc)))), '">',
  '<div class="g">', esc(genes), '</div>',
  '<div class="m">', ifelse(is.na(acc), "", esc(acc)),
  ifelse(is.na(len), "", paste0(" &middot; ", len, " aa")),
  ' &middot; ', mb,
  ifelse(nzchar(note), paste0('</div><div class="n">', esc(note)), ''),
  '</div></a>', collapse = "\n")

html <- sprintf('<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>%s</title>
<style>
 :root{--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--grid:#e6e5e1;--accent:#197EC0}
 *{box-sizing:border-box}
 body{margin:0;background:var(--surface);color:var(--ink);
      font:14px/1.55 ui-sans-serif,system-ui,-apple-system,"Helvetica Neue",Arial,sans-serif}
 header{padding:26px 22px 10px;max-width:1100px;margin:0 auto}
 h1{margin:0 0 6px;font-size:20px;letter-spacing:-0.01em}
 .sub{color:var(--ink2);font-size:12.5px;max-width:760px}
 .sub code{background:#f1f0ec;padding:1px 4px;border-radius:3px;font-size:11.5px}
 .wrap{max-width:1100px;margin:0 auto;padding:0 22px 40px}
 input[type=search]{width:100%%;max-width:340px;margin:16px 0 14px;padding:7px 10px;
   border:1px solid var(--grid);border-radius:6px;font:13px inherit;background:#fff}
 .grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(190px,1fr));gap:10px}
 .card{display:block;padding:11px 13px;border:1px solid var(--grid);border-radius:7px;
   background:#fff;text-decoration:none;color:inherit;transition:border-color .12s,transform .12s}
 .card:hover{border-color:var(--accent);transform:translateY(-1px)}
 .g{font-weight:650;font-size:14.5px;letter-spacing:-0.01em}
 .m{color:var(--ink2);font-size:11.5px;margin-top:2px}
 .n{font-size:11px;margin-top:3px;font-weight:600}
 .card[data-note*="all"] .n{color:#1B7837}
 .card:not([data-note*="all"]) .n{color:#C8690B}
 .none{color:var(--ink2);font-size:13px;padding:14px 0;display:none}
 footer{max-width:1100px;margin:0 auto;padding:0 22px 40px;color:var(--ink2);font-size:11.5px}
 @media (max-width:560px){ .grid{grid-template-columns:repeat(auto-fill,minmax(150px,1fr))} }
</style></head>
<body>
<header>
  <h1>%s</h1>
  <div class="sub">%d genes. Each page carries the per-residue class tracks
    (<b>peptide (inserting)</b>, <b>peptide (non-inserting)</b>, CT- and NT-context),
    the CNN attention head and the xgboost window head for both termini, and the
    window-free MLP and XGB residue models.
    Clicking a point on the <code>xgb window head</code> track opens the 36-residue
    window anchored at that residue &mdash; N points open the N window, C points the C
    window &mdash; drawn on the same residue axis, with the spread across the 20 seeds.
    Hovering marks that residue on every track.</div>
</header>
<div class="wrap">
  <input type="search" id="q" placeholder="filter by gene or accession…" autocomplete="off">
  <div class="grid" id="grid">
%s
  </div>
  <div class="none" id="none">no match</div>
</div>
<footer>Trunks: <code>lf_pepend_run_N_near</code> / <code>lf_pepend_run_C_near</code>,
  20 members each. Built %s.</footer>
<script>
 var q = document.getElementById("q"), cards = [].slice.call(document.querySelectorAll(".card"));
 q.addEventListener("input", function () {
   var v = q.value.trim().toUpperCase(), n = 0;
   cards.forEach(function (c) {
     var hit = !v || c.getAttribute("data-gene").indexOf(v) !== -1;
     c.style.display = hit ? "" : "none"; if (hit) n++;
   });
   document.getElementById("none").style.display = n ? "none" : "block";
 });
 q.focus();
</script>
</body></html>
', esc(ttl), esc(ttl), length(files), rows, format(Sys.Date(), "%e %B %Y"))

writeLines(html, out_p)
message(sprintf("%d page(s) -> %s (%.1f KB)", length(files), out_p, file.size(out_p) / 1024))
