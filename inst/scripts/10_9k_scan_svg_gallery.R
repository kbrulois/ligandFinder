#!/usr/bin/env Rscript
## ---- one HTML holding every scan-track SVG, switchable ------------------------
## Inlines the SVGs written by 10_9i_protein_scan_plot.R into a single
## self-contained page, grouped by where the gene's known end sits (training
## split, validation split, or neither), with click and arrow-key switching.
##
## ID COLLISION IS THE REASON THIS IS NOT A CAT. svglite emits internal
## <clipPath id="cpXXXX"> and references them as url(#cpXXXX). Those ids are only
## unique within one file, and several of these figures reuse the same ones -- so
## concatenated into one document, the later figures' clip paths resolve to the
## earlier figures' rectangles and panels render blank or clipped to the wrong
## box. Every id is therefore namespaced with the gene before inlining.
##
## The SVGs are up to 42 inches wide, so the page defaults to fitting them to the
## window (the track shape reads, the sequence letters do not) with a toggle to
## 1:1, where the panel scrolls horizontally and the letters are legible.
##
##   Rscript inst/scripts/10_9k_scan_svg_gallery.R
##
## Output: ~/AF2_analysis/lf_scan_gallery_<term>.html
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term   <- toupper(.opt("--term", "C"))
meta_p <- path.expand(.opt("--meta", "~/AF2_analysis/lf_scan_protein_meta.rds"))
out_p  <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_scan_gallery_%s.html", term)))
## The page used to hardcode its own description, which silently lied the moment
## the meta held a different kind of figure. Both are overridable.
ttl_s  <- .opt("--title", sprintf("Step-1 scan tracks \u2014 %s terminus", term))
sub_s  <- .opt("--subtitle", paste0(
  "One window per mature residue, scored by the 20-seed peptide-end model. ",
  "x is the window&rsquo;s <b>anchor</b>: the residue a cleavage there would leave as the ",
  "peptide&rsquo;s last. Band is &plusmn;1 sd across seeds."))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "viz_theme.R"))

meta <- readRDS(meta_p)
meta$group[is.na(meta$group)] <- "other"
## The headings used to be a fixed train/val/other, and `filter(group == g)` over
## that fixed set SILENTLY DROPPED every gene in any other group -- a meta with
## 71 rows rendered 59, because 12 genes have knowns in both splits and are
## grouped "train/val". Start from the known labels, then append whatever else
## the meta actually contains, so a gene cannot vanish from the page.
GRP <- c(train = "Training split", `train/val` = "Both splits",
         val = "Validation split", other = "Not in the training set")
.extra <- setdiff(unique(meta$group), names(GRP))
if (length(.extra)) {
  GRP <- c(GRP, stats::setNames(.extra, .extra))
  message("note: group(s) not in the standard set, kept as their own heading: ",
          paste(.extra, collapse = ", "))
}
GRP <- GRP[names(GRP) %in% unique(meta$group)]
meta <- meta %>% mutate(group = factor(group, levels = names(GRP))) %>%
  arrange(group, desc(best_score))
stopifnot(!any(is.na(meta$group)))
message(nrow(meta), " figures: ", paste(sprintf("%s %d", names(GRP),
        as.integer(table(factor(meta$group, levels = names(GRP))))), collapse = ", "))

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

## ---- inline one SVG, with its ids namespaced -------------------------------------
inline_svg <- function(path, pfx) {
  if (!file.exists(path)) return(NULL)
  s <- paste(readLines(path, warn = FALSE), collapse = "\n")
  s <- sub("(?s)^.*?(<svg)", "\\1", s, perl = TRUE)     # drop the XML prolog
  ## svglite quotes attributes with SINGLE quotes, and names its clip paths after
  ## the base64 of the clip rectangle -- cpMC4wMHw3OTIuMDB8MC4wMHw0MjQuODA= is
  ## "0.00|792.00|0.00|424.80", the whole canvas. Every figure of the same size
  ## therefore emits the SAME id, so this must catch both quote styles or the
  ## namespacing silently does nothing and the panels clip to each other.
  ids <- unique(c(
    gsub("^id='|'$", "", regmatches(s, gregexpr("id='[^']+'", s))[[1]]),
    gsub('^id="|"$', '', regmatches(s, gregexpr('id="[^"]+"', s))[[1]])))
  ## the ids are long and distinctive, and appear only as id='X' and url(#X), so
  ## a plain textual replace covers every reference form at once
  for (i in ids) s <- gsub(i, paste0(pfx, "_", i), s, fixed = TRUE)
  ## keep the viewBox, drop the fixed width/height so CSS can scale it
  hd <- regmatches(s, regexpr("<svg[^>]*>", s))
  if (length(hd)) {
    nh <- gsub(" (width|height)='[^']*'", "", gsub(' (width|height)="[^"]*"', '', hd))
    s <- sub(hd, nh, s, fixed = TRUE)
  }
  list(svg = s, n_ids = length(ids))
}

## ---- per-gene metadata line -------------------------------------------------------
meta_html <- function(r) {
  k <- r$knowns[[1]]
  kn <- if (!nrow(k)) "<span class='muted'>no known end in this set</span>" else
    paste(sprintf("<b>%s</b> ends at %d &mdash; score %.3f, rank %s of 2,980,852 <span class='muted'>(%s, %s)</span>",
                  esc(k$pep_name), k$anchor, k$score, format(k$rank_all, big.mark = ","),
                  esc(k$split), esc(k$motif)), collapse = "<br>")
  sprintf(paste0("<div class='m'><span class='k'>%s</span> <span class='muted'>%s</span>",
                 " &middot; %d residues, %s scanned &middot; %d local peaks",
                 " &middot; best <b>%.3f</b> at %d</div><div class='m'>%s</div>"),
          esc(r$gene), esc(r$accession), r$n_res, format(r$n_scanned, big.mark = ","),
          r$n_peaks, r$best_score, r$best_anchor, kn)
}

## ---- build ------------------------------------------------------------------------
btns <- list(); figs <- list(); metas <- list(); total_ids <- 0L
for (g in names(GRP)) {
  sub <- meta %>% filter(group == g)
  if (!nrow(sub)) next
  b <- sprintf(paste0("<button class='g' data-gene='%s'><span class='bn'>%s</span>",
                      "<span class='bs'>%.2f</span></button>"),
               esc(sub$gene), esc(sub$gene), sub$best_score)
  btns[[g]] <- sprintf("<div class='grp'><h2>%s <span class='muted'>(%d)</span></h2><div class='row'>%s</div></div>",
                       GRP[[g]], nrow(sub), paste(b, collapse = ""))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    got <- inline_svg(r$svg, r$gene)
    if (is.null(got)) { message("missing svg: ", r$svg); next }
    total_ids <- total_ids + got$n_ids
    figs[[r$gene]]  <- sprintf("<figure id='fig-%s' class='fig' hidden>%s</figure>", esc(r$gene), got$svg)
    metas[[r$gene]] <- sprintf("<div id='meta-%s' class='mw' hidden>%s</div>", esc(r$gene), meta_html(r))
  }
}
first <- meta$gene[1]

html <- sprintf('<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>%s</title>
<style>
 :root{--surface:%s;--ink:%s;--ink2:%s;--grid:%s;--s1:%s;--s2:%s}
 *{box-sizing:border-box}
 body{margin:0;background:var(--surface);color:var(--ink);
      font:14px/1.5 ui-sans-serif,system-ui,-apple-system,"Helvetica Neue",Arial,sans-serif}
 header{padding:18px 20px 10px}
 h1{margin:0 0 4px;font-size:18px}
 .sub{color:var(--ink2);font-size:12.5px;max-width:110ch}
 nav{padding:4px 20px 12px;border-bottom:1px solid var(--grid)}
 .grp{margin:10px 0 0}
 .grp h2{margin:0 0 6px;font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--ink2);font-weight:600}
 .row{display:flex;flex-wrap:wrap;gap:6px}
 button.g{display:flex;align-items:baseline;gap:7px;padding:5px 10px;border:1px solid var(--grid);
   border-radius:6px;background:#fff;color:var(--ink);cursor:pointer;font:inherit;font-size:12.5px}
 button.g:hover{border-color:var(--ink2)}
 button.g[aria-current=true]{background:var(--s1);border-color:var(--s1);color:#fff}
 button.g[aria-current=true] .bs{color:#fff;opacity:.85}
 .bn{font-weight:600}
 .bs{font-variant-numeric:tabular-nums;color:var(--ink2);font-size:11.5px}
 main{padding:14px 20px 40px}
 .bar{display:flex;align-items:center;gap:14px;flex-wrap:wrap;margin-bottom:10px}
 .m{font-size:12.5px;color:var(--ink);margin:2px 0}
 .k{font-weight:700;font-size:15px}
 .muted{color:var(--ink2)}
 .tog{margin-left:auto;display:flex;gap:6px;align-items:center;color:var(--ink2);font-size:12px}
 .tog button{padding:4px 9px;border:1px solid var(--grid);border-radius:6px;background:#fff;
   cursor:pointer;font:inherit;font-size:12px;color:var(--ink)}
 .tog button[aria-pressed=true]{background:var(--ink);color:#fff;border-color:var(--ink)}
 .viewer{border:1px solid var(--grid);border-radius:8px;background:#fff;padding:6px;overflow-x:auto}
 .fig{margin:0}
 .fig svg{display:block;height:auto}
 body[data-fit=fit] .fig svg{width:100%%}
 body[data-fit=actual] .fig svg{width:var(--w)}
 kbd{border:1px solid var(--grid);border-bottom-width:2px;border-radius:4px;padding:0 4px;font-size:11px}
</style></head>
<body data-fit="fit">
<header>
 <h1>%s</h1>
 <div class="sub">%s %d proteins.
  <kbd>&larr;</kbd><kbd>&rarr;</kbd> to step between them.</div>
</header>
<nav>%s</nav>
<main>
 <div class="bar"><div id="metabox">%s</div>
  <div class="tog">width
   <button id="bfit"  aria-pressed="true">fit</button>
   <button id="bact"  aria-pressed="false">1:1</button></div></div>
 <div class="viewer" id="viewer">%s</div>
</main>
<script>
(function(){
 var W = %s;
 var genes = %s;
 var cur = null;
 function show(g){
   if(!document.getElementById("fig-"+g)) return;
   document.querySelectorAll(".fig").forEach(function(f){f.hidden = f.id !== "fig-"+g;});
   document.querySelectorAll(".mw").forEach(function(m){m.hidden = m.id !== "meta-"+g;});
   document.querySelectorAll("button.g").forEach(function(b){
     b.setAttribute("aria-current", String(b.dataset.gene === g));});
   document.body.style.setProperty("--w", (W[g]*96)+"px");
   document.getElementById("viewer").scrollLeft = 0;
   cur = g;
   try{ history.replaceState(null,"","#"+g); }catch(e){}
 }
 document.querySelectorAll("button.g").forEach(function(b){
   b.addEventListener("click", function(){ show(b.dataset.gene); });});
 function step(d){ var i = genes.indexOf(cur); if(i<0) return;
   show(genes[(i + d + genes.length) %% genes.length]); }
 document.addEventListener("keydown", function(e){
   if(e.target.tagName === "INPUT") return;
   if(e.key === "ArrowRight"){ step(1); e.preventDefault(); }
   if(e.key === "ArrowLeft"){ step(-1); e.preventDefault(); }});
 function fit(on){
   document.body.dataset.fit = on ? "fit" : "actual";
   document.getElementById("bfit").setAttribute("aria-pressed", String(on));
   document.getElementById("bact").setAttribute("aria-pressed", String(!on)); }
 document.getElementById("bfit").addEventListener("click", function(){ fit(true); });
 document.getElementById("bact").addEventListener("click", function(){ fit(false); });
 var h = (location.hash||"").replace("#","");
 show(genes.indexOf(h) >= 0 ? h : %s);
})();
</script></body></html>',
  ttl_s, LF_VIZ$surface, LF_VIZ$ink, LF_VIZ$ink2, LF_VIZ$grid, LF_VIZ$s1, LF_VIZ$s2,
  ttl_s, sub_s, nrow(meta),
  paste(unlist(btns), collapse = "\n"),
  paste(unlist(metas), collapse = "\n"),
  paste(unlist(figs), collapse = "\n"),
  paste0("{", paste(sprintf('"%s":%.2f', meta$gene, meta$width_in), collapse = ","), "}"),
  paste0("[", paste(sprintf('"%s"', meta$gene), collapse = ","), "]"),
  sprintf('"%s"', first))

writeLines(html, out_p)
message(sprintf("\n%d figures, %s ids namespaced -> %s (%.1f MB)",
                length(figs), format(total_ids, big.mark = ","), out_p,
                file.info(out_p)$size / 1e6))
