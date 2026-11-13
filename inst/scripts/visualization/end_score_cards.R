## Per-end score chips on a two-model gallery index (Kevin, 2026-11-12: "redo
## pepend_pages_curated_AAswap_msacons_padg/index.html to have end score colors
## for each model type", then "make the color scale turbo").
##
## Each card lists every known end on its gene -- never-seen ends first (magenta:
## peptide ends; purple: large cleavage products), then the training ends with
## the fold that held each one out -- with one chip per model and window head:
## the curated baseline (lf_pepend_run_*_nomask_flatten, the pages' "curated"
## track) and ARM, each with its CNN linear head and its xgb window head. A chip
## is the best score within +/-TOL residues of the end, read from the per-gene
## scans the pages draw: the all-50 ensemble for a never-seen end (the heavy
## squares), the fold that held the end out for a training end (the triangles).
## Coloured on turbo, 0 -> 1.
##
## Each model is read with its OWN known set: the training ends of its own
## windows (INPUT_SFX), each scored by the fold that held it out in THAT model.
## Inside the +/-TOL window, an anchor that is another known end of the same
## model which the scoring members TRAINED on is skipped: GLP-1(7-37) ends one
## residue past GLP-1(7-36)amide, so without that rule the baseline's fold-5
## chip for 7-36 would be its score on a window it was trained on.
##
## Raw scores are not calibrated between models -- the _padg window fix lowers an
## arm's scores everywhere -- so the blurb says to compare chips down a column.
##
## Writes <GALLERY>/end_scores.csv and rebuilds <GALLERY>/index.html through
## inst/scripts/visualization/make_gallery_index.R (--card-html / --blurb / --footer), so the
## filter box and the cards are the stock ones.
##
##   ARM=esmonly36 LC_ALL=en_US.UTF-8 Rscript end_score_cards.R
##   env: ARM (required), GALLERY (~/AF2_analysis/pepend_pages_<ARM>), TOL (4),
##        TITLE (index title), OUT_CSV (<GALLERY>/end_scores.csv),
##        INDEX (<GALLERY>/index.html)
suppressMessages({library(data.table); library(reticulate)})
## package source tree: $LF_SRC, or the checkout these scripts sit in (lf_paths.R)
.lf_here <- c(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), getwd())[1]
source(file.path(.lf_here, "..", "lf_paths.R"))
## arm_var(): the arm's settings, read out of arms/<arm>.sh
source(lf_src("inst/scripts/modeling/arm_vars.R"))

A2  <- path.expand("~/AF2_analysis"); S <- path.expand(Sys.getenv("LF_RUN", file.path(A2, "lf_run")))
ARM <- sub("^_", "", Sys.getenv("ARM", ""))
if (!nzchar(ARM) || !file.exists(file.path(.lf_scripts, "arms", paste0(ARM, ".sh"))))
  stop("ARM must name an arm in arms/, e.g. ARM=esmonly36", call. = FALSE)
GAL     <- path.expand(Sys.getenv("GALLERY", file.path(A2, paste0("pepend_pages_", ARM))))
TOL     <- as.integer(Sys.getenv("TOL", "4"))
OUT_CSV <- path.expand(Sys.getenv("OUT_CSV", file.path(GAL, "end_scores.csv")))
INDEX   <- path.expand(Sys.getenv("INDEX", file.path(GAL, "index.html")))
SFX     <- arm_var(ARM, "INPUT_SFX")        # the arm's training windows; "" = the original ones
LAB     <- arm_var(ARM, "ARM_LABEL")
TITLE   <- Sys.getenv("TITLE", sprintf("LigandFinder &mdash; known-peptide pages, 5-fold CV: the original curated model vs %s", LAB))
np <- import("numpy", convert = TRUE)

genes <- sub("\\.html$", "", list.files(GAL, pattern = "\\.html$"))
genes <- setdiff(genes, sub("\\.html$", "", basename(INDEX)))
if (!length(genes)) stop("no <GENE>.html pages in ", GAL, call. = FALSE)

## ---- the two models -------------------------------------------------------------
## `npz` the per-gene scan suffix; `sfx` the training windows (INPUT_SFX)
MODELS <- list(cur = list(label = "curated", npz = "",               sfx = ""),
               arm = list(label = LAB,       npz = paste0("_", ARM), sfx = SFX))
HEADS  <- c(cnn = "CNN linear head", xgb_pm = "xgb window head")

## ---- known ends, per model -------------------------------------------------------
## the fold that holds each known out, in that model's own windows (the shared
## table for the original ones; arm_knowns_cache.R writes the _<sfx> tables)
knowns_of <- function(sfx) rbindlist(lapply(c("N", "C"), function(tm) {
  f <- file.path(A2, sprintf("lf_pepend_cv5_knowns_%s%s.csv", tm, sfx))
  if (!file.exists(f)) stop("no ", basename(f), " -- run INPUT_SFX=", sfx, " Rscript arm_knowns_cache.R",
                            call. = FALSE)
  d <- fread(f)
  data.table(accession = d$accession, term = tm, anchor = as.integer(d$anchor), fold = as.integer(d$fold))
}))
## a name for each end, from the window set's own knowns (one window per end,
## labelled from the shortest peptide that ends there)
names_of <- function(sfx) {
  k <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_nn_input_cv5_f1%s.rds", sfx)))$knowns)
  data.table(accession = k$accession, term = k$term, anchor = as.integer(k$anchor), gene = k$gene,
             name = ifelse(is.na(k$pep_name), sub("^[^_]+_", "", k$pep_id), k$pep_name))
}
kn <- lapply(MODELS, function(m) knowns_of(m$sfx))
nm <- unique(rbindlist(lapply(unique(vapply(MODELS, `[[`, "", "sfx")), names_of)),
             by = c("accession", "term", "anchor"))

## training ends: one row per end either model knows, with each model's fold
## (NA where it is not a known end of that model)
tr <- merge(kn$cur[, .(accession, term, anchor, fold_cur = fold)],
            kn$arm[, .(accession, term, anchor, fold_arm = fold)],
            by = c("accession", "term", "anchor"), all = TRUE)
tr <- nm[tr, on = c("accession", "term", "anchor")]
tr[, set := "trained"]
## never-seen ends: the arm's never-seen set (its negatives' overlaps removed),
## scored by both models' all-50 ensembles
un <- fread(file.path(S, sprintf("cv5_unseen_ends%s.csv", SFX)))
un <- un[, .(accession, term, anchor = as.integer(anchor), gene, name,
             set = ifelse(grepl("large", set), "never seen: large product", "never seen: peptide"))]
un[, `:=`(fold_cur = NA_integer_, fold_arm = NA_integer_)]
ends <- rbind(un, tr, use.names = TRUE)
## gene symbols as the pages and the scans key them
gs <- fread(file.path(A2, "lf_plm", "precursor_seqs.tsv"), select = c("accession", "gene"))
ends[is.na(gene), gene := gs$gene[match(accession, gs$accession)]]
ends <- ends[gene %in% genes]

## ---- scores ---------------------------------------------------------------------
.npz <- list()
scan <- function(tm, run, sfx) {
  f <- file.path(A2, sprintf("lf_winxgb_scan_%s_%s_nomask_flatten%s.npz", tm, run, sfx))
  if (is.null(.npz[[f]])) {
    if (!file.exists(f)) stop("no per-gene scan ", basename(f), call. = FALSE)
    .npz[[f]] <<- np$load(f, allow_pickle = TRUE)
  }
  .npz[[f]]
}
score_end <- function(r, k, head) {
  m <- MODELS[[k]]; a <- r$anchor; fo <- r[[paste0("fold_", k)]]
  never <- r$set != "trained"
  if (is.na(a) || (!never && is.na(fo))) return(list(v = NA_real_, at = NA_integer_))
  z <- scan(r$term, if (never) "cv5all50" else sprintf("cv5_f%d", fo), m$npz)
  key <- paste0(r$gene, "__anchor")
  if (!key %in% z$files) stop("scan has no ", r$gene, call. = FALSE)
  an <- as.integer(z[[key]]); s <- as.numeric(z[[paste0(r$gene, "__", head)]])
  ## the members trained on every known end of theirs outside the scoring fold:
  ## never score one of those as this end
  own <- kn[[k]][accession == r$accession & term == r$term]
  seen <- if (never) own$anchor else own[fold != fo]$anchor
  w <- abs(an - a) <= TOL & !(an %in% setdiff(seen, a))
  if (!any(w)) return(list(v = NA_real_, at = NA_integer_))
  i <- which(w)[which.max(s[w])]
  list(v = s[i], at = an[i])
}
for (k in names(MODELS)) for (h in names(HEADS)) {
  sc <- lapply(seq_len(nrow(ends)), function(i) score_end(ends[i], k, h))
  set(ends, j = sprintf("%s_%s", k, h),    value = vapply(sc, `[[`, 0, "v"))
  set(ends, j = sprintf("%s_%s_at", k, h), value = vapply(sc, `[[`, 0L, "at"))
}
fwrite(ends, OUT_CSV)

## ---- chips ----------------------------------------------------------------------
TURBO <- viridisLite::turbo(17)
ramp  <- grDevices::colorRamp(TURBO)
chip <- function(v, at, anc, k, h, src) {
  if (is.na(v)) return('<span class="c na" title="not a known end of this model">&ndash;</span>')
  rgb <- ramp(min(max(v, 0), 1))
  ## turbo is dark at both ends and bright in the middle: text by luminance
  fg  <- if (sum(c(0.2126, 0.7152, 0.0722) * rgb / 255) < 0.45) "#fff" else "#0b0b0b"
  sprintf('<span class="c" style="background:%s;color:%s" title="%s %s: %.3f at anchor %d (end at %d, best within &plusmn;%d); %s">%.2f</span>',
          grDevices::rgb(rgb, maxColorValue = 255), fg, MODELS[[k]]$label, HEADS[[h]], v, at, anc, TOL, src, v)
}
esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
label <- function(r) {
  nm_t <- sprintf("%s %s", esc(r$name), r$term)
  if (r$set != "trained") return(sprintf('<span class="%s">%s</span>', if (grepl("peptide", r$set)) "ns" else "nl", nm_t))
  f <- unique(na.omit(c(r$fold_cur, r$fold_arm)))
  sprintf('<span class="tr">%s</span> <span class="lp">f%s</span>', nm_t, paste(f, collapse = "/"))
}
ord <- c("never seen: peptide" = 0, "never seen: large product" = 1, trained = 2)
ends <- ends[order(gene, ord[set], term, anchor)]
cards <- ends[, {
  rows <- vapply(seq_len(.N), function(i) {
    r <- .SD[i]
    src <- if (r$set != "trained") "all-50 ensemble" else "the fold that held this end out"
    ch <- unlist(lapply(names(MODELS), function(k) lapply(names(HEADS), function(h)
      chip(r[[sprintf("%s_%s", k, h)]], r[[sprintf("%s_%s_at", k, h)]], r$anchor, k, h, src))))
    paste0('<div class="nsr"><span class="nsp">', label(r), "</span>", paste(ch, collapse = ""), "</div>")
  }, "")
  .(html = paste0('<div class="ns"><div class="nsh"><span>end</span><span>cur<br>CNN</span><span>cur<br>xgb</span>',
                  '<span>arm<br>CNN</span><span>arm<br>xgb</span></div>', paste(rows, collapse = ""), "</div>"))
}, by = gene]

## ---- blurb, footer, index -----------------------------------------------------------
nsp  <- ends[set == "never seen: peptide"]
summ <- unlist(lapply(names(MODELS), function(k) lapply(names(HEADS), function(h) {
  x <- na.omit(nsp[[sprintf("%s_%s", k, h)]])
  sprintf("<b>%s %s</b>: median %.2f, &ge; 0.5 at %d of %d", esc(MODELS[[k]]$label), HEADS[[h]], median(x), sum(x >= 0.5), length(x))
})))
plm  <- regmatches(arm_var(ARM, "PLM"), regexpr("[^/]+\\.parquet", arm_var(ARM, "PLM")))
drop <- arm_var(ARM, "DROP"); flags <- arm_var(ARM, "INPUT_FLAGS")
## where stage 5 writes this arm's metrics (pepend_arm_pipeline.sh's HEADCMP):
## lf_headcmp/unseen<sfx> for an arm on its own windows, else lf_headcmp
hc_rel <- if (nzchar(SFX)) paste0("lf_headcmp/unseen", SFX) else "lf_headcmp"
css <- paste0('<style>.ns{margin-top:7px;padding-top:6px;border-top:1px solid var(--grid);font-size:11px}',
  '.nsh,.nsr{display:grid;grid-template-columns:1fr 34px 34px 34px 34px;align-items:center;column-gap:3px}',
  '.nsh{color:var(--ink2);font-size:9.5px;line-height:1.15;margin-bottom:2px}.nsh span+span{text-align:center}',
  '.nsp{font-weight:600;overflow-wrap:anywhere;line-height:1.25}.nsr{margin-top:2px}',
  '.ns .ns,span.ns{color:#B0177E;font-weight:600}.nl,span.nl{color:#6A3D9A;font-weight:600}.tr{color:var(--ink)}',
  '.lp{color:var(--ink2);font-weight:400}',
  '.c{display:inline-block;text-align:center;border-radius:3px;padding:0 2px;font-size:10.5px;font-weight:650;line-height:16px;font-variant-numeric:tabular-nums}',
  '.c.na{background:transparent;color:var(--ink2);font-weight:400}',
  '.grid{grid-template-columns:repeat(auto-fill,minmax(250px,1fr))!important}</style>')
blurb <- paste0(css,
  sprintf('<p>%d of these pages carry known ends. Each page draws <b>two models</b>, each with both window heads (CNN linear head and xgb window head; 5 folds &times; 10 seeds per terminus, each fold a thin series, the all-50 ensemble the heavy line with squares):</p>', uniqueN(ends$gene)),
  '<ul style="margin:4px 0 8px 18px;padding:0">',
  '<li><b>curated</b> &mdash; the original 26-channel curated model on the original training windows.</li>',
  sprintf('<li><b>%s</b> (<code>%s</code>): the curated channels%s%s; trained on the <code>%s</code> windows%s.</li></ul>',
          esc(LAB), ARM, if (nzchar(drop)) sprintf(" minus <code>%s</code>", gsub(",", ", ", drop)) else "",
          if (length(plm)) sprintf(", plus the per-residue channels of <code>%s</code>", plm) else "",
          if (nzchar(SFX)) SFX else "original", if (nzchar(flags)) sprintf(" (<code>%s</code>)", flags) else ""),
  '<p>So the two tracks can differ in channels AND in training windows at once; an arm-vs-arm comparison on the same windows is in <code>', hc_rel, '/final_metrics.csv</code>.</p>',
  sprintf('<p><b>Chips</b>: every known end on the gene, with the score each model gives it &mdash; the highest score within &plusmn;%d residues of the end. <span class="ns">Magenta</span>: never-seen peptide ends; <span class="nl">purple</span>: never-seen large cleavage products &mdash; both scored by the all-50 ensemble (the heavy squares). Black with <span class="lp">f<i>k</i></span>: training ends, scored by the fold that held each one out in that model (the triangles); an anchor in the window that the scoring members trained on is skipped. &ndash; = not a known end of that model. Hover a chip for the exact value and anchor.</p>', TOL),
  sprintf('<div style="display:flex;align-items:center;gap:8px;margin:6px 0 4px;font-size:11px">score <span>0</span><span style="display:inline-block;width:220px;height:12px;border-radius:3px;background:linear-gradient(90deg,%s)"></span><span>1</span> <span style="color:var(--ink2)">(turbo)</span></div>', paste(TURBO, collapse = ",")),
  sprintf('<p style="font-size:11.5px">Never-seen peptide ends on these pages (%d): %s.</p>', nrow(nsp), paste(summ, collapse = "; ")),
  '<p style="font-size:11.5px"><b>Compare chips down a column, not across.</b> Raw scores are not calibrated between models: a window fix that stops padding counting as evidence lowers an arm\'s scores everywhere, so a lower chip is not by itself a worse call. Ranking quality over the whole proteome: <code>', hc_rel, '/final_metrics.csv</code>.</p>')
footer <- sprintf(paste0(
  'Window models: curated = <code>lf_pepend_run_{N,C}_cv5_f{1..5}_nomask_flatten</code>; %s = ',
  '<code>lf_pepend_run_{N,C}_cv5_f{1..5}_nomask_flatten_%s</code> (10 seeds each, pooled as cv5all50). ',
  'Chip scores from <code>lf_winxgb_scan_{N,C}_{cv5all50,cv5_f&lt;k&gt;}_nomask_flatten[_%s].npz</code>; ',
  'ends from <code>cv5_unseen_ends%s.csv</code> and <code>lf_pepend_cv5_knowns_{N,C}[%s].csv</code>; per-end table <code>%s</code>. ',
  'Chips by <code>end_score_cards.R</code>.'), esc(LAB), ARM, ARM, SFX, SFX, basename(OUT_CSV))
tmp <- tempfile("end_score_cards_"); dir.create(tmp)
saveRDS(as.data.frame(cards), file.path(tmp, "cards.rds"))
writeLines(blurb, file.path(tmp, "blurb.html")); writeLines(footer, file.path(tmp, "footer.html"))
rc <- system2("/usr/local/bin/Rscript", shQuote(c(lf_src("inst/scripts/visualization/make_gallery_index.R"),
  "--dir", GAL, "--out", INDEX, "--title", TITLE, "--no-stock-desc",
  "--blurb", file.path(tmp, "blurb.html"), "--footer", file.path(tmp, "footer.html"),
  "--card-html", file.path(tmp, "cards.rds"))), env = "LC_ALL=en_US.UTF-8")
if (rc != 0) stop("make_gallery_index.R failed", call. = FALSE)
message(sprintf("%d ends on %d of %d pages (%d never-seen peptide ends); chips in %s, table %s",
                nrow(ends), uniqueN(ends$gene), length(genes), nrow(nsp), INDEX, OUT_CSV))
for (s in summ) message("  ", gsub("<[^>]+>|&ge;", "", s))
