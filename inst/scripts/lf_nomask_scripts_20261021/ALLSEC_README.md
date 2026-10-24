# All-secretome pages — overnight run (set up 2026-10-22 evening)

## What is running

```
nohup ./allsec_overnight.sh > allsec_overnight.log 2>&1 &
```

Two stages, **both resumable**. If it dies, the machine sleeps, or you kill it,
just run the same command again: finished scans and existing pages are skipped.

1. **`allsec_scan.sh`** — per-gene CNN + xgb window scans, 5 folds × 2 termini,
   4 at a time, then the all-50 pool. Output `lf_winxgb_scan_*_nomask_flatten_all.npz`.
   This had to run even though `lf_winxgb_fullscan_*.npz` already covers every
   anchor: the fullscan carries `xgb_pm_sd` only, and the pages' spread band
   needs `cnn__sd`, which only `scan_genes` emits.
2. **`allsec_pages.sh`** — one page per gene into **`~/AF2_analysis/pepend_pages_allsec/`**
   (a new directory; `pepend_pages_flatten_cv5/` is left alone as the reference),
   shortest precursor first, then the index.

## Verified identical to the existing gallery

- Re-scanning genes the 104-gene npz already had reproduced **every value
  bit-for-bit** — anchor, cnn, cnn__sd, xgb_pm, xgb_pm__sd, both termini, all folds.
- Re-rendering DEFB4A, TGFA and BRINP3 gave pages **content-identical** to the
  gallery. The only difference is the 18 random `htmlwidget-<hex>` / `svg_<uuid>`
  element IDs ggiraph mints on every render; normalise those and the files match
  exactly. Same `final_pages.R`, same options, same `LF_POPUP_MAX` (2000, and 600
  for ANO8/C3/C5/EGF/COPA as before).

## Scope: 5,176 pages, not 5,616

| | genes | why |
|---|---|---|
| precursors in the residue cache | 5,616 | |
| − duplicate gene names | 2 | COL18A1, OTOA have two accessions each and would collide on `<gene>.html` |
| − gene names containing `/` | 5 | e.g. `UNQ6190/PRO20217` — cannot be a filename |
| − no entry in `plot_bundle.rds` | 431 | `the_input` has 5,180 genes; `final_pages.R` stops on `all(want %in% gs)` |
| **renderable** | **5,176** | `renderable_genes.txt` |

The 7 held back are in `skipped_genes.txt`; the 431 would need the plot bundle
rebuilt, which is its own job.

## This will NOT finish in one night

Measured on 5 genes (1,316 aa) under load: ≈ **20 s per gene + 0.17 s per residue**.
Over 5,176 genes / 2,969,477 residues that is ≈ 170 core-hours, so

- scan ≈ 3.5–5 h, then
- render ≈ **40+ h at 4 processes**, ~28 h at 6.

Overnight you should expect the scan plus roughly **1,800–2,500 pages** — all of
the short precursors and a good part of the 200–500 aa band. Re-run the same
command to continue; nothing is recomputed.

## Disk

Pages run ≈ **60 KB per residue**, so the full set is ≈ **170 GB** (524 GB free
at launch). `allsec_overnight.sh` watches free space every 5 min and stops the
render cleanly if it drops under **60 GB**, rather than leaving truncated HTML.

Shortest-first ordering is deliberate: if it stops early the casualties are the
giant precursors, whose pages are too big to open anyway — **MUC16 (14,507 aa)
projects to ~850 MB of HTML**, and 130 genes are over 2,000 aa. Worth deciding
whether those should be rendered at a reduced popup cap, or skipped.

## The 19 giants — answering that question (2026-10-24)

**A reduced popup cap does not work.** Measured on ADAMTS10 (1,103 aa), same
code path, only `LF_POPUP_MAX` changed: 85.2 MB at cap 2000 (not binding, all
1,103 anchors kept) against **55.0 MB at cap 50**. That is 28.7 KB per
suppressed anchor — real, but it leaves a 55 MB floor, 49.9 MB per 1,000
residues, which is ~720 MB for MUC16 before anything else.

**Most of a page is the species-alignment panel.** Counted in a loaded browser,
CP (1,065 aa, 87.6 MB) has 315,878 SVG elements and **273,373 of them — 87% —
are that one panel**: a `<rect>` and a `<text>` per species per residue over the
63 species in `species_dat`. 70,326 of its 70,597 text elements are single
amino-acid letters, i.e. the sequence drawn 66 times over. By bytes it is a
steady ~30% of a page in every length band. At 14,507 residues a residue is far
under a pixel wide, so none of it is readable at any zoom the page offers.

So the giants are rendered by **`allsec_giants.sh`** with three levers, all in
`final_pages_giant.R` and all off by default:

| env var | what it does |
|---|---|
| `LF_NO_ALIGN_MIN` | above N residues, drop the species-alignment panel |
| `LF_ENS_ONLY_MIN` | above N residues, keep only the two pooled all-50 ensembles, not the ten per-fold series |
| `LF_POPUP_MAX` | bound the popup payload (existing option) |

Dropping the alignment takes the same route as the 390 genes that have no
ortholog MSA: `plot_proteins.R` wraps it in `tryCatch` and logs "species
alignment SKIPPED". Verified on CCN2 that such a page is otherwise complete and
its nine widgets still share one width and one left edge, so the residue axis
still lines up.

`final_pages_giant.R` also **overrides `lf_plot_win_legend`** when the ens-only
mode is on. That function picks its keys from the `validation` flag, not from
which series exist, so an ens-only page would otherwise still advertise
"individual fold (mean of 10 seeds)" with no per-fold series drawn. The
override is a copy of the original with one key instead of two, installed after
`R/pepend_tracks.R` is sourced — deliberately here and not in the worktree,
which the overnight render is reading. Verified on MUSTN1 forced into ens-only
mode: the stale key is gone, "mean of all folds" remains.

**`final_pages_giant.R` is a copy of `final_pages.R`, and has to be.** R reads a
script incrementally as it runs, so editing `final_pages.R` while the overnight
render is using it kills that run at the end. The copy differs only by the two
guarded blocks above; with both env vars unset it behaves identically.

The 19 giants were chunk `c_hy`. That chunk was **emptied while still pending**
and its gene list saved to `giant_genes.txt` — otherwise the main run would have
rendered them at the default settings, and being one chunk it would have run
single-threaded on one worker for ~2 h while the other four idled.

**Status, 2026-10-24: the outcome is not in yet.** The 19 giants are rendering
as this is committed, so nothing above reports what they actually weigh or
whether they open — the largest two or three may still be too heavy, sizes
pending. `allsec_giants.sh` and `final_pages_giant.R` themselves are not in the
repository yet either; they live in `~/AF2_analysis/lf_nomask_scripts_20261021`,
which is not version controlled, which is why this section was committed ahead
of them. The measured per-gene table and the verdict replace this paragraph when
those renders land.

## Checking on it

```bash
tail -20 ~/AF2_analysis/lf_nomask_scripts_20261021/allsec_overnight.log
ls ~/AF2_analysis/pepend_pages_allsec/*.html | wc -l
du -sh ~/AF2_analysis/pepend_pages_allsec
```

## Where the page code comes from

`final_pages.R` sources `R/pepend_tracks.R` and `R/plot_proteins_win_new.R`, and
sets the python path to `inst/python`, out of a ligandFinder **source checkout** —
not out of the installed package, whose `inst/python` has neither `lf_winxgb` nor
`lf_resid` and whose plotting code is weeks behind. `lf_paths.R` / `lf_paths.sh`
resolve that checkout in one place for every script here:

1. `$LF_SRC`, if set
2. otherwise the checkout this script directory sits in — these scripts are
   committed at `inst/scripts/lf_nomask_scripts_20261021/`, so a copy that
   travels with the repository needs no configuration

A working copy kept **outside** a checkout (the one in `~/AF2_analysis`, say)
therefore needs `LF_SRC`, and says so if it is missing:

```bash
export LF_SRC=~/R_projects/ligandFinder
```

### historical note

Until 2026-10-24 these scripts each hardcoded the absolute path of the
`claude/class-position-mask-diagram-b398c7` worktree, so that worktree had to be
kept alive for the duration of a render. That is no longer true, and the branch
and worktree can be removed once nothing is running from a copy that still
points at them:

```bash
git push origin --delete claude/class-position-mask-diagram-b398c7
```
