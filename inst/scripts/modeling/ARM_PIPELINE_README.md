# The arm pipeline, and how to add an arm

Updated 2026-10-25; arms on their own training windows (2026-11-10/12) are the
last section. **The defaults are current practice**, so the bare command is the
right command. This supersedes the "how to run it" half of
`ESMONLY36_README.md`, which stays as the record of the ESM-C arm's results and
its five traps.

## What an arm is

An **arm** is a channel set. Everything else — model names, residue table, scan
npz, tracks, result CSVs, plots, page directory, attribution CSVs — is derived
from `$ARM`, so a new arm collides with nothing already on disk.

An arm lives in **one file**, `arms/<ARM>.sh`:

| variable | what it is |
|---|---|
| `ARM_LABEL` | the human name, used in banners and figure legends |
| `PLM` | the extra per-residue parquet (`accession, index, AA, plm_*`) |
| `DROP` | curated channels to drop, comma-separated; `""` to keep all 26 |
| `N_CHAN_EXPECT` | channels the trunk must end up with — asserted after training |
| `THRESH_DEFAULT` | the arm's segment threshold (a property of the arm, not of the script) |
| `mk_plm()` | how to build `$PLM` if it is missing |
| `mk_resid()` | how to build the arm's residue table |
| `INPUT_SFX` | *optional*: train on `lf_pepend_nn_input_cv5_f<k>${INPUT_SFX}.rds` instead of the shared windows (last section) |
| `INPUT_FLAGS` | *optional*: `key=value ...` that every one of those inputs' `params` must hold — asserted before training |
| `CONT_EXPECT` | *optional*: channels the trunks must noise-augment — checked in the first trunk's `config.json` |

### The arms that exist

| arm | channels | what it is |
|---|---|---|
| *(the curated model)* | 26 | the baseline; not an arm file, it predates the mechanism |
| `esmonly36` | 37 | `padding` + the 36-component PCA of ESM-C; all 25 curated channels dropped |
| `curated_wAAid` | 47 | all 26 curated **plus** a 21-letter one-hot of the residue's own identity |
| `curated_AAswap` | 43 | the four AA property channels **replaced** by that one-hot |
| `motifind` | 29 | all 26 curated plus three identity indicators: G, basic (K/R), P |
| `msa` | 43 | `curated_AAswap` with `cons_rs_n` replaced by `plm_cons_msa`, the only conservation channel; on the `_msa` windows (all four motif top-ups) and the v2 features (below). Kevin, 2026-11-12 |
| `msa_gap` | 43 | `msa` exactly, plus a `gap` class (7 classes); its `_msagap` windows are `_msa`'s with gap labels, nothing else. Evaluate it on `_msa`'s never-seen set (`UNSEEN_SFX=_msa`) |

`curated_wAAid` is the file-name token for what Kevin called **curated-wAAid**;
the display label is `curated + AA one-hot`. Underscore rather than hyphen
because the token also lands in R env vars and `paste0`-built file names.

## Running it

```bash
export LF_SRC=~/R_projects/ligandFinder      # the code: a checkout
export LF_RUN=~/AF2_analysis/lf_run          # its lists, blurbs and logs (this is the default)
P=$LF_SRC/inst/scripts/modeling/pepend_arm_pipeline.sh
cd $LF_RUN

zsh $P                                       # the default arm, all stages
ARM=esmonly36 zsh $P                         # another arm
STAGES=5,6 ARM=curated_wAAid zsh $P          # re-evaluate and re-plot
STAGES=shap zsh $P                           # just the attribution
```

Stages take a number **or** a name:

| # | name | what it does | resumable |
|---|---|---|---|
| 0 | `channels` | build the arm's extra per-residue parquet | skipped if it exists |
| 1 | `trunks` | train 5 folds × 2 termini, 10 seeds each | skipped per trunk |
| 2 | `residues` | the arm's residue table **+ the order check** | skipped if it exists |
| 2b | `groups` | an `INPUT_SFX` arm's own window groups **+ the row check**; other arms skip it | skipped per npz |
| 3 | `scans` | full-proteome scans, 10 of them | skipped per npz |
| 4 | `tracks` | assemble the tracks rds | skipped if it exists |
| 5 | `eval` | segment confusion/ROC, by motif, anchor ranking metrics | always re-runs |
| 6 | `plots` | threshold scan, motif table, metrics | always re-runs |
| 7 | `pages` | per-gene scans, known-peptide pages, gallery index | skipped per page |
| 8 | `shap` | per-residue-position attribution + its figures | always re-runs |

Stages 5, 6 and 8 are cheap and always re-run, which is what you want when only
`THRESH` or `TOL` changed.

## To add an arm

1. `cp arms/curated_wAAid.sh arms/<new>.sh`
2. Change `PLM`, `DROP`, `N_CHAN_EXPECT`, `THRESH_DEFAULT` and the two `mk_`
   functions.
3. `ARM=<new> ./pepend_arm_pipeline.sh`

Nothing else in this directory changes. The generalised runners are
`train_arm.sh` and `fullscan_arm.sh`; the older `train_esmonly36.sh` /
`fullscan_esmonly36.sh` remain on disk and still work, and `arms/esmonly36.sh`
reproduces what they did.

## Adding extra per-residue channels: use the `--plm` parquet

`10_8_pepend_train.R --plm <parquet>` is the project's existing, tested route for
"extra per-residue channels". `lf_plm_attach` slots them in beside the 26, and
`lf_winxgb.plm_residues` appends the same columns to the residue table the
exhaustive scan slices. Writing a new channel set **in that shape** means
neither the training script, nor the scan, nor the window builder needs an edit
— and `10_7_pepend_windows.R`, which builds the *shared* residue cache, is never
touched.

Requirements on the parquet:

* columns `accession`, `index`, `AA`, then the channels;
* channel names **must** start with `plm_`, because `lf_plm_read` selects them
  with `grep("^plm_")`. `plm_AA_R` is fine and self-documenting;
* one row per (accession, residue index). Residues absent from the table get
  zeros, which is the policy 10_8 already applies at training time.

`lf_plm_attach` checks each window's residue letters against the table's `AA`
column and errors on a mismatch, so a table cut from a different sequence
version cannot silently misalign. That check is worth having: it is a second,
independent confirmation that a new channel set is on the right residues.

## The five traps, still true

1. **The residue table is positional.** `lf_dcnn.scan.build_windows` does
   `x[ok] = res.feat[rows[ok]]` into a `(k, seq_len, n_channels)` array. A
   wrong-width table raises; a right-width table in the wrong order scans
   silently against the wrong channels and produces plausible nonsense.
   **`check_resid_order.py` is the guard**, and it checks *data*, not names:
   every curated channel must equal the curated table's own column of that
   name, and every `plm_*` channel must equal the parquet's column at that
   residue, re-joined independently on a sample of rows. It passes for both
   constructions in use ("26 curated + plm" and "padding + plm"). Stage 2 and
   stage 3 both run it.
2. **`10_8` caches by `--plm`, not by `--drop-channels`.** Two arms that differ
   only in dropped channels share a cache, and `lf_dcnn_run_isolated` reuses any
   `member_*.npz` already present on its presence alone — so one trains and
   reports as both. Always pass an explicit, distinct `--out`; `train_arm.sh`
   puts `$ARM` in the name for exactly this reason.
3. **Held-out triangles match a series BY NAME.** `lf_track_heldout_marks` emits
   `"C fold 1"`, which only matches the page's own series. Extra series need
   `options(lf.heldout_extra = ...)` with their own names. A held-out fold draws
   triangles along its *whole* series, not just at the end residue.
4. **ggiraph builds its SVG in JavaScript at load time.** Grepping the HTML for
   markers finds nothing — the raw file has zero `<polygon>` elements. Measure
   in the live DOM, or you will conclude a working change did nothing.
5. **Never install into `~/.virtualenvs/r-tensorflow`.** Installing matplotlib
   there once upgraded numpy and broke a running render. matplotlib lives in
   `~/.virtualenvs/lf-plm` (3.8.3); use that for every plot script.

A sixth, for `train_arm.sh` and `fullscan_arm.sh`: they re-enter themselves
through `xargs`, so **`ARM` is exported**. A child that fell back to the default
would train the wrong arm into the right file name.

## Which venv

| job | interpreter |
|---|---|
| embedding, PCA, every `*_plot.py`, `check_resid_order.py`, `aa_onehot_build.py` | `~/.virtualenvs/lf-plm/bin/python` |
| training, scanning, `shap_positions.py`, `shap_xgb_head.py` | `~/.virtualenvs/r-tensorflow/bin/python` |
| R | `/usr/local/bin/Rscript` with `LC_ALL=en_US.UTF-8` |

`metrics_final.R` takes the scripts directory as its **first positional
argument** — not optional; without it the run dies on `NA/cv5_unseen_ends.csv`.

## Attribution (stage 8)

Two scripts, because the two heads see different things and only one of them
sees residue positions at all.

* **`shap_positions.py`** — the flatten head reads the `(36 × n_channels)`
  window directly, so attribution is per residue position in the input.
  Estimator: **expected gradients** on the *pre-sigmoid* score, which is the
  Aumann–Shapley value of each `(position, channel)` cell with the arm's own
  candidate pool as the background distribution. This is what the `shap`
  package's `GradientExplainer` computes; it is implemented directly because
  `shap` is in neither venv and nothing is ever installed into r-tensorflow.
  The sigmoid is skipped deliberately — its gradient vanishes exactly where the
  model is confident. Completeness (`sum(phi)` vs `logit(x) − mean logit(bg)`)
  is reported per fold, not assumed: ~0.998 correlation at 64 background draws.
  Outputs `shap_position*.csv`, `shap_channel*.csv`, `shap_pos_family*.csv`,
  `shap_pos_letter*.csv`, `shap_completeness*.csv`.
* **`shap_xgb_head.py`** — the xgb head does **not** see the window; it sees 48
  pooled summaries of the trunk's per-index surface. So per-residue-position
  attribution is undefined for its inputs, and four of its eight statistics are
  positional instead (`at_anchor`, `near_anchor`, `pre_anchor`, `post_anchor`),
  which is the honest per-region answer. **Exact TreeSHAP**
  (`pred_contribs=True`), with the `sum(phi) == margin` identity asserted. The
  head is refit with `full_scan.py`'s own defaults because it is never saved to
  disk.

The per-position axis is an **offset from the peptide end**, not a 0–35 column
index: `PEPEND_ANCHOR` is 28 for C and 8 for N (1-based), so offset 0 is the end
residue, negative is inside the peptide and positive is past the end.

The explain set is **held out**: for fold *k* the `val` split is that fold's own
fifth of the knowns (16–18 ends per fold per terminus, the rest of `val` being
~400 deliberately hard negatives), so every positive is explained by ten members
that never trained on it, and the five folds together cover every known end once.

## The gallery index is package code now (2026-10-26)

The TP/FN index used to be a script here that had to be rewritten for each new
arm. It now lives in the package:

| file | what it is |
|---|---|
| `R/gallery_tpfn.R` | `lf_gallery_tpfn_index()`, `lf_tpfn_arm()`, `lf_tpfn_default_arms()`, `lf_tpfn_pick_head()`, `lf_tpfn_cards()`, `lf_tpfn_blurb()`, `lf_tpfn_fold_notes()` |
| `inst/scripts/visualization/make_tpfn_index.R` | its command-line face |

**An arm is data.** To add one, add an `lf_tpfn_arm()` entry to
`lf_tpfn_default_arms()` — key, label, file-name suffix, threshold — and
nothing else changes. Stage 7 of this pipeline calls the package script, so a
new arm appears in the index without touching the pipeline either.

```bash
# the usual: three arms, the defaults
Rscript $W/inst/scripts/visualization/make_tpfn_index.R --dir ~/AF2_analysis/pepend_pages_arms \
  --blurb final_blurb.html --footer final_footer.html

# or an explicit arm set, key:label:suffix:threshold ("-" = the empty suffix)
Rscript $W/inst/scripts/visualization/make_tpfn_index.R --dir <gallery> \
  --arms 'cur:curated:-:0.2,waa:curated+AAid:_curated_wAAid:0.5'
```

`index_tpfn_arms.R` in this directory is now a **shim** that forwards to it, so
the old command still works. `index_tpfn_esm.R` stays as the record of what
produced the two-arm ESM-C gallery.

### What it does not do

It does not compute the calls. `end_calls_T<thr><arm>_tol<tol>.csv` and
`segment_confusion<arm>_tol<tol>.csv` must already exist under `--headcmp`, one
pair per arm; a missing one is an error naming the command that makes it.

### Two traps it now encodes

1. **The head is chosen on `auc_T`, not on TP count.** Both heads of an arm sit
   at the same threshold, so TP alone rewards whichever head calls *more
   segments*: at T = 0.2 the curated xgb head finds exactly the same 26
   never-seen peptide ends as the linear head while calling 9,459
   false-positive segments against 6,991. `auc_T` = (1 + TPR − FPR)/2 charges
   for both.
2. **The confusion table has FOUR rows per (set, threshold)** — two termini ×
   two heads — and they are averaged. Reading one row picks a head on half the
   evidence, and the termini do disagree; `lf_tpfn_pick_head()` asserts the
   shape so that cannot happen silently.

A third, for anyone editing `R/gallery_tpfn.R`: inside `[.data.table` a bare
symbol resolves to a **column** of that name before a local variable, so an
argument named `set` silently compares the column with itself and matches every
row. The deciding set's argument is called `eval_set`, and masks are computed
outside the `[`.

## Arms on their own windows (2026-11-10/12)

**Kevin, 2026-11-10,** suspected the negatives lacked windows anchored at a
precursor terminus or at a monobasic site, so that the model learned "padding
present = peptide". It had: ~45% of positive windows carry padding against ~20%
of negatives, and the AAswap trunks scored 17% of unlabelled extracellular
anchors AT a precursor terminus > 0.5, against 0.1–0.2% of interior ones. He
asked, the same day, to replace Aminode conservation (`cons_rs_n`) with
conservation scored from each precursor's MSA (`inst/scripts/generate_residue_db/msa_cons/`;
AFDB's MSAs are 403, so ColabFold's).

The first two such arms, `curated_AAswap_padg` (the windows alone) and
`curated_AAswap_msacons_padg` (the windows and MSA conservation), were trained
on the `_padg` windows on 2026-11-10/11; their arm files are gone (History,
below), their trunks, scans, tables and gallery remain in `~/AF2_analysis`.

### The windows: `10_7` flags

All off by default; a build without them reproduces the old windows, their data
and their labels exactly (checked on fold 1).

| flag | `params` | what it does |
|---|---|---|
| `--pad-match` | `pad_match` | negatives topped up, per split, until each padding-profile bin is as common among them as among the positives; a negative's padded positions are labelled `padding`, as a positive's are |
| `--mono-match` | `mono_match` | the same for monobasic-anchored windows |
| `--db-match` | `db_match` | the same for dibasic-anchored windows (Kevin, 2026-11-12) |
| `--amid-match` | `amid_match` | the same for amidation (`G + dibasic`) windows; the pool is small (~1,100 non-ligand sites), so it is split by precursor between train and val |
| `--gap-class` | `gap_class` | the residues between a known end and its nearest dibasic site (`lf_pepend_gap()`: usually the amidation G) are `gap`, not cleavage context, and the vocabulary is `LF_PEPEND_CLASSES_GAP`; labels only, windows unchanged |
| `--feature-cache <rds>` | `dssp_index` | the per-residue feature cache to read (default the original); with `--rebuild-features` it is rebuilt there. A rebuild places each `secretome_aa` row at its UniProt position, not at the AlphaFold model's `index` (2026-11-12; `--legacy-dssp-index` for the old join) |

The motif top-ups run three rounds, motifs then padding; matching is close, not
exact. About 92% of positives sit at a motif or a precursor terminus and their
categories overlap, while one category's top-ups do not; on fold 1 the train
negatives land within 2-4 points of the positives (amidation short where the pool
runs out). Exact matching would take stratified top-ups and either ~10k negatives
per split or far fewer generic ones.

**The v2 features** (`~/AF2_analysis/lf_pepend_residue_cache_v2.rds`,
`lf_pepend_residues_v2.npz`, `lf_plm/aa_onehot_msacons_v2.parquet`): the rebuild
above. 215 precursors whose AlphaFold model was built on another sequence are
remapped; 53 gain 15,062 residues (AGT 29 -> 476 of 476), every other precursor
is bit-identical. An arm on them names its curated table with `RESID_BASE` (stage 2
checks the arm's table against it); the MSA channel is re-masked on it at the
original pinned scale. `run_msa_arms.sh` builds all of it, then trains, scans and
evaluates `msa` and `msa_gap`.

Top-up negatives come only from anchors whose window overlaps no annotated end.
The five fold inputs of a window set are built from the checkout root, e.g.

```bash
for k in 1 2 3 4 5; do
  LC_ALL=en_US.UTF-8 Rscript inst/scripts/modeling/10_7_pepend_windows.R --n-folds 5 --fold $k \
    --near-n 100 --seed 42 --pad-match --mono-match \
    --out ~/AF2_analysis/lf_pepend_nn_input_cv5_f${k}<sfx>.rds
done
```

and the arm that trains on them says so with `INPUT_SFX=<sfx>` and the
`INPUT_FLAGS` its inputs must carry (`arms/msa.sh` lists its own build).

### What `INPUT_SFX` changes downstream

| name | set by | effect |
|---|---|---|
| `INPUT_SFX` | the arm file | trunks train on `lf_pepend_nn_input_cv5_f<k>${INPUT_SFX}.rds`; every scan, per-gene scan and xgb attribution reads `lf_pepend_window_groups_<T>_cv5_f<k>${INPUT_SFX}.npz` |
| `UNSEEN_SFX` | env, default `INPUT_SFX` | the never-seen set, `cv5_unseen_ends${UNSEEN_SFX}.csv`: the old set minus every end whose window a negative of those inputs overlaps, written by `INPUT_SFX=<sfx> Rscript cv5_unseen_sfx.R` (run from the scripts directory; `OUT` overrides the path). It stops unless the ends those inputs hold out, and in which fold, are `lf_pepend_cv5_knowns_<T>.csv`'s (`KNOWNS_DIFFER_OK=1` writes the set anyway; `_padg` needs it) |
| `HEADCMP` | the pipeline | where stages 5–6 write: `lf_headcmp`, or `lf_headcmp/unseen${UNSEEN_SFX}`, so tables on another end set never overwrite the shared ones |

`UNSEEN_SFX=<sfx>` on an old arm puts it on the same footing for a paired
comparison (`UNSEEN_SFX=_padg STAGES=5,6 ARM=curated_AAswap ./pepend_arm_pipeline.sh`).
The paired tables of 2026-11-11 are in `~/AF2_analysis/lf_headcmp/unseen_padg/`
(written before 2026-11-12: History, below).

**Stage 2b, `groups`.** `full_scan`, `scan_genes` and `shap_xgb_head.py` take the
xgb head's inner split from a groups npz that must be the run's own training
windows row for row; a same-length wrong file misaligns silently. An
`INPUT_SFX` arm exports its own with `10_9r_export_window_groups.R`, from the
inputs it trained on and with an explicit `--out` (the default would overwrite
the shared unsuffixed file), then checks all 30 splits against the trunks'
`arrays.npz` (`__known` == `y_global`). `full_scan.py` and `scan_genes.py` now
assert the same, and `fullscan_arm.sh` and `cv5_post_esmonly36.sh` will not
start without all ten.

**Every arm is scored at the anchors the shared end tables carry** — the cv5
knowns, `cv5_unseen_ends*.csv`, `ligand_list`, UniProt peptides — by
`segments_roc.R`, `segments_motif.R`, `end_calls.R` and `metrics_final.R`, and
its held-out marks are drawn there by `final_pages_arms.R`, each end held out by
the fold those tables give it. That is exact only for a window set whose
held-out ends and folds ARE the tables' — `cv5_unseen_sfx.R` checks it; the
`_msa` and `_msagap` sets pass, and `_padg` (30 C ends anchored one past the
mature end, GLP-1's two C ends merged) does not, so its arms are scored one
residue off on those ends, as accepted on 2026-11-12. Two helpers read an
arm's own settings and windows: `arm_vars.R` (`arm_var(arm, var)`: one setting
out of `arms/<arm>.sh`, read as text, so the arm's training guards do not run;
`end_score_cards.R`, `terminus_bias.R`), and `arm_knowns_cache.R`
(`INPUT_SFX=<sfx>`: that window set's known ends with the fold that holds each
out, `lf_pepend_cv5_knowns_<T><sfx>.csv`, for `end_score_cards.R`; the inputs
are 92 MB each, and rebuilt inputs need the cache rebuilt).

### Pages

`cv5_post_esmonly36.sh` sources the arm file for its groups and names its logs
per arm. `esm_pages.sh` READS `INPUT_SFX` out of the arm file with sed —
sourcing would run the arm's training guards in a render — and EXPORTS
`UNSEEN_SFX`, or the xargs children would mark never-seen ends from the full
set. An arm with its own windows is drawn alone beside the curated baseline;
`final_pages_arms.R` needs an `ARM_REG` row for it. A gallery for such an arm:

```bash
PLM=<the arm's PLM parquet> LC_ALL=en_US.UTF-8 Rscript generate_residue_db/plm_pagegenes.R  # -> *_pagegenes.rds
ARM=<arm> PROCS=3 modeling/cv5_post_esmonly36.sh
ARM=<arm> PROCS=3 TITLE="..." visualization/esm_pages.sh
```

The MSA-conservation gallery of 2026-11-11
(`~/AF2_analysis/pepend_pages_curated_AAswap_msacons_padg/`, 104 pages) was
built this way. `msacons_gallery_qc.py` then checks every page's live CNN
series, both arms, against the per-gene scans (not the held-out marks).

Not stage 7: after the pages it rewrites `index.html` as the TP/FN index of
`lf_tpfn_default_arms()` over the tables in `lf_headcmp`, and neither holds an
arm on its own windows.

### The source tree, and the capability guard

Every runner sources `inst/scripts/lf_paths.sh`, which finds the source tree
(`LF_SRC`, or the checkout the script sits in) and names the folders the stages
live in (`$LF_MOD`, `$LF_EVAL`, `$LF_VIZ`, `$LF_RDB`); each child (the xargs
re-entries, the stage scripts) resolves it again. The data the runners read,
`page_genes.txt`, the blurb and footer html and (stage 5 hands it to
`metrics_final.R`) `cv5_unseen_ends*.csv`, and the logs they write live in
`$LF_RUN` (default `~/AF2_analysis/lf_run`). Until 2026-11 both sat in a working
copy of the scripts in `~/AF2_analysis/lf_nomask_scripts_20261021`; nothing runs
from a copy now.

An arm file may refuse a source tree that predates what it needs, testing the
capability, not the checkout's name: `arms/msa.sh` and `arms/msa_gap.sh` want
`plm_cons_msa` inside the `DEFAULT_CONT_CHANNELS` of
`$W/inst/python/lf_dcnn/config.py`. Without it the channel is not
noise-augmented as `cons_rs_n` was, and the comparison with the baseline
differs in one more way than intended. The arm's `CONT_EXPECT` checks the same
thing again in the trunks' own `config.json`.

### The MSA channel's pinned scale

Kevin, 2026-11-10: "train both arms now, scan after the API finishes". The
`curated_AAswap_msacons_padg` trunks trained on a partial build of the channel:
176 precursors still had no MSA and an all-zero channel. None holds a known
peptide, but they supply 4.5–5.4% of each fold's training negatives, a weak
"zero conservation = negative" cue open to that arm only; retraining would not
reproduce them. Once all 388 API alignments were on disk the channel was
rescored and rebuilt with the raw-JSD min/max PINNED to the build those trunks
trained on (`build_msacons_parquet.py --scale-from`), and
`check_msacons_rebuild.py` confirmed that only the precursors that had had no
MSA changed. That pinned scale (`aa_onehot_msacons.parquet.scale.json`) is the
one later arms reuse (`arms/msa.sh`).

## History

- 2026-11-12: the G-anchored C end convention (--amid-g-anchor, arms curated_AAswap_padg / curated_AAswap_msacons_padg) was removed; their artifacts remain in ~/AF2_analysis and now evaluate at mature anchors.
