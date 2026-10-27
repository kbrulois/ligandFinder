# The arm pipeline, and how to add an arm

Updated 2026-10-25. **The defaults are current practice**, so the bare command is
the right command. This supersedes the "how to run it" half of
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

### The arms that exist

| arm | channels | what it is |
|---|---|---|
| *(the curated model)* | 26 | the baseline; not an arm file, it predates the mechanism |
| `esmonly36` | 37 | `padding` + the 36-component PCA of ESM-C; all 25 curated channels dropped |
| `curated_wAAid` | 47 | all 26 curated **plus** a 21-letter one-hot of the residue's own identity |

`curated_wAAid` is the file-name token for what Kevin called **curated-wAAid**;
the display label is `curated + AA one-hot`. Underscore rather than hyphen
because the token also lands in R env vars and `paste0`-built file names.

## Running it

```bash
cd ~/AF2_analysis/lf_nomask_scripts_20261021

./pepend_arm_pipeline.sh                              # the default arm, all stages
ARM=esmonly36 ./pepend_arm_pipeline.sh                # another arm
STAGES=5,6 ARM=curated_wAAid ./pepend_arm_pipeline.sh # re-evaluate and re-plot
STAGES=shap ./pepend_arm_pipeline.sh                  # just the attribution
```

Stages take a number **or** a name:

| # | name | what it does | resumable |
|---|---|---|---|
| 0 | `channels` | build the arm's extra per-residue parquet | skipped if it exists |
| 1 | `trunks` | train 5 folds × 2 termini, 10 seeds each | skipped per trunk |
| 2 | `residues` | the arm's residue table **+ the order check** | skipped if it exists |
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
| `inst/scripts/make_tpfn_index.R` | its command-line face |

**An arm is data.** To add one, add an `lf_tpfn_arm()` entry to
`lf_tpfn_default_arms()` — key, label, file-name suffix, threshold — and
nothing else changes. Stage 7 of this pipeline calls the package script, so a
new arm appears in the index without touching the pipeline either.

```bash
# the usual: three arms, the defaults
Rscript $W/inst/scripts/make_tpfn_index.R --dir ~/AF2_analysis/pepend_pages_arms \
  --blurb final_blurb.html --footer final_footer.html

# or an explicit arm set, key:label:suffix:threshold ("-" = the empty suffix)
Rscript $W/inst/scripts/make_tpfn_index.R --dir <gallery> \
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
