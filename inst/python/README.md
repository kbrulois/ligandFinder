# lf_dcnn

The 1D-CNN that scores 36-residue peptide-cleavage windows, ported from
`inst/scripts/10_1dcnn_new6.R`. Two models — terminus `N` and terminus `C` —
each carry a **global** ranking head and a supervised **per-index** class head.

Everything upstream (data prep in `9.2_add_contact_data.R`) and everything
downstream (ranks, the amidation-motif block, nearest-known retrieval, plotting)
stays in R. This package owns the modeling only.

## Install / run

No install step: the package is plain Python and lives in `inst/python`. It
targets the same virtualenv R's `keras3` already uses, which is what makes the
reticulate bridge seamless.

```bash
cd inst/python && PYTHONPATH=. /Users/kbrulois/.virtualenvs/r-tensorflow/bin/python -m lf_dcnn selftest
```

## Two entry points

**In-process, via reticulate.** `R/dcnn_bridge.R` does the whole handoff; numpy
arrays pass in memory, nothing touches disk.

```r
source("R/dcnn_bridge.R")
res <- lf_dcnn_run(nn_input, all_params3)      # after 9.2 has built nn_input
nn_input_comb$pred     <- res$pred
nn_input_comb$pred_raw <- res$pred_raw
```

**Standalone, from disk.** Runnable and testable with no R at all.

```bash
python -m lf_dcnn train --input-dir /path/to/in --output-dir /path/to/out
```

`--input-dir` holds `arrays.npz`, `config.json` and an optional `meta.parquet`;
the run writes `outputs.npz`, `predictions.parquet` and `history.json`. Write
that directory from R with `lf_dcnn_export()`, or from Python with
`lf_dcnn.io.save_inputs()`.

## Array contract

Per term (`"N"`, `"C"`) and per split (`train`, `val`, `all`):

| name | shape | contents |
|---|---|---|
| `x` | `(n, 36, C)` | float32 input channels |
| `y_global` | `(n, 1)` | float32 0/1 window label |
| `y_per_index_cat` | `(n, 36, K_pi)` | cols `0:K_cat` are the `[6 real, none]` one-hot, col `K_pi-1` the padding flag |

**Do the `aperm` on the R side.** R is column-major and numpy is row-major, so
handing over the pre-`aperm` `(36, C, n)` array transposes silently and still
trains, just badly. `TermArrays.validate()` catches the shape half of that.

Returned, with rows in `term_order` (= `names(nn_input)`, so it lines up with
`bind_rows(lapply(nn_input, function(x) x$all))`):

| name | shape | contents |
|---|---|---|
| `pred_raw` | `(n,)` | uncalibrated global score |
| `pred` | `(n,)` | pooled-Platt calibrated |
| `per_index` | `(n, 36, K_pi)` | per-position softmax |
| `emb` | `(n, 16)` | L2-normalised `embed` layer output |

## Class bookkeeping

Class order is `["CT_cleavage_context", "DB", "gap", "NT_cleavage_context",
"pep_other", "pep_pocket", "padding", "none"]`. The per-index softmax reorders
to `[6 real, none, padding]`, so `K_cat = 7`, `K_pi = 8`, `none` sits at index
`6` and `padding` **last** at index `7`.

The per-index **loss** masks `none` and *keeps* `padding` — padding is a real
trained class whose label is leaked by the input padding channel, so it is
learned for free. The per-index **metric** masks both.

## Two switches where the R script diverges from its own design

`Config` defaults implement the intended architecture. `Config.r_exact()` (or
`--r-exact`) flips both back to reproduce `10_1dcnn_new6.R` bit-for-bit — use it
when comparing validation PR-AUC against an R run.

| switch | default | `r_exact` | what the R script does |
|---|---|---|---|
| `db_mask_both_termini` | `True` | `False` | The NT and CT masks are assigned to overlapping column sets in sequence, so the second assignment overwrites `DB` with the CT mask alone. `DB` ends up allowed only at positions 31–36, not at both termini. |
| `resample_each_epoch` | `True` | `False` | The `on_epoch_begin` callback rebinds `train_ds`, but `fit` already holds the original `tf.data` object, so the redraw never reaches training. The negative sample and the augmentation noise are drawn once for the whole run. |

Neither is an error the R would report; both change what the model sees.

## Layout

| file | what |
|---|---|
| `config.py` | `Config`: geometry, class bookkeeping, position masks, per-class weights, hyper-parameters |
| `model.py` | `build_model` / `compile_model` / `embed_model` and the two custom layers |
| `losses.py` | weighted BCE, focal per-index CE + smoothness, masked categorical accuracy |
| `data.py` | the array contract and the oversampled, noise-augmented sampler |
| `calibrate.py` | pooled Platt calibration (IRLS mirroring R's `glm`; scipy cross-check) |
| `pipeline.py` | `train` / `predict_all` / `embed` / `run` |
| `io.py` | the `.npz` + parquet exchange format |
| `cli.py` | `python -m lf_dcnn train｜selftest` |
| `synthetic.py` | contract-shaped fake windows for the tests |
| `selftest.py` | the checks, runnable without pytest |

scikit-learn is deliberately not used — it is not installed in the target venv,
which is why calibration is a hand-rolled IRLS rather than
`LogisticRegression`.
