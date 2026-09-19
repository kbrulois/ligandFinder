# lf_plm

Per-residue protein-language-model embeddings as extra input channels for the
`lf_dcnn` window models. Two resumable steps, no R in the loop:

```bash
# in a venv with torch (+ MPS/CUDA if you have it), transformers, sentencepiece,
# scikit-learn, h5py, pyarrow, pandas, joblib -- e.g. ~/.virtualenvs/lf-plm
cd inst/python
PYTHONPATH=. python -m lf_plm embed  --sequences ~/AF2_analysis/lf_plm/sequences.parquet \
                                     --out ~/AF2_analysis/lf_plm/prot_t5_raw.h5 --device mps
PYTHONPATH=. python -m lf_plm reduce --h5  ~/AF2_analysis/lf_plm/prot_t5_raw.h5 \
                                     --sequences ~/AF2_analysis/lf_plm/sequences.parquet \
                                     --k 32 --out ~/AF2_analysis/lf_plm/prot_t5_pca32.parquet
```

**Input is a table, not FASTA.** `sequences.parquet` has `accession`,
`gene`, `sequence_uni` -- what R exports from `secretome` for the genes in
`nn_input`. Whole precursors are embedded, signal peptide included: the models
were trained on full UniProt entries and a dibasic-site residue needs the
domain context a 36-residue window would strip. Slicing to windows happens on
the R side.

## embed

`Rostlab/prot_t5_xl_half_uniref50-enc` (ProtT5-XL-U50 encoder, 1024-d,
fp16) by default; `--backend esm_c` for ESM C 600M (1152-d, needs the `esm`
package). One `(L, D)` float16 dataset per accession in the h5; accessions
already present are skipped, so an interrupted run resumes. Sequences over
`--max-len` (2000) are embedded in overlapping chunks with a linear blend
across the overlap (cosine ~0.98 to the whole-sequence embedding on a test
protein). On an Apple-silicon GPU (`--device mps`) ProtT5 does roughly
1,000 residues/s; the 5,182-protein secretome set is ~45 min. On CPU it is
hours.

## reduce

1024 channels beside the 26 hand-built ones would put ~50k parameters in the
first conv of a ~21k-parameter model, against a few dozen positives. So: one
`IncrementalPCA` fit over a 500k-residue sample of the whole set (global, not
per protein -- the same rule 9.2's min-max ranges follow), cached beside the
output as `<out>_pca<k>.joblib`, applied identically to every protein. Each
component is z-scored by its sd from the fit, clipped to ±3 sd and mapped to
[0, 1], the range every other channel is on. Output: one row per residue,
`accession, index, AA, plm_01..plm_k`, plus a `.json` with the explained
variance.

## R side

`R/plm_features.R`: `lf_plm_read()` loads the table as one matrix per
accession; `lf_plm_attach(nn_input, plm, sequences, all_params3)` appends the
`plm_*` columns to every window's `data` tibble by `(accession, index)`,
checking the residue letters against the table (a mismatch means a different
sequence version and is an error, not a silent misalignment), with 0 at
padding rows like every other channel. It returns the extended `nn_input`
and channel list, ready for `lf_dcnn_export()` / `lf_dcnn_run()`.
`inst/scripts/10_5_benchmark_window_model.R --preset t5` is the with/without
comparison.

The PCA channels are deliberately NOT in `Config.cont_channel_names`, so they
get no train-time jitter; add them there if the with-PLM model overfits faster
than the baseline.
