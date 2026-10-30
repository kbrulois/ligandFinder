# Candidate peptide table — how to rebuild it

One command regenerates everything, at any threshold / tolerance:

```bash
./candidate_pipeline.sh                       # the defaults: xgb_pm, T 0.5, tol 4, 8-60 aa
THRESH=0.3 TOL=2 ./candidate_pipeline.sh      # another cut
HEAD=cnn ./candidate_pipeline.sh              # the linear head instead
```

Output `lf_headcmp/candidate_peptides_<HEAD>_T<THRESH>_tol<TOL>.csv`, with a
`.manifest.txt` beside it recording parameters, input sizes and mtimes, tool
versions, the exact mmseqs command and the full column list. Stage 1 also
writes `lf_headcmp/end_segments_<HEAD>_T<THRESH>.csv` (no TOL: a segment does
not depend on it), one row per end segment with its run, centre and how many
candidates it is in -- the protein pages draw it, and `final_pages.R`'s
`CAND_SETS` picks which sets a page shows.

## Stages

| | script | adds |
|---|---|---|
| 1 | `candidate_peptides.R` | the candidates themselves: pair N and C end segments, fan out by ±TOL, `group_id`/`group_rank`/`pair_score`, `end_score_max`, `total_offset`, motifs; plus the `end_segments_*.csv` sidecar |
| 2 | `candidate_mmseqs_prep.R` + `mmseqs easy-search` | search against known peptides 8–60 aa (`KNOWN_LEN_MIN`/`MAX`); no k-mer prefilter, so every candidate, however short, gets its closest one (settings and why: `candidate_mmseqs_args.sh`) |
| 3 | `candidate_annotate.R` | `mmseqs_*` best hit and per-class bests, `length_bin`, `overlap_*` |
| 4 | `candidate_multipep.sh` | 20 `mp_*` MultiPep classes + `mp_signalling_max` |
| 5 | `candidate_logistic.R` | `lr_prob`, `lr_prob_oof` |
| 6 | `candidate_xgb.R` | `xgb_score` (out of fold, grouped by gene), `xgb_nomp_score` (the same model without any MultiPep feature), `xgb_end_score` (the same on the 12 held-out end-score features alone), `label`, `curated_relation` / `uniprot_relation` (the known peptides of 8–60 aa a candidate touches; reported, never features), `mp_best_class` / `mp_best_score` (the best of peptide hormone or neuropeptide only), `*_neglog10_evalue` (−log10 of the MMseqs e-values), and the tidy `candidate_peptides_<...>_scored.csv`; features are held-out end scores, MultiPep and MMseqs hits on OTHER precursors only (see its header for the three leaks it avoids) |
| 7 | `candidate_best_per_group.R` | `candidate_peptides_<...>_best_per_group.csv`: per group (one N es × C es pair) the best cut by `xgb_score`, `xgb_nomp_score`, `end_score_max`, `mp_best_score` (peptide hormone or neuropeptide) and `mmseqs_sim` — 1 to 5 rows, `picked_by` naming the criteria each won |

`WITH_PDEEPPP=1` adds a stage 4b. It is **off by default**: on our own ground
truth PDeepPP's neuro head ranks real ligands backwards (AUC 0.21), behaving
like a length cut near 29 aa. Kept only so that finding stays checkable.

## Pages for a selection of candidates

```bash
O=~/AF2_analysis/pepend_pages_cand_nocurated_mp0.6_nomp0.45
FILTER='curated_relation == "none" & mp_best_score > 0.6 & xgb_nomp_score > 0.45' \
  LF_EXTRA_BUNDLE=~/AF2_analysis/plot_bundle_extra.rds OUT=$O Rscript candidate_select.R
GENES_FILE=$O/scan_genes.txt SFX=_cand_nocur_mp06_nomp045 ./arm_gene_scans.sh
OUT=$O ./candidate_pages.sh                # FORCE=1 redraws pages that already exist
```

- **`candidate_select.R`** evaluates FILTER, any R expression over the stage-6
  scored table, and writes into OUT:
  - the passing candidates and one row per precursor;
  - `genes.txt` (the drawable genes) and `scan_genes.txt` (the same as
    `GENE@ACCESSION`, the accession each page draws, e.g. COL18A1@P39060);
  - `focus.csv` (each card's anchor);
  - `selection.txt` (key, filter, extra bundle).
  The anchor is the gene's best passing candidate on that accession, by the first
  of xgb_nomp_score, xgb_score, mp_best_score, end_score_max, mmseqs_sim that
  FILTER names, whatever order FILTER names them in (`FOCUS=` overrides).
- **`arm_gene_scans.sh`** gives the ESM-C only and AAswap tracks their sd bands
  and records its suffix in OUT.
- **`candidate_pages.sh`** draws the genes with `final_pages.R` exactly as the
  candidate gallery is drawn, then the index, whose cards open each page at its
  anchor, best first.
  - It reads the key, the extra bundle and the arm-scan suffix from OUT, and
    refuses a given value that disagrees.
  - It is resumable. Each chunk renders in its own staging directory, and a page
    is accepted only when it is complete and self-contained.
  - A forced pass interrupted part-way is resumed by a plain re-run.
  - 3 processes (~8.5 GB resident each) draw 108 pages in about 40 min.

Every window-head track draws a ±1 sd band for each fold (across its 10 seeds)
and one for the all-folds mean (across all 50 members, pooled from the folds).
- **Curated arm:** the bands come from its per-gene scans and its live CNN
  members.
- **ESM-C only and AAswap arms:** their full-proteome tracks hold means only, so
  they need `arm_gene_scans.sh`. That script keeps the per-gene CNN mean and sd
  but takes xgb from the full-proteome scan. The per-gene refit of each member's
  xgb head drifted up to 0.047 for ESM-C, and the segments, the candidates and
  the earlier pages all use the full scan's values.

Every point on the page has one opacity, outline and fill alike. The exception
is the candidate track's class symbols, which stay opaque on the turbo-coloured
bars. Tooltips keep their own width (`width:max-content`), so the ones at a
track's right end no longer collapse to a few words per line.

- **`curated_relation == "none"` means no known peptide of 8–60 aa at any of
  the candidate's residues.** That covers both a curated ligand (`exact` /
  `overlap`) and, failing one, a UniProt peptide (`uniprot`; details in
  `uniprot_relation` / `uniprot_*`).
  - It can still pass a UniProt propeptide, or a cut within ±2 of a known
    peptide outside 8–60 aa (neuromedin-N, 6 aa). `& label %in% 0L` drops those
    too.
  - Tables scored before this change (2026-10-28) compared with curated ligands
    only. Their "none" let 222 UniProt peptides into the first gallery
    (gonadoliberin, orexin-A, …). Re-run stage 6 to apply the new meaning.
  - Each card's note lists the known peptides among its gene's passing cuts.
- **Genes `plot_bundle.rds` lacks.** About 430 secretome genes have no dibasic
  site after the signal peptide, which the old window model needed. Build their
  entries with `inst/scripts/plot_bundle_extra.R`, then select with
  `LF_EXTRA_BUNDLE`. Those pages take make_protein_plot_win's no-window branch.
- **COL18A1** has two accessions under one name, so the all-secretome scan
  skipped it. It is drawn from its own scan files (`NPZ_SFX=_col18a1`, built by
  `own_scan_npz.py` from the full-proteome scans). `candidate_pages.sh` does
  this for every gene in its `OWN_NPZ` table.

## Caches — why a re-run is cheap

Both external scorers cache **by sequence**, not by candidate, so changing
THRESH or TOL only pays for sequences never seen before:

- `multipep_ext/multipep_scores.csv`
- `pdeeppp_run/pdeeppp_neuro_scores.csv`

## Environments

- `~/.virtualenvs/multipep` — numpy<2, TF, pandas, sklearn, matplotlib
- `~/.virtualenvs/pdeeppp` — torch + fair-esm (ESM-2 650M cached in ~/.cache/torch)
- **Never install into `~/.virtualenvs/r-tensorflow`** — the protein-page render
  uses it, and a numpy upgrade there breaks TensorFlow for every new render
  process.

## Two things the pipeline guards

- **Leakage.** `lr_prob`'s outcome is defined from the `mmseqs_*` columns, so
  every `mmseqs_*` / `overlap_*` column is barred from its predictors;
  `candidate_logistic.R` asserts this.
- **Correlated rows.** Up to (2·TOL+1)² near-identical candidates share one
  segment pair, so CV folds are grouped **by gene**. A random split reports a
  fantasy AUC.
