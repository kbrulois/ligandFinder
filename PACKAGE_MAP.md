# ligandFinder package map

This map says where each kind of code lives, where new code goes, what the
pipelines are, and where their data comes from and goes. **Read it before
adding or moving code, and update it in the same commit when you add a folder
or change what a folder holds.** The rules are in `CLAUDE.md`.

## Two packages

- **ligandFinder** (this repo) predicts peptide-hormone precursor ends from
  per-residue features with 1-D CNN "arms". It also builds the residue database
  the models read and renders the protein pages and the public site.
- **fap** (`~/R_projects/fap`, "fast AlphaPulldown") runs the AlphaPulldown and
  AlphaFold 3 complex predictions of GPCRs with ligands on Sherlock, and turns
  them into metrics and contacts. It depends on ligandFinder; **ligandFinder
  never depends on fap**.
  - Code that builds, submits or reads AlphaFold complex models belongs in fap,
    which has its own README map.
  - What crosses back from fap is data (the pocket labels, below), never code.

## Top level

| path | what |
|---|---|
| `R/` | package functions. R ignores subfolders here, so it stays flat; group by file-name prefix (next section) |
| `inst/python/` | the Python library: `lf_dcnn` (the CNN), `lf_winxgb` (xgboost window head, full scans), `lf_resid` (window-free per-residue models), `lf_plm` (protein language model embeddings, a residue-db channel), `tests/` |
| `inst/scripts/` | every script, by kind (below) |
| `inst/templates/` | HTML templates the page and index scripts fill |
| `inst/extdata/` | small versioned reference tables: known ligands and pairings, GPCR lists, cell-type and location tables |
| `data/` | `id_mapping.rds` (UniProt entry / name / gene; fap reads it too) |
| `terraform/` | S3 and CloudFront hosting of the public site |
| `man/` | 2 Rd files and the README figure (the roxygen docs are not generated) |

### `R/`, by topic
| topic | files |
|---|---|
| residue database | `residue_db.R` (db path, downloads), `residue_mapping.R` (sequence and position mapping), `conservation.R`, `plm_features.R`, `eec_expression.R` |
| peptide ends and training windows | `pepend_windows.R`, `pepend_cv_folds.R`, `known_peptide_ends.R`, `amidation_motif.R`, `cap_peptides.R`, `uniprot_terminus_hits.R` |
| the Python bridge | `dcnn_bridge.R` (`lf_dcnn_run()` and friends) |
| pages and figures | `pepend_tracks.R`, `plot_proteins.R`, `plot_proteins_win_new.R`, `html_plot.R`, `viz_theme.R`, `plot_msa_conservation.R` |
| gallery indexes | `gallery_tpfn.R`, `gallery_known.R`, `gallery_candidate_plot.R` |
| shared helpers | `utils.R`, `lists_utils.R` |

## `inst/scripts/`

Scripts find each other through `lf_paths`, never through their own directory.

- **Shell:** `. ${0:A:h:h}/lf_paths.sh` sets `W` (the source tree: `$LF_SRC`, or
  the checkout the script sits in). It also sets `LF_RDB`, `LF_MOD`, `LF_EVAL`,
  `LF_VIZ` and `LF_SCRIPTS`, and `LF_RUN` (the run directory, under "Data").
  A script in one folder calls another as `$LF_EVAL/end_calls.R`, and so on.
- **R:** `source(file.path(.lf_here, "..", "lf_paths.R"))`, then
  `lf_src("inst/scripts/modeling/arm_vars.R")` or `lf_src("R/pepend_tracks.R")`.
  Use one `".."` per folder level.
- **Python:** put `inst/scripts/` on `sys.path`, then
  `import lf_paths; lf_paths.add_python()`.

| folder | what goes here | entry points |
|---|---|---|
| `generate_residue_db/` | everything that makes per-residue data (below) | `build_residue_db.R` |
| `modeling/` | training windows, training, the arm pipeline, scans, the candidate classifier. R, shell and Python drivers; the library they call is `inst/python/` | `10_7_pepend_windows.R`, `10_8_pepend_train.R`, `pepend_arm_pipeline.sh`, `run_msa_arms.sh`, `candidate_pipeline.sh` |
| `evaluation/` | metrics on a trained arm. Includes `cv5_unseen_sfx.R`, which builds the never-seen set; diagnostics like `trunk_val_compare.R`, `terminus_bias.R`, `tied_score_blocks.R`, `asa_afm_depth.R`; attribution (`shap_*`). **These write tables into `~/AF2_analysis/lf_headcmp`** | `end_calls.R`, `metrics_final.R`, `segments_roc.R`, `segments_motif.R`, `peaks_final.R` |
| `visualization/` | protein pages, gallery indexes, figures, the public site | `final_pages_arms.R` (stage 7), `arms2_allsec.sh` / `allsec_*.sh` (the all-secretome gallery), `make_gallery_index.R`, `make_known_index.R`, `make_candidate_plot.R`, `upload_plots_s3.sh` |
| `analysis/<question>/` | one-off scripts, one folder per question, each with a README saying what it found. They are records: nothing calls them | `analysis/README.md` indexes them |
| `analysis/archive/` | retired code, kept as a record: the db-anchored window model, replaced residue-db steps | |

### `generate_residue_db/`
- **The pipeline, in run order:** `1_parse_uniprot` … `9_expand_by_residue`, then
  `10_export_parquet`, driven by `build_residue_db.R`. It writes `secretome` and
  `secretome_aa` (one row per residue) and the parquet residue db at
  `get_db_path()`. `cleavage_sites.R`, `topology_initial.R` and `dssp.py` are
  steps of it.
- **AlphaFold-complex pocket labels:** `9.2_add_contact_data.R` once added them.
  It is not run now, but it is the step that should read fap's pocket CSV
  (under "Data").
- **The putative-peptide steps `11_` … `14_`:** replaced by the window CNN, and
  kept only because `visualization/10_4_plot_proteins.R` reads what
  `12_extract_peptides.R` and `13_plot_proteins.R` made. Kevin decides whether
  to retire them.
- **Extra per-residue channels:**
  - `msa_cons/`: MSA conservation from ColabFold MSAs, the `plm_cons_msa`
    channel.
  - `aa_onehot_build.py`, `motif_ind_build.py`, `plm_pagegenes.R`.
  - The arm residue tables: `residues_arm.py`, `residues_esmonly36.py`, and
    their checks.
  - The protein language model embeddings themselves come from
    `inst/python/lf_plm`.
- **Annotations:**
  - enteroendocrine and HPA expression (`enteroendocrine_expression.R`,
    `hpa_top_celltype.R`, `add_eec_*`, `merge_eec_*`);
  - secreted-protein and known-ligand lists;
  - disulfides (`get_dsb.R`, `get_afdsb.R`);
  - dibasic and cleavage-site helpers.
- **Check:** `test_residue_db.R` checks a built db.

### `modeling/`
- **Training windows:** `10_7_pepend_windows.R` builds the 5-fold inputs
  (`lf_pepend_nn_input_cv5_f<k><SFX>.rds`) and the per-residue feature cache.
  `10_9r_export_window_groups.R` (stage 2b) and `10_9p_export_residue_labels.R`
  export window groups and residue labels.
- **Training:** `10_8_pepend_train.R` trains one trunk, `train_arm.sh` trains an
  arm's 10 trunks, and `10_8t_pepend_tune.sh` searches hyperparameters.
- **The arm pipeline:**
  - **What an arm is:** a channel set, defined in one file, `arms/<arm>.sh`.
  - **Driver:** `pepend_arm_pipeline.sh` (`ARM=…`, `STAGES=…`) runs the stages:
    channels → trunks → scans → tracks → evaluation → plots → pages →
    attribution.
  - **Scripts per stage:** `fullscan_arm.sh`, `groups_esmcov.sh`,
    `cv5_post*.sh`, `arm_post.sh`, `tracks_esmonly36.R` (generic, despite the
    name).
  - **Unattended runs:** `run_msa_arms.sh` runs two arms end to end.
  - **Docs:** `ARM_PIPELINE_README.md` and `ARM_COMPARE_README.md`.
  - **Records of the ESM-C arm:** `esmonly36*_driver.sh`, `ESMONLY36_README.md`.
- **Scans:** the 10_8b scan, `arm_gene_scans.sh` (per-gene scans for pages),
  `score_anchors.py` and `own_scan_npz.py`.
- **Candidate peptides (stage 6):** `candidate_pipeline.sh` drives the
  `candidate_*` scripts (MMseqs2, MultiPep, xgboost classifier, selection);
  `CANDIDATES_README.md`.

## Pipelines, in order
1. **Residue db:** `generate_residue_db/build_residue_db.R`, then the extra
   channels (`msa_cons/`, `lf_plm`, `aa_onehot_build.py`, …).
2. **Windows:** `modeling/10_7_pepend_windows.R` (one run per fold and window set).
3. **An arm:**
   `LF_SRC=<checkout> zsh $LF_SRC/inst/scripts/modeling/pepend_arm_pipeline.sh`,
   run from `$LF_RUN`.
   - Stages 1-4 (trunks, scans, tracks) call `inst/python`.
   - Stage 5 is evaluation: `evaluation/end_calls.R`, `metrics_final.R`,
     `segments_*.R`.
   - Stage 6 (candidates) and stage 7 (pages, `visualization/`).
   - Then attribution (`evaluation/shap_*`).
4. **The public site:** a gallery (`visualization/arms2_allsec.sh` and
   friends), its indexes, then `upload_plots_s3.sh` to S3 and CloudFront
   (`terraform/`).

## Data: where inputs come from, where outputs go

Data never lives in a code folder. Big inputs and every output are under
`~/AF2_analysis` (or the residue-db directory). Only small, versioned
reference tables go in `inst/extdata/` or `data/`.

| input | where | made by |
|---|---|---|
| UniProt proteome and annotations, AlphaFold monomer models, DSSP, AlphaMissense, aminode | residue-db `raw/` | downloads; `generate_residue_db/1_…5_` |
| `secretome`, `secretome_aa` (one row per residue), parquet db | `~/peptide_alg/build_residue_db/processed/`, `get_db_path()` | `generate_residue_db/` steps 1-10 |
| `secretome_latest.rds` | `~/AF2_analysis` | the residue db's `secretome` table; `10_7 --secretome` reads it |
| **`knowns.rds`: known peptides docked on their receptors** | `~/AF2_analysis/knowns.rds` | **fap**: AlphaPulldown models → `metrics/` extraction → contacts. See below |
| `nn_dat_cache.rds` | `~/AF2_analysis` | the retired db-anchored model's data prep (`9.2_add_contact_data.R` era). `10_7` reads its channel scaling (`all_params3`, `chan_range`) and the old known windows from it |
| per-residue feature cache | `~/AF2_analysis/lf_pepend_residue_cache_v2.rds` | `10_7 --rebuild-features` (v2 = the 2026-11-12 residue-join fix) |
| UniProt peptide / propeptide ends | `~/AF2_analysis/pepend_locus/uniprot_peptides.csv` | UniProt export |
| known ligands, receptor pairings, candidate peptides | `inst/extdata/ligand_list.rds`, `GPCRdb_known_pairings_human_plus2more_unique.csv`, `Candidate_peps_KB_…txt` | curated |
| PLM and extra channels | `~/AF2_analysis/lf_plm/*.parquet` | `lf_plm`, `msa_cons/`, `aa_onehot_build.py`, `motif_ind_build.py` |
| the run directory: gene lists, blurbs, never-seen tables, index cards, logs | `$LF_RUN` (default `~/AF2_analysis/lf_run`) | the base never-seen set `cv5_unseen_ends.csv` was written once by `analysis/pepend_cv5_review/cv5_sets.R` (2026-10-21), and its per-window-set variants `cv5_unseen_ends_<sfx>.csv` by `evaluation/cv5_unseen_sfx.R`. The blurbs come from `visualization/index_tpfn*.R`. No script in the tree writes `page_genes.txt` or `gene_lengths.txt`. The logs come from the drivers |

**The pocket labels:**
- **What they are:** which residues of a known peptide sit in the receptor
  pocket, and how deep (`IC_dist`, the residue's Cα projected on the receptor's
  extracellular-to-intracellular axis).
- **Where they come from:** the AlphaFold complex models in `knowns.rds`.
  `10_7` uses one model per peptide (`location == "relevant"`, rank 1), and
  `R/afpd_metric_processing.R` in fap computes the contacts.
- **Plan (Kevin, 2026-11-13):** fap writes these labels as a CSV of
  pocket and non-pocket residues per peptide and model. `generate_residue_db/`
  reads that CSV, through a revived `9.2_add_contact_data.R` or its
  replacement. When that CSV exists, list it here with its path and the fap
  script that writes it.

| output | where |
|---|---|
| fold inputs, trunks, members, scans, tracks | `~/AF2_analysis/lf_pepend_*`, `lf_winxgb_scan_*.npz`, `lf_pepend_cv5_tracks_*.rds` |
| evaluation tables | `~/AF2_analysis/lf_headcmp/` (shared: `final_metrics.csv` is keyed by arm) |
| figures | `~/AF2_analysis/lf_viz_<date>/` |
| protein pages and galleries | `~/AF2_analysis/pepend_pages_*` |
| the public site | S3 behind CloudFront (`terraform/`, `upload_plots_s3.sh`) |

## Where new code goes
- **Builds or adds a per-residue feature or annotation** → `generate_residue_db/`.
  A new channel also needs an arm file and the lf_dcnn channel lists.
- **A new arm (channel set)** → `modeling/arms/<arm>.sh`. Training, scanning or
  pipeline logic → `modeling/`.
- **Scores a trained model** → `evaluation/`.
- **Draws a page, index or figure** → `visualization/`. A figure for one question
  goes in `analysis/<question>/`.
- **Answers one question once** → `analysis/<question>/`, with a README.
- **Python used by more than one script** → a module in `inst/python/`.
- **R used by more than one script** → a function in `R/<topic>.R`.
- **Builds, runs or reads AlphaFold complex models** → fap, not here.
- **Nothing fits** → add a folder, and add it to this map in the same commit.

## Adding a species (more training data)
The pipeline is human-only today. To add species, parameterise the scripts by
species rather than copying them per species.
- **Proteome.**
  - `1_parse_uniprot.R` reads the human UniProt proteome; give it a proteome or
    taxon id.
  - A 2025 mouse variant, `1_parse_uniprot_mouse.R`, was deleted with
    `putative_peptide/old/`. Recover it with
    `git show 2b268a1^:inst/scripts/putative_peptide/old/1_parse_uniprot_mouse.R`.
- **Structures.** fap's `structures/af_db_bulk_download.sh` already takes a
  species. DSSP (`4_integrate_dssp`) and relASA follow the structures.
- **Channels without a non-human source.** AlphaMissense (`5_integrate_afm`,
  the `min_afm` / `mean_afm` channels) exists for human only. The
  arm needs a policy for it: drop the channel, or treat it as missing, which
  is not the same as the by-design zeros of intracellular and
  transmembrane residues.
- **Channels that work for any species:** the MSA conservation (`msa_cons/`,
  ColabFold) and the PLM embeddings (`lf_plm`).
- **Known peptides.** `ligand_list.rds` and the UniProt peptide annotations per
  species.
- **Folds.** Orthologs of a held-out human peptide must sit in the same fold.
  `10_7` splits by sequence-similarity group, so build those groups across all
  species together, or the never-seen set leaks.
- **Outputs:** one residue-db directory per species (under `get_db_path()`),
  with the species in every output name.

## Outside the package, on purpose
Snapshots kept as records:
- `~/AF2_analysis/lf_src_arms2_20261103`, the code that rendered the public site;
- `~/AF2_analysis/_archived_code_20261112`;
- `~/AF2_analysis/lf_tune_C/code`, the tuner's copy of `lf_dcnn`.

Do not edit them, and do not add code beside them.
