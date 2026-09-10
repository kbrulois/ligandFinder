#!/bin/bash
## ---------------------------------------------------------------------------
## Build AlphaPulldown input features for sequences the feature database
## cannot serve.
##
##     sbatch --array=1-<N> afpd_create_features.sh <sequences.fasta> [out_dir]
##
## WHEN THIS IS NEEDED -- AND WHEN IT ISN'T
## ----------------------------------------
## Most peptides in this pipeline are windows of a protein that already has
## features: "CXL14,95-102" slices residues out of the full-length CXL14
## pickle at prediction time and needs nothing built.  This script is for the
## other case -- a sequence with no entry in the feature directory at all:
## designed peptides, non-human orthologs, mutants.
##
## afpd_check_features() in R tells the two apart before you submit anything.
##
## THE HEADER IS THE JOB-FILE NAME
## -------------------------------
## AlphaPulldown names each pickle after the first token of the FASTA header,
## and job files address proteins by that same name.  So ">KNG1x381x389"
## produces "KNG1x381x389.pkl.xz" and the job line reads
## "BKRB2;KNG1x381x389".  Follow the directory's convention -- UniProt entry
## name without the _HUMAN suffix.
##
## ARRAY RANGE MUST MATCH THE RECORD COUNT
## ---------------------------------------
## --seq_index is 1-based and picks one record out of the FASTA, so:
##     sbatch --array=1-$(grep -c '^>' seqs.fasta) afpd_create_features.sh seqs.fasta
## A short array leaves the tail unbuilt and says nothing about it; the gap
## resurfaces hours later as a FileNotFoundError from
## run_structure_prediction.py.  The task checks its own index below.
##
## NO GPU
## ------
## Feature building is MSA and template search -- CPU and memory, no GPU.
## It is submitted to a different partition than the prediction jobs.
## ---------------------------------------------------------------------------

#SBATCH --job-name=apd_features
#SBATCH --time=24:00:00
#SBATCH -e apd_features_%A_%a.err
#SBATCH -o apd_features_%A_%a.out
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

FASTA="${1:?usage: sbatch --array=1-N afpd_create_features.sh <fasta> [out_dir]}"
OUT_DIR="${2:-/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/peptides}"

DATA_DIR=/oak/stanford/groups/ebutcher/catherine/alphafold_db
MAX_TEMPLATE_DATE=2024-01-01

AP_ROOT=$GROUP_HOME/programs/anaconda3
AP_ENV=$AP_ROOT/envs/AlphaPulldown

## Two words, not one.  "source <activate> <env>" quoted into a single
## argument makes bash look for a file whose name contains a space, and the
## error names both paths joined together.
source "$AP_ROOT/bin/activate" "$AP_ENV"

CIF=$AP_ENV/bin/create_individual_features.py
[[ -f "$CIF" ]] || CIF=$(command -v create_individual_features.py || true)
[[ -n "$CIF" ]] || { echo "ERROR: create_individual_features.py not found in $AP_ENV" >&2; exit 1; }

[[ -f "$FASTA" ]] || { echo "ERROR: no such fasta: $FASTA" >&2; exit 1; }

mkdir -p "$OUT_DIR"
umask 0002

IDX="${SLURM_ARRAY_TASK_ID:-1}"
N_SEQ=$(grep -c '^>' "$FASTA")

if (( IDX > N_SEQ )); then
    echo "ERROR: array index $IDX exceeds the $N_SEQ record(s) in $FASTA" >&2
    exit 1
fi

NAME=$(grep '^>' "$FASTA" | sed -n "${IDX}p" | sed 's/^>//; s/[[:space:]].*//')

echo "fasta      : $FASTA ($N_SEQ records)"
echo "seq_index  : $IDX -> $NAME"
echo "out_dir    : $OUT_DIR"
echo "data_dir   : $DATA_DIR"

ARGS=(
    --fasta_paths "$FASTA"
    --data_dir "$DATA_DIR"
    --output_dir "$OUT_DIR"
    --max_template_date "$MAX_TEMPLATE_DATE"
    --seq_index "$IDX"
    --skip_existing
)

## The curated database is .pkl.xz.  --compress_features is not in every
## AlphaPulldown build, so ask this one rather than assume, and fall back to
## compressing the pickle afterwards -- the extension is what the feature
## lookup matches on.
COMPRESSES=0
if python "$CIF" --helpfull 2>&1 | grep -q -- '--compress_features'; then
    ARGS+=(--compress_features)
    COMPRESSES=1
else
    echo "note   : no --compress_features in this build; compressing afterwards"
fi

## Opt-in.  MMseqs2 runs the MSA on ColabFold's public server -- minutes
## instead of hours, but it uploads the sequence.  Off by default: unpublished
## designed peptides should not leave the cluster.
if [[ "${AP_MMSEQS2:-0}" == "1" ]]; then
    echo "msa    : ColabFold MMseqs2 server (sequence leaves the cluster)"
    ARGS+=(--use_mmseqs2)
else
    echo "msa    : local databases under $DATA_DIR"
fi

python "$CIF" "${ARGS[@]}"

if [[ "$COMPRESSES" -eq 0 && -f "$OUT_DIR/$NAME.pkl" ]]; then
    xz -T0 "$OUT_DIR/$NAME.pkl"
fi

## Verify rather than trust the exit status: --skip_existing makes a no-op
## look identical to a successful build, and a header the pickle was not
## named after would leave nothing behind at all.
echo
if [[ -f "$OUT_DIR/$NAME.pkl.xz" ]]; then
    echo "OK: $OUT_DIR/$NAME.pkl.xz ($(du -h "$OUT_DIR/$NAME.pkl.xz" | cut -f1))"
elif [[ -f "$OUT_DIR/$NAME.pkl" ]]; then
    echo "OK: $OUT_DIR/$NAME.pkl (uncompressed)"
else
    echo "ERROR: nothing written for '$NAME' -- check the header in $FASTA" >&2
    ls -la "$OUT_DIR" >&2
    exit 1
fi
