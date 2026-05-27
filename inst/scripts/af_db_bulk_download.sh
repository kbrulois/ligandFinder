#!/bin/bash
## ---------------------------------------------------------------------------
## Bulk-download an EBI AlphaFold-DB proteome tarball and lay it out for use
## by ligandFinder's cleave_peptide_from_af() / af_structure_to_pdb().
##
## EBI ships per-proteome tarballs containing one CIF + one PAE JSON per
## entry, all gzipped.  This script downloads the tarball, extracts it,
## decompresses, and sorts files into:
##   <out>/cif/AF-<UNI>-F1-model_v4.cif
##   <out>/pae/AF-<UNI>-F1-predicted_aligned_error_v4.json
##
## Usage:
##   bash inst/scripts/af_db_bulk_download.sh <out_dir> [species] [afver]
##
##   <out_dir>   Cache root (required).
##   species     human (default) | mouse | rat | zebrafish | yeast | ecoli
##   afver       AlphaFold DB version tag (default: v4)
##
## Example:
##   bash inst/scripts/af_db_bulk_download.sh \
##        /oak/stanford/groups/ebutcher/deorphan-AI-ze/af_db human v4
##
## Wall time on Sherlock: ~10 min download + ~5 min extract for human
## (~5 GB tarball, ~20k entries).  Final on-disk size ~25 GB (uncompressed).
## ---------------------------------------------------------------------------

set -euo pipefail

OUT="${1:-}"
SPECIES="${2:-human}"
AFVER="${3:-v4}"

if [[ -z "$OUT" ]]; then
  echo "Usage: $0 <out_dir> [species=human] [afver=v4]" >&2
  exit 1
fi

case "$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')" in
  human)     PROTEOME="UP000005640_9606_HUMAN" ;;
  mouse)     PROTEOME="UP000000589_10090_MOUSE" ;;
  rat)       PROTEOME="UP000002494_10116_RAT" ;;
  zebrafish) PROTEOME="UP000000437_7955_DANRE" ;;
  yeast)     PROTEOME="UP000002311_559292_YEAST" ;;
  ecoli)     PROTEOME="UP000000625_83333_ECOLI" ;;
  *) echo "Unsupported species: $SPECIES" >&2; exit 1 ;;
esac

URL="https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/${PROTEOME}_${AFVER}.tar"
TARBALL="${OUT}/$(basename "$URL")"
STAGE="${OUT}/_stage_$$"

mkdir -p "$OUT/cif" "$OUT/pae" "$STAGE"

echo "Source: $URL"
echo "Target: $OUT"

if [[ ! -f "$TARBALL" ]]; then
  echo "[1/4] Downloading tarball..."
  # --show-progress is GNU wget >= 1.16; fall back to default verbose output.
  if wget --help 2>&1 | grep -q -- --show-progress; then
    wget -q --show-progress -O "$TARBALL" "$URL"
  else
    wget -O "$TARBALL" "$URL"
  fi
else
  echo "[1/4] Tarball already present, skipping download."
fi

echo "[2/4] Extracting..."
tar -xf "$TARBALL" -C "$STAGE"

echo "[3/4] Sorting CIF + PAE files into place..."
find "$STAGE" -name '*model_'"${AFVER}"'.cif*'                       -exec mv -t "$OUT/cif/" {} +
find "$STAGE" -name '*predicted_aligned_error_'"${AFVER}"'.json*'    -exec mv -t "$OUT/pae/" {} +

echo "[4/4] Decompressing (parallel)..."
find "$OUT/cif" -name '*.gz' -print0 | xargs -0 -P 8 -n 50 gunzip
find "$OUT/pae" -name '*.gz' -print0 | xargs -0 -P 8 -n 50 gunzip

rm -rf "$STAGE"
rm -f  "$TARBALL"

CIF_N=$(find "$OUT/cif" -name 'AF-*.cif'  | wc -l)
PAE_N=$(find "$OUT/pae" -name 'AF-*.json' | wc -l)
echo ""
echo "Done."
echo "  cif/: $CIF_N entries"
echo "  pae/: $PAE_N entries"
echo ""
echo "Use cleave_peptide_from_af(..., cache_dir=\"$OUT\") to consume."
