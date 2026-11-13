#!/bin/bash
## Fetch the ColabFold/MMseqs2 MSAs (UniRef30 2023-02 + ColabFold env DB 2021-08)
## of the LiteFold human-proteome set for every precursor in the peptide-end
## residue cache whose sequence the set holds EXACTLY. AFDB's own MSAs would be
## preferred but have returned 403 to anonymous callers since Aug 2026
## (google-deepmind/alphafold#1137); the scorer reads any .a3m, so they can be
## swapped in later.
##
## Input:  $D/hf_needed.tsv (acc, a3m_bytes; written by build step)
## Output: $D/hf_colabfold/<acc>.a3m.gz   (resumable: finished files are skipped)
set -u
set -o pipefail
D=${D:-$HOME/AF2_analysis/lf_msa}
BASE=https://huggingface.co/datasets/LiteFold/human-proteome-wide-msa/resolve/main
P=${P:-2}
fetch() {
  acc=$1; out=$D/hf_colabfold/$acc.a3m.gz
  [[ -s $out ]] && return 0
  read path want < <(awk -F'\t' -v a="$acc" '$1==a {print $6, $8; exit}' $D/manifest.tsv)
  for try in 1 2 3 4 5 6; do
    ## raw file first, size-checked against the manifest's a3m_bytes, THEN gzip:
    ## piping curl into gzip tests gzip's status, so a transfer cut mid-stream
    ## would be stored as a valid, truncated alignment and skipped forever after
    if curl -sfL --retry 3 -o $out.raw "$BASE/$path" \
       && [[ $(wc -c < $out.raw | tr -d ' ') == "$want" ]] \
       && [[ $(head -c1 $out.raw) == ">" ]] && gzip -c $out.raw > $out.tmp; then
      rm -f $out.raw; mv $out.tmp $out; return 0
    fi
    sleep $((try * 60))     # HF rate-limits anonymous IPs (429) for minutes at a time
  done
  rm -f $out.raw $out.tmp; echo "FAILED $acc" >&2; return 1
}
export -f fetch; export D BASE
tail -n +2 $D/hf_needed.tsv | cut -f1 | xargs -P $P -I{} bash -c 'fetch {}'
echo "done: $(ls $D/hf_colabfold/*.a3m.gz | wc -l) files"
