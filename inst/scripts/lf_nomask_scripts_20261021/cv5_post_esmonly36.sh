#!/bin/zsh
# The ESM-C-only counterpart of cv5_post.sh (Kevin, 2026-10-24): all-50 ensemble
# dirs, per-gene scans of the 104 known-peptide page genes, and the pooled all-50
# npz -- for the TWO NEW PAGE TRACKS (the C set and the N set of ESM-C-only
# models). A copy rather than an edit, because cv5_post.sh is the record of what
# produced the curated pages and PR #14 is also changing it.
#
# TWO DIFFERENCES THAT MATTER beyond the TAG:
#  * --residues-npz MUST be this arm's 37-channel table. build_windows copies
#    res.feat into a (k, seq_len, n_channels) array POSITIONALLY, so the default
#    26-channel table raises, and a wrong-order table of the right width would
#    scan silently against the wrong channels.
#  * PROCS defaults to 4. Drop it if a large render is running: this machine is
#    memory-bound, so concurrent heavy jobs halve each other's throughput.
set -eu
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
ARM=${ARM:-esmonly36}
TAG=_nomask_flatten_${ARM}
RESID=${RESID:-$A2/lf_pepend_residues_${ARM}.npz}
PROCS=${PROCS:-4}
## GENES is overridable so a single gene can be scanned as a test
## (GENES=CXCL14 ...); unset it means all 104 page genes. A one-gene run
## REPLACES the npz, so re-run with the full list before rendering a gallery.
export GENES=${GENES:-$(paste -sd, $S/page_genes.txt)}

if [[ ${1:-} == scan ]]; then          # one fold/terminus scan (what xargs calls)
  k=$2 T=$3
  set +e                               # report a failed scan, don't abort the queue
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  python -m lf_winxgb.scan_genes --genes $GENES --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds \
    --residues-npz $RESID \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked \
    --out $A2/lf_winxgb_scan_${T}_cv5_f${k}${TAG}.npz > $S/cv5scan_${T}_f$k.log 2>&1
  echo "scan $T f$k rc=$?"
  exit 0
fi

for T in C N; do for k in 1 2 3 4 5; do
  c=$A2/lf_pepend_run_${T}_cv5_f${k}${TAG}_isolated/in/config.json
  grep -q '"class_position_mask": false' $c || { echo "MASK ON in $c"; exit 1; }
  grep -q '"global_head": "flatten"' $c || { echo "NOT FLATTEN: $c"; exit 1; }
  grep -q '"head_lr_mult": 1.0' $c || { echo "HEAD LR NOT 1: $c"; exit 1; }
  n=$(ls $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}_isolated/out/members/member_*_${T}.weights.h5 | wc -l)
  [[ $n -eq 10 ]] || { echo "$T f$k has $n members, want 10"; exit 1; }
done; done
echo "configs OK: flatten head, x1 lr, mask off, 10 members x 5 folds x 2 termini"

for T in C N; do
  E=$A2/lf_pepend_run_${T}_cv5all50${TAG}_isolated
  rm -rf $E; mkdir -p $E/in $E/out/members
  cp $A2/lf_pepend_run_${T}_cv5_f1${TAG}_isolated/in/config.json $E/in/
  i=0
  for k in 1 2 3 4 5; do
    for f in $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}_isolated/out/members/member_*_${T}.weights.h5; do
      ln -s $f $E/out/members/$(printf "member_%03d_%s.weights.h5" $i $T); i=$((i+1))
    done
  done
  echo "$T: $i members -> $E"
done

for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P $PROCS -n 2 ${0:A} scan

source ~/.virtualenvs/r-tensorflow/bin/activate
TAG=$TAG python - <<'EOF'
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis"); TAG = os.environ["TAG"]; GENES = os.environ["GENES"].split(",")
n, K = 10, 5
for T in ("C", "N"):
    z = [np.load(f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}.npz", allow_pickle=True) for k in range(1, K + 1)]
    out, miss = {}, []
    for g in GENES:
        if f"{g}__anchor" not in z[0].files: miss.append(g); continue
        an = z[0][f"{g}__anchor"]
        assert all(np.array_equal(an, zz[f"{g}__anchor"]) for zz in z), g
        out[f"{g}__anchor"] = an
        for h in ("cnn", "xgb_pm"):
            M = np.stack([zz[f"{g}__{h}"].astype("f8") for zz in z])
            Sd = np.stack([zz[f"{g}__{h}__sd"].astype("f8") for zz in z])
            m = M.mean(0)
            var = ((n - 1) * (Sd ** 2).sum(0) + n * ((M - m) ** 2).sum(0)) / (n * K - 1)
            out[f"{g}__{h}"] = m.astype("f4"); out[f"{g}__{h}__sd"] = np.sqrt(var).astype("f4")
    np.savez_compressed(f"{A2}/lf_winxgb_scan_{T}_cv5all50{TAG}.npz", **out)
    print(f"{T}: {len(out)//5} genes pooled over 5 folds x 10 members" +
          (f"; not in the scan: {', '.join(miss)}" if miss else ""))
EOF
echo "POST DONE"
