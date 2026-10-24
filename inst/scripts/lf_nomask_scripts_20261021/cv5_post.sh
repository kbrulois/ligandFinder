#!/bin/zsh
# After train_final.sh, for the 5-fold flatten models:
#   1. every config is the flatten head, mask off, 10 members per fold/terminus
#   2. the all-50 ensemble per terminus: one members/ dir of symlinks (fold k's
#      members 0-9 as 10(k-1)..10k-1), which is all load_member_models globs for
#   3. xgb window heads per fold on the page genes (lf_winxgb.scan_genes), 4 at
#      a time -- these also carry each fold's CNN mean, for the page's provenance check
#   4. the all-50 xgb/CNN tracks: mean over all 50 members, sd recovered exactly
#      from each fold's 10-member mean and sd (ddof = 1)
set -eu
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
TAG=_nomask_flatten
export GENES=$(paste -sd, $S/page_genes.txt)

if [[ ${1:-} == scan ]]; then          # one fold/terminus scan (what xargs calls)
  k=$2 T=$3
  set +e                               # report a failed scan, don't abort the queue
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  python -m lf_winxgb.scan_genes --genes $GENES --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds \
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

for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P 4 -n 2 ${0:A} scan

source ~/.virtualenvs/r-tensorflow/bin/activate
TAG=$TAG python -I - <<'EOF'
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
