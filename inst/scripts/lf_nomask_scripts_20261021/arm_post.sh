#!/bin/zsh
# For one trained arm (tag $1, e.g. _nomask_flatten): xgb window heads on each fold for the four page
# genes, then both folds pooled into one 10-member ensemble per terminus.
set -eu
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
GENES=CXCL14,CXCL17,ANO8,NPPA
export TAG=$1

# every member must have been trained WITHOUT the mask -- a config that lost the
# field would load as the masked model and nobody would notice
for T in C N; do for F in A B; do
  c=$A2/lf_pepend_run_${T}_fold${F}${TAG}_isolated/in/config.json
  grep -q '"class_position_mask": false' $c || { echo "MASK STILL ON in $c"; exit 1; }
  n=$(ls $A2/lf_pepend_run_${T}_fold${F}${TAG}_isolated/out/members/member_*_${T}.weights.h5 | wc -l)
  [[ $n -eq 5 ]] || { echo "$T fold$F has $n members, want 5"; exit 1; }
done; done
echo "configs OK: class_position_mask false, 5 members each"

source ~/.virtualenvs/r-tensorflow/bin/activate
cd $W/inst/python
for T in C N; do for F in A B; do
  ( python -m lf_winxgb.scan_genes --genes $GENES --term $T \
      --run $A2/lf_pepend_run_${T}_fold${F}${TAG}.rds \
      --groups $A2/lf_pepend_window_groups_${T}_fold${F}.npz \
      --pm-rounds 100 --skip-stacked \
      --out $A2/lf_winxgb_scan_${T}_fold${F}${TAG}.npz \
      > $S/scan_${T}_fold${F}${TAG}.log 2>&1
    echo "scan $T fold$F rc=$?" ) &
done; done
wait

# one 10-member ensemble per terminus: fold A's members 0-4, fold B's as 5-9,
# symlinked under one members/ dir -- all load_member_models globs for
for T in C N; do
  E=$A2/lf_pepend_run_${T}_2fold${TAG}_isolated
  rm -rf $E; mkdir -p $E/in $E/out/members
  cp $A2/lf_pepend_run_${T}_foldA${TAG}_isolated/in/config.json $E/in/
  i=0
  for F in A B; do
    for f in $A2/lf_pepend_run_${T}_fold${F}${TAG}_isolated/out/members/member_*_${T}.weights.h5; do
      ln -s $f $E/out/members/$(printf "member_%03d_%s.weights.h5" $i $T); i=$((i+1))
    done
  done
  echo "$T: $i members -> $E"
done

python -I - <<'EOF'
# pool the two folds' per-gene npz: mean over all 10 members, and the sd over
# all 10 recovered exactly from each fold's 5-member mean and sd (ddof=1)
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis"); TAG = os.environ["TAG"]
for T in ("C", "N"):
    a = np.load(f"{A2}/lf_winxgb_scan_{T}_foldA{TAG}.npz", allow_pickle=True)
    b = np.load(f"{A2}/lf_winxgb_scan_{T}_foldB{TAG}.npz", allow_pickle=True)
    out = {}
    for g in ("CXCL14", "CXCL17", "ANO8", "NPPA"):
        assert np.array_equal(a[f"{g}__anchor"], b[f"{g}__anchor"]), g
        out[f"{g}__anchor"] = a[f"{g}__anchor"]
        for h in ("cnn", "xgb_pm"):
            ma, mb = a[f"{g}__{h}"].astype("f8"), b[f"{g}__{h}"].astype("f8")
            sa, sb = a[f"{g}__{h}__sd"].astype("f8"), b[f"{g}__{h}__sd"].astype("f8")
            m = (ma + mb) / 2
            var = (4 * sa**2 + 4 * sb**2 + 5 * (ma - m)**2 + 5 * (mb - m)**2) / 9
            out[f"{g}__{h}"] = m.astype("f4")
            out[f"{g}__{h}__sd"] = np.sqrt(var).astype("f4")
        print(f"{T} {g:6s} {out[g+'__anchor'].size:5d} anchors | cnn peak "
              f"{out[g+'__cnn'].max():.3f} @{out[g+'__anchor'][out[g+'__cnn'].argmax()]} | "
              f"xgb_pm peak {out[g+'__xgb_pm'].max():.3f} @{out[g+'__anchor'][out[g+'__xgb_pm'].argmax()]}")
    np.savez_compressed(f"{A2}/lf_winxgb_scan_{T}_2fold{TAG}.npz", **out)
EOF
echo "POST DONE"
