#!/bin/zsh
# Window-identity tables for the _esmonly36x arm (Kevin, 2026-10-24).
#
# WHY THIS IS NEEDED, and why reusing the curated ones was wrong.
# fullscan_esmonly36.sh reuses lf_pepend_window_groups_<T>_cv5_f<k>.npz with the
# comment "they carry the fold/member structure, not the channels, so the same
# ten files serve both arms". That is true for an arm that differs only in its
# CHANNELS -- which _esmonly36 did. It is false for this arm, which differs in
# its TRAINING SET: full_scan.py reads <T>__train__accession and indexes it with
# the model's own y_global, so the two must agree row for row. They did not:
#
#   IndexError: boolean index did not match indexed array along dimension 0;
#   dimension is 469 but corresponding boolean dimension is 461
#
# 10_9r_export_window_groups.R's own header warns about exactly this class of
# mismatch for --positives ("pass the same --positives you trained with"); the
# same applies to --input. So: same inputs as training, distinct --out.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
for k in 1 2 3 4 5; do for T in C N; do
  out=$A2/lf_pepend_window_groups_${T}_cv5_f${k}_esmcov.npz
  if [[ -s $out ]]; then echo "$T f$k cached"; continue; fi
  echo "--- $T f$k"
  cd $W && LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_9r_export_window_groups.R \
    --term $T \
    --input $A2/lf_pepend_nn_input_cv5_f${k}_esmcov.rds \
    --out $out 2>&1 | grep -v "built under R version"
done; done
echo
echo "=== shapes, which must equal the trained arrays.npz ==="
source ~/.virtualenvs/r-tensorflow/bin/activate
python3 - <<'PY'
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis")
bad = 0
for k in range(1, 6):
    for T in ("C", "N"):
        g = np.load(f"{A2}/lf_pepend_window_groups_{T}_cv5_f{k}_esmcov.npz", allow_pickle=True)
        z = np.load(f"{A2}/lf_pepend_run_{T}_cv5_f{k}_nomask_flatten_esmonly36x_isolated/in/arrays.npz")
        ok = True
        for s in ("train", "val", "all"):
            a, b = g[f"{T}__{s}__accession"].shape[0], z[f"{T}__{s}__x"].shape[0]
            if a != b: ok = False; print(f"  MISMATCH {T} f{k} {s}: groups {a} vs arrays {b}")
        print(f"  {T} f{k} {'OK' if ok else 'BAD'}  train "
              f"{g[f'{T}__train__accession'].shape[0]}, val {g[f'{T}__val__accession'].shape[0]}")
        bad += (not ok)
print("ALL ALIGNED" if bad == 0 else f"{bad} MISALIGNED -- do not scan")
PY
