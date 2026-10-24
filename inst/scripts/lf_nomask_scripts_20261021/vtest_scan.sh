#!/bin/zsh
# Tiny end-to-end rehearsal of the all-secretome pipeline on 4 genes: 2 that the
# verified 104-gene gallery already has (HCRT, NPPA) and 2 it does not, scanned
# and pooled exactly as allsec_scan.sh does but into a _vtest npz set.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis; S=${0:A:h}; TAG=_nomask_flatten; SFX=_vtest
G=HCRT,NPPA,MUSTN1,DEFB113
if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  source ~/.virtualenvs/r-tensorflow/bin/activate; cd $W/inst/python
  python -m lf_winxgb.scan_genes --genes $G --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked \
    --out $A2/lf_winxgb_scan_${T}_cv5_f${k}${TAG}${SFX}.npz > $S/vtest_${T}_f$k.log 2>&1
  echo "vtest $T f$k rc=$?"
  exit 0
fi
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P 4 -n 2 ${0:A} one
source ~/.virtualenvs/r-tensorflow/bin/activate
TAG=$TAG SFX=$SFX GENES=$G python -I - <<'PY'
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis"); TAG=os.environ["TAG"]; SFX=os.environ["SFX"]
GENES = os.environ["GENES"].split(","); n, K = 10, 5
for T in ("C","N"):
    z=[np.load(f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}{SFX}.npz", allow_pickle=True) for k in range(1,K+1)]
    out={}
    for g in GENES:
        if f"{g}__anchor" not in z[0].files: continue
        an=z[0][f"{g}__anchor"]; out[f"{g}__anchor"]=an
        for h in ("cnn","xgb_pm"):
            M=np.stack([zz[f"{g}__{h}"].astype("f8") for zz in z])
            Sd=np.stack([zz[f"{g}__{h}__sd"].astype("f8") for zz in z])
            m=M.mean(0); var=((n-1)*(Sd**2).sum(0)+n*((M-m)**2).sum(0))/(n*K-1)
            out[f"{g}__{h}"]=m.astype("f4"); out[f"{g}__{h}__sd"]=np.sqrt(var).astype("f4")
    np.savez_compressed(f"{A2}/lf_winxgb_scan_{T}_cv5all50{TAG}{SFX}.npz", **out)
    print(f"{T}: {len(out)//5} genes pooled")
PY
echo "VTEST SCAN DONE"
