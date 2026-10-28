#!/bin/zsh
# Per-gene window-head scans of the OTHER arms (ESM-C only, AAswap) for a list
# of genes (Kevin, 2026-10-28): CNN and xgb mean AND sd per fold, then the
# all-50 pool. The pages' arm tracks need these for their sd bands -- the arms'
# full-proteome tracks (lf_pepend_cv5_tracks_*_flatten_<arm>.rds) carry means
# only, and no full-proteome file has a CNN sd.
#
# The call is cv5_post_esmonly36.sh's scan step (and fullscan_arm.sh's): its CNN
# mean equals the arm tracks rds to 1e-7 and adds each fold's CNN sd. The xgb
# keys are then replaced by the full-proteome scan's mean and per-fold sd (see
# below), so every mean the pages draw is still the rds's. The files carry a
# suffix of their own, so the arms' 104-gene scans are never replaced:
#   lf_winxgb_scan_<T>_cv5_f<k>_nomask_flatten_<arm><SFX>.npz    (k = 1..5)
#   lf_winxgb_scan_<T>_cv5all50_nomask_flatten_<arm><SFX>.npz    (pooled)
# final_pages.R draws them with ARM_SCAN_SFX=<SFX>.
#
#   GENES_FILE=<one GENE or GENE@ACCESSION per line> SFX=_cand_x ./arm_gene_scans.sh
#   ARMS="esmonly36 curated_AAswap" PROCS=4
# GENE@ACCESSION picks one of several accessions under a name (COL18A1@P39060).
# Resumable: a job whose npz holds every listed gene is skipped (the xgb swap
# and the pool always re-run). Logs go to LOGDIR (default .render/armscan beside
# GENES_FILE); on success the suffix is written to arm_scan_sfx.txt beside
# GENES_FILE, where candidate_pages.sh reads it.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
SFX=${SFX:?set SFX -- never empty, the un-suffixed files are the 104-gene scans of each arm}
GENES_FILE=${GENES_FILE:?one GENE or GENE@ACCESSION per line}
GENES_FILE=${~GENES_FILE}
ARMS=(${=ARMS:-esmonly36 curated_AAswap})
PROCS=${PROCS:-4}
LOGDIR=${LOGDIR:-${GENES_FILE:h}/.render/armscan}
export SFX GENES_FILE LOGDIR

if [[ ${1:-} == scan ]]; then          # one arm / fold / terminus (what xargs calls)
  arm=$2 k=$3 T=$4
  TAG=_nomask_flatten_${arm}
  out=$A2/lf_winxgb_scan_${T}_cv5_f${k}${TAG}${SFX}.npz
  source ~/.virtualenvs/r-tensorflow/bin/activate
  # cached only if the file holds EVERY listed gene: a longer list under the
  # same SFX is scanned again, not silently half-covered
  if [[ -s $out ]] && python -I -c '
import sys, zipfile
have = {n[:-len("__anchor.npy")] for n in zipfile.ZipFile(sys.argv[1]).namelist() if n.endswith("__anchor.npy")}
want = {l.strip().split("@")[0] for l in open(sys.argv[2]) if l.strip()}
sys.exit(0 if want <= have else 1)' $out $GENES_FILE; then
    echo "  $arm $T f$k cached"; exit 0
  fi
  cd $W/inst/python
  python -m lf_winxgb.scan_genes --genes $(grep . $GENES_FILE | paste -sd, -) --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds \
    --residues-npz $A2/lf_pepend_residues_${arm}.npz \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked \
    --out ${out:r}.part.npz > $LOGDIR/armscan_${arm}_${T}_f${k}${SFX}.log 2>&1
  rc=$?
  # only a finished scan takes the real name, so a killed job is redone
  [[ $rc -eq 0 && -s ${out:r}.part.npz ]] && mv ${out:r}.part.npz $out
  echo "  $(date +%H:%M:%S) $arm $T f$k rc=$rc"
  exit 0
fi

[[ -s $GENES_FILE ]] || { echo "no genes in $GENES_FILE"; exit 1; }
for a in $ARMS; do
  [[ -s $A2/lf_pepend_residues_${a}.npz ]] || { echo "no residue table for $a"; exit 1; }
  for T in C N; do for k in 1 2 3 4 5; do
    [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_${a}.rds ]] || { echo "no $a trunk $T f$k"; exit 1; }
    [[ -s $A2/lf_winxgb_fullscan_${T}_cv5_f${k}_nomask_flatten_${a}.npz ]] || { echo "no $a full scan $T f$k"; exit 1; }
  done; done
done
mkdir -p $LOGDIR
echo "ARM SCANS $(date '+%F %H:%M:%S')  arms: $ARMS; $(grep -c . $GENES_FILE) genes; suffix $SFX; $PROCS at a time"
for a in $ARMS; do for k in 1 2 3 4 5; do for T in C N; do print -r -- "$a $k $T"; done; done; done |
  xargs -P $PROCS -n 3 ${0:A} scan

missing=0
for a in $ARMS; do for T in C N; do for k in 1 2 3 4 5; do
  [[ -s $A2/lf_winxgb_scan_${T}_cv5_f${k}_nomask_flatten_${a}${SFX}.npz ]] || { echo "MISSING $a $T f$k"; missing=1; }
done; done; done
(( missing )) && { echo "SCANS INCOMPLETE -- see $LOGDIR/armscan_*${SFX}.log, then re-run to resume"; exit 1; }

# xgb from the full-proteome scans. A per-gene scan refits each member's xgb
# head, and for these arms the refit drifts from the full scan's (ESM-C up to
# 0.04 on 2026-10-28; the CNN agrees to 1e-7) -- while the arm tracks rds, the
# candidate segments and every page drawn so far use the full scan's. So each
# fold file's xgb_pm / xgb_pm__sd are replaced by the full scan's mean and its
# own per-fold sd (xgb_pm_sd), looked up by the gene's accession; the CNN keys
# stay the per-gene scan's. Idempotent.
source ~/.virtualenvs/r-tensorflow/bin/activate
ARMS_LIST="${ARMS[*]}" ORDER=${ORDER:-$A2/lf_nomask_scripts_20261021/residues_accession_order.txt} python -I - <<'EOF' \
  || { echo "XGB SWAP FAILED -- not pooling; re-run to redo it"; exit 1; }
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis"); SFX = os.environ["SFX"]
order = {a.strip(): i for i, a in enumerate(open(os.environ["ORDER"]))}
for arm in os.environ["ARMS_LIST"].split():
    TAG = f"_nomask_flatten_{arm}"
    for T in ("C", "N"):
        for k in range(1, 6):
            fs = np.load(f"{A2}/lf_winxgb_fullscan_{T}_cv5_f{k}{TAG}.npz", allow_pickle=True)
            pi, an_all = fs["prot_idx"], fs["anchor"]
            p = f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}{SFX}.npz"
            z = dict(np.load(p, allow_pickle=True))
            genes = [key[:-len("__anchor")] for key in z if key.endswith("__anchor")]
            drift = 0.0
            for g in genes:
                acc = str(z[f"{g}__accession"])
                w = pi == order[acc]
                assert np.array_equal(an_all[w], z[f"{g}__anchor"]), (arm, T, k, g)
                drift = max(drift, float(np.abs(fs["xgb_pm"][w].astype("f8") - z[f"{g}__xgb_pm"]).max()))
                z[f"{g}__xgb_pm"] = fs["xgb_pm"][w].astype("f4")
                z[f"{g}__xgb_pm__sd"] = fs["xgb_pm_sd"][w].astype("f4")
            tmp = p[:-4] + ".part.npz"
            np.savez_compressed(tmp, **z); os.replace(tmp, p)
            print(f"  {arm} {T} f{k}: xgb from the full scan for {len(genes)} genes (largest change {drift:.3g})")
EOF

# the all-50 pool, as allsec_scan.sh / cv5_post_esmonly36.sh: the mean over all
# 50 members, the sd recovered from each fold's 10-member mean and sd
ARMS_LIST="${ARMS[*]}" python -I - <<'EOF' || { echo "POOL FAILED -- re-run to redo it"; exit 1; }
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis"); SFX = os.environ["SFX"]
GENES = [l.strip().split("@")[0] for l in open(os.environ["GENES_FILE"]) if l.strip()]
n, K = 10, 5
for arm in os.environ["ARMS_LIST"].split():
    TAG = f"_nomask_flatten_{arm}"
    for T in ("C", "N"):
        z = [np.load(f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}{SFX}.npz", allow_pickle=True) for k in range(1, K + 1)]
        out, miss = {}, []
        for g in GENES:
            if f"{g}__anchor" not in z[0].files: miss.append(g); continue
            an = z[0][f"{g}__anchor"]
            assert all(np.array_equal(an, zz[f"{g}__anchor"]) for zz in z), g
            out[f"{g}__anchor"] = an
            out[f"{g}__accession"] = z[0][f"{g}__accession"]
            for h in ("cnn", "xgb_pm"):
                M = np.stack([zz[f"{g}__{h}"].astype("f8") for zz in z])
                Sd = np.stack([zz[f"{g}__{h}__sd"].astype("f8") for zz in z])
                m = M.mean(0)
                var = ((n - 1) * (Sd ** 2).sum(0) + n * ((M - m) ** 2).sum(0)) / (n * K - 1)
                out[f"{g}__{h}"] = m.astype("f4"); out[f"{g}__{h}__sd"] = np.sqrt(var).astype("f4")
        p = f"{A2}/lf_winxgb_scan_{T}_cv5all50{TAG}{SFX}.npz"
        np.savez_compressed(p, **out)
        print(f"  {arm} {T}: {len(GENES) - len(miss)} genes pooled over 5 folds x 10 members -> {os.path.basename(p)}"
              + (f"; not in the scan: {', '.join(miss)}" if miss else ""))
EOF
print -r -- $SFX > ${GENES_FILE:h}/arm_scan_sfx.txt
echo "ARM SCANS DONE $(date '+%F %H:%M:%S')  (suffix recorded in ${GENES_FILE:h}/arm_scan_sfx.txt)"
