"""CNN window score at EVERY mature anchor of the given precursors, one run.

The member mean of the global head, as on the page's CNN track, via the same
window builder lf_winxgb.scan_genes uses -- but keyed by accession (a gene name
can map to two precursors) and with no xgb fit.

    python score_anchors.py --run RUN.rds --term C --accessions A,B,... --out X.npz
"""
import argparse, os, sys
from pathlib import Path
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True)
ap.add_argument("--term", required=True)
ap.add_argument("--accessions", required=True, help="comma list, or @file with one per line")
ap.add_argument("--out", required=True)
ap.add_argument("--residues-npz", default="~/AF2_analysis/lf_pepend_residues.npz")
ap.add_argument("--labels-npz", default="~/AF2_analysis/lf_resid_labels.npz")
a = ap.parse_args()

## lf_dcnn from the checkout this script sits in, or $LF_SRC (lf_paths.py)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths; lf_paths.add_python()
from lf_dcnn.config import Config
from lf_dcnn.pipeline import predict_all
from lf_dcnn.scan import PEPEND_ANCHOR, Residues, build_windows, load_member_models

run = Path(os.path.expanduser(a.run))
iso = run.parent / (run.stem + "_isolated")
cfg = Config.from_json(iso / "in" / "config.json")
models = load_member_models(iso / "out", cfg, a.term)
print(f"{run.name}: {len(models)} members, head={cfg.global_head} x{cfg.head_lr_mult:g}, "
      f"mask={'on' if cfg.class_position_mask else 'off'}", flush=True)

accs = (Path(a.accessions[1:]).read_text().split() if a.accessions.startswith("@")
        else a.accessions.split(","))
res = Residues.from_npz(os.path.expanduser(a.residues_npz))
lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
acc_all = np.asarray([str(v) for v in lz["accession"]])
anchor = PEPEND_ANCHOR[a.term]
pad_col = cfg.channel_names.index("padding")

out, n = {}, 0
for acc in accs:
    hit = np.flatnonzero(acc_all == acc)
    if not hit.size:
        print(f"  {acc}: not in the residue cache -- skipped", flush=True); continue
    i = int(hit[0])
    anchors = np.arange(int(res.n_prot[i]), int(res.c_prot[i]) + 1)
    x = build_windows(res, i, anchors, cfg, anchor, pad_col)
    s = np.mean([np.asarray(predict_all(m, x)["global"], dtype="float32").reshape(-1)
                 for m in models], axis=0)
    out[f"{acc}__anchor"] = anchors.astype("int32")
    out[f"{acc}__cnn"] = s.astype("float32")
    n += anchors.size
np.savez_compressed(os.path.expanduser(a.out), **out)
print(f"{len(out)//2} precursors, {n} anchors -> {a.out}", flush=True)
