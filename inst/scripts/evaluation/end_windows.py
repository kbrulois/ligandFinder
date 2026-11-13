"""The 36 x 26 input window of every known peptide end, for a UMAP of the
flattened window features (Kevin, 2026-10-22).

Seen = the 5-fold CV knowns (lf_pepend_cv5_knowns_{C,N}.csv); never seen = the
two sets of cv5_unseen_ends.csv. Every window is built the way the scanner
builds them (lf_dcnn.scan.build_windows on lf_pepend_residues.npz), so seen and
never-seen ends go through one pipeline -- not the training arrays for one and
the scanner for the other. Each end carries its T = 0.2 call from end_calls.R
(linear head: held-out fold for seen ends, all-50 for never-seen).

    ~/.virtualenvs/r-tensorflow/bin/python end_windows.py
Writes end_windows.npz (X_C, X_N, channel names) and end_windows_meta.csv.
"""
import csv, os, sys
import numpy as np

## lf_dcnn from the checkout this script sits in, or $LF_SRC (lf_paths.py)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths; lf_paths.add_python()  # noqa: E402
from lf_dcnn.config import Config, PEPEND_ANCHOR  # noqa: E402
from lf_dcnn.scan import Residues, build_windows  # noqa: E402

A2 = os.path.expanduser("~/AF2_analysis"); S = A2 + "/lf_nomask_scripts_20261021"
acc = [l.strip() for l in open(S + "/residues_accession_order.txt")]
res = Residues.from_npz(A2 + "/lf_pepend_residues.npz")
assert len(acc) == len(res), (len(acc), len(res))
idx = {a: i for i, a in enumerate(acc)}

calls = {}
for r in csv.DictReader(open(A2 + "/lf_headcmp/end_calls_T0.2.csv")):
    if r["head"] == "cnn":
        calls[(r["set"], r["term"], r["accession"], int(r["anchor"]))] = (r["call"], r["reason"], r["score"])

ends, seen_keys = [], set()
def add(**e):
    k = (e["set"], e["term"], e["accession"], e["anchor"])
    if k in seen_keys: return
    seen_keys.add(k); ends.append(e)
for T in "CN":
    for r in csv.DictReader(open(f"{A2}/lf_pepend_cv5_knowns_{T}.csv")):
        add(set="seen", key="cross-validated", term=T, accession=r["accession"], anchor=int(r["anchor"]),
            gene=r["gene.x"], name=r["peps"], fold=r["fold"])
for r in csv.DictReader(open(S + "/cv5_unseen_ends.csv")):
    st = "never seen: peptides" if r["set"].startswith("never seen") else "never seen: large products"
    add(set=st, key=st, term=r["term"], accession=r["accession"], anchor=int(r["anchor"]),
        gene=r["gene"], name=r["name"], fold="")

out, meta = {}, []
for T in "CN":
    iso = f"{A2}/lf_pepend_run_{T}_cv5_f1_nomask_flatten_isolated"
    cfg = Config.from_json(iso + "/in/config.json")
    res.validate(cfg)
    pad = cfg.channel_names.index("padding")
    rows = [e for e in ends if e["term"] == T]
    X = np.stack([build_windows(res, idx[e["accession"]], np.array([e["anchor"]]), cfg,
                                PEPEND_ANCHOR[T], pad)[0] for e in rows])
    out[f"X_{T}"] = X.reshape(len(rows), -1)
    ## the same end's neighbours on its own precursor, anchor -3..+3: the windows
    ## the model has to tell the end apart from (35 of 36 residues shared at +/-1)
    SH = [-3, -2, -1, 1, 2, 3]
    Xs = np.stack([build_windows(res, idx[e["accession"]], np.array([e["anchor"] + d for d in SH]), cfg,
                                 PEPEND_ANCHOR[T], pad) for e in rows])
    out[f"Xshift_{T}"] = Xs.reshape(len(rows), len(SH), -1)
    out["shifts"] = np.array(SH)
    out["channels"] = np.array(cfg.channel_names)
    for e in rows:
        c = calls.get((e["key"], T, e["accession"], e["anchor"]), ("", "", ""))
        i_ = idx[e["accession"]]
        meta.append({**{k: e[k] for k in ("term", "set", "gene", "name", "accession", "anchor", "fold")},
                     "n_prot": int(res.n_prot[i_]), "c_prot": int(res.c_prot[i_]),
                     "call": c[0], "reason": c[1], "score": c[2]})
    print(T, X.shape, "windows;", sum(1 for e in rows if (e["key"], T, e["accession"], e["anchor"]) not in calls), "without a call")
np.savez_compressed(S + "/end_windows.npz", **out)
with open(S + "/end_windows_meta.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(meta[0])); w.writeheader(); w.writerows(meta)
print("channels:", list(out["channels"]))
