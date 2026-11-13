#!/usr/bin/env python3
"""Interactive version of end_window_umap.png (Kevin, 2026-10-22): hover any
end for what it is, both heads' calls at T = 0.2 and why, its 36-residue
window, and the nearest SEEN end in window space (PCA-50) with that end's call.
Coordinates come from end_umap_coords.csv, so the page matches the PNG.

    python3 end_umap_html.py      -> ~/AF2_analysis/end_window_umap.html
"""
import csv, json, os
import numpy as np
from sklearn.decomposition import PCA

A2 = os.path.expanduser("~/AF2_analysis"); S = A2 + "/lf_nomask_scripts_20261021"
PAGES = "pepend_pages_flatten_cv5"
pages = {f[:-5] for f in os.listdir(f"{A2}/{PAGES}") if f.endswith(".html") and f != "index.html"}
seqs = {r["accession"]: r["seq"] for r in csv.DictReader(open(S + "/end_window_seqs.csv"))}
co = list(csv.DictReader(open(S + "/end_umap_coords.csv")))
meta = list(csv.DictReader(open(S + "/end_windows_meta.csv")))
assert [(r["term"], r["accession"], r["anchor"], r["set"]) for r in co] == \
       [(r["term"], r["accession"], r["anchor"], r["set"]) for r in meta]
lens = {(r["accession"], r["term"], r["anchor"]): r["len"] for r in csv.DictReader(open(S + "/cv5_unseen_ends.csv"))}

## both heads' calls; seen ends on their held-out fold
calls = {}
for r in csv.DictReader(open(A2 + "/lf_headcmp/end_calls_T0.2.csv")):
    st = "seen" if r["set"] == "cross-validated" else r["set"]
    calls[(r["head"], st, r["term"], r["accession"], r["anchor"])] = r

def window(acc, anchor, term, n_prot, c_prot):
    s = seqs[acc]; a = int(anchor); pos = 28 if term == "C" else 8
    w = "".join(s[c - 1] if n_prot <= c <= c_prot else "·" for c in range(a - pos + 1, a - pos + 37))
    cut = pos if term == "C" else pos - 1          # C: after the anchor; N: before it
    return w[:cut], w[cut:]

z = np.load(S + "/end_windows.npz")
pts = []
for T in "CN":
    X = z[f"X_{T}"].astype("float64"); ii = [i for i, r in enumerate(co) if r["term"] == T]
    sd = X.std(0); keep = sd > 1e-9
    P = PCA(n_components=50, random_state=42).fit_transform((X[:, keep] - X[:, keep].mean(0)) / sd[keep])
    seen = np.array([co[i]["set"] == "seen" for i in ii])
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
    for j, i in enumerate(ii):
        r = co[i]; m = meta[i]
        dj = np.where(seen, D[j], np.inf); k = int(np.argmin(dj))
        nb = co[ii[k]]
        left, right = window(r["accession"], r["anchor"], T, int(m["n_prot"]), int(m["c_prot"]))
        c1 = calls[("cnn", r["set"], T, r["accession"], r["anchor"])]
        c2 = calls[("xgb_pm", r["set"], T, r["accession"], r["anchor"])]
        f = lambda c: dict(call=c["call"], reason=c["reason"],
                           score=round(float(c["score"]), 3) if c["score"] not in ("", "NA") else "n/a",
                           pc=c["peak_centre"], ph=round(float(c["peak_height"]), 2) if c["peak_height"] not in ("", "NA") else None)
        pts.append(dict(id=len(pts), t=T, set=r["set"], x=round(float(r["umap1"]), 3), y=round(float(r["umap2"]), 3),
                        gene=r["gene"], name=r["name"], acc=r["accession"], a=int(r["anchor"]), fold=r["fold"],
                        len=lens.get((r["accession"], T, r["anchor"]), ""), padded=int(r["padded"]),
                        lin=f(c1), xgb=f(c2), wl=left, wr=right,
                        nn=dict(name=nb["name"], gene=nb["gene"], a=int(nb["anchor"]), call=nb["call"],
                                d=round(float(dj[k]), 1), id=None, key=ii[k]),
                        page=f"{PAGES}/{r['gene']}.html" if r["gene"] in pages else ""))
key2id = {}
for p_, i in zip(pts, [i for T in "CN" for i, r in enumerate(co) if r["term"] == T]): key2id[i] = p_["id"]
for p_ in pts: p_["nn"]["id"] = key2id[p_["nn"].pop("key")]

html = open(S + "/end_umap_template.html").read().replace("/*DATA*/[]", json.dumps(pts, separators=(",", ":")))
out = A2 + "/end_window_umap.html"; open(out, "w").write(html)
print(f"{len(pts)} points -> {out} ({os.path.getsize(out) / 1024:.0f} KB)")
