"""Every known-peptide page of the MSA-conservation gallery against the
independent per-gene scans: the popup payload's CNN surface, scored LIVE on the
page from the members, must equal the scan's value at the same anchors.

Checks per page: page present, payload present, file not truncated, the expected
22 series, and every per-fold series of BOTH arms equal to its own scan npz. It
does NOT check the held-out marks (each arm's on its own anchors, the G for an
amidated C end): those were checked by hand, on ADCYAP1, 2026-11-11.

    ~/.virtualenvs/lf-plm/bin/python msacons_gallery_qc.py
"""
import csv
import json
import os
import sys

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")
ARM = "curated_AAswap_msacons_padg"
OUT = f"{A2}/pepend_pages_{ARM}"
S = f"{A2}/lf_nomask_scripts_20261021"
ARMS = {"curated": "", "MSA conservation": f"_{ARM}"}   # label -> npz suffix
TOL = 1.5e-4

Z = {}
for lb, sfx in ARMS.items():
    for T in "NC":
        for k in range(1, 6):
            Z[(lb, T, k)] = np.load(
                f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}_nomask_flatten{sfx}.npz", allow_pickle=True)

genes = [g.strip() for g in open(f"{S}/page_genes.txt") if g.strip()]
want = {f"{T}" for T in "NC"} | {f"{T} fold {k}" for T in "NC" for k in range(1, 6)} \
       | {f"MSA conservation {T} fold {k}" for T in "NC" for k in range(1, 6)}

bad, worst, checked, anchors = [], 0.0, 0, 0
for g in genes:
    p = f"{OUT}/{g}.html"
    if not os.path.exists(p):
        bad.append((g, "page missing")); continue
    try:
        with open(p, "rb") as fh:
            pay = None
            for line in fh:
                if line.startswith(b"<script>window.LF_PEP = "):
                    pay = line[len(b"<script>window.LF_PEP = "):].rstrip().rstrip(b";"); break
            fh.seek(-16, 2); tail = fh.read()
        if pay is None:
            bad.append((g, "no popup payload")); continue
        if b"</html>" not in tail:
            bad.append((g, "TRUNCATED")); continue
        d = json.loads(pay)
    except Exception as e:
        bad.append((g, f"read error {e}")); continue

    missing = want - set(d["series"])
    if missing:
        bad.append((g, "series missing: " + ", ".join(sorted(missing)))); continue

    for lb in ARMS:
        for T in "NC":
            for k in range(1, 6):
                key = f"{T} fold {k}" if lb == "curated" else f"{lb} {T} fold {k}"
                v = d["series"][key]
                z = Z[(lb, T, k)]
                if f"{g}__anchor" not in z.files:
                    bad.append((g, f"{key}: gene not in scan")); continue
                an, cnn = z[f"{g}__anchor"], z[f"{g}__cnn"]
                a = np.asarray(v["anchor"])
                idx = np.searchsorted(an, a)
                if (idx >= len(an)).any() or not np.array_equal(an[idx], a):
                    bad.append((g, f"{key}: anchors do not match the scan")); continue
                e = float(np.abs(cnn[idx] - np.asarray(v["g"], dtype="f8")).max())
                worst = max(worst, e); checked += 1; anchors += len(a)
                if e > TOL:
                    bad.append((g, f"{key}: max|page-scan|={e:.4g}"))

print(f"pages: {len(genes)}; series compared: {checked:,}; anchors: {anchors:,}")
print(f"worst |page - scan| over every series of every page: {worst:.3g}  (tolerance {TOL})")
if bad:
    print(f"\n{len(bad)} PROBLEM(S):")
    for g, m in bad[:40]:
        print("  ", g, m)
    sys.exit(1)
print("every page's drawn surface matches its own scan")
