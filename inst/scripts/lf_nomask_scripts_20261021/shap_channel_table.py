#!/usr/bin/env python3
"""Single-channel vs channel-family attribution, from shap_positions.py's
shap_channel_<arm>.csv -- the table behind "is relASA still the top channel?".

Kevin, 2026-11-10, on shap_positions_curated_AAswap_vs_esmonly36.png: "surprising
because previously relASA dominated as the overall most important metrics".
It still does. By SINGLE channel relASA is rank 1 in both models at both
termini (AAswap: 1 of 42, 13.5% of sum|phi| on C, 12.0% on N; curated: 18.4% /
18.5%), and its own attribution held up (AAswap / curated 1.38x on C, 0.94x on
N). The figure's family panel ADDS |phi| over a family's channels, which
lifts the 21-channel AA one-hot to the top (27-30% at ~1.3% per channel, net
signed near zero); relASA's share fell because the 21 new channels grew the
total (1.88x C, 1.46x N), not because the model leans on it less.

Not the Oct-20 grouped Shapley (pepend_shapley_importance.png, relASA 43-45%):
that credited channel GROUPS with the RANKING performance of the old curated
attention trunk by permutation. This is expected-gradients attribution of each
held-out positive's score -- another model and another question.

For each arm: the channel means over the five folds (set == "positive", the
held-out known ends; `padding` left out), the top channels by mean sum|phi|
with their share of the arm's total and net signed value, CHANNEL's rank, and
the families by total share, per-channel share and net signed. Then CHANNEL's
absolute attribution in every later arm against the FIRST arm, with each
arm's total and CHANNEL's share of the channels all the arms share.

    ~/.virtualenvs/lf-plm/bin/python shap_channel_table.py                # curated vs curated_AAswap
    ~/.virtualenvs/lf-plm/bin/python shap_channel_table.py ~/AF2_analysis/lf_headcmp curated curated_AAswap
    CHANNEL=mean_afm TOP=10 ~/.virtualenvs/lf-plm/bin/python shap_channel_table.py "" curated_AAswap esmonly36

Reads <dir>/shap_channel[_<arm>].csv ("curated" = the unsuffixed baseline file);
<dir> defaults to ~/AF2_analysis/lf_headcmp, shap_positions.py's --out-dir. Writes
$OUT_DIR (default <dir>)/shap_channel_table_<arm>_vs_<arm>.csv: arm, term, family,
channel, signed, sum_of_abs (fold means), share, rank. Standard library only, so
any venv runs it.
"""
import csv
import os
import pathlib
import sys
from collections import defaultdict

A2 = pathlib.Path.home() / "AF2_analysis"
## where shap_positions.py writes (its --out-dir), as shap_positions_plot.py reads
D = pathlib.Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 and sys.argv[1] else A2 / "lf_headcmp"
ARMS = sys.argv[2:] or ["curated", "curated_AAswap"]
OUT = pathlib.Path(os.environ.get("OUT_DIR", str(D))).expanduser()
CHANNEL = os.environ.get("CHANNEL", "relASA")
TOP = int(os.environ.get("TOP", "6"))
LABEL = {"curated": "curated", "curated_AAswap": "AAswap", "curated_wAAid": "curated + AA one-hot",
         "motifind": "curated + motif indicators", "esmonly36": "ESM-C only",
         "curated_AAswap_padg": "AAswap (padg windows)",
         "curated_AAswap_msacons_padg": "AAswap + MSA conservation (padg windows)"}


def sfx(arm):
    return "" if arm == "curated" else f"_{arm}"


def read_arm(arm):
    """(term, channel) -> {family, signed, sum_of_abs}: means over the folds of
    the held-out positives, padding left out; and the arm's channel count."""
    acc = defaultdict(lambda: {"signed": [], "sum_of_abs": [], "family": None})
    chans = set()
    with open(D / f"shap_channel{sfx(arm)}.csv") as f:
        for r in csv.DictReader(f):
            chans.add(r["channel"])
            if r["set"] != "positive" or r["channel"] == "padding":
                continue
            a = acc[(r["term"], r["channel"])]
            a["family"] = r["family"]
            a["signed"].append(float(r["signed"])); a["sum_of_abs"].append(float(r["sum_of_abs"]))
    mean = lambda x: sum(x) / len(x)
    return ({k: {"family": v["family"], "signed": mean(v["signed"]), "sum_of_abs": mean(v["sum_of_abs"])}
             for k, v in acc.items()}, len(chans))


def ranks(vals):
    """1 = largest; ties get their average rank, truncated (pandas rank().astype(int))."""
    out = {}
    for k, v in vals.items():
        above = sum(1 for w in vals.values() if w > v)
        same = sum(1 for w in vals.values() if w == v)
        out[k] = int(above + (same + 1) / 2)
    return out


tab = {arm: read_arm(arm) for arm in ARMS}
rows = []
for arm in ARMS:
    c, n_ch = tab[arm]
    fams = defaultdict(int)
    for (t, ch), v in c.items():
        if t == "C":
            fams[v["family"]] += 1
    print(f"\n######## {LABEL.get(arm, arm)} ({n_ch} ch): channels per family:", dict(sorted(fams.items())))
    for t in ("C", "N"):
        ct = {ch: v for (tt, ch), v in c.items() if tt == t}
        tot = sum(v["sum_of_abs"] for v in ct.values())
        rk = ranks({ch: v["sum_of_abs"] for ch, v in ct.items()})
        print(f"-- {t}: top {TOP} single channels by sum|phi| (share of all channels):")
        for ch, v in sorted(ct.items(), key=lambda kv: -kv[1]["sum_of_abs"])[:TOP]:
            print(f"     {ch:16} {v['family']:18} {v['sum_of_abs'] / tot:6.1%}   signed {v['signed']:+.2f}")
        if CHANNEL in ct:
            print(f"     {CHANNEL}: rank {rk[CHANNEL]} of {len(ct)}, {ct[CHANNEL]['sum_of_abs'] / tot:.1%}")
        fam = defaultdict(lambda: {"total": 0.0, "n": 0, "net": 0.0})
        for v in ct.values():
            g = fam[v["family"]]; g["total"] += v["sum_of_abs"]; g["n"] += 1; g["net"] += v["signed"]
        print("     by family:  share of sum|phi|   per-channel share   net signed")
        for fn, g in sorted(fam.items(), key=lambda kv: -kv[1]["total"]):
            print(f"       {fn:20} {g['total'] / tot:6.1%}  (n={g['n']:2d})   "
                  f"{g['total'] / g['n'] / tot:6.1%}         {g['net']:+.2f}")
        rows += [{"arm": arm, "term": t, "family": v["family"], "channel": ch, "signed": v["signed"],
                  "sum_of_abs": v["sum_of_abs"], "share": v["sum_of_abs"] / tot, "rank": rk[ch]}
                 for ch, v in sorted(ct.items(), key=lambda kv: rk[kv[0]])]

## CHANNEL's absolute attribution, every later arm against the first: did the
## model use it less, or did new channels only grow the total?
shared = set.intersection(*({ch for (_, ch) in tab[a][0]} for a in ARMS))
if len(ARMS) > 1 and CHANNEL in shared:
    lb = lambda a: LABEL.get(a, a)
    print(f"\n{CHANNEL} against {lb(ARMS[0])} ({len(shared)} channels shared by every arm):")
    for t in ("C", "N"):
        def stats(arm):
            ct = {ch: v["sum_of_abs"] for (tt, ch), v in tab[arm][0].items() if tt == t}
            return ct[CHANNEL], sum(ct.values()), sum(ct[ch] for ch in shared)
        r0, T0, s0 = stats(ARMS[0])
        for arm in ARMS[1:]:
            r1, T1, s1 = stats(arm)
            print(f"{t}: {CHANNEL} sum|phi| {lb(ARMS[0])} {r0:.2f} -> {lb(arm)} {r1:.2f} ({r1 / r0:.2f}x); "
                  f"all channels {T0:.2f} -> {T1:.2f} ({T1 / T0:.2f}x); "
                  f"{CHANNEL} share of the {len(shared)} shared channels {r0 / s0:.1%} -> {r1 / s1:.1%}")
elif len(ARMS) > 1:
    print(f"\n{CHANNEL} is not a channel of every arm; no comparison")

OUT.mkdir(parents=True, exist_ok=True)
out = OUT / f"shap_channel_table_{'_vs_'.join(ARMS)}.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print(f"\nwrote {out}")
