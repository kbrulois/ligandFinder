"""How well the per-index head calls each RESIDUE of a held-out window.

This is the head behind the per-residue class tracks on the protein pages, and
nothing in this session had measured it. Scored on the validation windows of
each trunk, 20 members averaged.

`none` is 70% of the residues -- every negative window is all `none` -- so
overall accuracy is dominated by it and is reported beside REAL-CLASS accuracy
(restricted to residues whose true label is one of the four real classes),
which is the number that says whether the head can label a peptide.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
## inst/python from the checkout this script sits in, or $LF_SRC (lf_paths.py).
## The installed package has no lf_winxgb and no lf_dcnn/scan.py, so it cannot
## stand in here; this used to be a hardcoded throwaway git worktree.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths; lf_paths.add_python()
from lf_winxgb.__main__ import _iso


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--term", default="C")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    from lf_dcnn.config import Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import PEPEND_ANCHOR, load_member_models

    iso = _iso(a.run)
    cfg = Config.from_json(iso / "in" / "config.json")
    models = load_member_models(iso / "out", cfg, a.term)
    pin = [str(c) for c in cfg.pi_names]
    anchor = PEPEND_ANCHOR[a.term]
    z = np.load(iso / "in" / "arrays.npz")
    X = np.ascontiguousarray(z[f"{a.term}__val__x"], dtype="float32")
    Y = z[f"{a.term}__val__y_per_index_cat"]
    known = z[f"{a.term}__val__y_global"].reshape(-1).astype(int)

    P = np.zeros((X.shape[0], cfg.seq_len, len(pin)), dtype="float64")
    for m in models:
        P += np.asarray(predict_all(m, X)["per_index_cat"], dtype="float64") / len(models)
    pred, true = P.argmax(-1), Y.argmax(-1)
    K = len(pin)
    pad_i, none_i = pin.index("padding"), pin.index("none")
    real = [i for i in range(K) if i not in (pad_i, none_i)]

    ok = pred == true
    not_pad = true != pad_i
    is_real = np.isin(true, real)
    out = dict(run=os.path.basename(a.run), term=a.term, classes=pin, anchor=int(anchor),
               n_windows=int(X.shape[0]), n_residues=int(true.size),
               acc_all=float(ok.mean()),
               acc_excl_padding=float(ok[not_pad].mean()),
               acc_real_only=float(ok[is_real].mean()),
               frac_real=float(is_real.mean()))
    ## per-class recall and precision
    out["per_class"] = {}
    for i, nm in enumerate(pin):
        t = true == i; p = pred == i
        out["per_class"][nm] = dict(
            n=int(t.sum()),
            recall=float(ok[t].mean()) if t.any() else None,
            precision=float((true[p] == i).mean()) if p.any() else None)
    ## accuracy along the window, over residues that are not padding
    by_pos = []
    for j in range(cfg.seq_len):
        m = not_pad[:, j]
        r = is_real[:, j]
        by_pos.append(dict(pos=j + 1,
                           acc=float(ok[:, j][m].mean()) if m.any() else None,
                           acc_real=float(ok[:, j][r].mean()) if r.any() else None,
                           n=int(m.sum()), n_real=int(r.sum())))
    out["by_position"] = by_pos
    ## accuracy on the windows that actually carry labels (knowns + near-end)
    lab_win = is_real.any(axis=1)
    out["acc_labelled_windows"] = float(ok[lab_win][not_pad[lab_win]].mean())
    out["n_labelled_windows"] = int(lab_win.sum())
    out["n_known_windows"] = int(known.sum())
    ## confusion, rows = true, cols = predicted
    cm = np.zeros((K, K), dtype=int)
    for t_, p_ in zip(true.ravel(), pred.ravel()): cm[t_, p_] += 1
    out["confusion"] = cm.tolist()
    json.dump(out, open(os.path.expanduser(a.out), "w"), indent=1)
    print(f"  {os.path.basename(a.run):42s} acc {out['acc_all']:.3f} | "
          f"excl padding {out['acc_excl_padding']:.3f} | real classes only {out['acc_real_only']:.3f} "
          f"({out['frac_real']*100:.0f}% of residues)", flush=True)
    for nm in pin:
        d = out["per_class"][nm]
        if d["n"]:
            print(f"      {nm:22s} n {d['n']:6,d}  recall {d['recall']:.3f}  "
                  f"precision {d['precision'] if d['precision'] is not None else float('nan'):.3f}", flush=True)


if __name__ == "__main__":
    main()
