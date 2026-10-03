"""Train and apply the per-residue classifiers.

    python -m lf_resid train --targets pep_pocket,ct_context --seeds 5
    python -m lf_resid train --monitor val_binary_accuracy   # the R's own monitor

Writes one prediction column per (target, model) over the WHOLE residue table,
in ``lf_pepend_residues.npz`` row order, so the plotting code can index straight
into it.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np

from . import data as D
from . import models as M


def _pythonpath() -> str:
    """The tree this module was imported from, so the child imports the same one."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cur = os.environ.get("PYTHONPATH", "")
    return here + (os.pathsep + cur if cur else "")


def _metrics(y, p):
    from sklearn.metrics import average_precision_score, roc_auc_score
    if y.sum() == 0 or y.sum() == y.size:
        return {"auc": float("nan"), "pr_auc": float("nan"), "prevalence": float(y.mean())}
    return {"auc": float(roc_auc_score(y, p)),
            "pr_auc": float(average_precision_score(y, p)),
            "prevalence": float(y.mean())}


def _member_path(out_dir, target, kind, seed):
    return os.path.join(os.path.expanduser(out_dir),
                        f"member_{target}_{kind}_{seed:03d}.npz")


def cmd_member(a):
    """Train ONE model in this process and write it. Used by the isolated path."""
    import json
    x = np.load(os.path.expanduser(a.matrix), mmap_mode="r")
    idx = np.load(os.path.expanduser(a.index), allow_pickle=True)
    tr, va, mature = idx["tr"], idx["va"], idx["mature"]
    y = idx[f"y_{a.target}"].astype("float32")
    names = [str(n) for n in idx["feat_names"]]
    xtr, xva = np.asarray(x[tr]), np.asarray(x[va])
    ytr, yva = y[tr], y[va]

    if a.model == "mlp":
        m, _ = M.fit_mlp(xtr, ytr, xva, yva, epochs=a.epochs, monitor=a.monitor,
                         seed=a.seed, verbose=1 if a.verbose else 0)
        pv = M.predict_mlp(m, xva)
        pa = np.concatenate([M.predict_mlp(m, np.asarray(x[mature[i:i + 200000]]))
                             for i in range(0, mature.size, 200000)])
    else:
        b, _ = M.fit_xgb(xtr, ytr, xva, yva, nrounds=a.nrounds, seed=a.seed,
                         feat_names=names, verbose=bool(a.verbose))
        pv = M.predict_xgb(b, xva, names)
        pa = np.concatenate([M.predict_xgb(b, np.asarray(x[mature[i:i + 200000]]), names)
                             for i in range(0, mature.size, 200000)])
    out = _member_path(a.out_dir, a.target, a.model, a.seed)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez_compressed(out, p_val=pv.astype("float32"), p_all=pa.astype("float32"),
                        metrics=np.asarray(json.dumps(_metrics(yva, pv))))
    print(f"wrote {out}  " + json.dumps({k: round(v, 4) for k, v in
                                         _metrics(yva, pv).items()}), flush=True)


def cmd_train(a):
    d = D.load(a.residues_npz, a.labels_npz,
               context=() if a.no_context else D.DEFAULT_CONTEXT)
    print(f"residues {len(d):,} | features {d.x.shape[1]} "
          f"({'no context' if a.no_context else 'with +/-2 context'})", flush=True)

    rows, lab_prot, ctrl = D.training_rows(d, n_ctrl=a.n_ctrl, seed=a.seed)
    tr, va = D.split_by_precursor(d, rows, val_frac=a.val_frac, seed=a.seed)
    print(f"training rows {len(rows):,} over {lab_prot.size} labelled + {ctrl.size} control "
          f"precursors -> train {len(tr):,}, val {len(va):,} "
          f"(split on whole precursors)", flush=True)

    mature = np.flatnonzero(d.mature == 1)

    # ---- isolate every model in its own process -------------------------------
    # Keras state accumulates across models built in one process: the 13th MLP of
    # a run dies inside the PR-AUC metric, and with clear_session() it dies inside
    # Adam instead. lf_dcnn hit this and answered it the same way. The matrix is
    # cached once and memory-mapped so a subprocess costs a model fit, not a
    # feature build.
    work = os.path.expanduser(a.work)
    os.makedirs(work, exist_ok=True)
    mat = os.path.join(work, "x.npy")
    if a.refresh or not os.path.exists(mat):
        D.cache_matrix(d, mat)
    idx_p = os.path.join(work, "index.npz")
    np.savez_compressed(idx_p, tr=tr, va=va, mature=mature,
                        y_pep_pocket=d.y["pep_pocket"], y_ct_context=d.y["ct_context"],
                        feat_names=np.asarray(d.feat_names))
    print(f"isolated work dir: {work}", flush=True)

    out = {}
    report = {}
    for target in a.targets.split(","):
        target = target.strip()
        y = d.y[target].astype("float32")
        ytr, yva = y[tr], y[va]
        print(f"\n=== {target} === train positives {int(ytr.sum())}/{len(tr):,} "
              f"({100*ytr.mean():.2f}%) | val {int(yva.sum())}/{len(va):,} "
              f"({100*yva.mean():.2f}%)", flush=True)
        if ytr.sum() == 0 or yva.sum() == 0:
            print("  no positives on one side of the split -- skipped", flush=True)
            continue

        xtr, xva = d.x[tr], d.x[va]
        for kind in a.models.split(","):
            kind = kind.strip()
            ps_val, ps_all = [], []
            for s in range(a.seeds):
                mp = _member_path(work, target, kind, a.seed + s)
                if not os.path.exists(mp) or a.refresh:
                    cmd = [sys.executable, "-m", "lf_resid", "member",
                           "--matrix", mat, "--index", idx_p, "--out-dir", work,
                           "--target", target, "--model", kind,
                           "--seed", str(a.seed + s), "--epochs", str(a.epochs),
                           "--nrounds", str(a.nrounds), "--monitor", a.monitor]
                    if a.verbose:
                        cmd.append("--verbose")
                    rc = subprocess.call(cmd, env={**os.environ,
                                                   "PYTHONPATH": _pythonpath()})
                    if rc != 0 or not os.path.exists(mp):
                        raise SystemExit(f"member {target}/{kind}/{a.seed+s} failed (rc={rc})")
                else:
                    print(f"  {kind} seed {s}: reusing {os.path.basename(mp)}", flush=True)
                z = np.load(mp, allow_pickle=True)
                ps_val.append(z["p_val"]); ps_all.append(z["p_all"])

            pv = np.mean(ps_val, axis=0)
            mm = _metrics(yva, pv)
            report[f"{target}__{kind}"] = mm
            print(f"  {kind} ENSEMBLE of {a.seeds}: " +
                  json.dumps({k: round(v, 4) for k, v in mm.items()}), flush=True)

            full = np.zeros(len(d), dtype="float32")
            full[mature] = np.mean(ps_all, axis=0).astype("float32")
            out[f"{target}__{kind}"] = full
            out[f"{target}__{kind}__sd"] = np.zeros(len(d), dtype="float32")
            if a.seeds > 1:
                out[f"{target}__{kind}__sd"][mature] = \
                    np.std(ps_all, axis=0, ddof=1).astype("float32")

    if not out:
        raise SystemExit("nothing trained")
    val_mask = np.zeros(len(d), dtype="int8"); val_mask[va] = 1
    tr_mask = np.zeros(len(d), dtype="int8"); tr_mask[tr] = 1
    np.savez_compressed(os.path.expanduser(a.out), **out,
                        prot_idx=d.prot_idx, resno=d.resno, mature=d.mature,
                        accession=d.accession, gene=d.gene,
                        y_pocket=d.y["pep_pocket"], y_ct=d.y["ct_context"],
                        in_train=tr_mask, in_val=val_mask,
                        report=np.asarray(json.dumps(report)))
    print(f"\npredictions -> {a.out}", flush=True)
    print(json.dumps(report, indent=1), flush=True)


def main(argv=None):
    p = argparse.ArgumentParser(prog="lf_resid")
    sub = p.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train")
    t.add_argument("--residues-npz", default="~/AF2_analysis/lf_pepend_residues.npz")
    t.add_argument("--labels-npz", default="~/AF2_analysis/lf_resid_labels.npz")
    t.add_argument("--out", default="~/AF2_analysis/lf_resid_preds.npz")
    t.add_argument("--targets", default="pep_pocket,ct_context")
    t.add_argument("--models", default="mlp,xgb")
    t.add_argument("--seeds", type=int, default=5)
    t.add_argument("--epochs", type=int, default=100)
    t.add_argument("--nrounds", type=int, default=500)
    t.add_argument("--n-ctrl", type=int, default=100)
    t.add_argument("--val-frac", type=float, default=0.4)
    t.add_argument("--seed", type=int, default=42)
    t.add_argument("--monitor", default=M.MONITOR_DEFAULT,
                   help="MLP early-stopping monitor; the R used val_binary_accuracy")
    t.add_argument("--no-context", action="store_true",
                   help="one residue's own features only, no +/-2 neighbours")
    t.add_argument("--importance", action="store_true")
    t.add_argument("--verbose", action="store_true")
    t.add_argument("--work", default="~/AF2_analysis/lf_resid_isolated",
                   help="where the cached matrix and per-model members live")
    t.add_argument("--refresh", action="store_true",
                   help="retrain members that are already on disk")
    t.set_defaults(func=cmd_train)

    m = sub.add_parser("member", help="train one model (used by the isolated path)")
    m.add_argument("--matrix", required=True)
    m.add_argument("--index", required=True)
    m.add_argument("--out-dir", required=True)
    m.add_argument("--target", required=True)
    m.add_argument("--model", required=True)
    m.add_argument("--seed", type=int, required=True)
    m.add_argument("--epochs", type=int, default=100)
    m.add_argument("--nrounds", type=int, default=500)
    m.add_argument("--monitor", default=M.MONITOR_DEFAULT)
    m.add_argument("--verbose", action="store_true")
    m.set_defaults(func=cmd_member)
    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
