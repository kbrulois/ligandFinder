"""Hyperparameter search for the peptide-end window CNN.

Tunes the learning rate, the U-Net's width (filters per level), the kernel size
and early stopping (patience, start epoch, epoch budget, monitor) over the
cross-validation folds 10_7_pepend_windows.R built, scoring every candidate on
known peptide ends that neither trained it nor chose its stopping epoch.

    python -m lf_dcnn tune plan   --root R --inputs R/in --n-configs 24
    python -m lf_dcnn tune run    --root R --workers 4
    python -m lf_dcnn tune report --root R
    python -m lf_dcnn tune extend --root R --top 5 --seeds 2,3     # then run again

``R/in/<fold>/`` is what ``10_8_pepend_train.R --export-only`` writes for one
fold: the array contract and the resolved Config that 10_8 would have trained,
plus window metadata. inst/scripts/10_8t_pepend_tune.sh drives all of it.

What is held out from what
--------------------------
Each fold trains on its own ``train`` split, unchanged -- the same windows and
positives the production trunk of that fold trains on. Its ``val`` split (the
fold's fifth of the known ends plus ~400 hard negatives) is cut into
``n_parts`` parts, and every rotation early-stops on all parts but one and is
scored on the one held out. Parts are drawn so that no precursor, and no pair
of known windows within ``sim_h`` sequence distance, straddles two parts, and a
near-end negative goes where its parent end goes (it shares 34 of 36 residues
with it). So every known end of the fold is scored exactly once, by a model
whose stopping epoch was chosen on other ends.

Why not just score on the windows early stopping watched: restoring the best
epoch takes a maximum over the monitor curve, and a longer patience takes it
over a longer stretch -- on the same windows it can only look better. That is
the trap the xgb head fell into (stopping on val inflated val PR AUC while
retrieval got worse).

Early stopping is replayed, not retrained
-----------------------------------------
Early stopping does not change the path training takes. Each trial trains once,
logs the monitor on every rotation's monitor windows after every epoch, and
replays every (monitor, start, patience, cap) setting online
(:class:`lf_dcnn.es_replay.ESGridTracker`, keras' own logic). When a setting's
answer is final, the weights of the epoch it restores -- kept in memory while
they could still be restored -- are scored. One run therefore answers every
setting, exactly, and settings are compared on the same trajectory.

The score
---------
Each held-out end is ranked among the fold's UNKNOWN candidate windows (the
``all`` split's ~41.8k motif candidates), as 10_8 ranks held-out ends:
``1 + #(candidates scoring higher)``. ``rank_fair`` also drops from that pool
every candidate sitting exactly on a real peptide end that is not a label here
(curated ligand or UniProt peptide/propeptide ends), which counting as false
positives would penalise (``10_9z_locus_pr.py::load_ends``). The objective is
the mean log10 rank (lower is better), reported as the geometric-mean rank.
AP and ROC AUC against the held-out part's own hard negatives -- the val-PR-AUC
view, which rewards boundary precision more than ranking -- are reported beside
it, never instead of it.

Caveats the report repeats
--------------------------
* The monitor sees (n_parts - 1)/n_parts of val; production watches all of it.
* Selection uses every fold's ends, so a model retrained with the winner and
  scored on these same folds is optimistic by the winner's curse. The report's
  leave-one-fold-out check estimates how much; the never-seen ends, untouched
  here, are the clean test.
* One trial = one model = one process (the second model in a keras process can
  die mid-fit), N processes in parallel. Trials are resumable: a trial with a
  ``result.json`` is never rerun.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .config import Config
from .es_replay import (
    MONITORS,
    ESGrid,
    ESGridTracker,
    average_precision,
    keras_pr_auc,
    rank_among,
    roc_auc,
)

PLAN_FILE = "plan.json"
TRIALS_DIR = "trials"
FOLDS_DIR = "folds"
REPORT_DIR = "report"
RESULT_FILE = "result.json"
CURVES_FILE = "curves.npz"
#: written by 10_8_pepend_train.R --export-only beside arrays.npz / config.json
VAL_META_FILE = "val_meta.parquet"
ALL_META_FILE = "meta.parquet"           # lf_dcnn_export's one row per `all` window

#: the default search space. ``unet_width`` is the first U-Net level's filter
#: count; every level below doubles it, so 16 is the production (16, 32, 64).
DEFAULT_SPACE: dict = {
    "learning_rate": {"type": "loguniform", "low": 1e-4, "high": 3e-3},
    "unet_width": {"type": "choice", "values": [8, 12, 16, 24, 32]},
    "conv_kernel": {"type": "choice", "values": [3, 5, 7, 9]},
}

#: log10-spaced edges for the learning-rate marginal in the report
LR_BIN_EDGES = (1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 1e-2)


# ------------------------------------------------------------ small helpers ---

def _write_json(path: Path, obj) -> None:
    """Write atomically: a crash mid-write must not leave a file that looks done."""
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def _sha1(path: Path) -> str:
    return hashlib.sha1(Path(path).read_bytes()).hexdigest()


def _fmt_time(sec: float) -> str:
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h{m:02d}m" if h else f"{m}m{s:02d}s"


def load_plan(root) -> dict:
    with open(Path(root).expanduser() / PLAN_FILE) as fh:
        return json.load(fh)


def trial_dir(root, config_id: str, fold: str, seed: int) -> Path:
    return Path(root).expanduser() / TRIALS_DIR / config_id / f"{fold}_s{int(seed)}"


def plan_trials(plan: Mapping) -> list[tuple[str, str, int]]:
    """Every (config, fold, seed), config-major with the baseline first, so an
    interrupted run leaves whole configs rather than a thin slice of each."""
    out = []
    for c in plan["configs"]:
        for seed in c["seeds"]:
            for fold in plan["folds"]:
                out.append((c["id"], fold, int(seed)))
    return out


def _config_entry(plan: Mapping, config_id: str) -> dict:
    for c in plan["configs"]:
        if c["id"] == config_id:
            return c
    raise KeyError(f"no config {config_id!r} in the plan")


# ----------------------------------------------------------- search space ---

def params_to_overrides(params: Mapping, base: Config) -> dict:
    """Search parameters -> Config fields. ``unet_width`` w becomes
    ``unet_filters = (w, 2w, 4w, ...)`` at the base trunk's depth."""
    out = {}
    for k, v in params.items():
        if k == "unet_width":
            out["unet_filters"] = [int(v) * 2**i for i in range(len(base.unet_filters))]
        elif k in Config.__dataclass_fields__:
            out[k] = v
        else:
            raise ValueError(f"unknown search parameter {k!r}: neither unet_width nor a Config field")
    return out


def base_params(base: Config, space: Mapping) -> dict:
    """The production values of the searched parameters."""
    out = {}
    for k in space:
        if k == "unet_width":
            w = int(base.unet_filters[0])
            if tuple(base.unet_filters) != tuple(w * 2**i for i in range(len(base.unet_filters))):
                raise ValueError(f"base unet_filters {base.unet_filters} do not double per level, "
                                 "so unet_width cannot describe them; search unet_filters directly")
            out[k] = w
        else:
            v = getattr(base, k)
            out[k] = list(v) if isinstance(v, tuple) else v
    return out


def config_tag(params: Mapping) -> str:
    bits = []
    for k, v in params.items():
        if k == "learning_rate":
            bits.append(f"lr{v:.2g}")
        elif k == "unet_width":
            bits.append(f"w{v}")
        elif k == "conv_kernel":
            bits.append(f"k{v}")
        else:
            bits.append(f"{k}={v}")
    return " ".join(bits)


def config_id(params: Mapping) -> str:
    return "c" + hashlib.sha1(json.dumps(params, sort_keys=True).encode()).hexdigest()[:7]


def sample_configs(space: Mapping, n: int, seed: int = 1, kind: str = "lhs") -> list[dict]:
    """``n`` parameter sets.

    ``lhs`` (default): a balanced random design -- every continuous dimension
    gets one draw from each of ``n`` equal strata (log-spaced for
    ``loguniform``), every categorical one each level ``n/levels`` times, the
    dimensions shuffled independently. Marginals are therefore balanced, which
    is what makes the per-parameter tables in the report readable at small n.
    ``grid``: the full factorial over each dimension's ``values`` (``n`` ignored).
    """
    rng = np.random.default_rng(seed)
    if kind == "grid":
        import itertools

        names = list(space)
        vals = [space[k].get("values") for k in names]
        if any(v is None for v in vals):
            raise ValueError("grid sampling needs explicit `values` for every dimension")
        return [dict(zip(names, combo)) for combo in itertools.product(*vals)]
    if kind != "lhs":
        raise ValueError(f"sampler must be 'lhs' or 'grid', got {kind!r}")
    cols = {}
    for name, spec in space.items():
        t = spec["type"]
        if t in ("loguniform", "uniform"):
            u = (rng.permutation(n) + rng.random(n)) / n
            lo, hi = float(spec["low"]), float(spec["high"])
            v = (np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo))) if t == "loguniform"
                 else lo + u * (hi - lo))
            cols[name] = [float(f"{x:.3g}") for x in v]
        elif t == "choice":
            levels = list(spec["values"])
            reps = -(-n // len(levels))
            idx = np.concatenate([rng.permutation(len(levels)) for _ in range(reps)])[:n]
            rng.shuffle(idx)
            cols[name] = [levels[i] for i in idx]
        else:
            raise ValueError(f"{name}: unknown type {t!r} (loguniform, uniform or choice)")
    return [{k: cols[k][i] for k in space} for i in range(n)]


# ------------------------------------------------------- folds and parts ---

def _seq_dist(a: str, b: str) -> float:
    """Proportion of differing residues over the positions both windows cover
    (``X`` = outside the precursor), as ``ape::dist.aa(scaled = TRUE)``."""
    n = d = 0
    for x, y in zip(a, b):
        if x == "X" or y == "X":
            continue
        n += 1
        d += x != y
    return d / n if n else 1.0


def val_groups(meta, y: np.ndarray, sim_h: float = 0.3) -> np.ndarray:
    """Group id per val window; windows of one group must share a part.

    Known ends join a group when they share a precursor or sit within
    ``sim_h`` sequence distance (single linkage, so no such pair can be split
    -- 10_7 clusters the folds at the same 0.3). A near-end negative takes its
    parent's group; every other negative gets -1 and is free.
    """
    n = len(y)
    g = np.full(n, -1, dtype="int64")
    pos = np.flatnonzero(y == 1)
    parent = list(range(pos.size))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    acc = meta["accession"].to_numpy()[pos]
    aa = meta["AA"].to_numpy()[pos] if "AA" in meta else None
    for i in range(pos.size):
        for j in range(i + 1, pos.size):
            if acc[i] == acc[j] or (aa is not None and _seq_dist(aa[i], aa[j]) <= sim_h):
                parent[find(i)] = find(j)
    roots = [find(i) for i in range(pos.size)]
    ids = {r: k for k, r in enumerate(dict.fromkeys(roots))}
    g[pos] = [ids[r] for r in roots]

    if "win_type" in meta and "parent" in meta:
        key = {(a, int(an)): g[i] for i, a, an in
               zip(pos, acc, meta["anchor"].to_numpy()[pos])}
        wt = meta["win_type"].to_numpy()
        par = meta["parent"].to_numpy()
        accs = meta["accession"].to_numpy()
        for i in np.flatnonzero(wt == "near_end"):
            if par[i] == par[i]:                      # not NaN
                g[i] = key.get((accs[i], int(par[i])), -1)
    return g


def assign_parts(meta, y: np.ndarray, groups: np.ndarray, n_parts: int, seed: int) -> np.ndarray:
    """Deal the val windows into ``n_parts`` parts.

    Groups of known ends go largest first to the part holding the fewest known
    ends; free negatives are then dealt per window type to the part holding the
    fewest negatives, so every part keeps roughly val's own negative mix.
    """
    rng = np.random.default_rng(seed)
    n = len(y)
    part = np.full(n, -1, dtype="int64")
    gids = np.unique(groups[groups >= 0])
    if gids.size < n_parts:
        raise ValueError(f"only {gids.size} independent groups of known ends in val, "
                         f"too few for {n_parts} parts; use fewer parts")
    npos = np.array([int(((groups == g) & (y == 1)).sum()) for g in gids])
    order = gids[rng.permutation(gids.size)]
    order = order[np.argsort(-npos[np.searchsorted(gids, order)], kind="stable")]
    load_pos = np.zeros(n_parts)
    load_all = np.zeros(n_parts)
    for g in order:
        members = groups == g
        j = min(range(n_parts), key=lambda k: (load_pos[k], load_all[k], k))
        part[members] = j
        load_pos[j] += int((members & (y == 1)).sum())
        load_all[j] += int(members.sum())

    free = np.flatnonzero(part < 0)
    wt = meta["win_type"].to_numpy() if "win_type" in meta else np.array(["neg"] * n)
    n_neg = np.array([int(((part == k) & (y == 0)).sum()) for k in range(n_parts)], dtype=float)
    for t in sorted(set(wt[free])):
        for i in rng.permutation(free[wt[free] == t]):
            j = min(range(n_parts), key=lambda k: (n_neg[k], k))
            part[i] = j
            n_neg[j] += 1
    for k in range(n_parts):
        sel = part == k
        if not (y[sel] == 1).any() or not (y[sel] == 0).any():
            raise ValueError(f"part {k} lacks positives or negatives; use fewer parts")
    return part


def other_end_mask(all_meta, term: str, other_ends: Path) -> np.ndarray:
    """Unknown candidates sitting exactly on a real peptide end that is not a
    label here -- curated ligand ends and UniProt peptide/propeptide ends, as in
    10_9z_locus_pr.py::load_ends. ``other_ends`` holds ``ligand_list.csv`` and
    ``uniprot_peptides.csv`` (~/AF2_analysis/pepend_locus/)."""
    import pandas as pd

    col = "end" if term == "C" else "start"
    ends = set()
    ll = pd.read_csv(Path(other_ends) / "ligand_list.csv")
    up = pd.read_csv(Path(other_ends) / "uniprot_peptides.csv")
    up = up[up["type"].isin(["peptide", "propeptide"])]
    for d in (ll, up):
        d = d[["accession", col]].dropna()
        ends.update(zip(d["accession"].astype(str), d[col].astype(int)))
    acc = all_meta["accession"].astype(str).to_numpy()
    anc = all_meta["anchor"].astype(int).to_numpy()
    hit = np.fromiter(((a, x) in ends for a, x in zip(acc, anc)), dtype=bool, count=acc.size)
    return hit & (all_meta["known"].to_numpy() == 0)


def prepare_fold(name: str, in_dir: Path, out_dir: Path, n_parts: int, part_seed: int,
                 other_ends: Path | None = None, sim_h: float = 0.3) -> dict:
    """Check one exported fold and fix its parts and candidate pools for the
    whole search (every trial of the fold must see the same parts)."""
    import pandas as pd

    cfg = Config.from_json(in_dir / "config.json")
    if len(cfg.term_order) != 1:
        raise ValueError(f"{in_dir}: one terminus per fold, got {cfg.term_order}")
    term = cfg.term_order[0]
    with np.load(in_dir / "arrays.npz") as z:
        val_x = z[f"{term}__val__x"]
        val_y = z[f"{term}__val__y_global"][:, 0]
        all_x = z[f"{term}__all__x"]
        all_y = z[f"{term}__all__y_global"][:, 0]
        train_y = z[f"{term}__train__y_global"][:, 0]
    vm = pd.read_parquet(in_dir / VAL_META_FILE)
    am = pd.read_parquet(in_dir / ALL_META_FILE)
    if len(vm) != val_y.size or len(am) != all_y.size:
        raise ValueError(f"{in_dir}: metadata rows ({len(vm)} val, {len(am)} all) do not match "
                         f"the arrays ({val_y.size}, {all_y.size})")
    if not (np.array_equal(vm["known"].to_numpy(), val_y) and np.array_equal(am["known"].to_numpy(), all_y)):
        raise ValueError(f"{in_dir}: metadata `known` disagrees with y_global -- row order differs")

    pos = np.flatnonzero(val_y == 1)
    row = {p: i for i, p in enumerate(am["peps"].to_numpy())}
    missing = [p for p in vm["peps"].to_numpy()[pos] if p not in row]
    if missing:
        raise ValueError(f"{in_dir}: val ends not in the `all` split: {missing[:5]}")
    pos_all = np.array([row[p] for p in vm["peps"].to_numpy()[pos]], dtype="int64")
    ## the rank of a held-out end is read from its `all` row; that is only the
    ## same thing as its val window if the two windows are identical
    if not np.allclose(val_x[pos], all_x[pos_all], atol=1e-6):
        raise ValueError(f"{in_dir}: a val end's window differs from its `all` copy")

    groups = val_groups(vm, val_y, sim_h=sim_h)
    parts = assign_parts(vm, val_y, groups, n_parts, part_seed)
    pool = all_y == 0
    fair = pool & ~other_end_mask(am, term, other_ends) if other_ends else None

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / f"{name}.npz", parts=parts, groups=groups, pos_rows=pos,
             pos_all_rows=pos_all, pool=pool,
             **({"pool_fair": fair} if fair is not None else {}))
    return {
        "dir": str(in_dir), "term": term,
        "fingerprint": {"config_sha1": _sha1(in_dir / "config.json"),
                        "arrays_bytes": (in_dir / "arrays.npz").stat().st_size},
        "n_val": int(val_y.size), "n_val_pos": int(pos.size),
        "parts_pos": [int(((parts == k) & (val_y == 1)).sum()) for k in range(n_parts)],
        "parts_neg": [int(((parts == k) & (val_y == 0)).sum()) for k in range(n_parts)],
        "n_groups": int(np.unique(groups[groups >= 0]).size),
        "n_pool": int(pool.sum()),
        "n_pool_fair": int(fair.sum()) if fair is not None else None,
        "n_train": int(train_y.size), "n_train_pos": int(train_y.sum()),
    }


# ------------------------------------------------------------------- plan ---

def make_plan(root, inputs: Sequence[Path] | Path, n_configs: int = 24, seeds=(1,),
              space: Mapping | None = None, sampler: str = "lhs", sampler_seed: int = 1,
              grid: ESGrid | None = None, n_parts: int = 4, part_seed: int = 20261102,
              other_ends: Path | None = None, sim_h: float = 0.3, force: bool = False,
              baseline_setting: Mapping | None = None, verbose: int = 1) -> dict:
    """Write ``<root>/plan.json``: the folds, their parts and pools, the
    early-stopping grid and the configs (the production baseline first).

    ``baseline_setting`` (monitor/start/patience/cap) is what every result is
    compared against; it defaults to the base Config's own early stopping,
    i.e. production's, and is always added to the grid."""
    root = Path(root).expanduser()
    if (root / PLAN_FILE).exists() and not force:
        raise SystemExit(f"{root / PLAN_FILE} exists; `tune extend` adds seeds, "
                         "a new search needs a new --root (or --force)")
    if (root / TRIALS_DIR).exists() and any((root / TRIALS_DIR).iterdir()) and force:
        raise SystemExit(f"{root / TRIALS_DIR} already holds results; refusing to re-plan "
                         "over them (their parts and grid would no longer match)")
    if isinstance(inputs, (str, Path)):
        d = Path(inputs).expanduser()
        fold_dirs = sorted(p for p in d.iterdir() if (p / "arrays.npz").exists())
    else:
        fold_dirs = [Path(p).expanduser() for p in inputs]
    if not fold_dirs:
        raise SystemExit(f"no fold inputs (directories with arrays.npz) under {inputs}")

    bases = [Config.from_json(d / "config.json") for d in fold_dirs]
    first = bases[0].to_dict()
    for d, b in zip(fold_dirs[1:], bases[1:]):
        diff = {k for k, v in b.to_dict().items() if v != first.get(k)}
        if diff:
            raise SystemExit(f"{d}: config differs from {fold_dirs[0]} in {sorted(diff)}; "
                             "every fold must be the same model")
    base = bases[0]
    space = dict(space or DEFAULT_SPACE)
    if baseline_setting is None:
        baseline_setting = {"monitor": "pr_auc", "start": base.start_from_epoch,
                            "patience": base.patience, "cap": base.epochs}
    baseline_setting = {"monitor": str(baseline_setting["monitor"]),
                        **{k: int(baseline_setting[k]) for k in ("start", "patience", "cap")}}
    grid = (grid or ESGrid()).with_setting(**baseline_setting)

    bp = base_params(base, space)
    configs = [{"id": "base", "tag": config_tag(bp) + " (production)", "params": bp,
                "overrides": params_to_overrides(bp, base), "seeds": [int(s) for s in seeds]}]
    seen = {json.dumps(bp, sort_keys=True)}
    for p in sample_configs(space, n_configs, seed=sampler_seed, kind=sampler):
        k = json.dumps(p, sort_keys=True)
        if k in seen:
            continue
        seen.add(k)
        configs.append({"id": config_id(p), "tag": config_tag(p), "params": p,
                        "overrides": params_to_overrides(p, base),
                        "seeds": [int(s) for s in seeds]})
    for c in configs:            # fail now, not in trial 37
        base.evolve(**c["overrides"])

    folds = {}
    for d in fold_dirs:
        folds[d.name] = prepare_fold(d.name, d, root / FOLDS_DIR, n_parts, part_seed,
                                     other_ends=other_ends, sim_h=sim_h)
    plan = {
        "version": 1,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "term": base.term_order[0],
        "folds": folds,
        "grid": grid.to_dict(),
        "n_parts": int(n_parts),
        "part_seed": int(part_seed),
        "sim_h": float(sim_h),
        "other_ends": str(other_ends) if other_ends else None,
        "space": space,
        "sampler": {"kind": sampler, "n": int(n_configs), "seed": int(sampler_seed)},
        "baseline": {"config": "base", "setting": baseline_setting},
        "configs": configs,
    }
    root.mkdir(parents=True, exist_ok=True)
    _write_json(root / PLAN_FILE, plan)
    if verbose:
        print(f"plan -> {root / PLAN_FILE}")
        print(f"  terminus {plan['term']}; {len(folds)} folds; {len(configs)} configs "
              f"(baseline + {len(configs) - 1}); seeds {list(seeds)}; "
              f"{len(plan_trials(plan))} trials")
        print(f"  early stopping replayed: monitors {list(grid.monitors)} x starts "
              f"{list(grid.starts)} x patiences {list(grid.patiences)} x caps {list(grid.caps)}")
        for f, s in folds.items():
            print(f"  {f}: val {s['n_val_pos']} ends in {s['n_groups']} groups -> parts "
                  f"{s['parts_pos']} ends / {s['parts_neg']} negatives; pool {s['n_pool']}"
                  + (f", fair {s['n_pool_fair']}" if s["n_pool_fair"] is not None else ""))
    return plan


def extend_plan(root, top: int = 5, seeds=(2, 3), objective: str | None = None,
                verbose: int = 1) -> dict:
    """Add ``seeds`` to the baseline and the ``top`` configs of the current
    report (each judged at its own best early-stopping setting)."""
    root = Path(root).expanduser()
    plan = load_plan(root)
    rep = summarise(root, objective=objective)
    best_per_cfg = (rep["table"].sort_values("objective")
                    .drop_duplicates("config"))
    keep = ["base"] + [c for c in best_per_cfg["config"] if c != "base"][:top]
    added = 0
    for c in plan["configs"]:
        if c["id"] in keep:
            new = [int(s) for s in seeds if int(s) not in c["seeds"]]
            c["seeds"] = c["seeds"] + new
            added += len(new) * len(plan["folds"])
    shutil.copy(root / PLAN_FILE, root / (PLAN_FILE + f".{time.strftime('%Y%m%d_%H%M%S')}.bak"))
    _write_json(root / PLAN_FILE, plan)
    if verbose:
        tags = {c["id"]: c["tag"] for c in plan["configs"]}
        print(f"extended {len(keep)} configs with seeds {list(seeds)}: {added} new trials")
        for k in keep:
            print(f"  {k}  {tags[k]}")
    return plan


# ------------------------------------------------------------------ trial ---

def _tune_callback_class():
    import keras

    class TuneCallback(keras.callbacks.Callback):
        """Monitor every rotation after every epoch, replay the early-stopping
        grid online, and score the restored epochs as their answers come in."""

        def __init__(self, grid: ESGrid, val_x, val_y, parts, pos_rows, pos_all_rows,
                     all_x, pool, pool_fair=None, keep_preds: bool = False, all_batch: int = 4096):
            super().__init__()
            self.grid = grid
            self.val_x = np.asarray(val_x, dtype="float32")
            self.val_y = np.asarray(val_y, dtype="float64")
            self.parts = np.asarray(parts)
            self.rotations = tuple(int(r) for r in np.unique(self.parts))
            self.monitor_mask = {r: self.parts != r for r in self.rotations}
            self.score_mask = {r: self.parts == r for r in self.rotations}
            self.pos_rows = np.asarray(pos_rows)
            self.pos_all_rows = np.asarray(pos_all_rows)
            self.all_x = all_x
            self.pool = np.asarray(pool, dtype=bool)
            self.pool_fair = None if pool_fair is None else np.asarray(pool_fair, dtype=bool)
            self.keep_preds = keep_preds
            self.all_batch = int(all_batch)
            self.tracker = ESGridTracker(grid, self.rotations)
            self.curves = {(m, r): [] for m in grid.monitors for r in self.rotations}
            self.loss, self.val_full = [], []
            self.val_preds, self.all_preds = [], []
            self.snap: dict[int, list] = {}
            self.evals: dict[int, dict] = {}
            self.finals: dict[tuple, tuple[int, int]] = {}
            self.eval_seconds = 0.0

        def on_train_begin(self, logs=None):
            # the single input tensor, not the inputs list: fed a bare array, a
            # list-input model warns about the structure on every process
            self._g = keras.Model(self.model.input, self.model.get_layer("global").output)

        def _predict_val(self):
            return np.asarray(self._g.predict_on_batch(self.val_x), dtype="float64").reshape(-1)

        def _predict_all(self):
            return np.asarray(self._g.predict(self.all_x, batch_size=self.all_batch, verbose=0),
                              dtype="float64").reshape(-1)

        def _evaluate(self, epoch: int, now: int, p_now: np.ndarray) -> dict:
            t0 = time.time()
            if epoch == now:
                pv, pa = p_now, self._predict_all()
            else:
                cur = self.model.get_weights()
                self.model.set_weights(self.snap[epoch])
                pv, pa = self._predict_val(), self._predict_all()
                self.model.set_weights(cur)
            pos_scores = pa[self.pos_all_rows]
            ev = {
                "rank": rank_among(np.sort(pa[self.pool]), pos_scores).tolist(),
                "score": pos_scores.tolist(),
                # the same windows scored as val rows: a cheap proof the all-row
                # mapping is right, kept so a mismatch would be visible later
                "val_all_diff": float(np.max(np.abs(pv[self.pos_rows] - pos_scores))),
                "ap": [], "roc": [], "pr_auc": [],
                "val_pr_auc_full": keras_pr_auc(self.val_y, pv),
            }
            if self.pool_fair is not None:
                ev["rank_fair"] = rank_among(np.sort(pa[self.pool_fair]), pos_scores).tolist()
            for r in self.rotations:
                m = self.score_mask[r]
                ev["ap"].append(average_precision(self.val_y[m], pv[m]))
                ev["roc"].append(roc_auc(self.val_y[m], pv[m]))
                ev["pr_auc"].append(keras_pr_auc(self.val_y[m], pv[m]))
            if self.keep_preds:
                ev["val_pred"] = pv.astype("float32").tolist()
            self.eval_seconds += time.time() - t0
            return ev

        def on_epoch_end(self, epoch, logs=None):
            p = self._predict_val()
            if self.keep_preds:
                self.val_preds.append(p.astype("float32"))
                self.all_preds.append(self._predict_all().astype("float32"))
            vals = {}
            for m in self.grid.monitors:
                f = MONITORS[m]
                for r in self.rotations:
                    mk = self.monitor_mask[r]
                    v = float(f(self.val_y[mk], p[mk]))
                    vals[m, r] = v
                    self.curves[m, r].append(v)
            self.val_full.append(keras_pr_auc(self.val_y, p))
            self.loss.append(float((logs or {}).get("loss", float("nan"))))

            snapshot, finals = self.tracker.step(epoch, vals)
            if snapshot:
                self.snap[epoch] = self.model.get_weights()
            for e in sorted({rest for _, _, rest in finals} - set(self.evals)):
                self.evals[e] = self._evaluate(e, epoch, p)
            for key, last, rest in finals:
                self.finals[key] = (last, rest)
            live = self.tracker.live_epochs()
            for e in [e for e in self.snap if e not in live]:
                del self.snap[e]
            if self.tracker.done:
                self.model.stop_training = True

        def result(self) -> dict:
            n_expected = (len(self.grid.monitors) * len(self.rotations) * len(self.grid.starts)
                          * len(self.grid.patiences) * len(self.grid.caps))
            if len(self.finals) != n_expected:
                raise RuntimeError(f"training ended with {len(self.finals)} of {n_expected} "
                                   "early-stopping settings unresolved")
            keys = sorted(self.finals)
            return {
                "evals": {str(e): v for e, v in sorted(self.evals.items())},
                "settings": {
                    "monitor": [k[0] for k in keys], "rotation": [k[1] for k in keys],
                    "start": [k[2] for k in keys], "patience": [k[3] for k in keys],
                    "cap": [k[4] for k in keys],
                    "last": [self.finals[k][0] for k in keys],
                    "restored": [self.finals[k][1] for k in keys],
                },
                "epochs_run": len(self.loss),
                "eval_seconds": round(self.eval_seconds, 2),
            }

    return TuneCallback


def _apply_threads() -> None:
    """Honour LF_TUNE_THREADS (set by `tune run --threads`) before keras builds
    anything; the TF_NUM_*_THREADS variables the runner also sets cover the
    case where the runtime was already initialised by the import."""
    n = os.environ.get("LF_TUNE_THREADS")
    if not n:
        return
    try:
        import tensorflow as tf

        tf.config.threading.set_intra_op_parallelism_threads(int(n))
        tf.config.threading.set_inter_op_parallelism_threads(2)
    except Exception:
        pass


def run_trial(root, config_id: str, fold: str, seed: int, verbose: int = 0,
              keep_preds: bool = False) -> Path:
    """Train and score ONE (config, fold, seed). Run it in its own process
    (`tune run` does): a second model trained in one keras process can die."""
    _apply_threads()
    from .data import TermArrays
    from .pipeline import set_seed, train

    root = Path(root).expanduser()
    plan = load_plan(root)
    entry = _config_entry(plan, config_id)
    finfo = plan["folds"][fold]
    in_dir = Path(finfo["dir"])
    fp = finfo["fingerprint"]
    if _sha1(in_dir / "config.json") != fp["config_sha1"] or \
            (in_dir / "arrays.npz").stat().st_size != fp["arrays_bytes"]:
        raise SystemExit(f"{in_dir} changed since the plan was made; its parts and "
                         "results would no longer match -- make a new plan")
    grid = ESGrid.from_dict(plan["grid"])
    base = Config.from_json(in_dir / "config.json")
    term = base.term_order[0]
    cfg = base.evolve(**entry["overrides"], epochs=grid.max_epochs, seed=int(seed))

    with np.load(in_dir / "arrays.npz") as z:
        tr = {k: z[f"{term}__train__{k}"] for k in ("x", "y_global", "y_per_index_cat")}
        if f"{term}__train__y_ins" in z.files:
            tr["y_ins"] = z[f"{term}__train__y_ins"]
        val_x = z[f"{term}__val__x"]
        val_y = z[f"{term}__val__y_global"][:, 0]
        all_x = z[f"{term}__all__x"]
    train_arr = TermArrays.from_mapping(tr).validate(cfg, where=f"{term}/train: ")
    with np.load(root / FOLDS_DIR / f"{fold}.npz") as fz:
        f = {k: fz[k] for k in fz.files}

    out = trial_dir(root, config_id, fold, seed)
    out.mkdir(parents=True, exist_ok=True)
    cb = _tune_callback_class()(grid, val_x, val_y, f["parts"], f["pos_rows"],
                                f["pos_all_rows"], all_x, f["pool"], f.get("pool_fair"),
                                keep_preds=keep_preds)
    # seeded as run_member seeds member k = seed - cfg.seed of a single-terminus
    # ensemble, so the baseline's seed s retraces production member s - 1
    set_seed(int(seed))
    rng = np.random.default_rng(int(seed))
    t0 = time.time()
    model, _ = train(term, {"train": train_arr}, cfg, verbose=verbose, extra_callbacks=[cb],
                     rng=rng, early_stopping=False, validate=False)
    res = cb.result()
    res.update({
        "version": 1, "config": config_id, "fold": fold, "seed": int(seed),
        "params": entry["params"], "overrides": entry["overrides"],
        "n_params": int(model.count_params()), "seconds": round(time.time() - t0, 1),
        "grid": plan["grid"], "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        "resolved_config": cfg.to_dict(),
    })
    curves = {"loss": np.asarray(cb.loss, "float32"), "val_pr_auc_full": np.asarray(cb.val_full, "float32")}
    for (m, r), v in cb.curves.items():
        curves[f"{m}__r{r}"] = np.asarray(v, "float32")
    if keep_preds:
        curves["val_preds"] = np.stack(cb.val_preds)
        curves["all_preds"] = np.stack(cb.all_preds)
    np.savez_compressed(out / CURVES_FILE, **curves)
    _write_json(out / RESULT_FILE, res)        # last: its presence marks the trial done
    if verbose:
        print(f"wrote {out / RESULT_FILE} ({res['epochs_run']} epochs, {res['seconds']}s, "
              f"{len(res['evals'])} restored epochs scored in {res['eval_seconds']}s)")
    return out / RESULT_FILE


# ----------------------------------------------------------------- runner ---

def run_plan(root, workers: int = 4, threads: int | None = None, retries: int = 3,
             python: str | None = None, configs: Sequence[str] | None = None,
             dry_run: bool = False) -> int:
    """Run every unfinished trial, ``workers`` processes at a time, one model
    per process. Re-running resumes."""
    root = Path(root).expanduser()
    plan = load_plan(root)
    trials = [t for t in plan_trials(plan) if configs is None or t[0] in configs]
    todo = [t for t in trials if not (trial_dir(root, *t) / RESULT_FILE).exists()]
    tags = {c["id"]: c["tag"] for c in plan["configs"]}
    print(f"{len(trials) - len(todo)}/{len(trials)} trials done; {len(todo)} to run, "
          f"{workers} at a time" + (f", {threads} threads each" if threads else ""), flush=True)
    if dry_run or not todo:
        for t in todo[:20]:
            print(f"  would run {t[0]} ({tags[t[0]]}) {t[1]} seed {t[2]}")
        return 0

    env = dict(os.environ)
    ## the children must import THIS tree, not whatever an installed
    ## ligandFinder would put first (lf-dcnn-path-shadowing)
    here = str(Path(__file__).resolve().parents[1])
    env["PYTHONPATH"] = here + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    if threads:
        env.update(LF_TUNE_THREADS=str(threads), TF_NUM_INTRAOP_THREADS=str(threads),
                   TF_NUM_INTEROP_THREADS="2", OMP_NUM_THREADS=str(threads))
    py = python or sys.executable

    def one(t):
        cid, fold, seed = t
        d = trial_dir(root, cid, fold, seed)
        d.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        for attempt in range(1, retries + 1):
            with open(d / "log.txt", "a") as log:
                log.write(f"### attempt {attempt} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                log.flush()
                rc = subprocess.call([py, "-u", "-m", "lf_dcnn", "tune", "trial", "--root", str(root),
                                      "--config", cid, "--fold", fold, "--seed", str(seed)],
                                     stdout=log, stderr=subprocess.STDOUT, env=env)
            if rc == 0 and (d / RESULT_FILE).exists():
                return t, True, attempt, time.time() - t0
        return t, False, retries, time.time() - t0

    t_start = time.time()
    n_ok = n_fail = 0
    durations = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(one, t) for t in todo]
        for k, fut in enumerate(as_completed(futs), 1):
            (cid, fold, seed), ok, att, sec = fut.result()
            n_ok += ok
            n_fail += not ok
            durations.append(sec)
            # per-trial time over the workers, not elapsed / done: the first
            # finishers ran in parallel with the rest of their batch
            eta = float(np.mean(durations)) * (len(todo) - k) / min(workers, len(todo))
            print(f"[{k}/{len(todo)}] {cid} ({tags[cid]}) {fold} s{seed} "
                  f"{'ok' if ok else 'FAILED'} in {_fmt_time(sec)}"
                  + (f" (attempt {att})" if att > 1 else "") + f"; ETA {_fmt_time(eta)}",
                  flush=True)
    print(f"done: {n_ok} ok, {n_fail} failed in {_fmt_time(time.time() - t_start)}"
          + ("; re-run to retry the failures (see their log.txt)" if n_fail else ""))
    return 1 if n_fail else 0


# ----------------------------------------------------------------- report ---

def _load_results(root: Path, plan: Mapping):
    """Long tables from every finished trial: one row per (trial, setting,
    held-out end) and one per (trial, setting)."""
    import pandas as pd

    pos_rows, trial_rows = [], []
    fold_cache = {}
    for (cid, fold, seed) in plan_trials(plan):
        p = trial_dir(root, cid, fold, seed) / RESULT_FILE
        if not p.exists():
            continue
        with open(p) as fh:
            res = json.load(fh)
        if fold not in fold_cache:
            with np.load(root / FOLDS_DIR / f"{fold}.npz") as fz:
                fold_cache[fold] = fz["parts"][fz["pos_rows"]]
        pos_part = fold_cache[fold]
        S = res["settings"]
        restored = {}
        last = {}
        for i in range(len(S["monitor"])):
            k = (S["monitor"][i], S["rotation"][i], S["start"][i], S["patience"][i], S["cap"][i])
            restored[k] = S["restored"][i]
            last[k] = S["last"][i]
        evals = res["evals"]
        rotations = sorted(set(S["rotation"]))
        settings = sorted({(k[0], k[2], k[3], k[4]) for k in restored})
        n_pos_rot = {r: int((pos_part == r).sum()) for r in rotations}
        for (m, s, pat, cap) in settings:
            ap = roc = prk = 0.0
            rest_eps, capped = [], 0
            for r in rotations:
                k = (m, r, s, pat, cap)
                ev = evals[str(restored[k])]
                w = n_pos_rot[r] / pos_part.size
                ap += w * ev["ap"][r]
                roc += w * ev["roc"][r]
                prk += w * ev["pr_auc"][r]
                rest_eps.append(restored[k])
                capped += last[k] == cap - 1
            trial_rows.append({
                "config": cid, "fold": fold, "seed": seed, "monitor": m, "start": s,
                "patience": pat, "cap": cap, "ap": ap, "roc": roc, "pr_auc": prk,
                "restored_median": float(np.median(rest_eps)),
                "capped_frac": capped / len(rotations), "epochs_run": res["epochs_run"],
                "n_params": res["n_params"], "seconds": res["seconds"],
            })
            for i, r in enumerate(pos_part):
                ev = evals[str(restored[(m, int(r), s, pat, cap)])]
                pos_rows.append((cid, fold, seed, m, s, pat, cap, i, ev["rank"][i],
                                 ev["rank_fair"][i] if "rank_fair" in ev else np.nan))
    pos = pd.DataFrame(pos_rows, columns=["config", "fold", "seed", "monitor", "start", "patience",
                                          "cap", "pos", "rank", "rank_fair"])
    return pos, pd.DataFrame(trial_rows)


SETTING = ["monitor", "start", "patience", "cap"]


def summarise(root, objective: str | None = None) -> dict:
    """Aggregate the finished trials. Returns the tables `report` writes.

    ``objective``: ``log_rank_fair`` (default when the plan had other-end
    files) or ``log_rank`` -- the mean log10 candidate rank of the held-out
    ends, each end first averaged over seeds.
    """
    import pandas as pd

    root = Path(root).expanduser()
    plan = load_plan(root)
    pos, trials = _load_results(root, plan)
    if pos.empty:
        raise SystemExit("no finished trials yet")
    if objective is None:
        objective = "log_rank_fair" if pos["rank_fair"].notna().all() else "log_rank"
    col = {"log_rank": "rank", "log_rank_fair": "rank_fair"}[objective]
    pos["lr"] = np.log10(pos[col].astype(float))
    for k in (100, 500, 1000):
        pos[f"top{k}"] = (pos[col] <= k).astype(float)

    n_folds = len(plan["folds"])
    key = ["config"] + SETTING
    # one value per held-out end: averaged over the seeds that config ran
    per_end = (pos.groupby(key + ["fold", "pos"], as_index=False)
               .agg(lr=("lr", "mean"), top100=("top100", "mean"), top500=("top500", "mean"),
                    top1000=("top1000", "mean"), n_seeds=("seed", "nunique")))
    agg = (per_end.groupby(key, as_index=False)
           .agg(objective=("lr", "mean"), lr_sd=("lr", "std"), n_ends=("lr", "size"),
                n_folds=("fold", "nunique"), n_seeds=("n_seeds", "min"),
                top100=("top100", "mean"), top500=("top500", "mean"), top1000=("top1000", "mean")))
    agg["gm_rank"] = 10 ** agg["objective"]
    med = (per_end.assign(r=10 ** per_end["lr"]).groupby(key)["r"].median()
           .rename("median_rank").reset_index())
    tagg = (trials.groupby(key, as_index=False)
            .agg(ap=("ap", "mean"), roc=("roc", "mean"), pr_auc=("pr_auc", "mean"),
                 restored_median=("restored_median", "median"),
                 capped_frac=("capped_frac", "mean"), epochs_run=("epochs_run", "median"),
                 n_params=("n_params", "first"), minutes=("seconds", lambda s: s.mean() / 60)))
    table = agg.merge(med, on=key).merge(tagg, on=key)
    table = table[table["n_folds"] == n_folds].copy()
    if table.empty:
        n_done = trials[["config", "fold", "seed"]].drop_duplicates().shape[0]
        raise SystemExit(f"no config has finished all {n_folds} folds yet "
                         f"({n_done} of {len(plan_trials(plan))} trials done)")
    tags = {c["id"]: c["tag"] for c in plan["configs"]}
    params = {c["id"]: c["params"] for c in plan["configs"]}
    table.insert(1, "tag", table["config"].map(tags))

    # ---- paired against the baseline: same held-out end, seed-averaged ----
    bl = plan["baseline"]
    bset = bl["setting"]
    bsel = (per_end["config"] == bl["config"])
    for k in SETTING:
        bsel &= per_end[k] == bset[k]
    base_end = per_end[bsel].set_index(["fold", "pos"])["lr"]
    if base_end.empty:
        raise SystemExit("the baseline has not finished on any fold yet")

    from scipy.stats import wilcoxon

    def paired(df):
        d = df.set_index(["fold", "pos"])["lr"] - base_end
        d = d.dropna()
        if d.empty:
            return pd.Series({"delta": np.nan, "delta_se": np.nan, "p_wilcoxon": np.nan,
                              "folds_better": np.nan})
        by_fold = d.groupby(level="fold").mean()
        p = np.nan
        if (d != 0).any():
            p = float(wilcoxon(d.to_numpy(), zero_method="wilcox").pvalue)
        return pd.Series({"delta": d.mean(), "delta_se": d.std(ddof=1) / math.sqrt(d.size),
                          "p_wilcoxon": p, "folds_better": int((by_fold < 0).sum())})

    pr = per_end.groupby(key).apply(paired, include_groups=False).reset_index()
    table = table.merge(pr, on=key, how="left")

    # A cap that no rotation of any trial reached gives exactly the answer of
    # the cap below it. Keep the smallest of such a run (the cheaper budget,
    # same model) and say how far it extends, so the table lists distinct
    # settings; the baseline's own row always stays.
    sig = (per_end.sort_values(["fold", "pos"]).groupby(key)["lr"]
           .apply(lambda s: hash(tuple(np.round(s.to_numpy(), 12)))).rename("_sig").reset_index())
    table = table.merge(sig, on=key).sort_values(key).reset_index(drop=True)
    ## the marginals select rows BY setting value (cap == 2000 ...), so they
    ## read the uncollapsed table; the collapsed one is for listing
    table_all = table.drop(columns="_sig").copy()
    grp = ["config", "monitor", "start", "patience"]
    run_id = (table["_sig"] != table.groupby(grp)["_sig"].shift(1)).cumsum()
    table["cap_upto"] = table.groupby(run_id)["cap"].transform("max")
    is_base = table["config"] == bl["config"]
    for k in SETTING:
        is_base &= table[k] == bset[k]
    first = run_id != run_id.shift(1)
    table = table[first | is_base].drop(columns="_sig")
    table = table.sort_values(["objective", "cap", "patience"]).reset_index(drop=True)

    # ---- leave-one-fold-out: the winner's curse, measured ----
    full = table[["config"] + SETTING]
    pe = per_end.merge(full, on=key)
    lofo_rows = []
    for f in sorted(plan["folds"]):
        other = pe[pe["fold"] != f].groupby(key)["lr"].mean()
        pick = other.idxmin()
        here = pe[pe["fold"] == f].set_index(key + ["pos"])["lr"]
        v_pick = here.xs(pick, level=key)
        v_base = base_end.xs(f, level="fold")
        best = tuple(table.loc[0, key])
        v_best = here.xs(best, level=key)
        lofo_rows.append({"fold": f, "pick_config": pick[0], "pick_tag": tags[pick[0]],
                          **dict(zip(SETTING, pick[1:])), "n_ends": int(v_pick.size),
                          "pick_lr": v_pick.mean(), "global_best_lr": v_best.mean(),
                          "baseline_lr": v_base.mean(),
                          "delta_vs_baseline": (v_pick - v_base).mean()})
    lofo = pd.DataFrame(lofo_rows)

    # ---- marginals ----
    marg = []
    for dim in SETTING:
        others = {k: bset[k] for k in SETTING if k != dim}
        sel = table_all.copy()
        for k, v in others.items():
            sel = sel[sel[k] == v]
        for lvl, g in sel.groupby(dim):
            marg.append({"dimension": dim, "level": lvl, "view": "other early-stopping at production",
                         "objective_mean": g["objective"].mean(), "gm_rank": 10 ** g["objective"].mean(),
                         "best_objective": g["objective"].min(), "n_configs": g["config"].nunique()})
    at_prod = table_all.copy()
    for k in SETTING:
        at_prod = at_prod[at_prod[k] == bset[k]]
    own_best = table.sort_values("objective").drop_duplicates("config")
    for view, src in (("production early-stopping", at_prod), ("each config's best early-stopping", own_best)):
        src = src.assign(**{k: src["config"].map(lambda c, k=k: params[c].get(k)) for k in plan["space"]})
        for dim in plan["space"]:
            vals = src[dim].map(lambda v: str(v) if isinstance(v, list) else v)
            if dim == "learning_rate":
                vals = pd.cut(src[dim], LR_BIN_EDGES, right=False).astype(str)
            for lvl, g in src.groupby(vals, sort=True):
                marg.append({"dimension": dim, "level": lvl, "view": view,
                             "objective_mean": g["objective"].mean(),
                             "gm_rank": 10 ** g["objective"].mean(),
                             "best_objective": g["objective"].min(), "n_configs": g["config"].nunique()})
    marginals = pd.DataFrame(marg)
    return {"plan": plan, "objective": objective, "table": table, "lofo": lofo,
            "marginals": marginals, "per_end": per_end, "trials": trials,
            "n_trials_done": int(trials[["config", "fold", "seed"]].drop_duplicates().shape[0]),
            "n_trials": len(plan_trials(plan))}


def _flags_10_8(row, plan) -> str:
    ov = _config_entry(plan, row["config"])["overrides"]
    f = [f"--learning-rate {ov['learning_rate']:g}"] if "learning_rate" in ov else []
    if "unet_filters" in ov:
        f.append("--unet-filters " + ",".join(str(int(v)) for v in ov["unet_filters"]))
    if "conv_kernel" in ov:
        f.append(f"--conv-kernel {int(ov['conv_kernel'])}")
    f += [f"--patience {int(row['patience'])}", f"--start-from-epoch {int(row['start'])}",
          f"--epochs {int(row['cap'])}"]
    if row["monitor"] != "pr_auc":
        f.append("--monitor val_global_auc")
    return " ".join(f)


def report(root, objective: str | None = None, top: int = 15, verbose: int = 1) -> dict:
    """Write ``<root>/report/``: the ranked table, the leave-one-fold-out
    check, the marginals, ``best.json`` and a ``summary.md``."""
    root = Path(root).expanduser()
    s = summarise(root, objective=objective)
    plan, table, lofo, marg = s["plan"], s["table"], s["lofo"], s["marginals"]
    out = root / REPORT_DIR
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "configs.csv", index=False)
    lofo.to_csv(out / "lofo.csv", index=False)
    marg.to_csv(out / "marginals.csv", index=False)
    s["trials"].to_csv(out / "trials.csv", index=False)

    bl = plan["baseline"]
    bsel = table["config"] == bl["config"]
    for k in SETTING:
        bsel &= table[k] == bl["setting"][k]
    if not bsel.any():
        raise SystemExit("the baseline has not finished every fold yet; nothing to compare against")
    brow = table[bsel].iloc[0]
    best = table.iloc[0]
    best_json = {
        "config": best["config"], "tag": best["tag"],
        "params": _config_entry(plan, best["config"])["params"],
        "overrides": _config_entry(plan, best["config"])["overrides"],
        "early_stopping": {k: (best[k].item() if hasattr(best[k], "item") else best[k]) for k in SETTING},
        "objective": s["objective"], "value": float(best["objective"]),
        "gm_rank": float(best["gm_rank"]), "baseline_gm_rank": float(brow["gm_rank"]),
        "delta_vs_baseline": float(best["delta"]), "p_wilcoxon": float(best["p_wilcoxon"]),
        "lofo_delta_vs_baseline": float(np.average(lofo["delta_vs_baseline"], weights=lofo["n_ends"])),
        "flags_10_8": _flags_10_8(best, plan),
        "trials_done": s["n_trials_done"], "trials_planned": s["n_trials"],
    }
    _write_json(out / "best.json", best_json)

    def fmt_row(r):
        p = "" if pd_isna(r["p_wilcoxon"]) else f"{r['p_wilcoxon']:.3g}"
        fb = "" if pd_isna(r["folds_better"]) else int(r["folds_better"])
        cap = (f"{int(r['cap'])}" if int(r["cap_upto"]) == int(r["cap"])
               else f"{int(r['cap'])}-{int(r['cap_upto'])}")
        return (f"| {r['tag']} | {r['monitor']} | {int(r['start'])} | {int(r['patience'])} | "
                f"{cap} | {r['gm_rank']:.0f} | {r['median_rank']:.0f} | "
                f"{r['top500']:.2f} | {r['ap']:.3f} | {r['delta']:+.3f} | {p} | "
                f"{fb}/{len(plan['folds'])} | {r['restored_median']:.0f} | "
                f"{r['capped_frac']:.2f} | {int(r['n_seeds'])} |")

    col = "rank_fair" if s["objective"] == "log_rank_fair" else "rank"
    lines = [
        f"# Hyperparameter search: {plan['term']} terminus",
        "",
        f"`{root}` -- {s['n_trials_done']} of {s['n_trials']} trials finished; "
        f"{len(plan['folds'])} folds, {int(table['n_ends'].iloc[0])} held-out ends; "
        f"val split in {plan['n_parts']} parts (early stopping on {plan['n_parts'] - 1}, "
        f"scored on the held-out one).",
        "",
        f"**Objective**: mean log10 `{col}` of each held-out end among the fold's unknown "
        "candidate windows (lower is better), shown as the geometric-mean rank (`gm rank`). "
        + ("`rank_fair` drops candidates sitting on real unlabelled peptide ends from the pool. "
           if col == "rank_fair" else "")
        + "`delta` is the paired difference in log10 rank against the baseline "
        "(production config and early stopping), averaged over the same ends; negative is better.",
        "",
        f"**Baseline** (production: {brow['tag']}, pr_auc / start {bl['setting']['start']} / "
        f"patience {bl['setting']['patience']} / cap {bl['setting']['cap']}): gm rank "
        f"{brow['gm_rank']:.0f}, median {brow['median_rank']:.0f}, top-500 {brow['top500']:.2f}.",
        "",
        f"**Best**: {best['tag']}, {best['monitor']} / start {int(best['start'])} / patience "
        f"{int(best['patience'])} / cap {int(best['cap'])}: gm rank {best['gm_rank']:.0f} "
        f"(delta {best['delta']:+.3f} log10, p = {best['p_wilcoxon']:.3g}, better on "
        f"{int(best['folds_better'])}/{len(plan['folds'])} folds).",
        "",
        "```",
        f"Rscript inst/scripts/10_8_pepend_train.R --term {plan['term']} ... {best_json['flags_10_8']}",
        "```",
        "",
        f"## Top {top}",
        "",
        "| config | monitor | start | patience | cap | gm rank | median | top-500 | AP (hard neg) "
        "| delta | p | folds better | restored ep | capped | seeds |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ] + [fmt_row(r) for _, r in table.head(top).iterrows()] + [
        "",
        "`cap` a-b: every cap from a to b gives this same answer (patience stopped first), "
        "so the smaller budget is listed. `capped`: share of rotations that ran out of "
        "epochs (stopped at the cap, not by patience) -- near 1 means a longer budget "
        "might still help.",
        "",
        "## Winner's curse: leave one fold out",
        "",
        "For each fold, the setting that is best on the OTHER folds, scored on this one. "
        "Its delta against the baseline is the honest estimate of what tuning buys; the "
        "best row's own delta above is selected on all folds and is optimistic.",
        "",
        "| fold | picked | monitor | start | patience | cap | picked log10 rank | global best | baseline | delta |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ] + [
        f"| {r['fold']} | {r['pick_tag']} | {r['monitor']} | {int(r['start'])} | {int(r['patience'])} "
        f"| {int(r['cap'])} | {r['pick_lr']:.3f} | {r['global_best_lr']:.3f} | {r['baseline_lr']:.3f} "
        f"| {r['delta_vs_baseline']:+.3f} |"
        for _, r in lofo.iterrows()
    ] + [
        "",
        f"Pooled leave-one-fold-out delta: {best_json['lofo_delta_vs_baseline']:+.3f} log10 "
        f"(x{10 ** best_json['lofo_delta_vs_baseline']:.2f} on the rank).",
        "",
        "## Marginals",
        "",
        "| dimension | level | view | gm rank (mean) | best | configs |",
        "|---|---|---|---|---|---|",
    ] + [
        f"| {r['dimension']} | {r['level']} | {r['view']} | {r['gm_rank']:.0f} | "
        f"{10 ** r['best_objective']:.0f} | {int(r['n_configs'])} |"
        for _, r in marg.iterrows()
    ] + [
        "",
        "## Read with care",
        "",
        f"* Early stopping here watches {plan['n_parts'] - 1}/{plan['n_parts']} of val; "
        "production watches all of it, so production's stopping is a little less noisy.",
        "* Every fold's ends took part in the selection. A trunk retrained with the winner "
        "and scored on these folds is optimistic -- use the leave-one-fold-out delta, and "
        "the never-seen ends, which this search never touched, as the clean test.",
        "* AP is against the held-out part's own hard negatives (near-end, motif and random "
        "windows): it rewards boundary precision more than genome-wide ranking. Rank decides.",
    ]
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    if verbose:
        print("\n".join(lines))
        print(f"\nreport -> {out}")
    return s


def pd_isna(v) -> bool:
    try:
        return v is None or (isinstance(v, float) and math.isnan(v))
    except TypeError:
        return False


# -------------------------------------------------------------------- CLI ---

def _csv_ints(s: str) -> tuple[int, ...]:
    return tuple(int(v) for v in s.split(",") if v.strip())


def add_parser(sub) -> None:
    """Register ``python -m lf_dcnn tune ...``."""
    t = sub.add_parser("tune", help="hyperparameter search (see lf_dcnn/tune.py)",
                       description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ts = t.add_subparsers(dest="tune_command", required=True)

    p = ts.add_parser("plan", help="write <root>/plan.json")
    p.add_argument("--root", required=True)
    p.add_argument("--inputs", required=True,
                   help="directory of fold inputs (one subdirectory per fold, as "
                        "10_8_pepend_train.R --export-only writes them)")
    p.add_argument("--n-configs", type=int, default=24)
    p.add_argument("--seeds", type=_csv_ints, default=(1,))
    p.add_argument("--space", help="JSON file with the search space (default: DEFAULT_SPACE)")
    p.add_argument("--sampler", choices=("lhs", "grid"), default="lhs")
    p.add_argument("--sampler-seed", type=int, default=1)
    p.add_argument("--monitors", default=",".join(ESGrid().monitors))
    p.add_argument("--starts", type=_csv_ints, default=ESGrid().starts)
    p.add_argument("--patiences", type=_csv_ints, default=ESGrid().patiences)
    p.add_argument("--caps", type=_csv_ints, default=ESGrid().caps)
    p.add_argument("--parts", type=int, default=4, help="val parts (early-stop on n-1, score 1)")
    p.add_argument("--part-seed", type=int, default=20261102)
    p.add_argument("--other-ends", help="directory with ligand_list.csv + uniprot_peptides.csv; "
                                        "enables rank_fair (~/AF2_analysis/pepend_locus)")
    p.add_argument("--force", action="store_true")

    r = ts.add_parser("run", help="run every unfinished trial")
    r.add_argument("--root", required=True)
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--threads", type=int, default=None, help="TF threads per trial process")
    r.add_argument("--retries", type=int, default=3)
    r.add_argument("--configs", default=None, help="comma-separated config ids (default all)")
    r.add_argument("--dry-run", action="store_true")

    tr = ts.add_parser("trial", help="one (config, fold, seed) in this process (used by `run`)")
    tr.add_argument("--root", required=True)
    tr.add_argument("--config", required=True)
    tr.add_argument("--fold", required=True)
    tr.add_argument("--seed", type=int, required=True)
    tr.add_argument("--verbose", type=int, default=0)
    tr.add_argument("--keep-preds", action="store_true", help="debug: log every epoch's predictions")

    rp = ts.add_parser("report", help="aggregate the finished trials into <root>/report/")
    rp.add_argument("--root", required=True)
    rp.add_argument("--objective", choices=("log_rank_fair", "log_rank"), default=None)
    rp.add_argument("--top", type=int, default=15)

    e = ts.add_parser("extend", help="add seeds to the baseline and the top configs")
    e.add_argument("--root", required=True)
    e.add_argument("--top", type=int, default=5)
    e.add_argument("--seeds", type=_csv_ints, default=(2, 3))
    e.add_argument("--objective", choices=("log_rank_fair", "log_rank"), default=None)

    t.set_defaults(func=main)


def main(args) -> int:
    c = args.tune_command
    if c == "plan":
        space = None
        if args.space:
            with open(args.space) as fh:
                space = json.load(fh)
        grid = ESGrid(tuple(m.strip() for m in args.monitors.split(",") if m.strip()),
                      args.starts, args.patiences, args.caps)
        make_plan(args.root, args.inputs, n_configs=args.n_configs, seeds=args.seeds,
                  space=space, sampler=args.sampler, sampler_seed=args.sampler_seed, grid=grid,
                  n_parts=args.parts, part_seed=args.part_seed,
                  other_ends=Path(args.other_ends).expanduser() if args.other_ends else None,
                  force=args.force)
        return 0
    if c == "run":
        return run_plan(args.root, workers=args.workers, threads=args.threads, retries=args.retries,
                        configs=args.configs.split(",") if args.configs else None,
                        dry_run=args.dry_run)
    if c == "trial":
        run_trial(args.root, args.config, args.fold, args.seed, verbose=args.verbose,
                  keep_preds=args.keep_preds)
        return 0
    if c == "report":
        report(args.root, objective=args.objective, top=args.top)
        return 0
    if c == "extend":
        extend_plan(args.root, top=args.top, seeds=args.seeds, objective=args.objective)
        return 0
    raise SystemExit(f"unknown tune command {c!r}")
