"""On-disk exchange format for the standalone (no-R) entry point.

An input directory holds::

    arrays.npz     the model arrays, keyed "<term>__<split>__<name>"
    config.json    the Config used to build them
    meta.parquet   (optional) one row per window of the concatenated `all`
                   splits, in term order -- carried through to the predictions

and a run writes back::

    outputs.npz        pred, pred_raw, pred_sd, per_index, per_index_sd, emb,
                       val_scores, val_labels, val_scores_members, pred_raw_members
    predictions.parquet  pred / pred_raw / term, plus any meta columns
    history.json       per-term training curves and the calibrator fit

An isolated run (one process per ensemble member) additionally leaves the
members behind, so a run can be resumed or re-combined without retraining::

    members/member_000_N.npz   "<term>__global|per_index|val|val_labels|emb"
    members/member_000_N.json  k, seed, terms, params, histories

one file pair per member AND terminus (each is one model, trained in its own
process); ``load_members`` merges a member's termini back together.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from .config import Config
from .data import ARRAY_NAMES

ARRAYS_FILE = "arrays.npz"
CONFIG_FILE = "config.json"
META_FILE = "meta.parquet"
OUTPUTS_FILE = "outputs.npz"
PREDICTIONS_FILE = "predictions.parquet"
HISTORY_FILE = "history.json"
MEMBERS_DIR = "members"

SEP = "__"


def _key(term: str, split: str, name: str) -> str:
    return f"{term}{SEP}{split}{SEP}{name}"


def save_inputs(
    path,
    data: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    cfg: Config,
    meta: "pd.DataFrame | None" = None,
    compress: bool = True,
) -> Path:
    """Write ``arrays.npz`` + ``config.json`` (+ ``meta.parquet``)."""
    out = Path(path).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    flat: dict[str, np.ndarray] = {}
    for term, splits in data.items():
        for split, arrays in splits.items():
            for name in ARRAY_NAMES:
                if name in arrays:
                    flat[_key(term, split, name)] = np.asarray(
                        arrays[name], dtype="float32"
                    )
    saver = np.savez_compressed if compress else np.savez
    saver(out / ARRAYS_FILE, **flat)
    cfg.to_json(out / CONFIG_FILE)
    if meta is not None:
        meta.to_parquet(out / META_FILE, index=False)
    return out


def load_inputs(path) -> tuple[dict, Config, "pd.DataFrame | None"]:
    """Read back what :func:`save_inputs` wrote."""
    src = Path(path).expanduser()
    cfg = Config.from_json(src / CONFIG_FILE)

    data: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    with np.load(src / ARRAYS_FILE) as z:
        for key in z.files:
            term, split, name = key.split(SEP)
            data.setdefault(term, {}).setdefault(split, {})[name] = z[key]

    meta_path = src / META_FILE
    meta = pd.read_parquet(meta_path) if meta_path.exists() else None
    return data, cfg, meta


def save_outputs(path, result, meta: "pd.DataFrame | None" = None) -> Path:
    """Write ``outputs.npz`` + ``predictions.parquet`` + ``history.json``."""
    out = Path(path).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        out / OUTPUTS_FILE,
        pred=result.pred,
        pred_raw=result.pred_raw,
        pred_sd=result.pred_sd,
        per_index=result.per_index,
        per_index_sd=result.per_index_sd,
        emb=result.emb,
        val_scores=result.val_scores,
        val_labels=result.val_labels,
        val_scores_members=result.val_scores_members,
        pred_raw_members=result.pred_raw_members,
    )

    term = np.concatenate(
        [np.repeat(t, result.n_by_term[t]) for t in result.term_order]
    )
    preds = pd.DataFrame(
        {"term": term, "pred": result.pred, "pred_raw": result.pred_raw}
    )
    if meta is not None:
        if len(meta) != len(preds):
            raise ValueError(
                f"meta has {len(meta)} rows but {len(preds)} windows were scored"
            )
        preds = pd.concat([meta.reset_index(drop=True), preds], axis=1)
    preds.to_parquet(out / PREDICTIONS_FILE, index=False)

    with open(out / HISTORY_FILE, "w") as fh:
        json.dump(
            {
                "term_order": list(result.term_order),
                "n_by_term": dict(result.n_by_term),
                "n_seeds": int(result.n_seeds),
                "params": dict(result.params) or
                          {t: int(m.count_params()) for t, m in result.models.items()},
                "pi_names": list(result.pi_names),
                "calibrator": result.calibrator.to_dict(),
                "histories": result.histories,
            },
            fh,
            indent=2,
        )
    return out


def load_outputs(path) -> dict:
    """Read back what :func:`save_outputs` wrote."""
    src = Path(path).expanduser()
    with np.load(src / OUTPUTS_FILE) as z:
        arrays = {k: z[k] for k in z.files}
    with open(src / HISTORY_FILE) as fh:
        arrays.update(json.load(fh))
    pred_path = src / PREDICTIONS_FILE
    if pred_path.exists():
        arrays["predictions"] = pd.read_parquet(pred_path)
    return arrays


# ---- ensemble members on disk ------------------------------------------------

def member_paths(path, k: int, term: str) -> tuple[Path, Path]:
    d = Path(path).expanduser() / MEMBERS_DIR
    stem = f"member_{int(k):03d}_{term}"
    return d / f"{stem}.npz", d / f"{stem}.json"


def save_member(path, member) -> list[Path]:
    """Write a :class:`lf_dcnn.pipeline.Member` under ``<path>/members/``, one
    npz + json pair per terminus it holds."""
    written = []
    for term in member.terms:
        npz, js = member_paths(path, member.k, term)
        npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            npz,
            **{
                _key(term, "global", "x"): np.asarray(member.global_all[term], dtype="float64"),
                _key(term, "per_index", "x"): np.asarray(member.per_index_all[term], dtype="float32"),
                _key(term, "val", "x"): np.asarray(member.val[term], dtype="float64"),
                _key(term, "val_labels", "x"): np.asarray(member.val_labels[term], dtype="float64"),
                _key(term, "emb", "x"): np.asarray(member.emb[term], dtype="float64"),
            },
        )
        with open(js, "w") as fh:
            json.dump({"k": int(member.k), "seed": member.seed, "terms": [term],
                       "params": {term: int(member.params[term])},
                       "histories": {term: member.histories[term]}}, fh)
        written.append(npz)
    return written


def load_members(path, term_order=None) -> list:
    """Read back what :func:`save_member` wrote under ``<path>/members/``,
    merging each member's per-terminus files into one Member.

    ``term_order`` (e.g. ``Config.term_order``) fixes the order of a merged
    member's ``terms``; otherwise it is the order the files are found in.
    Every member must end up with the same termini.
    """
    from .pipeline import Member

    d = Path(path).expanduser() / MEMBERS_DIR
    by_k: dict[int, Member] = {}
    for js in sorted(d.glob("member_*.json")):
        with open(js) as fh:
            meta = json.load(fh)
        npz = js.with_suffix(".npz")
        if not npz.exists():
            raise FileNotFoundError(f"{js} has no matching {npz.name}")
        k = int(meta["k"])
        m = by_k.setdefault(k, Member(k=k, seed=meta["seed"], terms=(), global_all={},
                                      per_index_all={}, val={}, val_labels={}, emb={},
                                      histories={}, params={}))
        with np.load(npz) as z:
            for t in meta["terms"]:
                if t in m.terms:
                    raise ValueError(f"member {k}: terminus {t!r} appears twice under {d}")
                m.terms = m.terms + (t,)
                m.global_all[t] = z[_key(t, "global", "x")]
                m.per_index_all[t] = z[_key(t, "per_index", "x")]
                m.val[t] = z[_key(t, "val", "x")]
                m.val_labels[t] = z[_key(t, "val_labels", "x")]
                m.emb[t] = z[_key(t, "emb", "x")]
                m.histories[t] = meta["histories"][t]
                m.params[t] = int(meta["params"][t])

    out = []
    for k in sorted(by_k):
        m = by_k[k]
        if term_order is not None:
            m.terms = tuple(t for t in term_order if t in m.terms) + \
                      tuple(t for t in m.terms if t not in term_order)
        out.append(m)
    have = {m.terms for m in out}
    if len(have) > 1:
        raise ValueError(f"members hold different termini: {sorted(have)}")
    return out
