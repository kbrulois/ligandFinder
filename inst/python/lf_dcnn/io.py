"""On-disk exchange format for the standalone (no-R) entry point.

An input directory holds::

    arrays.npz     the model arrays, keyed "<term>__<split>__<name>"
    config.json    the Config used to build them
    meta.parquet   (optional) one row per window of the concatenated `all`
                   splits, in term order -- carried through to the predictions

and a run writes back::

    outputs.npz        per_index, emb, pred, pred_raw, val_scores, val_labels
    predictions.parquet  pred / pred_raw / term, plus any meta columns
    history.json       per-term training curves and the calibrator fit
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
        per_index=result.per_index,
        emb=result.emb,
        val_scores=result.val_scores,
        val_labels=result.val_labels,
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
