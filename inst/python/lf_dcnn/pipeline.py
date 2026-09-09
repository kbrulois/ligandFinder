"""End-to-end training, scoring and calibration for both terminus models."""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Mapping, Sequence

import keras
import numpy as np

from .calibrate import PlattCalibrator
from .config import Config
from .data import OversampledWindows, TermArrays, as_term_data
from .model import build_model, compile_model, embed_model


def set_seed(seed: int | None) -> None:
    """Seed python, numpy and the Keras backend."""
    if seed is None:
        return
    keras.utils.set_random_seed(int(seed))
    random.seed(int(seed))
    np.random.seed(int(seed) % (2**32))


def l2norm(m: np.ndarray) -> np.ndarray:
    """Row-wise L2 normalisation, floored as in the R ``l2norm``."""
    m = np.asarray(m, dtype="float64")
    return m / np.sqrt(np.maximum((m**2).sum(axis=1), 1e-12))[:, None]


def train(
    term: str,
    splits: Mapping[str, TermArrays],
    cfg: Config | None = None,
    verbose: int = 1,
    extra_callbacks: Sequence[keras.callbacks.Callback] | None = None,
    rng: np.random.Generator | None = None,
) -> tuple[keras.Model, dict]:
    """Build, compile and fit one terminus model.

    Returns ``(model, history)`` where ``history`` is a plain dict of lists.
    """
    cfg = cfg or Config()
    # keras3::clear_session() equivalent: without it the two per-terminus models
    # share graph state.
    model = build_model(cfg, clear_session=True)
    compile_model(model, cfg, term)

    train_ds = OversampledWindows(splits["train"], cfg, rng=rng)

    callbacks = [
        keras.callbacks.EarlyStopping(
            monitor=cfg.monitor,
            mode="max",
            patience=cfg.patience,
            start_from_epoch=cfg.start_from_epoch,
            restore_best_weights=True,
        )
    ]
    if cfg.tensorboard_dir:
        import os

        log_dir = os.path.join(os.path.expanduser(cfg.tensorboard_dir), term)
        callbacks.append(
            keras.callbacks.TensorBoard(
                log_dir=log_dir, update_freq="epoch", histogram_freq=0, profile_batch=0
            )
        )
        if verbose:
            print(f"[{term}] tensorboard logs -> {log_dir}")
    if extra_callbacks:
        callbacks.extend(extra_callbacks)

    val = splits["val"]
    history = model.fit(
        train_ds,
        validation_data=(val.x, val.targets()),
        epochs=cfg.epochs,
        callbacks=callbacks,
        verbose=verbose,
    )
    return model, {k: [float(v) for v in vals] for k, vals in history.history.items()}


def predict_all(model: keras.Model, x: np.ndarray, batch_size: int = 256) -> dict[str, np.ndarray]:
    """Predict both heads.  Returns ``{"global": (n,), "per_index_cat": (n, seq, K_pi)}``."""
    out = model.predict(x, verbose=0, batch_size=batch_size)
    if not isinstance(out, dict):  # single-output fallback
        out = {"global": out}
    return {
        "global": np.asarray(out["global"], dtype="float64").reshape(-1),
        "per_index_cat": np.asarray(out["per_index_cat"], dtype="float32"),
    }


def embed(model: keras.Model, x: np.ndarray, batch_size: int = 256) -> np.ndarray:
    """L2-normalised output of the ``embed`` dense layer, ``(n, embed_units)``.

    No new loss: the embedding is whatever the global classifier already learned
    to sit on.  Used downstream for nearest-known-peptide retrieval.
    """
    em = embed_model(model)
    return l2norm(em.predict(x, verbose=0, batch_size=batch_size))


@dataclass
class Result:
    """Everything the R side needs back, in ``term_order`` row order."""

    term_order: tuple[str, ...]
    pred: np.ndarray           # (n,)  pooled-Platt calibrated global score
    pred_raw: np.ndarray       # (n,)  uncalibrated global score
    per_index: np.ndarray      # (n, seq_len, K_pi)
    emb: np.ndarray            # (n, embed_units), L2-normalised
    val_scores: np.ndarray     # (m,)  pooled validation raw scores
    val_labels: np.ndarray     # (m,)  pooled validation labels
    n_by_term: dict[str, int]
    pi_names: tuple[str, ...]
    calibrator: PlattCalibrator
    models: dict[str, keras.Model] = field(default_factory=dict, repr=False)
    histories: dict[str, dict] = field(default_factory=dict, repr=False)

    def to_dict(self) -> dict:
        """Plain arrays / lists -- the shape the R bridge and the CLI consume."""
        return {
            "term_order": list(self.term_order),
            "pred": self.pred,
            "pred_raw": self.pred_raw,
            "per_index": self.per_index,
            "emb": self.emb,
            "val_scores": self.val_scores,
            "val_labels": self.val_labels,
            "n_by_term": dict(self.n_by_term),
            "pi_names": list(self.pi_names),
            "calibrator": self.calibrator.to_dict(),
            "histories": self.histories,
        }

    def term_index(self, term: str) -> np.ndarray:
        """Row positions belonging to ``term`` in the concatenated outputs."""
        start = 0
        for t in self.term_order:
            n = self.n_by_term[t]
            if t == term:
                return np.arange(start, start + n)
            start += n
        raise KeyError(term)


def _order_terms(data: Mapping, cfg: Config) -> tuple[str, ...]:
    known = [t for t in cfg.term_order if t in data]
    extra = [t for t in data if t not in cfg.term_order]
    return tuple(known + extra)


def run(
    data: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    cfg: Config | None = None,
    verbose: int = 1,
    calibration_method: str = "irls",
) -> Result:
    """Train every terminus model, score the ``all`` split, calibrate, embed.

    ``data`` is the nested ``data[term][split][array_name]`` mapping described
    in :mod:`lf_dcnn.data`.  Rows of every returned array are the ``all`` splits
    concatenated in ``cfg.term_order``, which is the order R's
    ``bind_rows(lapply(nn_input, function(x) x$all))`` produces.
    """
    cfg = cfg or Config()
    set_seed(cfg.seed)
    terms = _order_terms(data, cfg)
    if not terms:
        raise ValueError("no terminus data supplied")
    td = as_term_data(data, cfg)

    models: dict[str, keras.Model] = {}
    histories: dict[str, dict] = {}
    raw_parts, per_index_parts, emb_parts = [], [], []
    val_score_parts, val_label_parts = [], []
    n_by_term: dict[str, int] = {}

    for term in terms:
        splits = td[term]
        for needed in ("train", "val", "all"):
            if needed not in splits:
                raise KeyError(f"term {term!r} is missing the {needed!r} split")

        if verbose:
            print(f"[{term}] training ({len(splits['train'])} train, "
                  f"{int(splits['train'].y_global.sum())} positive)")
        # a per-term generator keeps the two models' sampling independent of
        # each other's epoch counts
        rng = np.random.default_rng(None if cfg.seed is None else cfg.seed + terms.index(term))
        model, history = train(term, splits, cfg, verbose=verbose, rng=rng)
        if verbose:
            print(f"[{term}] trainable params: {model.count_params():,}")
        models[term] = model
        histories[term] = history

        preds = predict_all(model, splits["all"].x)
        raw_parts.append(preds["global"])
        per_index_parts.append(preds["per_index_cat"])
        emb_parts.append(embed(model, splits["all"].x))
        n_by_term[term] = len(splits["all"])

        val_score_parts.append(predict_all(model, splits["val"].x)["global"])
        val_label_parts.append(splits["val"].y_global[:, 0].astype("float64"))

    pred_raw = np.concatenate(raw_parts)
    val_scores = np.concatenate(val_score_parts)
    val_labels = np.concatenate(val_label_parts)

    # ONE pooled calibrator across both models (see lf_dcnn.calibrate)
    calibrator = PlattCalibrator.fit(val_scores, val_labels, method=calibration_method)
    pred = calibrator.predict(pred_raw)

    return Result(
        term_order=terms,
        pred=pred,
        pred_raw=pred_raw,
        per_index=np.concatenate(per_index_parts, axis=0),
        emb=np.concatenate(emb_parts, axis=0),
        val_scores=val_scores,
        val_labels=val_labels,
        n_by_term=n_by_term,
        pi_names=cfg.pi_names,
        calibrator=calibrator,
        models=models,
        histories=histories,
    )
