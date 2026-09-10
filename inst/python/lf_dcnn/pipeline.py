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
    pred_raw: np.ndarray       # (n,)  uncalibrated global score, ensemble MEAN
    pred_sd: np.ndarray        # (n,)  sd of the raw score across ensemble members
    per_index: np.ndarray      # (n, seq_len, K_pi)  ensemble MEAN
    per_index_sd: np.ndarray   # (n, seq_len, K_pi)  sd across ensemble members
    n_seeds: int
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
            "pred_sd": self.pred_sd,
            "per_index": self.per_index,
            "per_index_sd": self.per_index_sd,
            "n_seeds": self.n_seeds,
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


def _mean_sd(total, total_sq, n):
    """Mean and sample sd from running sums -- avoids holding every member."""
    mean = total / n
    if n < 2:
        return mean, np.zeros_like(mean)
    var = (total_sq - n * mean**2) / (n - 1)
    return mean, np.sqrt(np.maximum(var, 0.0))


def run(
    data: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    cfg: Config | None = None,
    verbose: int = 1,
    calibration_method: str = "irls",
    n_seeds: int = 1,
) -> Result:
    """Train every terminus model, score the ``all`` split, calibrate, embed.

    ``data`` is the nested ``data[term][split][array_name]`` mapping described
    in :mod:`lf_dcnn.data`.  Rows of every returned array are the ``all`` splits
    concatenated in ``cfg.term_order``, which is the order R's
    ``bind_rows(lapply(nn_input, function(x) x$all))`` produces.

    With ``n_seeds > 1`` each terminus trains that many members (seeds
    ``cfg.seed`` .. ``cfg.seed + n_seeds - 1``) and the global score and the
    per-residue softmax are reduced to their mean and sd across members. The sd
    is the useful half: a single model's softmax cannot distinguish "confidently
    0.5" from "the members disagree violently", and those mean opposite things
    when reading a per-residue profile.

    The embedding is taken from the FIRST member only. Averaging embeddings
    across independently initialised models is meaningless -- each learns its
    own basis -- so nearest-known retrieval uses one member's space.
    """
    cfg = cfg or Config()
    set_seed(cfg.seed)
    terms = _order_terms(data, cfg)
    if not terms:
        raise ValueError("no terminus data supplied")
    td = as_term_data(data, cfg)

    models: dict[str, keras.Model] = {}
    histories: dict[str, dict] = {}
    raw_parts, raw_sd_parts, emb_parts = [], [], []
    per_index_parts, per_index_sd_parts = [], []
    val_score_parts, val_label_parts = [], []
    n_by_term: dict[str, int] = {}

    for term in terms:
        splits = td[term]
        for needed in ("train", "val", "all"):
            if needed not in splits:
                raise KeyError(f"term {term!r} is missing the {needed!r} split")

        if verbose:
            print(f"[{term}] training {n_seeds} member(s) "
                  f"({len(splits['train'])} train, "
                  f"{int(splits['train'].y_global.sum())} positive)")

        g_sum = g_sq = pi_sum = pi_sq = v_sum = None
        base = 0 if cfg.seed is None else cfg.seed
        for k in range(n_seeds):
            cfg_k = cfg if n_seeds == 1 else cfg.evolve(seed=base + k)
            # Re-seed per member so the ensemble members differ. With a single
            # member, leave the stream alone: run() already seeded once, and
            # re-seeding here would reset it per TERM, changing the second
            # term's model relative to the pre-ensemble behaviour.
            if n_seeds > 1:
                set_seed(cfg_k.seed)
            # a per-term generator keeps the termini's sampling independent
            rng = np.random.default_rng(None if cfg_k.seed is None
                                        else cfg_k.seed + terms.index(term))
            model, history = train(term, splits, cfg_k, verbose=verbose, rng=rng)
            if verbose:
                print(f"[{term}] member {k + 1}/{n_seeds} (seed {cfg_k.seed}), "
                      f"{model.count_params():,} params")

            p = predict_all(model, splits["all"].x)
            v = predict_all(model, splits["val"].x)["global"]
            g, pi = p["global"], p["per_index_cat"].astype("float64")
            if g_sum is None:
                g_sum, g_sq = g.copy(), g**2
                pi_sum, pi_sq = pi.copy(), pi**2
                v_sum = v.copy()
                models[term] = model              # first member is the kept one
                histories[term] = history
                emb_parts.append(embed(model, splits["all"].x))
            else:
                g_sum += g; g_sq += g**2
                pi_sum += pi; pi_sq += pi**2
                v_sum += v

        g_mean, g_sd = _mean_sd(g_sum, g_sq, n_seeds)
        pi_mean, pi_sd = _mean_sd(pi_sum, pi_sq, n_seeds)
        raw_parts.append(g_mean)
        raw_sd_parts.append(g_sd)
        per_index_parts.append(pi_mean.astype("float32"))
        per_index_sd_parts.append(pi_sd.astype("float32"))
        n_by_term[term] = len(splits["all"])

        val_score_parts.append(v_sum / n_seeds)
        val_label_parts.append(splits["val"].y_global[:, 0].astype("float64"))

    pred_raw = np.concatenate(raw_parts)
    pred_sd = np.concatenate(raw_sd_parts)
    val_scores = np.concatenate(val_score_parts)
    val_labels = np.concatenate(val_label_parts)

    # ONE pooled calibrator across both models (see lf_dcnn.calibrate)
    calibrator = PlattCalibrator.fit(val_scores, val_labels, method=calibration_method)
    pred = calibrator.predict(pred_raw)

    return Result(
        term_order=terms,
        pred=pred,
        pred_raw=pred_raw,
        pred_sd=pred_sd,
        per_index=np.concatenate(per_index_parts, axis=0),
        per_index_sd=np.concatenate(per_index_sd_parts, axis=0),
        n_seeds=int(n_seeds),
        emb=np.concatenate(emb_parts, axis=0),
        val_scores=val_scores,
        val_labels=val_labels,
        n_by_term=n_by_term,
        pi_names=cfg.pi_names,
        calibrator=calibrator,
        models=models,
        histories=histories,
    )
