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
    val_scores: np.ndarray     # (m,)  pooled validation raw scores, ensemble MEAN
    val_labels: np.ndarray     # (m,)  pooled validation labels
    #: per-member raw scores, ``(n_seeds, m)`` / ``(n_seeds, n)``: what a
    #: benchmark needs to quote an AUC as mean +/- sd ACROSS SEEDS rather than
    #: the single AUC of the averaged score (which is a different, usually
    #: higher, number). Row k is member k, in seed order.
    val_scores_members: np.ndarray
    pred_raw_members: np.ndarray
    n_by_term: dict[str, int]
    pi_names: tuple[str, ...]
    calibrator: PlattCalibrator
    #: trainable parameters per terminus; kept apart from ``models`` because a
    #: run combined from on-disk members has no live models to count
    params: dict[str, int] = field(default_factory=dict)
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
            "val_scores_members": self.val_scores_members,
            "pred_raw_members": self.pred_raw_members,
            "n_by_term": dict(self.n_by_term),
            "pi_names": list(self.pi_names),
            "calibrator": self.calibrator.to_dict(),
            "params": dict(self.params),
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


@dataclass
class Member:
    """One ensemble member's outputs over every terminus, before combining.

    Everything :func:`combine` needs, and nothing that cannot be written to disk
    (``models`` is the exception and is only populated in-process).  ``k`` is the
    member's index, ``seed`` the Config seed it trained with.
    """

    k: int
    seed: int | None
    terms: tuple[str, ...]
    global_all: dict[str, np.ndarray]      # term -> (n,)          raw global score
    per_index_all: dict[str, np.ndarray]   # term -> (n, seq, K)   per-index softmax
    val: dict[str, np.ndarray]             # term -> (m,)          raw val scores
    val_labels: dict[str, np.ndarray]      # term -> (m,)
    emb: dict[str, np.ndarray]             # term -> (n, embed)    L2-normalised
    histories: dict[str, dict]
    params: dict[str, int]
    models: dict[str, keras.Model] = field(default_factory=dict, repr=False)


def member_seed(cfg: Config, k: int, n_seeds: int) -> int | None:
    """Config seed of member ``k``: ``cfg.seed + k``, or ``cfg.seed`` untouched
    for a single-member run (see the note in :func:`run_member`)."""
    if n_seeds == 1 or cfg.seed is None:
        return cfg.seed
    return int(cfg.seed) + k


def run_member(
    data: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    cfg: Config | None = None,
    k: int = 0,
    n_seeds: int = 1,
    verbose: int = 1,
    keep_models: bool = True,
    terms: Sequence[str] | None = None,
) -> Member:
    """Train ONE ensemble member (seed ``cfg.seed + k``) on every terminus,
    or only on ``terms``.

    The SECOND model trained in a Keras/TF process can die mid-fit in a
    retraced tf.function -- an optimizer slot, an AUC-metric variable or the
    regularization-loss sum suddenly reads as shape ``[0]``, as if it captured
    resources the previous model freed.  The standalone CLI's ``--isolated``
    mode therefore runs ONE model per process: this function with a single
    ``terms`` entry, and :func:`combine` afterwards.  ``Member.terms`` records
    what was trained; per-terminus members merge back in ``io.load_members``.
    """
    cfg = cfg or Config()
    all_terms = _order_terms(data, cfg)
    if not all_terms:
        raise ValueError("no terminus data supplied")
    if terms is None:
        terms = all_terms
    else:
        terms = tuple(terms)
        unknown = [t for t in terms if t not in all_terms]
        if unknown:
            raise KeyError(f"terms {unknown} not in the data ({all_terms})")
    td = as_term_data(data, cfg)
    cfg_k = cfg if n_seeds == 1 else cfg.evolve(seed=member_seed(cfg, k, n_seeds))

    m = Member(k=int(k), seed=cfg_k.seed, terms=tuple(terms), global_all={}, per_index_all={},
               val={}, val_labels={}, emb={}, histories={}, params={})
    for term in terms:
        splits = td[term]
        for needed in ("train", "val", "all"):
            if needed not in splits:
                raise KeyError(f"term {term!r} is missing the {needed!r} split")
        # Re-seed per member so the ensemble members differ. With a single
        # member, leave the stream alone: run() already seeded once, and
        # re-seeding here would reset it per TERM, changing the second term's
        # model relative to the pre-ensemble behaviour.
        if n_seeds > 1:
            set_seed(cfg_k.seed)
        # a per-term generator keeps the termini's sampling independent
        rng = np.random.default_rng(None if cfg_k.seed is None
                                    else cfg_k.seed + all_terms.index(term))
        if verbose:
            print(f"[{term}] member {k + 1}/{n_seeds} (seed {cfg_k.seed}): "
                  f"{len(splits['train'])} train, "
                  f"{int(splits['train'].y_global.sum())} positive")
        model, history = train(term, splits, cfg_k, verbose=verbose, rng=rng)
        if verbose:
            print(f"[{term}] member {k + 1}/{n_seeds} (seed {cfg_k.seed}), "
                  f"{model.count_params():,} params")

        p = predict_all(model, splits["all"].x)
        m.global_all[term] = p["global"]
        m.per_index_all[term] = p["per_index_cat"]
        m.val[term] = predict_all(model, splits["val"].x)["global"]
        m.val_labels[term] = splits["val"].y_global[:, 0].astype("float64")
        m.emb[term] = embed(model, splits["all"].x)
        m.histories[term] = history
        m.params[term] = int(model.count_params())
        if keep_models:
            m.models[term] = model
    return m


def combine(
    members: Sequence[Member],
    cfg: Config | None = None,
    calibration_method: str = "irls",
) -> Result:
    """Reduce members to the ensemble :class:`Result`.

    Global score and per-residue softmax become their mean and sd across
    members; the embedding, models and histories are the FIRST member's (see
    :func:`run`); one pooled calibrator is fit on the mean validation scores.
    """
    cfg = cfg or Config()
    if not members:
        raise ValueError("no members to combine")
    members = sorted(members, key=lambda m: m.k)
    first = members[0]
    terms = first.terms
    n_seeds = len(members)
    for m in members:
        if m.terms != terms:
            raise ValueError(f"member {m.k} has terms {m.terms}, expected {terms}")

    raw_parts, raw_sd_parts, emb_parts = [], [], []
    per_index_parts, per_index_sd_parts = [], []
    val_score_parts, val_label_parts = [], []
    raw_member_parts, val_member_parts = [], []   # per term: (n_seeds, n) / (n_seeds, m)
    n_by_term: dict[str, int] = {}

    for term in terms:
        # running sums, as before: avoids stacking every member's per-index
        # softmax in float64 at once
        g_sum = g_sq = pi_sum = pi_sq = v_sum = None
        for m in members:
            g = np.asarray(m.global_all[term], dtype="float64")
            pi = np.asarray(m.per_index_all[term], dtype="float64")
            v = np.asarray(m.val[term], dtype="float64")
            if g_sum is None:
                g_sum, g_sq = g.copy(), g**2
                pi_sum, pi_sq = pi.copy(), pi**2
                v_sum = v.copy()
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
        n_by_term[term] = int(g_mean.size)
        emb_parts.append(first.emb[term])
        val_score_parts.append(v_sum / n_seeds)
        val_label_parts.append(np.asarray(first.val_labels[term], dtype="float64"))
        raw_member_parts.append(np.stack([np.asarray(m.global_all[term], dtype="float64")
                                          for m in members]))
        val_member_parts.append(np.stack([np.asarray(m.val[term], dtype="float64")
                                          for m in members]))

    pred_raw = np.concatenate(raw_parts)
    pred_sd = np.concatenate(raw_sd_parts)
    val_scores = np.concatenate(val_score_parts)
    val_labels = np.concatenate(val_label_parts)
    # members stay on the row axis; windows concatenate across terms like the means
    pred_raw_members = np.concatenate(raw_member_parts, axis=1)
    val_scores_members = np.concatenate(val_member_parts, axis=1)

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
        val_scores_members=val_scores_members,
        pred_raw_members=pred_raw_members,
        n_by_term=n_by_term,
        pi_names=cfg.pi_names,
        calibrator=calibrator,
        params=dict(first.params),
        models=dict(first.models),
        histories=dict(first.histories),
    )


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

    This is :func:`run_member` for every member followed by :func:`combine`,
    all in one process.  For more than a few members prefer the CLI's
    ``--isolated`` mode, which runs each member in its own process.
    """
    cfg = cfg or Config()
    set_seed(cfg.seed)
    members = [
        run_member(data, cfg, k=k, n_seeds=n_seeds, verbose=verbose, keep_models=(k == 0))
        for k in range(int(n_seeds))
    ]
    return combine(members, cfg, calibration_method=calibration_method)
