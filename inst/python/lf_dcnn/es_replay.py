"""Replay keras EarlyStopping offline, and the metrics the tuner scores with.

Early stopping never changes the path training takes -- it only decides when to
stop and which epoch's weights to hand back. So ONE run that keeps going until
every setting of interest would have stopped answers every (monitor, start,
patience, epoch cap) combination at once: replay each setting's decision on the
logged monitor curve and read the model at the epoch it would have restored.

That is exact, not an approximation, provided the replay is EarlyStopping's own
logic. :class:`KerasEarlyStopping` reproduces ``on_epoch_end`` /
``on_train_end`` of keras 3.6 line for line for the configuration
``lf_dcnn.pipeline.train`` uses (``mode="max"``, ``restore_best_weights=True``,
``min_delta=0``, no baseline), and the selftest checks it against the real
callback on random curves with ties and plateaus.

Two details of that logic matter and are easy to get wrong:

* nothing before ``start_from_epoch`` counts -- a peak at epoch 250 under
  ``start_from_epoch=300`` is never restored (see the member histories trap:
  ``max(g[300:])``, not ``max(g)``);
* weights are restored at the end of training EVEN WITHOUT an early stop, so a
  run that hits its epoch cap still hands back its best epoch since the start.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


# ---------------------------------------------------------------- metrics ---

def keras_pr_auc(y_true, y_score, num_thresholds: int = 200) -> float:
    """``keras.metrics.AUC(curve="PR")`` -- the production monitor -- in numpy.

    Same 200 evenly spaced thresholds with the epsilon end points, the same
    bucketing (a score counts as above threshold ``t`` only when strictly
    greater) and the same Davis & Goadrich interpolation, so it agrees with the
    keras metric to float precision. It is NOT average precision: the
    interpolation and the bucketing both differ, and on ~17 positives the two
    can differ in the second decimal.
    """
    y = np.asarray(y_true).reshape(-1) > 0.5
    s = np.clip(np.asarray(y_score, dtype="float32").reshape(-1), 0.0, 1.0)
    T = int(num_thresholds)
    # keras: bucket = relu(ceil(p * (T - 1)) - 1), computed in float32
    b = np.maximum(np.ceil(s * np.float32(T - 1)) - 1, 0).astype(np.int64)
    tp = np.cumsum(np.bincount(b[y], minlength=T)[::-1])[::-1].astype("float64")
    fp = np.cumsum(np.bincount(b[~y], minlength=T)[::-1])[::-1].astype("float64")
    fn = float(y.sum()) - tp

    def _div(a, d):          # ops.divide_no_nan: 0 wherever the divisor is 0
        out = np.zeros_like(a)
        np.divide(a, d, out=out, where=d != 0)
        return out

    dtp = tp[:-1] - tp[1:]
    p = tp + fp
    dp = p[:-1] - p[1:]
    prec_slope = _div(dtp, np.maximum(dp, 0))
    intercept = tp[1:] - prec_slope * p[1:]
    both = (p[:-1] > 0) & (p[1:] > 0)
    safe_p_ratio = np.where(both, _div(p[:-1], np.maximum(p[1:], 0)), 1.0)
    inc = _div(prec_slope * (dtp + intercept * np.log(np.where(both, safe_p_ratio, 1.0))),
               np.maximum(tp[1:] + fn[1:], 0))
    return float(inc.sum())


def _average_ranks(x: np.ndarray) -> np.ndarray:
    """1-based ranks, ties sharing their average rank (scipy's ``rankdata``)."""
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    new = np.r_[True, xs[1:] != xs[:-1]]
    grp = np.cumsum(new) - 1
    counts = np.bincount(grp)
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    r = np.empty(x.size, dtype="float64")
    r[order] = (starts + (counts + 1) / 2.0)[grp]
    return r


def roc_auc(y_true, y_score) -> float:
    """ROC AUC as the Mann-Whitney statistic, ties counted half. NaN without
    both classes."""
    y = np.asarray(y_true).reshape(-1) > 0.5
    s = np.asarray(y_score, dtype="float64").reshape(-1)
    n1 = int(y.sum())
    n0 = y.size - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = _average_ranks(s)
    return float((r[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def average_precision(y_true, y_score) -> float:
    """Average precision, scikit-learn's definition: sum over the distinct
    thresholds of (recall step) x precision. NaN without positives."""
    y = np.asarray(y_true).reshape(-1) > 0.5
    s = np.asarray(y_score, dtype="float64").reshape(-1)
    n1 = int(y.sum())
    if n1 == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    ys, ss = y[order], s[order]
    last = np.r_[ss[1:] != ss[:-1], True]          # last row of each tie group
    tps = np.cumsum(ys)[last].astype("float64")
    fps = np.cumsum(~ys)[last].astype("float64")
    recall = tps / n1
    precision = tps / (tps + fps)
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def rank_among(pool_sorted: np.ndarray, scores) -> np.ndarray:
    """``1 + #(pool > score)`` for each score -- the candidate-rank convention
    of 10_8_pepend_train.R (``sum(unk > v) + 1``): a tie with the pool does not
    push a known down. ``pool_sorted`` must be ascending."""
    s = np.asarray(scores, dtype=pool_sorted.dtype)
    return 1 + (pool_sorted.size - np.searchsorted(pool_sorted, s, side="right"))


#: monitors the replay can use, each a function of (labels, global scores) on
#: the monitor windows. ``pr_auc`` is what production early-stops on.
MONITORS = {
    "pr_auc": keras_pr_auc,
    "roc_auc": roc_auc,
}


# ----------------------------------------------------- early stopping replay ---

class KerasEarlyStopping:
    """The state of one ``keras.callbacks.EarlyStopping(mode="max",
    restore_best_weights=True, patience=..., start_from_epoch=...)``.

    :meth:`update` is ``on_epoch_end``. ``best_epoch`` is the epoch whose
    weights keras would restore -- ``None`` until the first epoch at or after
    ``start_from_epoch``, in which case keras restores nothing and the last
    epoch's weights stand.
    """

    __slots__ = ("start", "patience", "best", "best_epoch", "wait", "stopped_epoch")

    def __init__(self, start_from_epoch: int, patience: int):
        self.start = int(start_from_epoch)
        self.patience = int(patience)
        self.best = -math.inf
        self.best_epoch: int | None = None
        self.wait = 0
        self.stopped_epoch: int | None = None

    def update(self, epoch: int, value: float) -> str | None:
        """One epoch. Returns ``"best"`` when ``epoch`` becomes the epoch to
        restore, ``"stop"`` when training stops after it, else ``None``."""
        if self.stopped_epoch is not None or value is None or epoch < self.start:
            return None
        event = None
        if self.best_epoch is None:
            # keras: "if best weights were never set, the current weights are
            # the best" -- before the improvement test, so it holds for a NaN too
            self.best_epoch = epoch
            event = "best"
        self.wait += 1
        if value - 0.0 > self.best:                 # min_delta = 0, strict
            self.best = value
            self.best_epoch = epoch
            self.wait = 0
            return "best"
        if self.wait >= self.patience and epoch > 0:
            self.stopped_epoch = epoch
            return "stop"
        return event


def replay_early_stopping(values: Sequence[float], start_from_epoch: int,
                          patience: int, epochs: int) -> tuple[int, int]:
    """What ``fit(epochs=epochs)`` with that EarlyStopping would have done on
    the monitor curve ``values`` (``values[e]`` = the monitor after epoch ``e``,
    0-based).

    Returns ``(last_epoch_trained, restored_epoch)``. A curve shorter than
    ``epochs`` is accepted only if the setting stops inside it; otherwise the
    run that logged it ended before this setting would have, and its answer is
    not in the curve.
    """
    es = KerasEarlyStopping(start_from_epoch, patience)
    n = min(len(values), int(epochs))
    for e in range(n):
        if es.update(e, float(values[e])) == "stop":
            return e, es.best_epoch
    if n < int(epochs):
        raise ValueError(
            f"curve of {len(values)} epochs ends before start={start_from_epoch}, "
            f"patience={patience}, epochs={epochs} would have stopped")
    last = n - 1
    return last, (es.best_epoch if es.best_epoch is not None else last)


@dataclass(frozen=True)
class ESGrid:
    """The early-stopping settings one tuning run replays: every combination of
    ``monitors`` x ``starts`` x ``patiences`` x ``caps`` (epoch budgets).

    The run itself lasts ``max(caps)`` epochs, or less if every setting has
    stopped before that.
    """

    monitors: tuple[str, ...] = ("pr_auc", "roc_auc")
    starts: tuple[int, ...] = (0, 100, 300, 500)
    patiences: tuple[int, ...] = (50, 100, 200, 300, 500)
    caps: tuple[int, ...] = (1000, 2000, 3000)

    def __post_init__(self):
        for f in ("monitors", "starts", "patiences", "caps"):
            v = tuple(getattr(self, f))
            if not v:
                raise ValueError(f"ESGrid.{f} is empty")
            if f != "monitors":
                v = tuple(sorted({int(x) for x in v}))
            object.__setattr__(self, f, v)
        bad = [m for m in self.monitors if m not in MONITORS]
        if bad:
            raise ValueError(f"unknown monitor(s) {bad}; known: {sorted(MONITORS)}")
        if min(self.caps) < 1 or min(self.patiences) < 0 or min(self.starts) < 0:
            raise ValueError("caps must be >= 1, patiences and starts >= 0")

    @property
    def max_epochs(self) -> int:
        return max(self.caps)

    def with_setting(self, monitor: str, start: int, patience: int, cap: int) -> "ESGrid":
        """This grid, widened (if need be) to include one setting."""
        return ESGrid(tuple(dict.fromkeys(self.monitors + (monitor,))),
                      self.starts + (int(start),), self.patiences + (int(patience),),
                      self.caps + (int(cap),))

    def to_dict(self) -> dict:
        return {"monitors": list(self.monitors), "starts": list(self.starts),
                "patiences": list(self.patiences), "caps": list(self.caps)}

    @classmethod
    def from_dict(cls, d: Mapping) -> "ESGrid":
        return cls(tuple(d["monitors"]), tuple(d["starts"]), tuple(d["patiences"]),
                   tuple(d["caps"]))


class ESGridTracker:
    """Every setting of an :class:`ESGrid`, for every rotation, run online.

    One :class:`KerasEarlyStopping` per ``(monitor, rotation, start,
    patience)``; caps need no state of their own, because a cap only decides
    whether a setting ran out of epochs before its patience did. :meth:`step`
    returns, for one epoch, whether some setting's best moved to this epoch (so
    its weights must be kept) and the settings whose answer became final, each
    as ``((monitor, rotation, start, patience, cap), last_epoch, restored)``.

    All settings of one ``(monitor, rotation, start)`` share the same best epoch
    until they stop, so at most ``monitors x rotations x starts`` snapshots are
    ever live (:meth:`live_epochs`).
    """

    def __init__(self, grid: ESGrid, rotations: Iterable[int]):
        self.grid = grid
        self.rotations = tuple(rotations)
        self.states = {
            (m, r, s, p): KerasEarlyStopping(s, p)
            for m in grid.monitors for r in self.rotations
            for s in grid.starts for p in grid.patiences
        }
        self.finished = False

    def step(self, epoch: int, values: Mapping[tuple[str, int], float]):
        snapshot = False
        finals = []
        for key, st in self.states.items():
            if st.stopped_epoch is not None:
                continue
            ev = st.update(epoch, values[key[0], key[1]])
            if ev == "best":
                snapshot = True
            elif ev == "stop":
                # patience ran out at `epoch`: the answer for every cap that
                # would still have been running here
                for cap in self.grid.caps:
                    if cap - 1 >= epoch:
                        finals.append((key + (cap,), epoch, st.best_epoch))
        for cap in self.grid.caps:
            if epoch == cap - 1:
                # the cap ends training: keras still restores the best epoch
                for key, st in self.states.items():
                    if st.stopped_epoch is None:
                        finals.append((key + (cap,), epoch,
                                       st.best_epoch if st.best_epoch is not None else epoch))
        if epoch >= self.grid.max_epochs - 1:
            self.finished = True
        return snapshot, finals

    @property
    def done(self) -> bool:
        """Every setting's answer is final, so training can stop."""
        return self.finished or all(st.stopped_epoch is not None for st in self.states.values())

    def live_epochs(self) -> set[int]:
        """Epochs some still-running setting could yet restore."""
        if self.finished:
            return set()
        return {st.best_epoch for st in self.states.values()
                if st.stopped_epoch is None and st.best_epoch is not None}
