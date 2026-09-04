"""The array contract and the oversampled, noise-augmented training sampler.

Per terminus (``"N"``, ``"C"``) and per split (``train``, ``val``, ``all``):

==================  =========================  ==============================
name                shape                      contents
==================  =========================  ==============================
``x``               ``(n, seq_len, C)``        float32 input channels
``y_global``        ``(n, 1)``                 float32 0/1 window label
``y_per_index_cat`` ``(n, seq_len, K_pi)``     float32; cols ``0:K_cat`` are
                                               the ``[6 real, none]`` one-hot,
                                               col ``K_pi - 1`` the padding flag
==================  =========================  ==============================

R builds these by ``aperm``-ing a ``(seq_len, C, n)`` array to ``(n, seq_len,
C)``.  Do the ``aperm`` on the R side: R is column-major and numpy is
row-major, so handing over the pre-``aperm`` array transposes silently and
still trains, just badly.  :func:`TermArrays.validate` catches the shape half
of that mistake.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

import keras
import numpy as np

from .config import Config

SPLITS = ("train", "val", "all")
ARRAY_NAMES = ("x", "y_global", "y_per_index_cat")


@dataclass
class TermArrays:
    """One split's three arrays."""

    x: np.ndarray
    y_global: np.ndarray
    y_per_index_cat: np.ndarray

    @classmethod
    def from_mapping(cls, m: Mapping[str, np.ndarray]) -> "TermArrays":
        return cls(**{k: np.asarray(m[k], dtype="float32") for k in ARRAY_NAMES})

    def __len__(self) -> int:
        return int(self.x.shape[0])

    def validate(self, cfg: Config, where: str = "") -> "TermArrays":
        n = len(self)
        want = {
            "x": (n, cfg.seq_len, cfg.n_channels),
            "y_global": (n, 1),
            "y_per_index_cat": (n, cfg.seq_len, cfg.K_pi),
        }
        for name, shape in want.items():
            got = tuple(getattr(self, name).shape)
            if got != shape:
                raise ValueError(
                    f"{where}{name}: expected shape {shape}, got {got}. "
                    "If the sample axis is not first, the array was handed over "
                    "before its aperm/transpose."
                )
        return self

    def targets(self) -> dict[str, np.ndarray]:
        return {"global": self.y_global, "per_index_cat": self.y_per_index_cat}


def as_term_data(
    data: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]], cfg: Config
) -> dict[str, dict[str, TermArrays]]:
    """Coerce and validate the nested ``data[term][split][name]`` mapping."""
    out: dict[str, dict[str, TermArrays]] = {}
    for term in data:
        out[term] = {}
        for split in data[term]:
            out[term][split] = TermArrays.from_mapping(data[term][split]).validate(
                cfg, where=f"{term}/{split}: "
            )
    return out


class OversampledWindows(keras.utils.PyDataset):
    """All positives + ``k x n_pos`` negatives, reshuffled and re-noised.

    The window set is tiny and heavily imbalanced, so every epoch draws a fresh
    negative sample and fresh augmentation noise.  Continuous channels are
    jittered at REAL (non-padding) positions only, by
    ``N(0, noise_frac * sd_of_that_channel)`` where the sd is taken over the
    real positions of the drawn subset.  One-hot / flag channels are left
    untouched so augmented inputs stay on the manifold; validation data is
    never touched.

    Set ``cfg.resample_each_epoch = False`` to draw once and reuse, which is
    what ``10_1dcnn_new6.R`` actually does -- its ``on_epoch_begin`` callback
    rebinds a variable that ``fit`` no longer reads.
    """

    def __init__(self, arrays: TermArrays, cfg: Config, rng: np.random.Generator | None = None):
        super().__init__(workers=1, use_multiprocessing=False, max_queue_size=10)
        self.arrays = arrays
        self.cfg = cfg
        self.rng = rng if rng is not None else np.random.default_rng(cfg.seed)

        y = arrays.y_global[:, 0]
        self.pos_idx = np.flatnonzero(y == 1)
        self.neg_idx = np.flatnonzero(y == 0)
        if self.pos_idx.size == 0:
            raise ValueError("training split has no positive windows")
        if self.neg_idx.size == 0:
            raise ValueError("training split has no negative windows")

        self._resample()

    # --- sampling -----------------------------------------------------------
    @property
    def n_selected(self) -> int:
        return int(self.pos_idx.size * (1 + self.cfg.n_negatives_per_positive))

    @property
    def steps_per_epoch(self) -> int:
        return int(math.ceil(self.n_selected / self.cfg.batch_size))

    def _shuffle(self) -> None:
        """Re-order the current draw.

        The R pipeline puts a full-buffer ``dataset_shuffle`` before ``batch``,
        and tf.data reshuffles it on every iteration -- so batch composition
        varies each epoch even when the underlying draw does not.  Mirror that,
        or the fixed-draw path would train on identical batches every epoch.
        """
        order = self.rng.permutation(self._x.shape[0])
        self._x, self._yg, self._yc = self._x[order], self._yg[order], self._yc[order]

    def _resample(self) -> None:
        cfg = self.cfg
        n_neg_target = int(self.pos_idx.size * cfg.n_negatives_per_positive)
        # with replacement only when there are not enough distinct negatives
        replace = self.neg_idx.size < n_neg_target
        sampled_neg = self.rng.choice(self.neg_idx, size=n_neg_target, replace=replace)
        idx = np.concatenate([self.pos_idx, sampled_neg])
        self.rng.shuffle(idx)

        x_sel = self.arrays.x[idx].copy()
        if cfg.noise_frac > 0 and cfg.cont_channel_ids:
            # TRUE at real positions; padding positions stay exactly 0
            real_mask = x_sel[:, :, cfg.pad_channel_id] == 0
            if real_mask.any():
                for ch in cfg.cont_channel_ids:
                    sl = x_sel[:, :, ch]
                    sd_ch = np.nanstd(sl[real_mask], ddof=1)
                    if np.isfinite(sd_ch) and sd_ch > 0:
                        noise = self.rng.normal(0.0, cfg.noise_frac * sd_ch, size=sl.shape)
                        x_sel[:, :, ch] = sl + noise * real_mask

        self._x = x_sel.astype("float32")
        self._yg = self.arrays.y_global[idx].astype("float32")
        self._yc = self.arrays.y_per_index_cat[idx].astype("float32")

    # --- keras.utils.PyDataset ---------------------------------------------
    def __len__(self) -> int:
        return self.steps_per_epoch

    def __getitem__(self, i):
        lo = i * self.cfg.batch_size
        hi = min(lo + self.cfg.batch_size, self._x.shape[0])
        return (
            self._x[lo:hi],
            {"global": self._yg[lo:hi], "per_index_cat": self._yc[lo:hi]},
        )

    def on_epoch_end(self) -> None:
        # Keras calls this after every epoch, so the next epoch trains on a
        # fresh draw -- equivalent to the intended on_epoch_begin rebuild.
        # Without the redraw we still reshuffle, as the R tf.data pipeline does.
        if self.cfg.resample_each_epoch:
            self._resample()
        else:
            self._shuffle()


def build_dataset(
    arrays: TermArrays, cfg: Config, rng: np.random.Generator | None = None
) -> OversampledWindows:
    """Convenience wrapper matching the R ``make_oversampled_dataset`` name."""
    return OversampledWindows(arrays, cfg, rng=rng)
