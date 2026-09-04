"""Synthetic windows matching the array contract, for tests and demos.

Deliberately mildly learnable: positives carry a ``pep_pocket`` stretch in the
middle plus a lift in one continuous channel, so a short run produces a
non-degenerate PR-AUC and a calibrator that actually has two classes to fit.
"""

from __future__ import annotations

import numpy as np

from .config import Config


def make_term(
    cfg: Config,
    n: int,
    n_pos: int,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """One split's ``x`` / ``y_global`` / ``y_per_index_cat``."""
    seq, C, K_pi = cfg.seq_len, cfg.n_channels, cfg.K_pi
    pi_names = list(cfg.pi_names)
    i_pad = pi_names.index("padding")
    i_none = pi_names.index("none")
    i_pocket = pi_names.index("pep_pocket")
    i_db = pi_names.index("DB")
    lo, hi = cfg.mid_span

    x = rng.random((n, seq, C)).astype("float32")
    y_global = np.zeros((n, 1), dtype="float32")
    y_global[:n_pos, 0] = 1.0

    cls = np.full((n, seq), i_none, dtype="int64")
    # a couple of padded positions at the front of every window
    n_pad = rng.integers(0, 4, size=n)
    for i in range(n):
        cls[i, : n_pad[i]] = i_pad
        cls[i, cfg.ct_span[0] - 1 : cfg.ct_span[0] + 1] = i_db
        if y_global[i, 0] == 1:
            a = rng.integers(lo - 1, hi - 6)
            cls[i, a : a + 5] = i_pocket
            x[i, a : a + 5, 1] += 0.6            # the learnable lift

    x[:, :, cfg.pad_channel_id] = (cls == i_pad).astype("float32")
    y_per_index_cat = np.eye(K_pi, dtype="float32")[cls]

    order = rng.permutation(n)
    return {
        "x": x[order],
        "y_global": y_global[order],
        "y_per_index_cat": y_per_index_cat[order],
    }


def make_data(
    cfg: Config | None = None,
    n_train: int = 96,
    n_pos_train: int = 12,
    n_val: int = 40,
    n_pos_val: int = 6,
    n_all: int = 60,
    n_pos_all: int = 8,
    seed: int = 0,
) -> dict[str, dict[str, dict[str, np.ndarray]]]:
    """A full ``data[term][split][name]`` mapping for both termini."""
    cfg = cfg or Config()
    rng = np.random.default_rng(seed)
    return {
        term: {
            "train": make_term(cfg, n_train, n_pos_train, rng),
            "val": make_term(cfg, n_val, n_pos_val, rng),
            "all": make_term(cfg, n_all, n_pos_all, rng),
        }
        for term in cfg.term_order
    }
