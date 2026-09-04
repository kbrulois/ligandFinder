"""Configuration for the lf_dcnn peptide-cleavage-window models.

Everything the model, the losses and the sampler need to agree on lives here:
class order, the position masks, the per-index class weights and the training
hyper-parameters.  The R original (``inst/scripts/10_1dcnn_new6.R``) states all
of these in 1-based, column-major terms; the translation to 0-based numpy
indices is done once, here, and nowhere else.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from typing import Mapping

import numpy as np

# --- class vocabulary --------------------------------------------------------
# Order is load-bearing: it is the order R's `classes` vector uses, and the
# one-hot columns of `y_per_index_cat` are derived from it.
CLASS_NAMES: tuple[str, ...] = (
    "CT_cleavage_context",
    "DB",
    "gap",
    "NT_cleavage_context",
    "pep_other",
    "pep_pocket",
    "padding",
    "none",
)

#: default input channel order, mirroring `all_params3` in 9.2_add_contact_data.R
DEFAULT_CHANNEL_NAMES: tuple[str, ...] = (
    "cons_rs_n", "min_afm", "mean_afm", "relASA",
    "SS_P", "SS_S", "SS_E", "SS_-", "SS_T", "SS_G", "SS_B", "SS_H", "SS_I",
    "Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin",
    "NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy",
    "AA_hydro", "AA_charge", "AA_mw", "AA_pI",
    "padding",
)

#: channels that receive train-time noise augmentation. One-hot / flag channels
#: (SS_*, padding) are left intact so augmented inputs stay on the manifold, and
#: the AA_* property channels are excluded too: they are a deterministic lookup
#: on the residue letter, so perturbing them yields amino acids that do not exist.
DEFAULT_CONT_CHANNELS: tuple[str, ...] = (
    "cons_rs", "cons_rs_n", "min_afm", "mean_afm", "relASA",
)


@dataclass(frozen=True)
class Config:
    """Immutable model/training configuration.

    Use :meth:`r_exact` for a configuration that reproduces
    ``10_1dcnn_new6.R`` bit-for-bit, including the two places where that script
    diverges from its own documented intent (see ``db_mask_both_termini`` and
    ``resample_each_epoch``).
    """

    # --- geometry ---
    seq_len: int = 36
    n_channels: int = len(DEFAULT_CHANNEL_NAMES)
    channel_names: tuple[str, ...] = DEFAULT_CHANNEL_NAMES
    class_names: tuple[str, ...] = CLASS_NAMES
    term_order: tuple[str, ...] = ("N", "C")

    # --- position masks, 1-based inclusive as in the R source ---
    nt_span: tuple[int, int] = (1, 5)
    ct_span: tuple[int, int] = (31, 36)
    mid_span: tuple[int, int] = (6, 30)

    # --- architecture ---
    l2: float = 1e-3
    conv_filters: tuple[int, ...] = (16, 8)
    conv_kernel: int = 3
    conv_dropout: tuple[float, ...] = (0.2, 0.3)
    attention_heads: int = 2
    attention_key_dim: int = 8
    embed_units: int = 16
    embed_dropout: float = 0.3

    # --- losses ---
    gamma: float = 2.0
    smoothness_weight: float = 0.01
    loss_weight_global: float = 0.05
    loss_weight_per_index: float = 1.0
    # The oversampler already rebalances each batch to ~1:3, so ALSO upweighting
    # positives here would double-correct the imbalance (the old ~20x weight).
    bce_weight_0: float = 1.0
    bce_weight_1: float = 1.0
    # per-index class weights: everything 1, DB doubled, and the anchor-side
    # cleavage context tripled (CT for the C model, NT for the N model).
    pi_weight_db: float = 2.0
    pi_weight_anchor: float = 3.0

    # --- optimisation ---
    learning_rate: float = 3e-4
    clipnorm: float = 1.0
    epochs: int = 2000
    batch_size: int = 32
    n_negatives_per_positive: int = 3
    noise_frac: float = 0.1
    monitor: str = "val_global_pr_auc"
    patience: int = 300
    start_from_epoch: int = 300

    # --- channel roles ---
    pad_channel_name: str = "padding"
    cont_channel_names: tuple[str, ...] = DEFAULT_CONT_CHANNELS

    # --- documented-intent switches (see class docstring) ---
    #: ``True`` (default, and what the design calls for) allows the DB class at
    #: BOTH termini in the global head's position mask.  The R script assigns the
    #: NT and CT masks to overlapping column sets in sequence, so its second
    #: assignment silently overwrites DB with the CT mask alone; ``False``
    #: reproduces that.
    db_mask_both_termini: bool = True
    #: ``True`` (default, and what the design calls for) redraws the negative
    #: sample and the augmentation noise every epoch.  The R script's
    #: ``on_epoch_begin`` callback rebinds its ``train_ds`` variable, but ``fit``
    #: already holds the original dataset, so the redraw never reaches training;
    #: ``False`` reproduces that (one fixed draw for the whole run).
    resample_each_epoch: bool = True

    seed: int | None = 1

    # ------------------------------------------------------------------ derived
    #: fields R may hand over as doubles; Keras needs real ints
    _INT_FIELDS = (
        "seq_len", "n_channels", "conv_kernel", "attention_heads",
        "attention_key_dim", "embed_units", "epochs", "batch_size",
        "n_negatives_per_positive", "patience", "start_from_epoch",
    )

    def __post_init__(self) -> None:
        # reticulate hands R numerics over as floats; coerce so Keras gets ints
        for f in self._INT_FIELDS:
            object.__setattr__(self, f, int(getattr(self, f)))
        if self.seed is not None:
            object.__setattr__(self, "seed", int(self.seed))
        for f in ("conv_dropout", "nt_span", "ct_span", "mid_span"):
            object.__setattr__(self, f, tuple(getattr(self, f)))
        object.__setattr__(self, "conv_filters", tuple(int(v) for v in self.conv_filters))
        object.__setattr__(self, "channel_names", tuple(self.channel_names))
        object.__setattr__(self, "class_names", tuple(self.class_names))
        object.__setattr__(self, "term_order", tuple(self.term_order))
        object.__setattr__(self, "cont_channel_names", tuple(self.cont_channel_names))
        if self.channel_names and len(self.channel_names) != self.n_channels:
            raise ValueError(
                f"n_channels={self.n_channels} but got {len(self.channel_names)} channel names"
            )
        for name in ("padding", "none"):
            if name not in self.class_names:
                raise ValueError(f"class_names must contain {name!r}")

    @classmethod
    def r_exact(cls, **kwargs) -> "Config":
        """Config reproducing ``10_1dcnn_new6.R`` including its two quirks."""
        kwargs.setdefault("db_mask_both_termini", False)
        kwargs.setdefault("resample_each_epoch", False)
        return cls(**kwargs)

    def evolve(self, **kwargs) -> "Config":
        """Return a copy with ``kwargs`` overridden."""
        return replace(self, **kwargs)

    # --- class column bookkeeping (all 0-based) ---
    @property
    def class_index(self) -> dict[str, int]:
        return {name: i for i, name in enumerate(self.class_names)}

    @property
    def real_cols(self) -> tuple[int, ...]:
        """Columns of the real classes: everything but ``none`` and ``padding``."""
        return tuple(
            i for i, n in enumerate(self.class_names) if n not in ("none", "padding")
        )

    @property
    def none_col(self) -> int:
        return self.class_index["none"]

    @property
    def padding_col(self) -> int:
        return self.class_index["padding"]

    @property
    def cat_cols(self) -> tuple[int, ...]:
        """The 6 real classes plus an explicit ``none`` (last)."""
        return self.real_cols + (self.none_col,)

    @property
    def K_cat(self) -> int:
        return len(self.cat_cols)

    @property
    def pi_cols(self) -> tuple[int, ...]:
        """Per-index softmax column order: ``[6 real, none, padding]``."""
        return self.cat_cols + (self.padding_col,)

    @property
    def K_pi(self) -> int:
        return len(self.pi_cols)

    @property
    def pi_names(self) -> tuple[str, ...]:
        return tuple(self.class_names[c] for c in self.pi_cols)

    @property
    def none_index(self) -> int:
        """Index of ``none`` within the per-index softmax (== ``K_cat - 1``)."""
        return self.K_cat - 1

    @property
    def padding_index(self) -> int:
        """Index of ``padding`` within the per-index softmax (== ``K_pi - 1``)."""
        return self.K_pi - 1

    # --- position masks ---
    def _span_mask(self, span: tuple[int, int]) -> np.ndarray:
        """1-based inclusive ``span`` -> a 0/1 vector of length ``seq_len``."""
        lo, hi = span
        m = np.zeros(self.seq_len, dtype="float32")
        m[lo - 1 : hi] = 1.0
        return m

    @property
    def nt_mask(self) -> np.ndarray:
        return self._span_mask(self.nt_span)

    @property
    def ct_mask(self) -> np.ndarray:
        return self._span_mask(self.ct_span)

    @property
    def mid_mask(self) -> np.ndarray:
        return self._span_mask(self.mid_span)

    @property
    def mask_matrix_cat(self) -> np.ndarray:
        """``(seq_len, K_cat)`` 0/1 mask over ``[6 real, none]``.

        ``none`` is allowed everywhere; the terminus contexts only at their own
        terminus; ``gap``/``pep_other``/``pep_pocket`` only in the middle.
        """
        by_class = {
            "NT_cleavage_context": self.nt_mask,
            "CT_cleavage_context": self.ct_mask,
            "gap": self.mid_mask,
            "pep_other": self.mid_mask,
            "pep_pocket": self.mid_mask,
            "DB": (
                np.maximum(self.nt_mask, self.ct_mask)
                if self.db_mask_both_termini
                else self.ct_mask
            ),
            "none": np.ones(self.seq_len, dtype="float32"),
        }
        cols = [
            by_class.get(self.class_names[c], np.zeros(self.seq_len, dtype="float32"))
            for c in self.cat_cols
        ]
        return np.stack(cols, axis=1).astype("float32")

    # --- per-index class weights ---
    def pi_weights(self, term: str) -> np.ndarray:
        """Per-model per-index class weights, normalised to mean 1.

        A hand-set prior, not learned: emphasise the anchor-side cleavage
        context (CT for the C model, NT for the N model) plus a moderate DB
        boost.  Mean-1 keeps the per-index loss magnitude stable.
        """
        w = np.ones(self.K_pi, dtype="float32")
        names = self.pi_names
        anchor = "CT_cleavage_context" if term == "C" else "NT_cleavage_context"
        for i, n in enumerate(names):
            if n == "DB":
                w[i] = self.pi_weight_db
            if n == anchor:
                w[i] = self.pi_weight_anchor
        return (w / w.mean()).astype("float32")

    # --- channel roles ---
    @property
    def pad_channel_id(self) -> int:
        """0-based index of the padding input channel."""
        if self.channel_names and self.pad_channel_name in self.channel_names:
            return self.channel_names.index(self.pad_channel_name)
        return self.n_channels - 1

    @property
    def cont_channel_ids(self) -> tuple[int, ...]:
        """0-based indices of the continuous channels eligible for noise."""
        if not self.channel_names:
            return ()
        wanted = set(self.cont_channel_names)
        return tuple(i for i, n in enumerate(self.channel_names) if n in wanted)

    # --- serialisation ---
    def to_dict(self) -> dict:
        from dataclasses import asdict

        return asdict(self)

    @classmethod
    def from_dict(cls, d: Mapping) -> "Config":
        fields = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: _retuple(v) for k, v in d.items() if k in fields})

    def to_json(self, path) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @classmethod
    def from_json(cls, path) -> "Config":
        with open(path) as fh:
            return cls.from_dict(json.load(fh))


def _retuple(v):
    return tuple(v) if isinstance(v, list) else v
