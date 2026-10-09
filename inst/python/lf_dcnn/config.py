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
# Order is load-bearing: the one-hot columns of `y_per_index_cat` are derived
# from it.
#
# `DB` and `gap` were removed once the peptide-end set became the production
# one. Both were defined RELATIVE TO THE DIBASIC SLOT that set does not anchor
# on -- `DB` was the pair itself and `gap` the space between it and the peptide
# -- so neither has a referent here: whatever precedes the peptide is NT
# context and whatever follows it CT context, dibasic or not.
CLASS_NAMES: tuple[str, ...] = (
    "CT_cleavage_context",
    "NT_cleavage_context",
    "pep_other",
    "pep_pocket",
    "padding",
    "none",
)

#: class vocabulary of the peptide-END set (inst/scripts/10_7_pepend_windows.R).
#: No ``DB`` and no ``gap``: both were defined relative to the dibasic slot that
#: set does not anchor on, so whatever precedes the peptide is NT context and
#: whatever follows it CT context, dibasic or not.  Mirrors
#: ``LF_PEPEND_CLASSES`` in R/pepend_windows.R; order is load-bearing, as for
#: :data:`CLASS_NAMES`.
PEPEND_CLASS_NAMES: tuple[str, ...] = (
    "CT_cleavage_context",
    "NT_cleavage_context",
    "pep_other",
    "pep_pocket",
    "padding",
    "none",
)

#: the insertion classes of the peptide-end head, in softmax column order.
#: Decided by where the `pep_pocket` residues sit relative to the end a window
#: scores -- see lf_pepend_insertion_class() in R/pepend_windows.R.
PEPEND_INS_CLASSES: tuple[str, ...] = ("inserting", "loop", "non_inserting")

#: where the peptide's own terminal residue sits in a peptide-end window,
#: 1-based (``LF_PEPEND$anchor`` in R/pepend_windows.R).
PEPEND_ANCHOR: dict[str, int] = {"N": 8, "C": 28}

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

    :meth:`r_exact` is the flat trunk with the position ramp. It NO LONGER
    reproduces ``10_1dcnn_new6.R``: that model's per-index head carried the
    ``DB`` and ``gap`` classes, and those went with the dibasic-anchored set, so
    its 2,438 parameters cannot be rebuilt from this vocabulary. The preset
    survives as an architecture choice, not as a reproduction.
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
    #: ``"unet"`` (default since 2026-09-18, see 10_5_benchmark_window_model.R)
    #: pools down and upsamples back, widening the channel count on the way
    #: down; ``"flat"`` keeps full 36-position resolution through the whole
    #: trunk (see :meth:`r_exact`).  See
    #: :func:`lf_dcnn.model.build_model`.
    trunk: str = "unet"
    #: where the global (window-score) head reads from.
    #: ``"attn"``       the masked class softmax -> attention -> pool (default)
    #: ``"gap"``        the same masked class softmax, average-pooled with NO
    #:                  attention. The ablation of ``attn``: identical input to
    #:                  the embedding, identical width, the pooling alone
    #:                  differs, and it carries no parameters of its own.
    #: ``"bottleneck"`` the U-Net bottleneck -> pool. Gives the score direct
    #:                  access to the pooled whole-window representation instead
    #:                  of forcing it through the 7-channel class softmax.
    #: ``"both"``       concatenate the two pooled vectors.
    global_head: str = "attn"
    #: a second window-level head: a 3-way softmax over the insertion classes
    #: (see :data:`PEPEND_INS_CLASSES`), reading the same pooled embedding the
    #: window score does.
    #:
    #: Off by default, because the label exists only for the peptide-end set and
    #: only for its KNOWN windows. A negative is not a peptide end, so it has no
    #: insertion class at all -- it carries an all-zero row and is masked out of
    #: this head's loss, exactly as an unlabelled position is masked out of the
    #: per-index loss. See :class:`lf_dcnn.losses.MaskedWindowCatLoss`.
    ins_head: bool = False
    ins_class_names: tuple[str, ...] = PEPEND_INS_CLASSES
    #: give the per-index softmax an explicit ``none`` column.
    #:
    #: ``none`` is never a training target either way -- positions labelled
    #: ``none`` are masked out of the per-index loss (see
    #: :class:`lf_dcnn.losses.PerIndexCatLoss`). The column only decides whether
    #: the model has somewhere to put probability mass at a position it is not
    #: scored on, and whether the global head's attention sees that column.
    #: ``False`` drops it: ``K_cat`` 7 -> 6, ``K_pi`` 8 -> 7, and a ``none``
    #: position becomes an all-zero label row, which is what the loss then masks
    #: on. See inst/scripts/10_5_benchmark_window_model.R --preset none_class.
    include_none: bool = True
    #: train on the ``none`` positions instead of masking them out of the
    #: per-index loss.
    #:
    #: ON by default since 2026-09-25. At 20 seeds on the C terminus
    #: (10_5_benchmark_window_model.R --preset none20) it recovers more
    #: held-out known peptide ends (median candidate rank 204 vs 318, 37 vs 31
    #: of 52 in the top 500) and collapses ensemble disagreement (member sd
    #: .05 vs .28), at 0.716 vs 0.845 per-residue real-class accuracy.
    #:
    #: The caveat the numbers do not carry: the evidence is C-TERMINUS ONLY.
    #: The 20-seed arms ran at pi_weight_none = 1.0, the default since
    #: 2026-09-25 (see that field).
    #:
    #: What masking costs: a NEGATIVE window is all ``none``, so it contributes
    #: nothing at all to the per-index loss, and the ``none`` column never
    #: receives a positive gradient. Measured on the trained default it is
    #: effectively dead -- mean probability 0.007 at scored positions, and
    #: never the argmax at any of 1,065,564 positions. Training on it saturates
    #: validation PR AUC, recovers more held-out known peptide ends and
    #: tightens the ensemble, at ~5 points of per-residue real-class accuracy.
    #:
    #: ``none`` is 96% of the raw training positions (75% after the 1:3 window
    #: oversampling), so it needs ``pi_weight_none`` to not swamp the six real
    #: classes.  Requires ``include_none``.
    none_in_loss: bool = True
    #: append the in-graph ``[0, 1]`` position ramp as an extra input channel
    #: (see :class:`lf_dcnn.model.PositionRamp`). Off by default since
    #: 2026-09-18: the trunk gets the raw channels only. The original R model
    #: had it on; :meth:`r_exact` still turns it on.
    position_ramp: bool = False
    conv_filters: tuple[int, ...] = (16, 8)
    conv_kernel: int = 3
    conv_dropout: tuple[float, ...] = (0.2, 0.3)
    #: U-Net only: channels per level, last entry is the bottleneck.
    #: ``(16, 32, 64)`` means 36->18->9 with 16, 32 then 64 filters.
    unet_filters: tuple[int, ...] = (16, 32, 64)
    unet_dropout: tuple[float, ...] = (0.2, 0.2, 0.2)
    pool_size: int = 2
    attention_heads: int = 2
    attention_key_dim: int = 8
    embed_units: int = 16
    embed_dropout: float = 0.3

    # --- losses ---
    gamma: float = 2.0
    smoothness_weight: float = 0.01
    loss_weight_ins: float = 1.0
    loss_weight_global: float = 0.05
    loss_weight_per_index: float = 1.0
    # The oversampler already rebalances each batch to ~1:3, so ALSO upweighting
    # positives here would double-correct the imbalance (the old ~20x weight).
    bce_weight_0: float = 1.0
    bce_weight_1: float = 1.0
    # per-index class weights: everything 1, DB doubled, and the anchor-side
    # cleavage context tripled (CT for the C model, NT for the N model).
    pi_weight_anchor: float = 3.0
    #: weight of the ``none`` class when ``none_in_loss``, relative to an
    #: ordinary (weight-1) class. The oversampler already brings the in-batch
    #: ``none``:real ratio down to ~3.2:1, so ~0.31 equalises their TOTAL loss
    #: mass; the focal term (``gamma``) suppresses easy ``none`` positions
    #: further. 1.0 lets ``none`` dominate, 0.0 reproduces masking it out.
    #:
    #: 1.0 by default since 2026-09-25 -- the weight the 20-seed evidence for
    #: ``none_in_loss`` was gathered at (--preset none20, the ``unet_none_w1``
    #: arm), and better than 0.1 on every aggregate there. It does not collapse
    #: onto ``none``: it takes 97.7% of background argmaxes but 3.8% of scored
    #: ones, so the oversampler and the focal term, not this weight, are what
    #: hold the majority back. It was 0.1 before, chosen on a 5-seed
    #: 0.1/0.3/0.6 sweep (--preset none_weight) where every weight was within
    #: noise; 0.1 costs less real-class accuracy (0.766 vs 0.716 at 0.6). The
    #: weights rank CANDIDATES quite differently despite matching on every
    #: summary metric (Spearman 0.27 between w=0.1 and w=0.6).
    pi_weight_none: float = 1.0

    # --- optimisation ---
    learning_rate: float = 3e-4
    clipnorm: float = 1.0
    epochs: int = 2000
    batch_size: int = 32
    n_negatives_per_positive: int = 3
    noise_frac: float = 0.1
    #: when set, each term writes TensorBoard logs under ``<dir>/<term>`` so
    #: training curves can be watched live (the Python stand-in for keras3's
    #: RStudio viewer, which only engages when R drives fit()).
    tensorboard_dir: str | None = None
    monitor: str = "val_global_pr_auc"
    patience: int = 300
    start_from_epoch: int = 300

    # --- channel roles ---
    pad_channel_name: str = "padding"
    cont_channel_names: tuple[str, ...] = DEFAULT_CONT_CHANNELS

    # --- documented-intent switches (see class docstring) ---
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
        "n_negatives_per_positive", "patience", "start_from_epoch", "pool_size",
    )

    def __post_init__(self) -> None:
        # reticulate hands R numerics over as floats; coerce so Keras gets ints
        for f in self._INT_FIELDS:
            object.__setattr__(self, f, int(getattr(self, f)))
        if self.seed is not None:
            object.__setattr__(self, "seed", int(self.seed))
        object.__setattr__(self, "position_ramp", bool(self.position_ramp))
        object.__setattr__(self, "include_none", bool(self.include_none))
        object.__setattr__(self, "none_in_loss", bool(self.none_in_loss))
        if self.none_in_loss and not self.include_none:
            raise ValueError("none_in_loss needs include_none: there is no `none` column to train")
        for f in ("conv_dropout", "unet_dropout", "nt_span", "ct_span", "mid_span"):
            object.__setattr__(self, f, tuple(getattr(self, f)))
        for f in ("conv_filters", "unet_filters"):
            object.__setattr__(self, f, tuple(int(v) for v in getattr(self, f)))
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
        if self.trunk not in ("flat", "unet"):
            raise ValueError(f"trunk must be 'flat' or 'unet', got {self.trunk!r}")
        if self.global_head not in ("attn", "gap", "bottleneck", "both"):
            raise ValueError(
                "global_head must be 'attn', 'gap', 'bottleneck' or 'both', got "
                f"{self.global_head!r}")
        ## `gap` reads the class softmax like `attn` does, so it needs no
        ## bottleneck and works on either trunk
        if self.global_head in ("bottleneck", "both") and self.trunk != "unet":
            raise ValueError(
                f"global_head={self.global_head!r} needs a bottleneck; it is only "
                "available with trunk='unet'")
        if self.trunk == "unet":
            # every pooling step must divide the sequence exactly, or the
            # upsampled decoder will not line back up with its skip connection
            n_levels = len(self.unet_filters) - 1
            if n_levels < 1:
                raise ValueError("unet_filters needs at least 2 entries")
            if len(self.unet_dropout) != len(self.unet_filters):
                raise ValueError("unet_dropout must match unet_filters in length")
            step = self.seq_len
            for lvl in range(n_levels):
                if step % self.pool_size:
                    raise ValueError(
                        f"seq_len={self.seq_len} is not divisible by pool_size="
                        f"{self.pool_size} at level {lvl + 1} (length {step}); "
                        "reduce the number of unet_filters levels or pad the window"
                    )
                step //= self.pool_size

    @classmethod
    def r_exact(cls, **kwargs) -> "Config":
        """The flat trunk with the position ramp and ``none`` masked out.

        It is named for the pre-port R model of ``10_1dcnn_new6.R`` and no
        longer reproduces it: that model's per-index head had ``DB`` and
        ``gap``, which were removed with the dibasic-anchored set. What is left
        is the architecture, not the reproduction.
        """
        kwargs.setdefault("trunk", "flat")
        kwargs.setdefault("position_ramp", True)
        kwargs.setdefault("none_in_loss", False)      # the R model masked `none`
        kwargs.setdefault("resample_each_epoch", False)
        return cls(**kwargs)

    @property
    def K_ins(self) -> int:
        """Width of the insertion-class softmax."""
        return len(self.ins_class_names)

    @classmethod
    def pepend(cls, term: str, **kwargs) -> "Config":
        """Config for the peptide-END set of ``inst/scripts/10_7_pepend_windows.R``.

        That set anchors on the peptide's own terminal residue -- its last
        residue at position 28 for a C window, its first at position 8 for an N
        window -- rather than on a dibasic site, and drops ``DB`` and ``gap``
        (see :data:`PEPEND_CLASS_NAMES`).

        The stock ``nt_span``/``ct_span``/``mid_span`` encode the dibasic layout
        (NT 1-5, CT 31-36) and are simply WRONG for this one, so they are
        re-derived from the anchor to match ``lf_pepend_labels``: everything
        strictly beyond the peptide terminus is that side's context, the peptide
        itself runs to the window edge, and the near-side context fills only
        what a peptide short enough to end inside the window leaves over.

        The spans are term-specific but :attr:`mask_matrix_cat` is not, so this
        pins ``term_order`` to the single terminus asked for; build a second
        Config for the other one rather than training both from this.
        """
        if term not in PEPEND_ANCHOR:
            raise ValueError(
                f"term must be one of {sorted(PEPEND_ANCHOR)}, got {term!r}")
        a = PEPEND_ANCHOR[term]
        seq_len = int(kwargs.get("seq_len", cls.seq_len))
        if term == "C":
            spans = dict(mid_span=(1, a), ct_span=(a + 1, seq_len), nt_span=(1, a - 1))
        else:
            spans = dict(mid_span=(a, seq_len), nt_span=(1, a - 1), ct_span=(a + 1, seq_len))
        for k, v in spans.items():
            kwargs.setdefault(k, v)
        kwargs.setdefault("class_names", PEPEND_CLASS_NAMES)
        kwargs.setdefault("term_order", (term,))
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
        """The 6 real classes, plus an explicit ``none`` last when included."""
        return self.real_cols + ((self.none_col,) if self.include_none else ())

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
        """Index of ``none`` within the per-index softmax (== ``K_cat - 1``),
        or ``-1`` when there is no ``none`` column.

        ``-1`` tells the loss and the metric to mask on an all-zero label row
        instead of on a ``none`` one-hot; both say "this position is background,
        do not score it".
        """
        return self.K_cat - 1 if self.include_none else -1

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
            "pep_other": self.mid_mask,
            "pep_pocket": self.mid_mask,
            "none": np.ones(self.seq_len, dtype="float32"),
        }
        ## A class with no rule used to fall back to an all-zero column, which
        ## silently forbids it everywhere. Be loud instead: that is how a stale
        ## vocabulary (one still carrying DB or gap) would otherwise get through.
        missing = [self.class_names[c] for c in self.cat_cols if self.class_names[c] not in by_class]
        if missing:
            raise ValueError(
                f"no position rule for class(es) {missing}. The vocabulary is "
                f"{list(self.class_names)}; DB and gap were removed with the "
                "dibasic-anchored set.")
        cols = [by_class[self.class_names[c]] for c in self.cat_cols]
        return np.stack(cols, axis=1).astype("float32")

    # --- per-index class weights ---
    def pi_weights(self, term: str) -> np.ndarray:
        """Per-model per-index class weights, normalised to mean 1.

        A hand-set prior, not learned: emphasise the anchor-side cleavage
        context (CT for the C model, NT for the N model) plus a moderate DB
        boost.  Mean-1 keeps the per-index loss magnitude stable.

        ``pi_weight_none`` is applied AFTER that normalisation, so turning
        ``none_in_loss`` on scales only the ``none`` entry and leaves every
        real class at exactly the weight a masked model trains with -- the two
        differ by the `none` positions alone, not by a shifted loss magnitude.
        The returned vector's mean is then not 1 (below it for any
        ``pi_weight_none`` < 1); that is the intended trade, since renormalising
        afterwards would push every real class up and reintroduce exactly the
        confound this avoids.
        """
        w = np.ones(self.K_pi, dtype="float32")
        names = self.pi_names
        anchor = "CT_cleavage_context" if term == "C" else "NT_cleavage_context"
        for i, n in enumerate(names):
            if n == anchor:
                w[i] = self.pi_weight_anchor
        w = (w / w.mean()).astype("float32")
        if self.include_none and self.none_in_loss:
            w[self.none_index] *= float(self.pi_weight_none)
        return w

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
