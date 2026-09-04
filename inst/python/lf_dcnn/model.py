"""The two-headed 1D-CNN over 36-residue peptide-cleavage windows."""

from __future__ import annotations

import keras
import numpy as np
from keras import layers, ops

from .config import Config
from .losses import PerIndexCatLoss, WeightedBinaryCrossentropy, make_masked_cat_accuracy


@keras.saving.register_keras_serializable(package="lf_dcnn")
class PositionRamp(layers.Layer):
    """Emit a normalised position ramp ``[0, 1]`` shaped ``(batch, seq, 1)``.

    The windows are anchored and aligned (fixed cleavage/DB layout), so absolute
    position is meaningful -- but conv is translation-equivariant and the
    attention/pooling above it are order-agnostic.  Appending the ramp as an
    extra input channel lets the trunk learn where along the window it is.
    Built in-graph so the raw inputs (and the padding-channel logic that keys
    off them) stay untouched.
    """

    def __init__(self, seq_len: int, **kw):
        super().__init__(**kw)
        self.seq_len = int(seq_len)

    def call(self, x):
        ones = ops.sum(x, axis=-1, keepdims=True) * 0.0 + 1.0   # (batch, seq, 1)
        idx = ops.cumsum(ones, axis=1) - 1.0                    # 0 .. seq-1
        return idx / float(self.seq_len - 1)

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.seq_len, 1)

    def get_config(self):
        return {**super().get_config(), "seq_len": self.seq_len}


@keras.saving.register_keras_serializable(package="lf_dcnn")
class ClassPositionMask(layers.Layer):
    """Drop the trailing padding logit, then zero class/position combinations
    that cannot occur.

    Fuses the two R lambdas (``op_take`` of the first ``K_cat`` columns, then a
    multiply by ``mask_matrix_cat``).  Note this multiplies LOGITS by 0/1, so a
    masked class is neutralised rather than driven to zero probability -- the
    same soft behaviour the R model has.
    """

    def __init__(self, mask_matrix, **kw):
        super().__init__(**kw)
        self.mask_matrix = np.asarray(mask_matrix, dtype="float32")

    def call(self, x):
        # A graph constant, not a weight: it must not show up in count_params().
        mask = ops.convert_to_tensor(self.mask_matrix, dtype=x.dtype)
        return x[..., : self.mask_matrix.shape[1]] * mask

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1], self.mask_matrix.shape[1])

    def get_config(self):
        return {**super().get_config(), "mask_matrix": self.mask_matrix.tolist()}


def build_model(cfg: Config | None = None, clear_session: bool = True) -> keras.Model:
    """Build the functional model.

    Outputs are a dict: ``global`` (scalar sigmoid ranking score) and
    ``per_index_cat`` (per-position softmax over ``K_pi`` classes).

    ``clear_session`` mirrors ``keras3::clear_session()`` in the R loop; without
    it the two per-terminus models share graph state and auto-generated layer
    names drift between runs.
    """
    cfg = cfg or Config()
    if clear_session:
        keras.backend.clear_session()

    inputs = layers.Input(shape=(cfg.seq_len, cfg.n_channels), name="window")
    reg = keras.regularizers.l2(cfg.l2)

    pos_channel = PositionRamp(cfg.seq_len, name="pos_ramp")(inputs)
    conv_in = layers.Concatenate(name="conv_in")([inputs, pos_channel])

    # Slim trunk for the tiny (<100 example) training sets: two small conv
    # layers rather than four wide ones, which held ~26k of the old ~33k params.
    x = conv_in
    for i, (filters, drop) in enumerate(zip(cfg.conv_filters, cfg.conv_dropout)):
        x = layers.Conv1D(
            filters,
            kernel_size=cfg.conv_kernel,
            activation="gelu",
            padding="same",
            kernel_regularizer=reg,
            name=f"conv{i + 1}",
        )(x)
        x = layers.Dropout(drop, name=f"drop{i + 1}")(x)
    shared = x

    per_index_logits = layers.Conv1D(
        filters=cfg.K_pi, kernel_size=1, kernel_regularizer=reg, name="per_index_logits"
    )(shared)

    # Per-index head: supervised softmax over [6 real, none, padding].
    per_index_cat = layers.Activation("softmax", name="per_index_cat")(per_index_logits)

    # Global head input: masked softmax over the [6 real + none] logits, so the
    # ranking head's class view includes background/none but not padding.
    masked_sum = ClassPositionMask(cfg.mask_matrix_cat, name="class_position_mask")(
        per_index_logits
    )
    masked_sum = layers.Activation("softmax", name="masked_sum")(masked_sum)

    # A SMALL multi-head attention over the masked class softmax, then pool. A
    # plain GAP of the softmax washed out all positional signal and tanked
    # global AUC; the attention restores it at ~1/6 the cost of the old
    # 4-head/key_dim-32 version. The dense named "embed" is the representation
    # reused (no new loss) for nearest-known retrieval.
    attn = layers.MultiHeadAttention(
        num_heads=cfg.attention_heads, key_dim=cfg.attention_key_dim, name="attn"
    )(masked_sum, masked_sum)
    g = layers.LayerNormalization(name="attn_norm")(attn)
    g = layers.GlobalAveragePooling1D(name="gap")(g)
    g = layers.Dense(cfg.embed_units, activation="gelu", name="embed")(g)
    g = layers.Dropout(cfg.embed_dropout, name="embed_drop")(g)
    global_output = layers.Dense(1, activation="sigmoid", name="global")(g)

    return keras.Model(
        inputs=inputs,
        outputs={"global": global_output, "per_index_cat": per_index_cat},
        name="lf_dcnn",
    )


def compile_model(model: keras.Model, cfg: Config, term: str) -> keras.Model:
    """Attach the optimiser, the two losses and the metrics."""
    model.compile(
        # The old cosine schedule barely moved over the ~600-1200 total steps
        # here (LR stayed ~0.98x initial), so it was a constant LR anyway.
        optimizer=keras.optimizers.Adam(
            learning_rate=cfg.learning_rate, clipnorm=cfg.clipnorm
        ),
        loss={
            "global": WeightedBinaryCrossentropy(
                weight_0=cfg.bce_weight_0, weight_1=cfg.bce_weight_1
            ),
            "per_index_cat": PerIndexCatLoss(
                class_weights=cfg.pi_weights(term),
                none_index=cfg.none_index,
                gamma=cfg.gamma,
                smoothness_weight=cfg.smoothness_weight,
            ),
        },
        loss_weights={
            "global": cfg.loss_weight_global,
            "per_index_cat": cfg.loss_weight_per_index,
        },
        metrics={
            "global": [
                keras.metrics.AUC(name="auc"),
                keras.metrics.AUC(name="pr_auc", curve="PR"),
            ],
            "per_index_cat": [
                make_masked_cat_accuracy(cfg.none_index, cfg.padding_index)
            ],
        },
    )
    return model


def embed_model(model: keras.Model) -> keras.Model:
    """A model that returns the 16-d ``embed`` representation."""
    return keras.Model(model.input, model.get_layer("embed").output)
