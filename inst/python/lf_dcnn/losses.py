"""Losses and metrics for the per-index and global heads.

Index conventions (0-based; the R source states them 1-based):

* ``y_true`` has ``K_pi`` channels, the 8-class one-hot
  ``[6 real, none, padding]``.
* ``none`` sits at ``cfg.none_index`` (== ``K_cat - 1``), ``padding`` last at
  ``cfg.padding_index`` (== ``K_pi - 1``).
* The per-index LOSS masks ``none`` and keeps ``padding`` -- padding is a real
  trained class whose label is leaked by the input padding channel, so it is
  learned for free.  The per-index METRIC masks both.
"""

from __future__ import annotations

import keras
from keras import ops

EPS = 1e-7


@keras.saving.register_keras_serializable(package="lf_dcnn")
class WeightedBinaryCrossentropy(keras.losses.Loss):
    """Class-weighted BCE on the global head.

    Mirrors the R ``weighted_binary_crossentropy`` closure.  The R version does
    not clip, which is exact only while the sigmoid stays off its saturation
    points; the clip here is the standard Keras epsilon and is numerically
    identical everywhere the R form is finite.
    """

    def __init__(self, weight_0=1.0, weight_1=1.0, name="weighted_binary_crossentropy", **kw):
        super().__init__(name=name, **kw)
        self.weight_0 = float(weight_0)
        self.weight_1 = float(weight_1)

    def call(self, y_true, y_pred):
        y_true = ops.cast(y_true, y_pred.dtype)
        p = ops.clip(y_pred, EPS, 1.0 - EPS)
        return -(
            self.weight_1 * y_true * ops.log(p)
            + self.weight_0 * (1.0 - y_true) * ops.log(1.0 - p)
        )

    def get_config(self):
        return {**super().get_config(), "weight_0": self.weight_0, "weight_1": self.weight_1}


@keras.saving.register_keras_serializable(package="lf_dcnn")
class PerIndexCatLoss(keras.losses.Loss):
    """Focal categorical CE over the ``K_pi`` per-index classes, plus smoothness.

    ``none`` positions are masked out; ``padding`` positions are trained.  The
    smoothness term penalises large changes between adjacent positions of the
    predicted softmax, skipping any adjacency that touches a masked position.
    """

    def __init__(
        self,
        class_weights,
        none_index,
        gamma=2.0,
        smoothness_weight=0.01,
        name="per_index_cat_loss",
        **kw,
    ):
        super().__init__(name=name, **kw)
        self.class_weights = [float(w) for w in class_weights]
        self.none_index = int(none_index)
        self.gamma = float(gamma)
        self.smoothness_weight = float(smoothness_weight)

    def call(self, y_true, y_pred):
        y_true = ops.cast(y_true, y_pred.dtype)
        labels = y_true                                  # (batch, seq, K_pi) one-hot
        none_flag = y_true[:, :, self.none_index]        # (batch, seq): 1 at background
        keep = 1.0 - none_flag                           # train real + padding; mask none

        p = ops.clip(y_pred, EPS, 1.0 - EPS)
        cw = ops.reshape(
            ops.convert_to_tensor(self.class_weights, dtype=y_pred.dtype),
            (1, 1, len(self.class_weights)),
        )

        # focal categorical CE per position:
        #   -sum_c w_c * y_c * (1 - p_c)^gamma * log(p_c)
        ce = -ops.sum(
            cw * labels * ops.power(1.0 - p, self.gamma) * ops.log(p), axis=-1
        )                                                # (batch, seq)
        main = ops.sum(ce * keep, axis=1) / (ops.sum(keep, axis=1) + EPS)   # (batch,)

        # smoothness over adjacent positions of the predicted softmax
        pred_curr = y_pred[:, 1:, :]
        pred_prev = y_pred[:, :-1, :]
        sm_mask = keep[:, 1:] * keep[:, :-1]                               # (batch, seq-1)
        sq_pos = ops.mean(ops.square(pred_curr - pred_prev), axis=2)       # (batch, seq-1)
        smoothness = ops.sum(sq_pos * sm_mask, axis=1) / (
            ops.sum(sm_mask, axis=1) + EPS
        )                                                                  # (batch,)

        return main + self.smoothness_weight * smoothness

    def get_config(self):
        return {
            **super().get_config(),
            "class_weights": self.class_weights,
            "none_index": self.none_index,
            "gamma": self.gamma,
            "smoothness_weight": self.smoothness_weight,
        }


def make_masked_cat_accuracy(none_index: int, padding_index: int):
    """Argmax accuracy over positions that are neither ``padding`` nor ``none``.

    ``padding`` is excluded because it is trivially predictable from the input
    padding channel, and ``none`` because it is masked out of the loss.
    """

    def masked_cat_accuracy(y_true, y_pred):
        y_true = ops.cast(y_true, y_pred.dtype)
        pad_flag = y_true[:, :, padding_index]
        none_flag = y_true[:, :, none_index]
        keep = (1.0 - pad_flag) * (1.0 - none_flag)
        correct = (
            ops.cast(
                ops.equal(ops.argmax(y_pred, axis=-1), ops.argmax(y_true, axis=-1)),
                y_pred.dtype,
            )
            * keep
        )
        return ops.sum(correct) / (ops.sum(keep) + EPS)

    masked_cat_accuracy.__name__ = "masked_cat_accuracy"
    return masked_cat_accuracy
