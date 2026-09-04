"""Pooled Platt calibration of the global ranking score.

Per-model calibration is unstable here (a handful of validation positives, so
near-perfect separation), and both models feed ONE shared global ranking.  So a
single logistic calibrator is fit on the validation scores and labels POOLED
across the models, then applied to every raw score.  The map is monotone, so it
leaves each model's ROC/PR-AUC unchanged; it only puts the two score scales on a
common footing.

scikit-learn is not available in the target environment, so the fit is a
hand-rolled IRLS that mirrors R's ``glm(family = binomial)`` -- same
initialisation, same convergence test, same 25-iteration cap -- with a
``scipy.optimize`` path available as a cross-check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: R's glm.control() defaults
GLM_EPSILON = 1e-8
GLM_MAXIT = 25


def _sigmoid(z: np.ndarray) -> np.ndarray:
    out = np.empty_like(z, dtype="float64")
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


@dataclass
class PlattCalibrator:
    """Univariate logistic calibrator: ``P(label) = sigmoid(a + b * score)``."""

    intercept: float = 0.0
    slope: float = 1.0
    n_obs: int = 0
    n_pos: int = 0
    iterations: int = 0
    converged: bool = False
    deviance: float = float("nan")
    method: str = "irls"

    # --- fitting ------------------------------------------------------------
    @classmethod
    def fit(
        cls,
        scores,
        labels,
        method: str = "irls",
        maxit: int = GLM_MAXIT,
        epsilon: float = GLM_EPSILON,
    ) -> "PlattCalibrator":
        s = np.asarray(scores, dtype="float64").ravel()
        y = np.asarray(labels, dtype="float64").ravel()
        if s.shape != y.shape:
            raise ValueError(f"scores {s.shape} and labels {y.shape} must match")
        if s.size == 0:
            raise ValueError("no observations to calibrate on")
        if len(np.unique(y)) < 2:
            raise ValueError(
                "pooled validation labels are all one class; cannot fit a calibrator"
            )

        if method == "irls":
            beta, it, converged, dev = _irls_logistic(s, y, maxit=maxit, epsilon=epsilon)
        elif method == "lbfgs":
            beta, it, converged, dev = _scipy_logistic(s, y)
        else:
            raise ValueError(f"unknown method {method!r} (use 'irls' or 'lbfgs')")

        return cls(
            intercept=float(beta[0]),
            slope=float(beta[1]),
            n_obs=int(s.size),
            n_pos=int(y.sum()),
            iterations=int(it),
            converged=bool(converged),
            deviance=float(dev),
            method=method,
        )

    # --- applying -----------------------------------------------------------
    def predict(self, scores) -> np.ndarray:
        s = np.asarray(scores, dtype="float64").ravel()
        return _sigmoid(self.intercept + self.slope * s)

    __call__ = predict

    @property
    def coef(self) -> np.ndarray:
        return np.array([self.intercept, self.slope], dtype="float64")

    def to_dict(self) -> dict:
        return {
            "intercept": self.intercept,
            "slope": self.slope,
            "n_obs": self.n_obs,
            "n_pos": self.n_pos,
            "iterations": self.iterations,
            "converged": self.converged,
            "deviance": self.deviance,
            "method": self.method,
        }


def _design(s: np.ndarray) -> np.ndarray:
    return np.column_stack([np.ones_like(s), s])


def _binomial_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    eps = np.finfo("float64").eps
    mu = np.clip(mu, eps, 1.0 - eps)
    return float(-2.0 * np.sum(y * np.log(mu) + (1.0 - y) * np.log(1.0 - mu)))


def _irls_logistic(s, y, maxit=GLM_MAXIT, epsilon=GLM_EPSILON):
    """IRLS exactly as R's ``glm.fit`` runs it for a binomial/logit model.

    Perfect separation makes the MLE diverge; R stops at ``maxit`` with a
    warning and returns the last iterate, and so does this.  The resulting map
    is still monotone, which is all the pooled ranking needs.
    """
    X = _design(s)
    # R's binomial()$initialize: mustart <- (weights * y + 0.5) / (weights + 1)
    mu = (y + 0.5) / 2.0
    eta = np.log(mu / (1.0 - mu))
    beta = np.zeros(X.shape[1])
    dev = _binomial_deviance(y, mu)
    converged = False
    it = 0

    for it in range(1, maxit + 1):
        mu_eta = mu * (1.0 - mu)                 # d(mu)/d(eta) == variance here
        mu_eta = np.maximum(mu_eta, np.finfo("float64").eps)
        z = eta + (y - mu) / mu_eta              # working response
        w = np.sqrt(mu_eta)                      # working weights
        beta, *_ = np.linalg.lstsq(X * w[:, None], z * w, rcond=None)
        eta = X @ beta
        mu = _sigmoid(eta)
        dev_new = _binomial_deviance(y, mu)
        if abs(dev_new - dev) / (abs(dev_new) + 0.1) < epsilon:
            dev = dev_new
            converged = True
            break
        dev = dev_new

    return beta, it, converged, dev


def _scipy_logistic(s, y):
    """Direct MLE via ``scipy.optimize`` -- the cross-check on the IRLS path."""
    from scipy.optimize import minimize

    X = _design(s)

    def nll(beta):
        eta = X @ beta
        # log(1 + exp(eta)) computed stably
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta))

    def grad(beta):
        return X.T @ (_sigmoid(X @ beta) - y)

    res = minimize(nll, np.zeros(X.shape[1]), jac=grad, method="L-BFGS-B")
    mu = _sigmoid(X @ res.x)
    return res.x, int(res.nit), bool(res.success), _binomial_deviance(y, mu)


def calibrate(val_scores, val_labels, raw_scores, method: str = "irls"):
    """Fit on pooled validation, apply to ``raw_scores``.

    Returns ``(calibrated, calibrator)``.
    """
    cal = PlattCalibrator.fit(val_scores, val_labels, method=method)
    return cal.predict(raw_scores), cal
