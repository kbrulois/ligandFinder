"""lf_dcnn -- the 1D-CNN that scores 36-residue peptide-cleavage windows.

Two models (terminus ``N`` and ``C``) each carry a global ranking head and a
supervised per-index class head.  Ported from
``inst/scripts/10_1dcnn_new6.R``; see ``README.md`` for the R <-> Python
boundary and the two documented-intent switches on :class:`Config`.
"""

from . import io, synthetic
from .config import CLASS_NAMES, DEFAULT_CHANNEL_NAMES, Config
from .calibrate import PlattCalibrator, calibrate
from .data import OversampledWindows, TermArrays, build_dataset
from .model import build_model, compile_model, embed_model
from .pipeline import Result, embed, predict_all, run, set_seed, train

__all__ = [
    "CLASS_NAMES",
    "DEFAULT_CHANNEL_NAMES",
    "Config",
    "OversampledWindows",
    "PlattCalibrator",
    "Result",
    "TermArrays",
    "build_dataset",
    "build_model",
    "calibrate",
    "compile_model",
    "embed",
    "embed_model",
    "io",
    "predict_all",
    "run",
    "set_seed",
    "synthetic",
    "train",
]

__version__ = "0.1.0"
