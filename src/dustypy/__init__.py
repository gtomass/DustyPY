__version__ = "1.0.0"
__author__ = "G. Tomassini"
__all__ = ["Model", "Runner", "Result", "Dataset", "Fitter", "Filter"]

import os
for var in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
    if var not in os.environ:
        os.environ[var] = "1"

from .core.model import Model
from .core.runner import Runner
from .core.result import Result
from .core.dataset import Dataset
from .core.fitter import Fitter
from .core.filter import Filter