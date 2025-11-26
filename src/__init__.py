<<<<<<< HEAD
from .polynomial import Polynomial
=======
from .polynomial import Polynomial
try:
    from . import polylib_ext as _ext
except Exception:
    _ext = None

__all__ = ["Polynomial"]
>>>>>>> ab47c69 (Implemented in C)
