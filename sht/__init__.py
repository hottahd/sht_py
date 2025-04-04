from .sht_py import sht_py
from .sht_ctype import sht, legendre_polynomial, legendre_transform
try:
    from ._version import version as __version__
except ImportError:
    __version__ = 'unknown'