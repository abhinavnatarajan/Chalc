"""
Chalc is a package to compute the persistent homology of chromatic complexes of points clouds in Euclidean space.
"""
from __future__ import annotations
from . import _version
from . import chromatic
from . import filtration
__all__ = ['chromatic', 'filtration']
__version__: str = '0.2.0'
