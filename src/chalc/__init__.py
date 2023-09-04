from __future__ import annotations

from ._version import __version__

import os as _os
from importlib.metadata import files as _files
_bindir = [x 
     for x in _files('chalc') 
     if len(x.parent.parts) > 0 
     and x.parent.parts[-1] == 'bin'][0].locate().parent
with _os.add_dll_directory(_bindir):
    from . import filtration
    from . import chromatic
    from . import sixpack

__doc__ = "Chalc is a package to compute the persistent homology of chromatic complexes of points clouds in Euclidean space."