from __future__ import annotations

from ._version import __version__
import platform

from . import filtration
if platform.system() == 'Windows':
    # need to load up the GMP and MPFR DLLs
    import os as _os
    from importlib.metadata import files as _files
    _bindir = [x 
        for x in _files('chalc') 
        if len(x.parent.parts) > 0 
        and x.parent.parts[-1] == 'bin'][0].locate().parent
    with _os.add_dll_directory(_bindir):
        from . import chromatic
else:
    from . import chromatic

from . import sixpack
from . import plotting

__all__ = ['chromatic', 'filtration', 'sixpack', 'plotting']

__doc__ = "Chalc is a package to compute the persistent homology of chromatic complexes of points clouds in Euclidean space."