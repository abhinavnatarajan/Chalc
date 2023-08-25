from .version import __version__
del globals()['version']

from . import filtration
from . import chromatic
from .chromatic import(
    chromatic_alpha_complex,
    chromatic_delrips_complex,
    chromatic_delcech_complex,
    delaunay_complex, 
    stratify
)