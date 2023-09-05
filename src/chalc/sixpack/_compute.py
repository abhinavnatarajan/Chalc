from __future__ import annotations

from .. import chromatic
from ._types import DiagramEnsemble, Diagram
import numpy as np
import numpy
from collections import defaultdict
from phimaker import compute_ensemble
import typing, re

_ChromaticMethod = defaultdict(lambda : "chromatic alpha", 
                               {
                                   "chromatic alpha" : chromatic.alpha,
                                   "chromatic delcech" : chromatic.delcech,
                                   "chromatic delrips" : chromatic.delrips
                               }
)

def format_docstring(func):
    docstring = func.__doc__.replace('{', '{{').replace('}', '}}')
    docstring = re.sub(r'\${{([^{}]*)}}', r'{\1}', docstring)
    func.__doc__ = eval('f"""' + docstring + '"""')
    return func

def _get_boundary_matrix(x : numpy.ndarray[numpy.float64[m, n]], 
                        colours : list[int], 
                        *, # rest are keyword only
                        dom : list[int] = None, 
                        k : int = None,
                        method : str = _ChromaticMethod.default_factory(), 
                        max_dgm_dim : int = 2) -> tuple[numpy.ndarray[numpy.float64[m, n]], list[float], list[int]] :
    if dom is not None and k is None:
        # new colours: 1 -> domain, 0 -> codomain
        colours = np.isin(colours, dom).astype(np.int64)
        def check_in_domain(bitmask : int) -> bool:
            return bitmask == 2 # bitmask == 2 means colour 1
        
    elif k is not None and dom is None:
        def check_in_domain(bitmask : int) -> bool:
            return _num_colours_in_bitmask(bitmask) <= k # k-chromatic simplex
    
    else:
        raise RuntimeError("Only one of k or dom is allowed")

    # Compute chromatic complex
    if method in _ChromaticMethod.keys():
        K = _ChromaticMethod[method](x, colours)
    else:
        raise RuntimeError(f"method must be one of {_ChromaticMethod.keys()}")

    # Build up the matrix
    matrix = []
    entrance_times = []
    dimensions = []
    for column in K.serialised():
        facet_idxs = column[0]
        dimension = max(0, len(facet_idxs) - 1)
        # Only need 2-skeleton
        if dimension > max_dgm_dim + 1:
            break
        dimensions.append(dimension)
        entrance_times.append(column[2])
        in_domain = check_in_domain(column[3])

        matrix.append((in_domain, dimension, facet_idxs))

    return matrix, entrance_times, dimensions

def _num_colours_in_bitmask(n : int) -> int :
    sum = 0
    while (n):
        sum += (n & 1)
        n >>= 1
    return sum

@format_docstring
def compute(x : numpy.ndarray[numpy.float64[m, n]], 
            colours : list[int], 
            *, # rest are keyword only
            dom : typing.Optional[list[int]] = None, 
            k : typing.Optional[int] = None, 
            method : str = _ChromaticMethod.default_factory(), 
            max_dgm_dim : int = 2) -> DiagramEnsemble :
    """
    Compute the 6-pack of persistence diagrams of a coloured point-cloud.

    This function constructs a chromatic complex :math:`L` from the point cloud, 
    and computes persistence diagrams associated with the inclusion 
    :math:`f : K \hookrightarrow L` of a filtered subcomplex into :math:`L`.

    Parameters
    ----------
    x : 
        Numpy matrix whose columns are points.
    colours :
        List of integers describing the colours of the points.

    Keyword Args
    ------------
    dom : 
        List of integers describing the colours of the points in the domain.
    k :
        If not ``None``, then the domain is taken to be the :math:`k`-chromatic 
        subcomplex of :math:`L`, i.e., the subcomplex of simplices of at 
        most :math:`k` colours.
    method:
        Filtration used to construct the chromatic complex. Must be one of 
        ``${str(list(_ChromaticMethod.keys()))}``.
    max_dgm_dim :
        Maximum homological dimension for which the persistence diagrams are computed. 
    
    Returns
    -------
    dgms : DiagramEnsemble
        Diagrams corresponding to the following persistence modules (where 
        :math:`H_*` is the persistent homology functor and :math:`f_*` is the 
        induced map on persistent homology):

        #. :math:`H_*(K)` (domain)
        #. :math:`H_*(L)` (codomain)
        #. :math:`\ker(f_*)` (kernel)
        #. :math:`\mathrm{coker}(f_*)` (cokernel)
        #. :math:`\mathrm{im}(f_*)` (image)
        #. :math:`H_*(L, K)` (relative homology)

        Each diagram is represented by sets of paired and unpaired simplices, 
        and contain simplices of all dimensions. ``dgms`` also contains the 
        entrance times of the simplices and their dimensions.

    """
    
    matrix, entrance_times, dimensions = _get_boundary_matrix(
        x, colours,
        dom=dom, 
        k=k, 
        method=method, 
        max_dgm_dim=max_dgm_dim)
    d =  compute_ensemble(matrix)
    dgms = DiagramEnsemble(
        Diagram(d.ker), 
        Diagram(d.cok), 
        Diagram(d.f), 
        Diagram(d.g), 
        Diagram(d.im), 
        Diagram(d.rel), 
        entrance_times, 
        dimensions)
    return dgms