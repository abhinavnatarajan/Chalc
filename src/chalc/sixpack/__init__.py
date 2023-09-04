from __future__ import annotations

import numpy as np
import numpy # for annotations
from enum import Enum
import functools
from .. import chromatic
from phimaker import compute_ensemble, WashboardServer

__all__ = ['get_boundary_matrix', 'compute', 'ChromaticMethod']

class Wrapper:
    def __init__(self, function):
        self.function = function
        functools.update_wrapper(self, function)
    def __call__(self,*args, **kwargs):
        return self.function(*args, **kwargs)
    def __repr__(self):
        return self.function.__repr__()
    
class ChromaticMethod(Enum):
    alpha = Wrapper(chromatic.alpha)
    delcech = Wrapper(chromatic.delcech)
    delrips = Wrapper(chromatic.delrips)
    def __call__(self, *args, **kwargs):
        return self.value(*args, **kwargs)

# domain specified by colours
def get_boundary_matrix(x : numpy.ndarray[numpy.float64[m, n]], 
                        colours : list[int], 
                        *, # rest are keyword only
                        dom : list[int] = None, 
                        k : int = None,
                        method : ChromaticMethod = ChromaticMethod.alpha, 
                        max_dgm_dim : int = 1) -> tuple[numpy.ndarray[numpy.float64[m, n]], list[float], list[int]] :
    
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
    if isinstance(method, ChromaticMethod):
        K = method(x, colours)
    else:
        raise RuntimeError("method must be an instance of ``SixPackMethod``.")

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

def compute(x : numpy.ndarray[numpy.float64[m, n]], 
            colours : list[int], 
            *, # rest are keyword only
            dom : list[int] = None, 
            k : int = None,
            method : ChromaticMethod = ChromaticMethod.alpha, 
            max_dgm_dim : int = 1) -> None :
    
    matrix, entrance_times, dimensions = get_boundary_matrix(
        x, colours,
        dom=dom, 
        k=k, 
        method=method, 
        max_dgm_dim=max_dgm_dim)
    
    truncation = truncation = max(entrance_times) * 1.01
    
    dgms =  compute_ensemble(matrix)
    
    # washboard = build_washboard_object(
    #         dgms, max_dgm_dim, truncation, dimensions, entrance_times
    #     )
    # Washboard
    washboard_server = WashboardServer.build(
        dgms,
        max_dgm_dim,
        truncation,
        dimensions,
        entrance_times,
    )
    washboard_server.open()


def _num_colours_in_bitmask(n : int) -> int :
    sum = 0
    while (n):
        sum += (n & 1)
        n >>= 1
    return sum

