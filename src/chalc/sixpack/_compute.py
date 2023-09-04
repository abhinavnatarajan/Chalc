from __future__ import annotations

import numpy # for annotations
from phimaker import compute_ensemble, WashboardServer
from ._get_boundary_matrix import get_boundary_matrix, ChromaticMethod
from ._plotting import plot_sixpack

def compute(x : numpy.ndarray[numpy.float64[m, n]], 
            colours : list[int], 
            *, # rest are keyword only
            dom : list[int] = None, 
            k : int = None,
            method : ChromaticMethod = ChromaticMethod.alpha, 
            max_dgm_dim : int = 2) -> None :
    
    matrix, entrance_times, dimensions = get_boundary_matrix(
        x, colours,
        dom=dom, 
        k=k, 
        method=method, 
        max_dgm_dim=max_dgm_dim)
    
    dgms =  compute_ensemble(matrix)
    return dgms, entrance_times, dimensions