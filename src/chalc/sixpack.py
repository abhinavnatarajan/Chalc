from __future__ import annotations

__doc__ = "Routines for computing six-packs of persistence diagrams."

from . import chromatic
import numpy as np
from collections import defaultdict
from phimaker import compute_ensemble
from dataclasses import dataclass, fields, field
from numbers import Real
import typing, re, h5py

_ChromaticMethod = defaultdict(lambda : "chromatic alpha", 
                               {
                                   "chromatic alpha" : chromatic.alpha,
                                   "chromatic delcech" : chromatic.delcech,
                                   "chromatic delrips" : chromatic.delrips
                               }
)

def _format_docstring(func):
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
        new_colours = np.isin(colours, dom).astype(np.int64)
        check_in_domain = lambda b : b == 2
        
    elif k is not None and dom is None:
        check_in_domain = lambda b: _num_colours_in_bitmask(b) <= k # k-chromatic simplex
    
    else:
        raise RuntimeError("Only one of k or dom is allowed")

    # Compute chromatic complex
    if method in _ChromaticMethod.keys():
        K = _ChromaticMethod[method](x, new_colours)
    else:
        raise RuntimeError(f"method must be one of {_ChromaticMethod.keys()}")

    # Build up the matrix
    matrix = []
    entrance_times = []
    dimensions = []
    for column in K.serialised():
        facet_idxs = column[0]
        dimension = max(0, len(facet_idxs) - 1)
        # Only need k-skeleton for k <= max_dgm_dim + 1
        if dimension > max_dgm_dim + 1:
            break
        dimensions.append(dimension)
        entrance_times.append(column[2])
        in_domain = check_in_domain(column[3])

        matrix.append((in_domain, dimension, facet_idxs))

    return matrix, entrance_times, dimensions

def _colours_to_bitmask(colours : list[int]) -> int :
    return sum(2**i for i in colours)

def _bitmask_to_colours(b : int) -> list[int] :
    i = 0
    res = []
    while (b) :
        if (b & 1):
            res.append(i)            
        i += 1
        b >>= 1
    return res

def _num_colours_in_bitmask(b : int) -> int :
    return len(_bitmask_to_colours(b))

# check if the colours specified by bitmask b1 are a subset of those specified by b2
def _colours_are_subset(b1 : int, b2 : int) -> bool :
    return not (~b2 & b1)

@_format_docstring
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
    d = compute_ensemble(matrix)
    dgms = DiagramEnsemble(
        Diagram._fromPhimaker(d.ker), 
        Diagram._fromPhimaker(d.cok), 
        Diagram._fromPhimaker(d.g), 
        Diagram._fromPhimaker(d.f), 
        Diagram._fromPhimaker(d.im), 
        Diagram._fromPhimaker(d.rel), 
        entrance_times, 
        dimensions)
    return dgms


@dataclass
class Diagram:
    """ 
    Persistence diagram object, represented by a 
    set of simplex pairings and a set of unpaired simplices.
    """
    #: Set of tuples of indices of paired simplices. 
    paired: set[tuple[int, int]] = field(default_factory=lambda:{})
    #: Set of indices of unpaired simplices.
    unpaired: set[int] = field(default_factory=lambda:{})

    @classmethod
    def _fromPhimaker(cls, obj):
        return cls(obj.paired, obj.unpaired)
    
    @classmethod
    def _from_matrices(cls, p, u):
        paired = set((p[i, 0], p[i, 1]) for i in range(p.shape[0]))
        unpaired = set(u)
        return cls(paired, unpaired)

    def paired_as_matrix(self):
        return np.concatenate(
            list(self.paired)).reshape((len(self.paired), 2))
    


@dataclass
class DiagramEnsemble:
    """
    Six pack of persistence diagrams.
    """
    #: Kernel
    ker: Diagram = field(default_factory = lambda : Diagram())
    #: Cokernel
    cok: Diagram = field(default_factory = lambda : Diagram())
    #: Domain
    dom: Diagram = field(default_factory = lambda : Diagram())
    #: Codomain
    cod: Diagram = field(default_factory = lambda : Diagram())
    #: Image
    im: Diagram = field(default_factory = lambda : Diagram())
    #: Relative
    rel: Diagram = field(default_factory = lambda : Diagram())
    #: Entrance times of the simplices
    entrance_times : list[Real] = field(default_factory = lambda : [])
    #: Dimensions of the simplices
    dimensions: list[int] = field(default_factory = lambda : [])

def save_diagrams(dgms : DiagramEnsemble, file : h5py.Group) -> None :
    """
    Save a six-pack of persistence diagrams to a HDF5 file or group.

    Parameters
    ----------
    dgms :
        6-pack of diagrams to save to file/group.
    file :
        A h5py file or group.
    """
    file.create_dataset(
        'entrance_times', data=np.array(dgms.entrance_times))
    file.create_dataset(
        'dimensions', data=np.array(dgms.dimensions))
    dgm_types = (f.name for f in fields(dgms) if f.name not in ['entrance_times', 'dimensions'])
    for dgm in dgm_types:
        grp = file.create_group(dgm)
        grp.create_dataset('paired', data = getattr(dgms, dgm).paired_as_matrix())
        grp.create_dataset('unpaired', data = np.array(list(getattr(dgms, dgm).unpaired)))

def load_diagrams(file : h5py.Group) -> DiagramEnsemble :
    """
    Load a six-pack of persistence diagrams from a HDF5 file or group.
    """
    dgms = DiagramEnsemble()
    dgms.entrance_times = list(file['entrance_times'])
    dgms.dimensions = list(file['dimensions'])
    dgm_types = (f.name for f in fields(DiagramEnsemble) 
                 if f.name not in ['entrance_times', 'dimensions'])
    for dgm in dgm_types:
        grp = file[dgm]
        setattr(dgms, dgm, Diagram._from_matrices(grp['paired'], grp['unpaired']))
        
    return dgms