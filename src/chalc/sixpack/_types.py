from dataclasses import dataclass, fields
import numpy as np
import h5py
from numbers import Real

@dataclass
class Diagram:
    """ 
    Persistence diagram object, represented by a 
    set of simplex pairings and a set of unpaired simplices.
    """
    #: Set of tuples of indices of paired simplices. 
    paired: set[tuple[int, int]]
    #: Set of indices of unpaired simplices.
    unpaired: set[int]

    def __init__(self, obj):
        self.paired = obj.paired
        self.unpaired = obj.unpaired

    def pair_matrix(self):
        return np.concatenate(
            list(self.paired)).reshape((len(self.paired), 2))
    


@dataclass
class DiagramEnsemble:
    """
    Six pack of persistence diagrams.
    """
    #: Kernel
    ker: Diagram
    #: Cokernel
    cok: Diagram
    #: Domain
    dom: Diagram
    #: Codomain
    cod: Diagram
    #: Image
    im: Diagram
    #: Relative
    rel: Diagram
    #: Entrance times of the simplices
    entrance_times : list[Real]
    #: Dimensions of the simplices
    dimensions: list[int]

def save_diagrams(dgms : DiagramEnsemble, file : h5py.File) -> None :
    file.create_dataset(
        'entrance_times', data=np.array(dgms.entrance_times))
    file.create_dataset(
        'dimensions', data=np.array(dgms.dimensions))
    dgm_types = (f for f in fields(dgms) if f.name not in ['entrance_times', 'dimensions'])
    for dgm in dgm_types:
        grp = file.create_group(dgm.name)
        grp.create_dataset('paired', data = getattr(dgms, dgm.name).pair_matrix())
        grp.create_dataset('unpaired', data = np.array(list(getattr(dgms, dgm.name).unpaired)))
        
        
