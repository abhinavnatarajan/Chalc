"""
Module containing utilities to store and manipulate abstract filtered simplicial complexes.
"""
from __future__ import annotations
import typing
__all__ = ['Filtration', 'Simplex', 'complete_complex', 'standard_simplex']
class Filtration:
    """
    Class representing a filtered simplicial complex.
    """
    def __contains__(self, vertices: list[int]) -> bool:
        """
        Check for membership of a simplex in the complex.
        
        Args:
        	vertices : Vertex labels of the simplex to check for.
        """
    def __init__(self, n: int, k: int) -> None:
        """
        Construct a discrete filtered simplicial complex with default filtration time of 0.
        
        Args:
        	n : Number of vertices. Cannot be changed after initialisation.
        	k : Maximum dimension of a simplex that the complex can have.
        		This parameter is required for memory efficiency,
        		and cannot be changed after initialisation.
        """
    def __iter__(self) -> typing.Iterator[Simplex]:
        """
        Iterate over the simplices in the complex, ordered by dimension.
        
        There are no guarantees on the relative order of simplices with the same dimension.
        """
    def __len__(self) -> int:
        """
        The total number of simplices in the complex.
        """
    def __repr__(self) -> str:
        ...
    def add_simplex(self, vertices: list[int], filt_value: float) -> bool:
        """
        Add a simplex to a filtered simplicial complex.
        
        Args:
        	vertices : List of vertex labels corresponding to existing vertices in the complex.
        	filt_value : Filtration value to associate to the new simplex.
        
        Note:
        	Faces of the added simplex that are already present
        	in the simplicial complex will have their filtration values reduced if necessary.
        """
    def boundary_matrix(self, max_dimension: int = -1) -> list[tuple[list[int], int, float, list[int]]]:
        """
        Compute the boundary matrix of the filtration.
        
        Args:
        	max_dimension: The maximum dimension of simplices to be considered.
        
        Returns:
        	A list `x` of simplices in the simplicial complex ordered by
        	filtration time, dimension, and label (in that order).
        	Each simplex :math:`\\sigma` has dimension at most ``max_dimension`` and is
        	represented by a tuple containing the following items.
        
        	1.  A list containing the indices in `x` of the facets of :math:`\\sigma`, sorted in ascending order.
        	2.  The lexicographic key of :math:`\\sigma` in the simplicial complex.
        	3.  The filtration time of :math:`\\sigma`.
        	4.  The set of colours of the vertices of :math:`\\sigma`.
        """
    def get_label_from_vertex_labels(self, vertices: list[int]) -> int:
        """
        Get the dictionary key of a simplex.
        
        The key of a simplex is its lexicographic index with respect to
        its vertex labels sorted in ascending order, counting all
        possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`,
        where :math:`N` is the number of vertices in the complex.
        The simplex need not be present in the simplicial complex.
        
        Args:
        	vertices : List of vertex labels of the simplex.
        """
    def is_filtration(self) -> bool:
        """
        Check if the filtration property is satisfied.
        
        Returns true if each simplex has a filtration value at least as large as each of its faces.
        """
    def propagate_colours(self) -> None:
        """
        Ensure that simplex colours are consistent
        with the colours of their vertices.
        
        You should call this whenever you change the colour of any vertices.
        """
    def propagate_filt_values(self, start_dim: int, upwards: bool = True) -> None:
        """
        Propagate filtration values upwards or downwards.
        
        If propagating upwards, the filtration value of each simplex
        will be set to the maximum filtration value of any of its faces.
        If propagating downwards, the filtration value of each simplex
        will be set to the minimum filtration value of any of its cofaces.
        For example, setting the filtration values in dimension 1
        to edge lengths and propagating upwards gives Rips-type filtration.
        
        Args:
        	start_dim : Dimension from which to start propagating (exclusive).
        	upwards : If true then values are propagated upwards, downwards otherwise. Defaults to true.
        """
    def skeleton(self, k: int) -> Filtration:
        """
        Get a copy of the k-skeleton of the filtration.
        
        Args:
        	k: Dimension of the skeleton to return.
        """
    @property
    def dimension(self) -> int:
        """
        Current maximum dimension of a maximal simplex in the complex.
        """
    @property
    def max_dimension(self) -> int:
        """
        Maximum dimension of simplex that this complex can store.
        
        Set during initialisation.
        """
    @property
    def num_vertices(self) -> int:
        """
        Number of vertices in the simplicial complex.
        
        Set during initialisation.
        """
    @property
    def simplices(self) -> list[dict[int, Simplex]]:
        """
        A list such that ``simplices[k]`` is a dictionary of handles
        to the :math:`k`-simplices in the complex.
        
        The key of a :math:`k`-simplex in ``simplices[k]`` is the lexicographic index
        of that simplex with respect to its vertex labels sorted in ascending order,
        counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.
        """
class Simplex:
    """
    Class representing a simplex in a filtered simplicial complex.
    """
    def __repr__(self) -> str:
        ...
    def set_colour(self, colour: int) -> None:
        """
        Change the colour of a vertex.
        
        Raises:
        	ValueError: If the simplex is not a vertex or if
        		`colour >=` :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>`.
        
        Tip:
        	It is recommended to call the member function
        	:meth:`propagate_colours() <chalc.filtration.Filtration.propagate_colours>`
        	from the parent simplicial complex after changing the colour of a vertex.
        """
    @property
    def colours(self) -> list[int]:
        """
        Set of colours of the vertices of the simplex.
        """
    @property
    def dimension(self) -> int:
        """
        			Dimension of the simplex.
        """
    @property
    def facets(self) -> list[Simplex]:
        """
        Read-only list of handles to the facets of the simplex.
        """
    @property
    def filtration_value(self) -> float:
        """
        Filtration value of the simplex.
        
        If you modify this value, you should call
        :meth:`propagate_filt_values() <chalc.filtration.Filtration.propagate_filt_values>`
        from the parent complex to ensure that filtration times remain monotonic.
        """
    @filtration_value.setter
    def filtration_value(self, arg1: float) -> None:
        ...
    @property
    def label(self) -> int:
        """
        Label of the simplex in its parent filtered complex.
        
        A :math:`k`-simplex :math:`\\sigma` is labelled by the lexicographic index of :math:`\\sigma`
        with respect to its vertex labels sorted in ascending order,
        counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.
        """
    @property
    def vertices(self) -> list[int]:
        """
        List of (sorted, ascending) vertex labels of the simplex.
        """
def complete_complex(n: int, k: int) -> Filtration:
    """
    Compute the :math:`k`-skeleton of the complete simplicial complex on :math:`n` vertices.
    
    Filtration values are initialised to zero and all vertices coloured with the colour 0.
    
    Raises:
    	RuntimeError:
    		If ``n<= 0`` or ``k >= n`` or ``k < 0``.
    """
def standard_simplex(n: int) -> Filtration:
    """
    Compute the filtered simplicial complex corresponding to the standard abstract :math:`n`-simplex.
    
    Filtration values are initialised to zero and all vertices coloured with the colour 0.
    """
