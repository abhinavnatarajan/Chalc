"""

            Module containing utilities to store and manipulate abstract filtered simplicial complexes. 
        
"""
from __future__ import annotations
import typing
__all__ = ['FilteredComplex', 'Simplex', 'clique_complex', 'standard_simplex']
class FilteredComplex:
    """
    
                Class representing a filtered simplicial complex.
            
    """
    @typing.overload
    def __init__(self, n: int, k: int) -> None:
        """
                    Construct a discrete filtered simplicial complex 
                    with default filtration time of 0.
                    
                    Parameters
                    ----------
                    n :
                        Number of vertices. Cannot be changed after initialisation.
                    k :
                        Maximum dimension of a simplex that the complex can have.
                        This parameter is required for memory efficiency, and cannot
                        be changed after initialisation.
        """
    @typing.overload
    def __init__(self, other: FilteredComplex, k: int) -> None:
        """
                    Returns a handle to the :math:`k`-skeleton 
                    of an existing (filtered) simplicial complex.
                    
                    Parameters
                    ----------
                    other :
                        Filtered simplicial complex.
                    k :
                        Skeleton dimension.
        """
    def __repr__(self) -> str:
        ...
    def add_simplex(self, vertices: typing.List[int], filt_value: float) -> bool:
        """
                        Add a simplex to a filtered simplicial complex.
                        
                        Parameters
                        ----------
                        vertices :
                            List of vertex labels corresponding to existing vertices 
                            in the complex.
                        filt_value :
                            Filtration value to associate to the new simplex. 
                        
                        Note
                        ----
                        Faces of the added simplex that are already present 
                        in the simplicial complex will have their filtration 
                        values reduced if necessary.                 
        """
    def count_simplices_in_dim(self, dimension: int) -> int:
        """
                        Count the number of simplices of a given dimension.
        """
    def get_label_from_vertex_labels(self, vertices: typing.List[int]) -> int:
        """
                        Returns the dictionary key of a simplex with respect to 
                        its vertex labels sorted in ascending order, counting all 
                        possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`, 
                        where :math:`N` is the number of vertices in the complex. 
                        The simplex need not be present in the simplicial complex.
                        
                        Parameters
                        ----------
                        vertices : 
                            List of vertex labels of the simplex.
        """
    def has_simplex(self, vertices: typing.List[int]) -> bool:
        """
                        Check for membership of a simplex in the complex.
                        
                        Parameters
                        ----------
                        vertices : 
                            Vertex labels of the simplex to check for. 
        """
    def propagate_colours(self) -> None:
        """
                        Function to make sure that simplex colours are consistent 
                        with the colours of their vertices. You should call this 
                        whenever you change the colour of any vertices. 
        """
    def propagate_filt_values(self, start_dim: int, upwards: bool = ...) -> None:
        """
                        Propagate filtration values upwards or downwards to ensure that
                        every simplex appears after its faces.
                        For example, setting the filtration values in dimension 1
                        and propagating upwards is akin to the Rips filtration.
                        
                        Parameters
                        ----------
                        start_dim :
                            Dimension from which to start propagating (exclusive). 
                        upwards :
                            If true then values are propagated upwards, downwards otherwise.
        """
    def serialised(self) -> typing.List[typing.Tuple[typing.List[int], int, float, int]]:
        """
                        Serialised representation of the simplicial complex in a format
                        suitable for persistent homology computations.
        
                        :return:
                            A list `x` of simplices in the simplicial complex ordered by dimension 
                            followed by filtration time. Each simplex :math:`\sigma` is 
                            represented by a tuple containing the following items.
        
                            1.  A list containing the indices in `x` of the facets of :math:`\sigma`.
                            2.  The lexicographic key of :math:`\sigma` in the simplicial complex.
                            3.  The filtration time of :math:`\sigma`.
                            4.  The colour bitmask of :math:`\sigma`.
        """
    @property
    def dimension(arg0: FilteredComplex) -> int:
        """
        
                        Current maximum dimension of a maximal simplex in the complex.
                    
        """
    @property
    def max_dimension(self) -> int:
        """
        
                        Maximum dimension of simplex that this complex can store. Set during initialisation. 
                    
        """
    @property
    def num_simplices(arg0: FilteredComplex) -> int:
        """
        
                        The total number of simplices in the complex.
                    
        """
    @property
    def num_vertices(self) -> int:
        """
        
        
                        Number of vertices in the simplicial complex. Set during initialisation.
                    
        """
    @property
    def simplices(arg0: FilteredComplex) -> typing.List[typing.Dict[int, Simplex]]:
        """
        
        
                        A list such that `simplices[k]` contains \
                        the dictionary of :math:`k`-simplices of the complex. \
                        The key of a :math:`k`-simplex in `simplices[k]` is \
                        the lexicographic index of that simplex with respect \
                        to its vertex labels sorted in ascending order, \
                        counting all possible sorted subsequences of :math:`(0, ..., N-1)` \
                        of length :math:`k`.
                    
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
                        
                        Raises
                        ------
                        ValueError
                            If the simplex is not a vertex or if `colour >= MAX_NUM_COLOURS`.
                        
                        Tip
                        ---
                        It is recommended to call the member function 
                        `propagate_colours` of the parent simplicial complex 
                        after changing the colour of a vertex.
        """
    @property
    def colours(arg0: Simplex) -> int:
        """
        
                        Bitmask of colours in the simplex, where the rightmost \
                        bit represents the smallest colour index. 
                    
        """
    @property
    def dimension(self) -> int:
        """
        
                    Dimension of the simplex.
                
        """
    @property
    def facets(arg0: Simplex) -> typing.List[Simplex]:
        """
        
                        Returns a read-only list of handles to the facets of the simplex.
                    
        """
    @property
    def filtration_value(self) -> float:
        """
        
                        Filtration value of the simplex. If you modify this value, \
                        you must call `propagate_filtration_value` on the parent \
                        complex to ensure that filtration times remain monotonic.
                    
        """
    @filtration_value.setter
    def filtration_value(self, arg0: float) -> None:
        ...
    @property
    def label(self) -> int:
        """
        
                        Label of the simplex in its parent filtered complex. \
                        A :math:`k`-simplex :math:`\sigma` is labelled 
                        by the lexicographic index of :math:`\sigma` with respect 
                        to its vertex labels sorted in ascending order, 
                        counting all possible sorted subsequences of :math:`(0, ..., N-1)`
                        of length :math:`k`.
                    
        """
    @property
    def vertices(arg0: Simplex) -> typing.List[int]:
        """
        
                        List of (sorted, ascending) vertex labels of the simplex.
                    
        """
def clique_complex(n: int, k: int) -> FilteredComplex:
    """
                Returns the :math:`k`-skeleton of the complete simplicial complex
                on :math:`n` vertices, with filtration values initialised to zero.
    """
def standard_simplex(n: int) -> FilteredComplex:
    """
                Returns the simplicial complex corresponding to 
                the standard abstract :math:`n`-simplex.
    """
