"""

            Module containing geometry routines to compute chromatic complexes. 
        
"""
from __future__ import annotations
import chalc.filtration
import numpy
import typing
__all__ = ['MAX_NUM_COLOURS', 'chromatic_alpha_complex', 'chromatic_delcech_complex', 'chromatic_delrips_complex', 'delaunay_complex']
def chromatic_alpha_complex(x: numpy.ndarray[numpy.float64[m, n]], colours: typing.List[int]) -> chalc.filtration.FilteredComplex:
    """
                Computes the chromatic alpha complex of a coloured point cloud. 
    
                Parameters
                ----------
                x : 
                    A numpy matrix whose columns are points in the point cloud.
                colours :
                    A list of integers describing the colours of the points.
                    Note that the actual colours of vertices in the output filtration
                    may not correspond to the input colours unless the set of values in
                    `colours` is contiguous and `colours[0] = 0`.
                
                Returns
                -------
                FilteredComplex
                    The chromatic alpha complex.
    
                Raises
                ------
                ValueError
                    If the number of unique values in `colours` is greater than `MAX_NUM_COLOURS`.
    
                See Also
                --------
                chromatic_delrips_complex, chromatic_delcech_complex 
    """
def chromatic_delcech_complex(x: numpy.ndarray[numpy.float64[m, n]], colours: typing.List[int]) -> chalc.filtration.FilteredComplex:
    """
                Returns the chromatic Delaunay-Cech complex of a coloured point cloud.
    
                Parameters
                ----------
                x : 
                    A numpy matrix whose columns are points in the point cloud.
                colours : 
                    A list of integers describing the colours of the points.
                    Note that the actual colours of vertices in the output filtration
                    may not correspond to the input colours unless the set of values in
                    `colours` is contiguous and `colours[0] = 0`.
    
                Returns
                -------
                FilteredComplex
                    The chromatic Delaunay-Cech complex.
    
                Raises
                ------
                ValueError
                    If the number of unique values in `colours` is greater than `MAX_NUM_COLOURS`.
    
                Notes
                -----
                The chromatic Delaunay-Cech complex of the point cloud
                has the same set of simplices as the chromatic alpha complex, 
                but with Cech filtration times.
    
                See Also
                --------
                chromatic_alpha_complex, chromatic_delrips_complex 
    """
def chromatic_delrips_complex(x: numpy.ndarray[numpy.float64[m, n]], colours: typing.List[int]) -> chalc.filtration.FilteredComplex:
    """
                Computes the chromatic Delaunay-Rips complex of a coloured point cloud.
    
                Parameters
                ----------
                x :
                    A numpy matrix whose columns are points in the point cloud.
                colours :
                    A list of integers describing the colours of the points.
                    Note that the actual colours of vertices in the output filtration
                    may not correspond to the input colours unless the set of values in
                    `colours` is contiguous and `colours[0] = 0`.
    
                Returns
                -------
                FilteredComplex
                    The chromatic Delaunay-Rips complex.
    
                Raises
                ------
                ValueError
                    If the number of unique values in `colours` is greater than `MAX_NUM_COLOURS`.
    
                Notes
                -----
                The chromatic Delaunay-Rips complex of the point cloud
                has the same set of simplices as the chromatic alpha complex, 
                but with Vietoris-Rips filtration times.
    
                See Also
                --------
                chromatic_alpha_complex, chromatic_delcech_complex 
    """
def delaunay_complex(x: numpy.ndarray[numpy.float64[m, n]]) -> chalc.filtration.FilteredComplex:
    """
                Returns the Delaunay triangulation of a point cloud in Euclidean space.
    
                Parameters
                ----------
                x :
                    A numpy matrix whose columns represent points.
    
                Returns
                -------
                FilteredComplex
                    The Delaunay triangulation.
    """
MAX_NUM_COLOURS: int = 4
