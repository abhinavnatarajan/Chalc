"""

            Module containing geometry routines to compute chromatic complexes.

            Attributes:
                MaxColoursChromatic (int): Maximum number of colours that can be handled by the methods in this module.
        
"""
from __future__ import annotations
import chalc.filtration
import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['MaxColoursChromatic', 'alpha', 'delaunay', 'delcech', 'delrips']
def alpha(x: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], colours: list[int]) -> tuple[chalc.filtration.FilteredComplex, bool]:
    """
                Computes the chromatic alpha complex of a coloured point cloud.
    
                Args:
                    x : Numpy matrix whose columns are points in the point cloud.
                    colours : List of integers describing the colours of the points.
    
                Returns:
                    The chromatic alpha complex and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.
    
                Raises:
                    ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.
    
                See Also:
                    :func:`delrips`, :func:`delcech`
    """
def delaunay(x: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> chalc.filtration.FilteredComplex:
    """
                Returns the Delaunay triangulation of a point cloud in Euclidean space.
    
                Args:
                    x : Numpy matrix whose columns are points in the point cloud.
    
                Returns:
                    The Delaunay triangulation.
    """
def delcech(x: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], colours: list[int]) -> tuple[chalc.filtration.FilteredComplex, bool]:
    """
                Returns the chromatic Delaunay-Cech complex of a coloured point cloud.
    
                Args:
                    x : Numpy matrix whose columns are points in the point cloud.
                    colours : List of integers describing the colours of the points.
    
                Returns:
                    The chromatic Delaunay-Cech complex and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.
    
                Raises:
                    ValueError : If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.
    
                Notes:
                    The chromatic Delaunay-Cech complex of the point cloud has the same set of simplices as the chromatic alpha complex, but with Cech filtration times.
    
                See Also:
                    :func:`alpha`, :func:`delrips`
    """
def delrips(x: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], colours: list[int]) -> tuple[chalc.filtration.FilteredComplex, bool]:
    """
                Computes the chromatic Delaunay-Rips complex of a coloured point cloud.
    
                Args:
                    x : Numpy matrix whose columns are points in the point cloud.
                    colours : List of integers describing the colours of the points.
    
                Returns:
                    The chromatic Delaunay-Rips complex and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.
    
                Raises:
                    ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.
    
                Notes:
                    The chromatic Delaunay-Rips complex of the point cloud has the same set of simplices as the chromatic alpha complex, but with Vietoris-Rips filtration times. The convention used is that the filtration time of a simplex is half the maximum edge length in that simplex. With this convention, the chromatic Delaunay-Rips complex and chromatic alpha complex have visually similar persistence diagrams.
    
                See Also:
                    :func:`alpha`, :func:`delcech`
    """
MaxColoursChromatic: int = 64
