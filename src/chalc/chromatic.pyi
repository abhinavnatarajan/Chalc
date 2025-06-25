"""
Module containing geometry routines to compute chromatic Delaunay filtrations.
"""
from __future__ import annotations
import chalc.filtration
import numpy
import typing
__all__ = ['MaxColoursChromatic', 'alpha', 'delaunay', 'delcech', 'delrips']
M = typing.TypeVar("M", bound=int)
N = typing.TypeVar("N", bound=int)
@typing.overload
def alpha(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: numpy.ndarray[tuple[M, typing.Literal[1]], numpy.dtype[numpy.uint16]], max_num_threads: int = 1) -> tuple[chalc.filtration.Filtration, bool]:
    ...
@typing.overload
def alpha(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: list[int], max_num_threads: int = 1) -> tuple[chalc.filtration.Filtration, bool]:
    """
    Compute the chromatic alpha filtration of a coloured point cloud.
    
    Args:
    	points : Numpy matrix whose columns are points in the point cloud.
    	colours : List or numpy array of integers describing the colours of the points.
    	max_num_threads: Maximum number of parallel threads to use.
    		If non-positive, the number of threads to use is automatically determined
    		by the threading library (Intel OneAPI TBB). Note that this may be less
    		than the number of available CPU cores depending on the number of points
    		and the system load. The default is 1, which means no parallelism.
    
    Returns:
    	The chromatic alpha filtration and a boolean flag to
    	indicate if numerical issues were encountered.
    	In case of numerical issues, a warning is also raised.
    
    Raises:
    	ValueError:
    		If any value in ``colours`` is
    		>= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0,
    		or if the length of ``colours`` does not match the number of points.
    	RuntimeError:
    		If the dimension of the point cloud + the number of colours is too large
    		for computations to run without overflowing.
    
    Notes:
    	This function is included for pedantic reasons.
    	For most purposes you should instead consider using :func:`chalc.chromatic.delcech`,
    	which is faster to compute, more numerically stable, and has the same persistent homology.
    
    See Also:
    	:func:`delrips`, :func:`delcech`
    """
@typing.overload
def delaunay(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: numpy.ndarray[tuple[M, typing.Literal[1]], numpy.dtype[numpy.uint16]]) -> chalc.filtration.Filtration:
    ...
@typing.overload
def delaunay(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: list[int]) -> chalc.filtration.Filtration:
    """
    Compute the chromatic Delaunay triangulation of a coloured point cloud in Euclidean space.
    
    Args:
    	points : Numpy matrix whose columns are points in the point cloud.
    	colours : List or numpy array of integers describing the colours of the points.
    
    Raises:
    	ValueError:
    		If any value in ``colours`` is
    		>= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0,
    		or if the length of ``colours`` does not match the number of points.
    	RuntimeError:
    		If the dimension of the point cloud + the number of colours is too large
    		for computations to run without overflowing.
    
    Returns:
    	The Delaunay triangulation.
    """
@typing.overload
def delcech(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: numpy.ndarray[tuple[M, typing.Literal[1]], numpy.dtype[numpy.uint16]], max_num_threads: int = 1) -> tuple[chalc.filtration.Filtration, bool]:
    ...
@typing.overload
def delcech(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: list[int], max_num_threads: int = 1) -> tuple[chalc.filtration.Filtration, bool]:
    """
    Compute the chromatic Delaunay--Čech filtration of a coloured point cloud.
    
    Args:
    	points : Numpy matrix whose columns are points in the point cloud.
    	colours : List or numpy array of integers describing the colours of the points.
    	max_num_threads: Maximum number of parallel threads to use.
    		If non-positive, the number of threads to use is automatically determined
    		by the threading library (Intel OneAPI TBB). Note that this may be less
    		than the number of available CPU cores depending on the number of points
    		and the system load. The default is 1, which means no parallelism.
    
    Returns:
    	The chromatic Delaunay--Čech filtration and a boolean flag to indicate
    	if numerical issues were encountered.
    	In case of numerical issues, a warning is also raised.
    
    Raises:
    	ValueError :
    		If any value in ``colours`` is
    		>= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0,
    		or if the length of ``colours`` does not match the number of points.
    	RuntimeError:
    		If the dimension of the point cloud + the number of colours is too large
    		for computations to run without overflowing.
    
    Notes:
    	The chromatic Delaunay--Čech filtration of the point cloud
    	has the same set of simplices as the chromatic alpha filtration,
    	but with Čech filtration times.
    
    See Also:
    	:func:`alpha`, :func:`delrips`
    """
@typing.overload
def delrips(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: numpy.ndarray[tuple[M, typing.Literal[1]], numpy.dtype[numpy.uint16]], max_num_threads: int = 1) -> tuple[chalc.filtration.Filtration, bool]:
    ...
@typing.overload
def delrips(points: numpy.ndarray[tuple[M, N], numpy.dtype[numpy.float64]], colours: list[int], max_num_threads: int = 1) -> tuple[chalc.filtration.Filtration, bool]:
    """
    Compute the chromatic Delaunay--Rips filtration of a coloured point cloud.
    
    Args:
    	points : Numpy matrix whose columns are points in the point cloud.
    	colours : List or numpy array of integers describing the colours of the points.
    	max_num_threads: Maximum number of parallel threads to use.
    		If non-positive, the number of threads to use is automatically determined
    		by the threading library (Intel OneAPI TBB). Note that this may be less
    		than the number of available CPU cores depending on the number of points
    		and the system load. The default is 1, which means no parallelism.
    
    Returns:
    	The chromatic Delaunay--Rips filtration and a boolean flag to indicate
    	if numerical issues were encountered.
    	In case of numerical issues, a warning is also raised.
    
    Raises:
    	ValueError:
    		If any value in ``colours`` is
    		>= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0,
    		or if the length of ``colours`` does not match the number of points.
    	RuntimeError:
    		If the dimension of the point cloud + the number of colours is too large
    		for computations to run without overflowing.
    
    Notes:
    	The chromatic Delaunay--Rips filtration of the point cloud
    	has the same set of simplices as the chromatic alpha filtration,
    	but with Vietoris--Rips filtration times.
    	The convention used is that the filtration time of a simplex
    	is half the maximum edge length in that simplex.
    	With this convention, the chromatic Delaunay--Rips filtration
    	and chromatic alpha filtration have the same persistence diagrams
    	in degree zero.
    
    See Also:
    	:func:`alpha`, :func:`delcech`
    """
MaxColoursChromatic: int = 16
