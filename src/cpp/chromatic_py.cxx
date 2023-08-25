#include "chromatic.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

PYBIND11_MODULE(chromatic, m)
{
    using namespace chalc;
    using namespace chalc::common;
    namespace py = pybind11;
    m.doc() =
        R"docstring(
            Module containing routines to compute geometric filtered simplicial complexes. 
        )docstring";
    m.def("delaunay_complex", &delaunay_complex,
        R"docstring(
            Returns the Delaunay triangulation of a point cloud in Euclidean space.

            Parameters
            ----------
            x : numpy.ndarray[numpy.float64[m, n]]
                A numpy matrix whose columns represent points.

            Returns
            -------
            K : FilteredComplex
                The Delaunay triangulation.
        )docstring",
        py::arg("x"));
    m.def("chromatic_delrips_complex", &chromatic_delrips_complex,
        R"docstring(
            Computes the chromatic Delaunay-Rips complex of a coloured point cloud.

            Parameters
            ----------
            x : numpy.ndarray[numpy.float64[m, n]]
                A numpy matrix whose columns are points in the point cloud.
            colours : List[int]
                A list of integers describing the colours of the points.
                Note that the actual colours of vertices in the output filtration
                may not correspond to the input colours unless the set of values in
                `colours` is contiguous and `colours[0] = 0`.

            Returns
            -------
            K : FilteredComplex
                The chromatic Delaunay-Rips complex.

            Notes
            -----
            The chromatic Delaunay-Rips complex of the point cloud
            has the same set of simplices as the chromatic alpha complex, 
            but with Vietoris-Rips filtration times.

            See Also
            --------
            chromatic_alpha_complex, chromatic_delcech_complex 
        )docstring",
        py::arg("x"), py::arg("colours"));
    m.def("chromatic_alpha_complex", &chromatic_alpha_complex,
        R"docstring(
            Computes the chromatic alpha complex of a coloured point cloud.
            This is a filtration of the Delaunay complex of the point cloud
            after stratification by the colours, and has the same 
            persistent homotopy type as the Cech complex of the point cloud. 

            Parameters
            ----------
            x : numpy.ndarray[numpy.float64[m, n]]
                A numpy matrix whose columns are points in the point cloud.
            colours : List[int]
                A list of integers describing the colours of the points.
                Note that the actual colours of vertices in the output filtration
                may not correspond to the input colours unless the set of values in
                `colours` is contiguous and `colours[0] = 0`.
            
            Returns
            -------
            K : FilteredComplex
                The chromatic alpha complex.

            See Also
            --------
            chromatic_delrips_complex, chromatic_delcech_complex 
        )docstring",
        py::arg("x"), py::arg("colours"));
    m.def("chromatic_delcech_complex", &chromatic_delcech_complex,
        R"docstring(
            Returns the chromatic Delaunay-Cech complex of a coloured point cloud.

            Parameters
            ----------
            x : numpy.ndarray[numpy.float64[m, n]]
                A numpy matrix whose columns are points in the point cloud.
            colours : List[int]
                A list of integers describing the colours of the points.
                Note that the actual colours of vertices in the output filtration
                may not correspond to the input colours unless the set of values in
                `colours` is contiguous and `colours[0] = 0`.

            Returns
            -------
            K : FilteredComplex
                The chromatic Delaunay-Cech complex.

            Notes
            -----
            The chromatic Delaunay-Cech complex of the point cloud
            has the same set of simplices as the chromatic alpha complex, 
            but with Cech filtration times.

            See Also
            --------
            chromatic_alpha_complex, chromatic_delrips_complex 
        )docstring",
        py::arg("x"), py::arg("colours"));

    py::add_ostream_redirect(m);
}