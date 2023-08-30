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
            Module containing geometry routines to compute chromatic complexes. 
        )docstring";
    m.attr("MAX_NUM_COLOURS") = MAX_NUM_COLOURS;
    m.def("delaunay_complex", &delaunay_complex,
        R"docstring(
            Returns the Delaunay triangulation of a point cloud in Euclidean space.

            Parameters
            ----------
            x :
                A numpy matrix whose columns represent points.

            Returns
            -------
            FilteredComplex
                The Delaunay triangulation.
        )docstring",
        py::arg("x"))
    .def("chromatic_delrips_complex", &chromatic_delrips_complex,
        R"docstring(
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
        )docstring",
        py::arg("x"), py::arg("colours"))
    .def("chromatic_alpha_complex", 
    [](const Eigen::MatrixXd& points, const vector<index_t>& colours) {
        std::ostringstream ostream;
        auto res = chromatic_alpha_complex(points, colours, ostream);
        if (ostream.tellp() > 0) {
            PyErr_WarnEx(PyExc_RuntimeWarning, 
                        ostream.str().c_str(), 
                        1);
        }
        return res;
    },
        R"docstring(
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
        )docstring",
        py::arg("x"), py::arg("colours"))
    .def("chromatic_delcech_complex", 
    [](const Eigen::MatrixXd& points, const vector<index_t>& colours) {
        std::ostringstream ostream;
        auto res = chromatic_delcech_complex(points, colours, ostream);
        if (ostream.tellp() > 0) {
            PyErr_WarnEx(PyExc_RuntimeWarning, 
                        ostream.str().c_str(), 
                        1);
        }
        return res;
    },
        R"docstring(
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
        )docstring",
        py::arg("x"), py::arg("colours"));
}