#include "chromatic.h"
#include "../common.h"
#include <pybind11/eigen.h>

PYBIND11_MODULE(chromatic, m)
{
    using namespace chalc::chromatic;
    using namespace chalc::stl;
    using namespace chalc;
    namespace py = pybind11;
    m.doc() =
        R"docstring(
            Module containing geometry routines to compute chromatic complexes.

            Attributes
            ----------
            MaxColoursChromatic : int
                Maximum number of colours that can be handled by the methods in this module.
        )docstring";
    m.attr("MaxColoursChromatic") = py::int_(MAX_NUM_COLOURS);
    m.def("delaunay", &delaunay,
        R"docstring(
            Returns the Delaunay triangulation of a point cloud in Euclidean space.

            Parameters
            ----------
            x :
                Numpy matrix whose columns are points in the point cloud.

            Returns
            -------
            FilteredComplex
                The Delaunay triangulation.
        )docstring",
        py::arg("x"))
    .def("delrips", &delrips,
        R"docstring(
            Computes the chromatic Delaunay-Rips complex of a coloured point cloud.

            Parameters
            ----------
            x :
                Numpy matrix whose columns are points in the point cloud.
            colours :
                List of integers describing the colours of the points.
            
            Returns
            -------
            FilteredComplex
                The chromatic Delaunay-Rips complex.

            Raises
            ------
            ValueError
                If any value in ``colours`` is >= 
                :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` 
                or < 0, or if the length of ``colours`` does not match the number of points.

            Notes
            -----
            The chromatic Delaunay-Rips complex of the point cloud
            has the same set of simplices as the chromatic alpha complex, 
            but with Vietoris-Rips filtration times.

            See Also
            --------
            alpha, delcech
        )docstring",
        py::arg("x"), py::arg("colours"))
    .def("alpha",
    [](const Eigen::MatrixXd& points, const vector<index_t>& colours) {
        std::ostringstream ostream;
        auto res = alpha(points, colours, ostream);
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
                Numpy matrix whose columns are points in the point cloud.
            colours :
                List of integers describing the colours of the points.

            Returns
            -------
            FilteredComplex
                The chromatic alpha complex.

            Raises
            ------
            ValueError
                If any value in ``colours`` is >= 
                :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` 
                or < 0, or if the length of ``colours`` does not match the number of points.

            See Also
            --------
            delrips, delcech
        )docstring",
        py::arg("x"), py::arg("colours"))
    .def("delcech",
    [](const Eigen::MatrixXd& points, const vector<index_t>& colours) {
        std::ostringstream ostream;
        auto res = delcech(points, colours, ostream);
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
                Numpy matrix whose columns are points in the point cloud.
            colours :
                List of integers describing the colours of the points.

            Returns
            -------
            FilteredComplex
                The chromatic Delaunay-Cech complex.

            Raises
            ------
            ValueError
                If any value in ``colours`` is >= 
                :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` 
                or < 0, or if the length of ``colours`` does not match the number of points.

            Notes
            -----
            The chromatic Delaunay-Cech complex of the point cloud
            has the same set of simplices as the chromatic alpha complex, 
            but with Cech filtration times.

            See Also
            --------
            alpha, delrips
        )docstring",
        py::arg("x"), py::arg("colours"));
}
