#include <chalc/chromatic/chromatic.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(chromatic, m) {
    using namespace chalc::chromatic;
    using namespace chalc::stl;
    using namespace chalc;
    namespace py = pybind11;
    py::object log_warn = py::module_::import("logging")
                              .attr("getLogger")("chalc.chromatic")
                              .attr("warning");
    m.doc() =
        R"docstring(
            Module containing geometry routines to compute chromatic Delaunay filtrations.

            Attributes:
                MaxColoursChromatic (int): Maximum number of colours that can be handled by the methods in this module.
        )docstring";
    m.attr("MaxColoursChromatic") = py::int_(MAX_NUM_COLOURS);
    m.def("delaunay", &delaunay,
          R"docstring(
            Returns the chromatic Delaunay triangulation of a coloured point cloud in Euclidean space.

            Args:
                x : Numpy matrix whose columns are points in the point cloud.
                colours : List of integers describing the colours of the points.

            Raises:
                ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

            Returns:
                The Delaunay triangulation.
        )docstring",
          py::arg("x"),py::arg("colours"))
        .def(
            "delrips",
            [](const Eigen::MatrixXd &points, const vector<index_t> &colours) {
                auto res = delrips(points, colours);
                return tuple{res, false};
            },
            R"docstring(
            Computes the chromatic Delaunay-Rips filtration of a coloured point cloud.

            Args:
                x : Numpy matrix whose columns are points in the point cloud.
                colours : List of integers describing the colours of the points.

            Returns:
                The chromatic Delaunay-Rips filtration and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.

            Raises:
                ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

            Notes:
                The chromatic Delaunay-Rips filtration of the point cloud has the same set of simplices as the chromatic alpha filtration, but with Vietoris-Rips filtration times. The convention used is that the filtration time of a simplex is half the maximum edge length in that simplex. With this convention, the chromatic Delaunay-Rips filtration and chromatic alpha filtration have visually similar persistence diagrams.

            See Also:
                :func:`alpha`, :func:`delcech`
        )docstring",
            py::arg("x"), py::arg("colours"))
        .def(
            "alpha",
            [log_warn](const Eigen::MatrixXd &points,
                       const vector<index_t> &colours) {
                std::ostringstream ostream;
                auto res = alpha(points, colours, ostream);
                bool issues = false;
                if (ostream.tellp() > 0) {
                    issues = true;
                    log_warn(ostream.str().c_str());
                }
                return tuple{res, issues};
            },
            R"docstring(
            Computes the chromatic alpha filtration of a coloured point cloud.

            Args:
                x : Numpy matrix whose columns are points in the point cloud.
                colours : List of integers describing the colours of the points.

            Returns:
                The chromatic alpha filtration and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.

            Raises:
                ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

            See Also:
                :func:`delrips`, :func:`delcech`
        )docstring",
            py::arg("x"), py::arg("colours"))
        .def(
            "delcech",
            [log_warn](const Eigen::MatrixXd &points,
                       const vector<index_t> &colours) {
                std::ostringstream ostream;
                auto res = delcech(points, colours, ostream);
                bool issues = false;
                if (ostream.tellp() > 0) {
                    issues = true;
                    log_warn(ostream.str().c_str());
                }
                return tuple{res, issues};
            },
            R"docstring(
            Returns the chromatic Delaunay-Cech filtration of a coloured point cloud.

            Args:
                x : Numpy matrix whose columns are points in the point cloud.
                colours : List of integers describing the colours of the points.

            Returns:
                The chromatic Delaunay-Cech filtration and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.

            Raises:
                ValueError : If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

            Notes:
                The chromatic Delaunay-Cech filtration of the point cloud has the same set of simplices as the chromatic alpha filtration, but with Cech filtration times.

            See Also:
                :func:`alpha`, :func:`delrips`
        )docstring",
            py::arg("x"), py::arg("colours"));
}
