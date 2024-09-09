#include <chalc/chromatic/chromatic.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(chromatic, m) {
	using namespace chalc::chromatic;
	using namespace chalc::stl;
	using namespace chalc;
	namespace py = pybind11;
	py::object log_warn =
		py::module_::import("logging").attr("getLogger")("chalc.chromatic").attr("warning");
	m.doc() =
		R"docstring(
Module containing geometry routines to compute chromatic Delaunay filtrations.
        )docstring";
	m.attr("MaxColoursChromatic") = py::int_(MAX_NUM_COLOURS);
	m.def(
		 "delaunay",
		 [](const Eigen::MatrixXd& points, const Eigen::VectorXi& colours) {
			 std::vector<long long int> colours_vec(colours.data(),
		                                            colours.data() + colours.size());
			 return delaunay(points, colours_vec);
		 },
		 py::arg("x"),
		 py::arg("colours"))
		.def("delaunay",
	         &delaunay,
	         R"docstring(
Returns the chromatic Delaunay triangulation of a coloured point cloud in Euclidean space.

Args:
	x : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.

Raises:
	ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

Returns:
	The Delaunay triangulation.
	      )docstring",
	         py::arg("x"),
	         py::arg("colours"))
		.def(
			"delrips",
			[](const Eigen::MatrixXd& points, const Eigen::VectorXi& colours) {
				std::vector<long long int> colours_vec(colours.data(),
		                                               colours.data() + colours.size());
				return tuple{delrips(points, colours_vec), false};
			},
			py::arg("x"),
			py::arg("colours"))
		.def(
			"delrips",
			[](const Eigen::MatrixXd& points, const vector<index_t>& colours) {
				return tuple{delrips(points, colours), false};
			},
			R"docstring(
Computes the chromatic Delaunay--Rips filtration of a coloured point cloud.

Args:
	x : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.

Returns:
	The chromatic Delaunay--Rips filtration and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.

Raises:
	ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

Notes:
	The chromatic Delaunay--Rips filtration of the point cloud has the same set of simplices as the chromatic alpha filtration, but with Vietoris--Rips filtration times. The convention used is that the filtration time of a simplex is half the maximum edge length in that simplex. With this convention, the chromatic Delaunay--Rips filtration and chromatic alpha filtration have the same persistence diagrams in degree zero.

See Also:
	:func:`alpha`, :func:`delcech`
			)docstring",
			py::arg("x"),
			py::arg("colours"))
		.def(
			"alpha",
			[](const Eigen::MatrixXd& points, const Eigen::VectorXi& colours) {
				std::vector<long long int> colours_vec(colours.data(),
		                                               colours.data() + colours.size());
				return alpha(points, colours_vec);
			},
			py::arg("x"),
			py::arg("colours"))
		.def("alpha",
	         &alpha,
	         R"docstring(
Computes the chromatic alpha filtration of a coloured point cloud.

Args:
	x : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.

Returns:
	The chromatic alpha filtration and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.

Raises:
	ValueError: If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

Notes:
	This function is included for pedantic reasons. For most purposes you should instead consider using :func:`chalc.chromatic.delcech`, which is faster to compute, more numerically stable, and has the same persistent homology.

See Also:
	:func:`delrips`, :func:`delcech`
	         )docstring",
	         py::arg("x"),
	         py::arg("colours"))
		.def(
			"delcech",
			[](const Eigen::MatrixXd& points, const Eigen::VectorXi& colours) {
				std::vector<long long int> colours_vec(colours.data(),
		                                               colours.data() + colours.size());
				return delcech(points, colours_vec);
			},
			py::arg("x"),
			py::arg("colours"))
		.def("delcech",
	         &delcech,
	         R"docstring(
Returns the chromatic Delaunay--Čech filtration of a coloured point cloud.

Args:
	x : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.

Returns:
	The chromatic Delaunay--Čech filtration and a boolean flag to indicate if numerical issues were encountered. In case of numerical issues, a warning is also raised.

Raises:
	ValueError : If any value in ``colours`` is >= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0, or if the length of ``colours`` does not match the number of points.

Notes:
	The chromatic Delaunay--Čech filtration of the point cloud has the same set of simplices as the chromatic alpha filtration, but with Čech filtration times.

See Also:
	:func:`alpha`, :func:`delrips`
	         )docstring",
	         py::arg("x"),
	         py::arg("colours"));
}
