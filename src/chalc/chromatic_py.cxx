#include <boost/type_traits/type_with_alignment.hpp>
#include <chalc/chromatic/chromatic.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(chromatic, m) {  // NOLINT
	using chalc::alpha;
	using chalc::alpha_parallel;
	using chalc::colour_t;
	using chalc::delaunay;
	using chalc::delcech;
	using chalc::delcech_parallel;
	using chalc::delrips;
	using chalc::delrips_parallel;
	using chalc::index_t;
	using chalc::label_t;
	using chalc::MAX_NUM_COLOURS;
	using std::tuple;
	using std::vector;
	namespace py = pybind11;

	py::object log_warn =
		py::module_::import("logging").attr("getLogger")("chalc.chromatic").attr("warning");
	m.doc() = "Module containing geometry routines to compute chromatic Delaunay filtrations.";
	m.attr("MaxColoursChromatic") = py::int_(MAX_NUM_COLOURS);
	m.def(
		 "delaunay",
		 [](const Eigen::MatrixXd&          points,
	        const Eigen::VectorX<colour_t>& colours_ndarray,
	        bool                            parallel) {
			 const vector<colour_t> colours(colours_ndarray.cbegin(), colours_ndarray.cend());
			 if (parallel) {
				 return delaunay<CGAL::Parallel_tag>(points, colours);
			 } else {
				 return delaunay<CGAL::Sequential_tag>(points, colours);
			 }
		 },
		 py::arg("points"),
		 py::arg("colours"),
		 py::arg("parallel") = true
	)
		.def(
			"delaunay",
			[](const Eigen::MatrixXd& points, const vector<colour_t>& colours, bool parallel) {
				if (parallel) {
					return delaunay<CGAL::Parallel_tag>(points, colours);
				} else {
					return delaunay<CGAL::Sequential_tag>(points, colours);
				}
			},
			R"docstring(Compute the chromatic Delaunay triangulation of a coloured point cloud in Euclidean space.

Args:
	points : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.
	parallel: If true, use parallel computation during the spatial sorting phase of the triangulation.

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

)docstring",
			py::arg("points"),
			py::arg("colours"),
			py::arg("parallel") = true
		)
		.def(
			"delrips",
			[](const Eigen::MatrixXd& points,
	           const Eigen::VectorX<colour_t>& colours_ndarray,
	           const int max_num_threads) {
				const vector<colour_t> colours(colours_ndarray.cbegin(), colours_ndarray.cend());
				if (max_num_threads == 1) {
					return delrips(points, colours);
				} else {
					return delrips_parallel(points, colours, max_num_threads);
				}
			},
			py::arg("points"),
			py::arg("colours"),
			py::arg("max_num_threads") = 0
		)
		.def(
			"delrips",
			[](const Eigen::MatrixXd& points,
	           const vector<colour_t>& colours,
	           const int max_num_threads) {
				if (max_num_threads == 1) {
					return delrips(points, colours);
				} else {
					return delrips_parallel(points, colours, max_num_threads);
				}
			},
			R"docstring(Compute the chromatic Delaunay--Rips filtration of a coloured point cloud.

Args:
	points : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.
	max_num_threads: Hint for maximum number of parallel threads to use.
		If non-positive, the number of threads to use is automatically determined
		by the threading library (Intel OneAPI TBB). Note that this may be less
		than the number of available CPU cores depending on the number of points
		and the system load.

Returns:
	The chromatic Delaunay--Rips filtration.

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

)docstring",
			py::arg("points"),
			py::arg("colours"),
			py::arg("max_num_threads") = 0
		)
		.def(
			"alpha",
			[](const Eigen::MatrixXd& points,
	           const Eigen::VectorX<colour_t>& colours_ndarray,
	           const int max_num_threads) {
				const vector<colour_t> colours(colours_ndarray.cbegin(), colours_ndarray.cend());
				if (max_num_threads == 1) {
					return alpha(points, colours);
				} else {
					return alpha_parallel(points, colours, max_num_threads);
				}
			},
			py::arg("points"),
			py::arg("colours"),
			py::arg("max_num_threads") = 0
		)
		.def(
			"alpha",
			[](const Eigen::MatrixXd& points,
	           const vector<colour_t>& colours,
	           const int max_num_threads) {
				if (max_num_threads == 1) {
					return alpha(points, colours);
				} else {
					return alpha_parallel(points, colours, max_num_threads);
				}
			},
			R"docstring(Compute the chromatic alpha filtration of a coloured point cloud.

Args:
	points : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.
	max_num_threads: Hint for maximum number of parallel threads to use.
		If non-positive, the number of threads to use is automatically determined
		by the threading library (Intel OneAPI TBB). Note that this may be less
		than the number of available CPU cores depending on the number of points
		and the system load.

Returns:
	The chromatic alpha filtration.

Raises:
	ValueError:
		If any value in ``colours`` is
		>= :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>` or < 0,
		or if the length of ``colours`` does not match the number of points.
	RuntimeError:
		If the dimension of the point cloud + the number of colours is too large
		for computations to run without overflowing.

Notes:
	:func:`chalc.chromatic.delcech` has the same 6-pack of persistent homology, and often
	has slightly better performance.

See Also:
	:func:`delrips`, :func:`delcech`

)docstring",
			py::arg("points"),
			py::arg("colours"),
			py::arg("max_num_threads") = 0
		)
		.def(
			"delcech",
			[](const Eigen::MatrixXd& points,
	           const Eigen::VectorX<colour_t>& colours_ndarray,
	           const int max_num_threads) {
				const vector<colour_t> colours(colours_ndarray.cbegin(), colours_ndarray.cend());
				if (max_num_threads == 1) {
					return delcech(points, colours);
				} else {
					return delcech_parallel(points, colours, max_num_threads);
				}
			},
			py::arg("points"),
			py::arg("colours"),
			py::arg("max_num_threads") = 0
		)
		.def(
			"delcech",
			[](const Eigen::MatrixXd& points,
	           const vector<colour_t>& colours,
	           const int max_num_threads) {
				if (max_num_threads == 1) {
					return delcech(points, colours);
				} else {
					return delcech_parallel(points, colours, max_num_threads);
				}
			},
			R"docstring(Compute the chromatic Delaunay--Čech filtration of a coloured point cloud.

Args:
	points : Numpy matrix whose columns are points in the point cloud.
	colours : List or numpy array of integers describing the colours of the points.
	max_num_threads: Hint for maximum number of parallel threads to use.
		If non-positive, the number of threads to use is automatically determined
		by the threading library (Intel OneAPI TBB). Note that this may be less
		than the number of available CPU cores depending on the number of points
		and the system load.

Returns:
	The chromatic Delaunay--Čech filtration.

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
	but with Čech filtration times. Despite the different filtration values,
	it has the same persistent homology as the chromatic alpha filtration.

See Also:
	:func:`alpha`, :func:`delrips`

)docstring",
			py::arg("points"),
			py::arg("colours"),
			py::arg("max_num_threads") = 0
		);
}
