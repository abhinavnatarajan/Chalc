/*
    This file is part of Chalc.

    Chalc: Chromatic Alpha Complexes.
    Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha
   Complexes”. Online preprint available at http://arxiv.org/abs/2212.03128.
    Accessed: 2023-02-28 22:07:23 UTC.
    DOI: 10.48550/arXiv.2212.03128.

    Project homepage:    http://abhinavnatarajan.github.io/Chalc

    Copyright (c) 2023 Abhinav Natarajan

    Contributors:
    Abhinav Natarajan

    Licensing:
    Chalc is released under the GNU General Public License ("GPL").

    GNU General Public License ("GPL") copyright permissions statement:
    **************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include "chromatic.h"
#include "chalc/filtration/filtration.h"
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Dimension.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Mpzf.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/tags.h>
#include <ConstrainedMiniball/cmb.hpp>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <array>
#include <boost/property_map/property_map.hpp>
#include <numeric>
#include <oneapi/tbb.h>
#include <ranges>
#include <stdexcept>
#include <string>

namespace {
using chalc::colour_t;
using chalc::Filtration;
using chalc::index_t;
using chalc::MAX_NUM_COLOURS;

using CGAL::Mpzf;
using CGAL::Quotient;
using CGAL::spatial_sort;

using cmb::constrained_miniball;
using cmb::miniball;
using cmb::SolutionPrecision;
using cmb::utility::equidistant_subspace;
using cmb::utility::TypeConverter;

using Eigen::all;
using Eigen::lastN;
using Eigen::MatrixX;
using Eigen::MatrixXd;
using Eigen::VectorX;

using oneapi::tbb::blocked_range;
using oneapi::tbb::enumerable_thread_specific;
using oneapi::tbb::parallel_for;
using oneapi::tbb::task_arena;

using std::array;
using std::domain_error;
using std::iota;
using std::map;
using std::min;
using std::numeric_limits;
using std::ranges::any_of;
using std::ranges::sort;
using std::ranges::unique;
using std::runtime_error;
using std::shared_ptr;
using std::tie;
using std::to_string;
using std::tuple;
using std::vector;

// Typedefs for the CGAL Delaunay triangulation.
using Kernel_d                     = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using Triangulation_data_structure = CGAL::Triangulation_data_structure<
	Kernel_d::Dimension,
	CGAL::Triangulation_vertex<Kernel_d, index_t>,
	CGAL::Triangulation_full_cell<Kernel_d>>;
using DelaunayTriangulation = CGAL::Delaunay_triangulation<Kernel_d, Triangulation_data_structure>;
using Point_d               = Kernel_d::Point_d;

// Typedefs for spatial sorting the points before the triangulation.
using PointVectorPropertyMap = boost::iterator_property_map<
	vector<Point_d>::const_iterator,
	boost::typed_identity_property_map<index_t>,
	const Point_d,
	const Point_d&>;
using SpatialSortingTraits_d =
	CGAL::Spatial_sort_traits_adapter_d<Kernel_d, PointVectorPropertyMap>;

// Convert a matrix of coordinate column vectors to a vector of CGAL points.
auto matrix_columns_to_points_vec(const MatrixXd& x_arr) -> vector<Point_d> {
	vector<Point_d> points(x_arr.cols());
	auto            cols_begin = x_arr.colwise().cbegin();
	auto            cols_end   = x_arr.colwise().cend();
	for (auto&& [i, column] = tuple<index_t, decltype(cols_begin)>{0, cols_begin};
	     column != cols_end;
	     column++, i++) {
		points[i] = Point_d(column->cbegin(), column->cend());
	}
	return points;
}

template <typename T> auto reorder(const vector<T>& v, const vector<index_t>& idx) -> vector<T> {
	vector<T> result(v.size());
	for (auto i = 0; i < v.size(); i++) {
		result[i] = v[idx[i]];
	}
	return result;
}

template <typename T> void remove_duplicates_inplace(vector<T>& vec) {
	sort(vec.begin(), vec.end());
	const auto [last, end] = unique(vec);
	vec.erase(last, end);
}

template <typename T> auto compare(const void* a, const void* b) -> int {
	const auto& arg1 = *(static_cast<const T*>(a));
	const auto& arg2 = *(static_cast<const T*>(b));
	return arg1 < arg2 ? -1 : arg1 > arg2 ? +1 : 0;
}

// Make a vector contiguous and start at zero.
// Returns new vector and number of distinct elements.
auto canonicalise(const vector<colour_t>& vec) -> tuple<vector<colour_t>, index_t> {
	vector<colour_t>        new_vec(vec.size());
	map<colour_t, colour_t> m;
	for (auto&& [c, i] = tuple{vec.begin(), static_cast<colour_t>(0)}; c != vec.end(); c++) {
		if (!m.contains(*c)) {
			m[*c] = i++;
		}
	}
	for (size_t i = 0; i < vec.size(); i++) {
		new_vec[i] = m[vec[i]];
	}
	return tuple{new_vec, m.size()};
}

// Stratify a coloured point set.
// Points are provided as columns of a matrix or matrix expression.
// Colours are provided as a vector.
auto chromatic_lift(const MatrixXd& points, const vector<colour_t>& colours) -> MatrixXd {
	// Make sure colours are contiguous and start at zero
	auto&& [new_colours, num_colours] = canonicalise(colours);
	// Check that the number of colours is valid.
	if (any_of(new_colours, [](const index_t& colour) {
			return (colour >= MAX_NUM_COLOURS);
		})) {
		throw domain_error(
			"Too many colours; at most " + to_string(MAX_NUM_COLOURS - 1) +
			" colours are supported."
		);
	}
	auto     dim = points.rows();
	MatrixXd result(dim + num_colours - 1, points.cols());
	result.topRows(dim) = points;
	if (num_colours != 1) {
		result.bottomRows(num_colours - 1).setZero();
		for (auto i = 0; i < result.cols(); i++) {
			if (new_colours[i] != 0) {
				result((dim + new_colours[i]) - 1, i) = 1.0;
			}
		}
	}
	return result;
}

}  // namespace

namespace chalc {

// Create a Delaunay triangulation from a collection of coordinate vectors.
template <typename Concurrency_tag>
auto delaunay(const MatrixXd& X, const vector<colour_t>& colours) -> Filtration {
	if (X.cols() > numeric_limits<index_t>::max()) {
		throw runtime_error("Number of points is too large.");
	}
	if (colours.size() != X.cols()) {
		throw domain_error("len(colours) must equal number of points.");
	}
	if (X.cols() == 0) {
		return Filtration{0, 0};
	}

	// Chromatic lift of the point cloud.
	MatrixXd Y(chromatic_lift(X, colours));
	if (Y.rows() > static_cast<Eigen::Index>(numeric_limits<int>::max())) {
		throw runtime_error("Dimension of stratified points is too large.");
	}
	auto points = matrix_columns_to_points_vec(Y);
	// Spatial sorting for the points.
	// This makes the insertion O(nlogn) instead of O(n^{ceil(d/2) + 1}).
	vector<index_t> indices(Y.cols());
	iota(indices.begin(), indices.end(), 0);
	spatial_sort<Concurrency_tag>(
		indices.begin(),
		indices.end(),
		SpatialSortingTraits_d(PointVectorPropertyMap(points.cbegin()))
	);
	// Insert the points into the Delaunay triangulation one by one
	// and add their indexing data.
	int                                     max_dim = Y.rows();  // NOLINT
	auto                                    delY    = DelaunayTriangulation(max_dim);
	DelaunayTriangulation::Full_cell_handle hint;
	DelaunayTriangulation::Vertex_handle    v;
	for (auto&& idx: indices) {
		v         = delY.insert(points[idx], hint);
		v->data() = idx;
		hint      = v->full_cell();
	}

	// Initialise the filtration with the actual dimension of the triangulation.
	// This way we have the correct dimension even if there are degenerate cases.
	auto dim = delY.current_dimension();
	// dim = -2 -> empty triangulation, not possible in our case
	// dim = -1 -> only one point, corner case
	// dim = 0 -> only two points but no 1-simplex, impossible in our case
	// dim > 0 -> at least one edge
	if (dim < 0) {
		dim = 0;
	}
	Filtration      result(Y.cols(), dim);
	vector<index_t> max_cell_vertex_labels(dim + 1);
	for (auto cell_it = delY.finite_full_cells_begin(); cell_it != delY.finite_full_cells_end();
	     cell_it++) {
		// Iterate over the vertices of the cell and get their labels.
		for (auto&& [label_it, vert_it] =
		         tuple{max_cell_vertex_labels.begin(), cell_it->vertices_begin()};
		     // vert_it != cell_it->vertices_end();
		     // vert_it++) { // BUG IN CGAL: vertices_end segfaults for cells with less vertices
		     // than maximal_dimension
		     vert_it != cell_it->vertices_end() && label_it != max_cell_vertex_labels.end();
		     label_it++, vert_it++) {
			*label_it = (*vert_it)->data();
		}
		result.add_simplex(max_cell_vertex_labels, 0.0);
	}
	// Modify the colours of the vertices.
	for (auto& [idx, vert]: result.get_simplices()[0]) {
		vert->set_colour(colours[idx]);
	}
	result.propagate_colours();
	return result;
}

// Create the chromatic Delaunay--Rips filtration.
auto delrips(const MatrixXd& points, const vector<colour_t>& colours) -> Filtration {
	// Get the delaunay triangulation
	Filtration delX = delaunay(points, colours);

	// Modify the filtration values.
	if (delX.dimension() >= 1) {
		auto points_exact = points.template cast<Mpzf>();
		for (auto& [idx, edge]: delX.get_simplices()[1]) {
			auto&& verts = edge->get_vertex_labels();
			edge->value() =
				((points_exact.col(verts[0]) - points_exact.col(verts[1])).norm() * Mpzf(0.5))
					.to_double();
		}
		delX.propagate_filt_values(1, true);
	}
	return delX;
}

// Create the chromatic Delaunay--Rips filtration with parallelisation.
auto delrips_parallel(
	const MatrixXd&         points,
	const vector<colour_t>& colours,
	const int               max_num_threads
) -> Filtration {
	// Get the delaunay triangulation.
	Filtration delX = delaunay<CGAL::Parallel_tag>(points, colours);

	// Modify the filtration values.
	if (delX.dimension() >= 1) {
		auto       points_exact = points.template cast<Mpzf>();
		task_arena arena(max_num_threads == 0 ? task_arena::automatic : max_num_threads);
		// Store all the simplex pointers in a vector,
		// which is convenient to iterate over using the TBB API.
		vector<const shared_ptr<Filtration::Simplex>*> edges;
		edges.reserve(delX.get_simplices()[1].size());
		for (auto&& [_, edge]: delX.get_simplices()[1]) {
			edges.push_back(&edge);
		}
		arena.execute([&] {
			parallel_for(
				blocked_range<size_t>(0, edges.size(), 10000),
				[&](const blocked_range<size_t>& r) {
					for (size_t idx = r.begin(); idx < r.end(); idx++) {
						auto&& edge  = *edges[idx];
						auto&& verts = edge->get_vertex_labels();
						edge->value() =
							((points_exact.col(verts[0]) - points_exact.col(verts[1])).norm() *
					         Mpzf(0.5))
								.to_double();
					}
				}
			);
		});
		delX.propagate_filt_values(1, true);
	}
	return delX;
}

// Compute the chromatic alpha complex.
auto alpha(const MatrixXd& points, const vector<colour_t>& colours) -> Filtration {
	// Get the delaunay triangulation.
	Filtration delX(delaunay(points, colours));

	// Partition the vertices by colour.
	// We will need this later to check if stacks are empty.
	array<vector<index_t>, MAX_NUM_COLOURS> verts_by_colour;
	for (auto&& [i, colour] = tuple{static_cast<index_t>(0), colours.begin()};
	     colour != colours.end();
	     colour++, i++) {
		verts_by_colour.at(*colour).push_back(i);
	}
	if (delX.dimension() >= 1) {
		auto to_double = TypeConverter<cmb::SolutionExactType, double>{};

		auto&& points_exact   = points.template cast<Mpzf>();
		auto&& points_exact_q = points.template cast<Quotient<Mpzf>>();

		// Start at the current dimension.
		for (auto&& p = delX.dimension(); p >= 1; p--) {
			// Iterate over p-simplices
			for (auto&& [_, simplex]: delX.get_simplices()[p]) {
				auto&& verts = simplex->get_vertex_labels();

				// Partition the vertices of this simplex by colour.
				array<vector<index_t>, MAX_NUM_COLOURS> verts_by_colour_in_simplex;
				for (auto&& v: verts) {
					verts_by_colour_in_simplex.at(colours[v]).push_back(v);
				}

				// Partition the vertices of all cofaces of the simplex by colour.
				// verts_by_colour_all_cofaces[j] will have duplicates but this
				// does not affect correctness.
				array<vector<index_t>, MAX_NUM_COLOURS> verts_by_colour_all_cofaces;
				for (auto&& cofacet: simplex->get_cofacets()) {
					for (auto&& v: cofacet->get_vertex_labels()) {
						verts_by_colour_all_cofaces.at(colours[v]).push_back(v);
					}
				}

				// For each colour j in the simplex, find the affine subspace
				// of points equidistant to the vertices of colour j in the simplex.
				// Suppose this subspace is E_j and is defined by a matrix equation:
				// E_j * x = b_j.
				// We append E_j and b_j to the bottom of a matrix E and a vector b.
				// Then the solution set E of the matrix equation
				// E * x = b
				// is the intersection over all j of the affine subspaces E_j.
				MatrixX<Mpzf> E(0, points.rows());
				VectorX<Mpzf> b;
				for (auto&& verts_j: verts_by_colour_in_simplex) {
					if (verts_j.empty()) {
						continue;
					}
					// Each of these has size at least 1, so the -1 in the following is
					// safe. The narrowing cast is safe because verts_j.size() is always
					// at most the number of points, which is less than
					// Eigen::Index::max().
					auto num_new_rows = static_cast<Eigen::Index>(verts_j.size() - 1);
					E.conservativeResize(E.rows() + num_new_rows, Eigen::NoChange);
					b.conservativeResize(b.rows() + num_new_rows);
					MatrixX<Mpzf>::BlockXpr           E_new_rows = E.bottomRows(num_new_rows);
					Eigen::VectorBlock<VectorX<Mpzf>> b_new_rows = b(lastN(num_new_rows));
					tie(E_new_rows, b_new_rows) = equidistant_subspace(points_exact(all, verts_j));
				}

				// Get the smallest bounding ball of the points
				// in the simplex, with the added constraint that
				// the centre x of the ball must satisfy the equation
				// E * x = b i.e., it lies in the affine subspace E.
				auto&& [centre, sqRadius, success] =
					constrained_miniball<SolutionPrecision::EXACT>(points_exact(all, verts), E, b);
				assert(success && "Constrained miniball failed.");
				bool stack_is_empty = true;
				if (p == delX.dimension()) {
					// For maximal simplices there is nothing more to do.
					simplex->value() = sqrt(to_double(sqRadius));
				} else {
					// If the simplex is not maximal, check if the stack is empty.
					for (auto&& [j, verts_j] = tuple{0, verts_by_colour_in_simplex.begin()};
					     verts_j != verts_by_colour_in_simplex.end();
					     verts_j++, j++) {
						if (verts_j->empty()) {
							continue;
						}
						// Get the radius of the j-coloured sphere in the stack.
						auto rj_squared = (points_exact_q(all, verts_j[0]) - centre).squaredNorm();
						// Get the distance of the nearest point of colour j to the centre,
						// among all vertices in cofaces of this simplex.
						auto squared_dist_to_nearest_pt_of_colour_j =
							(points_exact_q(all, verts_by_colour_all_cofaces.at(j)).colwise() -
						     centre)
								.colwise()
								.squaredNorm()
								.minCoeff();
						// Check if the nearest point is not in the interior of the sphere.
						stack_is_empty &= (squared_dist_to_nearest_pt_of_colour_j >= rj_squared);
					}
					// If the stack is empty, assign the filtration value.
					if (stack_is_empty) {
						simplex->value() = sqrt(to_double(sqRadius));
					} else {
						auto&& cofacets = simplex->get_cofacets();
						if (cofacets.size() == 0) {
							continue;
						}
						simplex->value() = cofacets[0]->value();
						for (auto&& cofacet: cofacets) {
							simplex->value() = min(simplex->value(), cofacet->value());
						}
					}
				}
			}
		}
	}
	return delX;
}

// Construct the chromatic alpha filtration with parallelisation.
auto alpha_parallel(
	const MatrixXd&         points,
	const vector<colour_t>& colours,
	const int               max_num_threads
) -> Filtration {
	// Get the delaunay triangulation.
	Filtration delX(delaunay<CGAL::Parallel_tag>(points, colours));

	// Partition the vertices by colour.
	// We will need this later to check if stacks are empty.
	array<vector<index_t>, MAX_NUM_COLOURS> verts_by_colour{};
	for (auto&& [i, colour] = tuple{static_cast<index_t>(0), colours.begin()};
	     colour != colours.end();
	     colour++, i++) {
		verts_by_colour.at(*colour).push_back(i);
	}
	if (delX.dimension() >= 1) {
		task_arena arena(max_num_threads == 0 ? task_arena::automatic : max_num_threads);
		// Modify the filtration values
		auto&& points_exact   = points.template cast<Mpzf>();
		auto&& points_exact_q = points.template cast<Quotient<Mpzf>>();

		// Each worker thread will get its own copy of the quotient-to-double approximator.
		enumerable_thread_specific<TypeConverter<cmb::SolutionExactType, double>>
			thread_local_to_double;
		// Start at the current dimension.
		for (auto&& p = delX.dimension(); p >= 1; p--) {
			// Store all the simplex pointers in a vector,
			// which is convenient to iterate over using the TBB API.
			vector<const shared_ptr<Filtration::Simplex>*> simplices;
			simplices.reserve(delX.get_simplices()[p].size());
			for (auto&& [_, simplex]: delX.get_simplices()[p]) {
				simplices.push_back(&simplex);
			}
			arena.execute([&] {
				parallel_for(
					blocked_range<size_t>(0, simplices.size(), 500),
					[&](const blocked_range<size_t>& r) {
						auto&& to_double = thread_local_to_double.local();
						for (size_t idx = r.begin(); idx < r.end(); idx++) {
							auto   simplex = *simplices[idx];
							auto&& verts   = simplex->get_vertex_labels();

							// Partition the vertices of this simplex by colour.
							array<vector<index_t>, MAX_NUM_COLOURS> verts_by_colour_in_simplex;
							for (auto&& v: verts) {
								verts_by_colour_in_simplex.at(colours[v]).push_back(v);
							}

							// Partition the vertices of all cofaces of the simplex by colour.
						    // verts_by_colour_all_cofaces[j] will have duplicates but this
						    // does not affect correctness.
							array<vector<index_t>, MAX_NUM_COLOURS> verts_by_colour_all_cofaces;
							for (auto&& cofacet: simplex->get_cofacets()) {
								for (auto&& v: cofacet->get_vertex_labels()) {
									verts_by_colour_all_cofaces.at(colours[v]).push_back(v);
								}
							}

							// For each colour j in the simplex, find the affine subspace
						    // of points equidistant to the vertices of colour j in the simplex.
						    // Suppose this subspace is E_j and is defined by a matrix equation:
						    // E_j * x = b_j.
						    // We append E_j and b_j to the bottom of a matrix E and a vector b.
						    // Then the solution set E of the matrix equation
						    // E * x = b
						    // is the intersection over all j of the affine subspaces E_j.
							MatrixX<Mpzf> E(0, points.rows());
							VectorX<Mpzf> b;
							for (auto&& verts_j: verts_by_colour_in_simplex) {
								if (verts_j.empty()) {
									continue;
								}
								// Each of these has size at least 1, so the -1 in the following is
							    // safe. The narrowing cast is safe because verts_j.size() is always
							    // at most the number of points, which is less than
							    // Eigen::Index::max().
								auto num_new_rows = static_cast<Eigen::Index>(verts_j.size() - 1);
								E.conservativeResize(E.rows() + num_new_rows, Eigen::NoChange);
								b.conservativeResize(b.rows() + num_new_rows);
								MatrixX<Mpzf>::BlockXpr E_new_rows = E.bottomRows(num_new_rows);
								Eigen::VectorBlock<VectorX<Mpzf>> b_new_rows =
									b(lastN(num_new_rows));
								tie(E_new_rows, b_new_rows) =
									equidistant_subspace(points_exact(all, verts_j));
							}

							// Get the smallest bounding ball of the points
						    // in the simplex, with the added constraint that
						    // the centre x of the ball must satisfy the equation
						    // E * x = b i.e., it lies in the affine subspace E.
							auto&& [centre, sqRadius, success] =
								constrained_miniball<SolutionPrecision::EXACT>(
									points_exact(all, verts),
									E,
									b
								);
							bool stack_is_empty = true;
							if (p == delX.dimension()) {
								// For maximal simplices there is nothing more to do
								simplex->value() = sqrt(to_double(sqRadius));
							} else {
								// If the simplex is not maximal, check if the stack is empty
								for (auto&& [j, verts_j] =
							             tuple{0, verts_by_colour_in_simplex.begin()};
							         verts_j != verts_by_colour_in_simplex.end();
							         j++, verts_j++) {
									if (verts_j->empty()) {
										continue;
									}
									// Get the radius of the j-coloured sphere in the stack.
									auto rj_squared =
										(points_exact_q(all, verts_j[0]) - centre).squaredNorm();
									// Get the distance of the nearest point of colour j to the
								    // centre, among all vertices in cofaces of this simplex.
									auto squared_dist_to_nearest_pt_of_colour_j =
										(points_exact_q(all, verts_by_colour_all_cofaces.at(j))
								             .colwise() -
								         centre)
											.colwise()
											.squaredNorm()
											.minCoeff();
									// Check if the nearest point is not in the interior of the
								    // sphere.
									stack_is_empty &=
										(squared_dist_to_nearest_pt_of_colour_j >= rj_squared);
								}
								// If the stack is empty, assign the filtration value.
								if (stack_is_empty) {
									simplex->value() = sqrt(to_double(sqRadius));
								} else {
									auto&& cofacets = simplex->get_cofacets();
									if (cofacets.size() == 0) {
										continue;
									}
									simplex->value() = cofacets[0]->value();
									for (auto&& cofacet: cofacets) {
										simplex->value() = min(simplex->value(), cofacet->value());
									}
								}
							}
						}
					}
				);
			});
		}
		thread_local_to_double.clear();  // Clear the thread-local storage.
	}

	return delX;
}

// Create the chromatic Delaunay-Cech complex.
auto delcech(const MatrixXd& points, const vector<colour_t>& colours) -> Filtration {
	// Get the delaunay triangulation.
	Filtration delX(delaunay(points, colours));
	// Modify the filtration values.
	auto       points_exact_q = points.template cast<Quotient<Mpzf>>();
	const auto one_by_four    = Quotient<Mpzf>(Mpzf(0.25));  // Used for the edge lengths.
	auto       to_double      = TypeConverter<cmb::SolutionExactType, double>{};
	if (delX.dimension() >= 1) {
		for (auto&& p = delX.dimension(); p > 1; p--) {
			for (auto&& [_, simplex]: delX.get_simplices()[p]) {
				auto&& verts = simplex->get_vertex_labels();
				auto&& [centre, sqRadius, success] =
					miniball<SolutionPrecision::EXACT>(points(all, verts));
				assert(success && "Miniball failed.");
				simplex->value() = sqrt(to_double(sqRadius));
			}
		}
		// Fast version for dimension 1.
		for (auto&& [idx, edge]: delX.get_simplices()[1]) {
			auto&& verts = edge->get_vertex_labels();
			// We use exact types here for consistency with the calculations in higher
			// dimensions.
			edge->value() = sqrt(to_double(
				(points_exact_q.col(verts[0]) - points_exact_q.col(verts[1])).squaredNorm() *
				one_by_four
			));
		}
	}
	return delX;
}

// Create the chromatic Delaunay--Cech complex.
auto delcech_parallel(
	const MatrixXd&         points,
	const vector<colour_t>& colours,
	const int               max_num_threads
) -> Filtration {
	// Start
	// Get the delaunay triangulation.
	Filtration delX(delaunay<CGAL::Parallel_tag>(points, colours));
	// Modify the filtration values.
	auto       points_exact_q = points.template cast<Quotient<Mpzf>>();
	const auto one_by_four = Quotient<Mpzf>(Mpzf(0.25), Mpzf(1.0));  // Used for the edge lengths.
	if (delX.dimension() >= 1) {
		task_arena arena(max_num_threads == 0 ? task_arena::automatic : max_num_threads);
		// Each worker thread will get its own copy of the quotient-to-double approximator.
		enumerable_thread_specific<TypeConverter<cmb::SolutionExactType, double>>
			thread_local_to_double;
		for (auto&& p = delX.dimension(); p > 1; p--) {
			// Store all the simplex pointers in a vector,
			// which is convenient to iterate over using the TBB API.
			vector<const shared_ptr<Filtration::Simplex>*> simplices;
			simplices.reserve(delX.get_simplices()[p].size());
			for (auto&& [_, simplex]: delX.get_simplices()[p]) {
				simplices.push_back(&simplex);
			}
			arena.execute([&] {
				parallel_for(
					blocked_range<size_t>(0, simplices.size(), 500),
					[&](const blocked_range<size_t>& r) {
						auto&& to_double = thread_local_to_double.local();
						for (size_t idx = r.begin(); idx < r.end(); idx++) {
							auto   simplex = *simplices[idx];
							auto&& verts   = simplex->get_vertex_labels();
							auto&& [centre, sqRadius, success] =
								miniball<SolutionPrecision::EXACT>(points(all, verts));
							simplex->value() = sqrt(to_double(sqRadius));
						}
					}
				);
			});
		}
		// Fast version for dimension 1.
		vector<const shared_ptr<Filtration::Simplex>*> edges;
		edges.reserve(delX.get_simplices()[1].size());
		for (auto&& [_, edge]: delX.get_simplices()[1]) {
			edges.push_back(&edge);
		}
		arena.execute([&] {
			parallel_for(
				blocked_range<size_t>(0, edges.size(), 500),
				[&](const blocked_range<size_t>& r) {
					auto&& to_double = thread_local_to_double.local();
					for (size_t idx = r.begin(); idx < r.end(); idx++) {
						auto&& edge  = *edges[idx];
						auto&& verts = edge->get_vertex_labels();
						// We use exact types here for consistency with the calculations in higher
					    // dimensions.
						edge->value() = sqrt(to_double(
							(points_exact_q.col(verts[0]) - points_exact_q.col(verts[1]))
								.squaredNorm() *
							one_by_four
						));
					}
				}
			);
		});
		thread_local_to_double.clear();  // Clear the thread-local storage.
	}
	return delX;
}

}  // namespace chalc
