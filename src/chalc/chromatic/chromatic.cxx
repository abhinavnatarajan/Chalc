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
using chalc::Colour;
using chalc::Filtration;
using chalc::Index;
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
auto& automatic   = oneapi::tbb::task_arena::automatic;
using constraints = task_arena::constraints;

using std::array;
using std::domain_error;
using std::iota;
using std::min;
using std::numeric_limits;
using std::ranges::any_of;
using std::ranges::sort;
using std::ranges::unique;
using std::runtime_error;
using std::tie;
using std::to_string;
using std::tuple;
using std::unordered_map;
using std::vector;

// Typedefs for the CGAL Delaunay triangulation.
using Kernel                     = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
using TriangulationDataStructure = CGAL::Triangulation_data_structure<
	Kernel::Dimension,
	CGAL::Triangulation_vertex<Kernel, Index>,
	CGAL::Triangulation_full_cell<Kernel>>;
using DelaunayTriangulation = CGAL::Delaunay_triangulation<Kernel, TriangulationDataStructure>;
using Point                 = Kernel::Point_d;

// Typedefs for spatial sorting the points before the triangulation.
using PointVectorPropertyMap = boost::iterator_property_map<
	vector<Point>::const_iterator,
	boost::typed_identity_property_map<Index>,
	const Point,
	const Point&>;
using SpatialSortingTraits = CGAL::Spatial_sort_traits_adapter_d<Kernel, PointVectorPropertyMap>;

// Convert a matrix of coordinate column vectors to a vector of CGAL points.
auto matrix_columns_to_points_vec(const MatrixXd& x_arr) -> vector<Point> {
	vector<Point> points(x_arr.cols());
	auto          cols_begin = x_arr.colwise().cbegin();
	auto          cols_end   = x_arr.colwise().cend();
	for (auto&& [i, column] = tuple<Index, decltype(cols_begin)>{0, cols_begin}; column != cols_end;
	     column++, i++) {
		points[i] = Point(column->cbegin(), column->cend());
	}
	return points;
}

template <typename T> auto reorder(const vector<T>& v, const vector<Index>& idx) -> vector<T> {
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
auto canonicalise(const vector<Colour>& vec) -> tuple<vector<Colour>, Index> {
	vector<Colour>                new_vec(vec.size());
	unordered_map<Colour, Colour> m;
	for (auto&& [c, i] = tuple{vec.cbegin(), static_cast<Colour>(0)}; c != vec.cend(); c++) {
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
auto chromatic_lift(const MatrixXd& points, const vector<Colour>& colours) -> MatrixXd {
	// Make sure colours are contiguous and start at zero
	auto&& [new_colours, num_colours] = canonicalise(colours);
	// Check that the number of colours is valid.
	if (any_of(new_colours, [](const Index& colour) {
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
auto delaunay(const MatrixXd& X, const vector<Colour>& colours) -> Filtration {
	if (X.cols() > numeric_limits<Index>::max()) {
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
	vector<Index> indices(Y.cols());
	iota(indices.begin(), indices.end(), 0);  // NOLINT
	spatial_sort<Concurrency_tag>(
		indices.begin(),
		indices.end(),
		SpatialSortingTraits(PointVectorPropertyMap(points.cbegin()))
	);
	// Insert the points into the Delaunay triangulation one by one
	// and add their indexing data.
	auto                                    max_dim = static_cast<int>(Y.rows());
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
	Filtration    result(Y.cols(), dim);
	vector<Index> max_cell_vertex_labels(dim + 1);
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
	for (auto& [idx, vert]: result.simplices()[0]) {
		vert->set_colour(colours[idx]);
	}
	result.propagate_colours();
	return result;
}

// Create the chromatic Delaunay--Rips filtration.
auto delaunay_rips(const MatrixXd& points, const vector<Colour>& colours) -> Filtration {
	Filtration delX = delaunay(points, colours);

	if (delX.dimension() >= 1) {
		auto       to_double      = TypeConverter<cmb::SolutionExactType, double>{};
		const auto one_by_four    = Quotient<Mpzf>(1, 4);  // Used for the edge lengths.
		auto       points_exact_q = points.template cast<Quotient<Mpzf>>();
		for (auto& [idx, edge]: delX.simplices()[1]) {
			auto&& verts  = edge->vertex_labels();
			edge->value() = sqrt(to_double(
				(points_exact_q.col(verts[0]) - points_exact_q.col(verts[1])).squaredNorm() *
				one_by_four
			));
		}
		delX.propagate_filt_values(1, true);
	}
	return delX;
}

// Create the chromatic Delaunay--Rips filtration with parallelisation.
auto delaunay_rips_parallel(
	const MatrixXd&       points,
	const vector<Colour>& colours,
	const int             max_num_threads
) -> Filtration {
	Filtration delX = delaunay<CGAL::Parallel_tag>(points, colours);

	if (delX.dimension() >= 1) {
		auto       points_exact_q = points.template cast<Quotient<Mpzf>>();
		const auto one_by_four    = Quotient<Mpzf>(1, 4);  // Used for the edge lengths.

		task_arena arena(
			constraints{
				automatic,                                                       // numa_id
				max_num_threads == 0 ? task_arena::automatic : max_num_threads,  // max_concurrency
				automatic,                                                       // core_type
				1  // max_threads_per_core
			}
		);

		enumerable_thread_specific<TypeConverter<cmb::SolutionExactType, double>>
			thread_local_to_double;

		vector<Filtration::Simplex*> edges;
		edges.reserve(delX.simplices()[1].size());
		for (auto&& [_, edge_ptr]: delX.simplices()[1]) {
			edges.push_back(edge_ptr);
		}
		arena.execute([&] {
			parallel_for(
				blocked_range<size_t>(0, edges.size(), 10000),
				[&](const blocked_range<size_t>& r) {
					for (size_t idx = r.begin(); idx < r.end(); idx++) {
						auto&& to_double = thread_local_to_double.local();
						auto&& edge      = edges[idx];
						auto&& verts     = edge->vertex_labels();
						edge->value()    = sqrt(to_double(
                            (points_exact_q.col(verts[0]) - points_exact_q.col(verts[1]))
                                .squaredNorm() *
                            one_by_four
                        ));
					}
				}
			);
		});
		delX.propagate_filt_values(1, true);
	}
	return delX;
}

// Compute the chromatic alpha complex.
auto alpha(const MatrixXd& points, const vector<Colour>& colours) -> Filtration {
	Filtration delX(delaunay(points, colours));

	if (delX.dimension() >= 1) {
		auto to_double = TypeConverter<cmb::SolutionExactType, double>{};

		auto&& points_exact   = points.template cast<Mpzf>();
		auto&& points_exact_q = points.template cast<Quotient<Mpzf>>();

		// Partition the vertices by colour.
		// We will need this later to check if stacks are empty.
		array<vector<Index>, MAX_NUM_COLOURS> verts_by_colour;
		for (auto&& [i, colour] = tuple{static_cast<Index>(0), colours.cbegin()};
		     colour != colours.cend();
		     colour++, i++) {
			verts_by_colour.at(*colour).push_back(i);
		}
		for (auto&& p = delX.dimension(); p >= 1; p--) {
			for (auto&& [_, simplex]: delX.simplices()[p]) {
				auto&& verts = simplex->vertex_labels();

				// Partition the vertices of this simplex by colour.
				array<vector<Index>, MAX_NUM_COLOURS> verts_by_colour_in_simplex;
				for (auto&& v: verts) {
					verts_by_colour_in_simplex.at(colours[v]).push_back(v);
				}

				// Partition the vertices of all cofaces of the simplex by colour.
				// verts_by_colour_all_cofaces[j] will have duplicates but this
				// does not affect correctness.
				array<vector<Index>, MAX_NUM_COLOURS> verts_by_colour_all_cofaces;
				for (auto&& cofacet: simplex->cofacets()) {
					for (auto&& v: cofacet->vertex_labels()) {
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
					for (auto&& [j, verts_j] = tuple{0, verts_by_colour_in_simplex.cbegin()};
					     verts_j != verts_by_colour_in_simplex.cend();
					     verts_j++, j++) {
						if (verts_j->empty()) {
							continue;
						}
						// Get the radius of the j-coloured sphere in the stack.
						auto rj_squared =
							(points_exact_q(all, (*verts_j)[0]) - centre).squaredNorm();
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
						auto&& cofacets = simplex->cofacets();
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
	const MatrixXd&       points,
	const vector<Colour>& colours,
	const int             max_num_threads
) -> Filtration {
	Filtration delX(delaunay<CGAL::Parallel_tag>(points, colours));

	if (delX.dimension() >= 1) {
		task_arena arena(
			constraints{
				automatic,                                                       // numa_id
				max_num_threads == 0 ? task_arena::automatic : max_num_threads,  // max_concurrency
				automatic,                                                       // core_type
				1  // max_threads_per_core
			}
		);

		enumerable_thread_specific<TypeConverter<cmb::SolutionExactType, double>>
			thread_local_to_double;

		auto&& points_exact   = points.template cast<Mpzf>();
		auto&& points_exact_q = points.template cast<Quotient<Mpzf>>();

		// Partition the vertices by colour.
		// We will need this later to check if stacks are empty.
		array<vector<Index>, MAX_NUM_COLOURS> verts_by_colour{};
		for (auto&& [i, colour] = tuple{static_cast<Index>(0), colours.cbegin()};
		     colour != colours.cend();
		     colour++, i++) {
			verts_by_colour.at(*colour).push_back(i);
		}

		for (auto&& p = delX.dimension(); p >= 1; p--) {
			vector<Filtration::Simplex*> simplices;
			simplices.reserve(delX.simplices()[p].size());
			for (auto&& [_, simplex_ptr]: delX.simplices()[p]) {
				simplices.push_back(simplex_ptr);
			}
			arena.execute([&] {
				parallel_for(
					blocked_range<size_t>(0, simplices.size(), 500),
					[&](const blocked_range<size_t>& r) {
						auto&& to_double = thread_local_to_double.local();
						for (size_t idx = r.begin(); idx < r.end(); idx++) {
							auto&& simplex = simplices[idx];
							auto&& verts   = simplex->vertex_labels();

							// Partition the vertices of this simplex by colour.
							array<vector<Index>, MAX_NUM_COLOURS> verts_by_colour_in_simplex;
							for (auto&& v: verts) {
								verts_by_colour_in_simplex.at(colours[v]).push_back(v);
							}

							// Partition the vertices of all cofaces of the simplex by colour.
						    // verts_by_colour_all_cofaces[j] will have duplicates but this
						    // does not affect correctness.
							array<vector<Index>, MAX_NUM_COLOURS> verts_by_colour_all_cofaces;
							for (auto&& cofacet: simplex->cofacets()) {
								for (auto&& v: cofacet->vertex_labels()) {
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
							             tuple{0, verts_by_colour_in_simplex.cbegin()};
							         verts_j != verts_by_colour_in_simplex.cend();
							         j++, verts_j++) {
									if (verts_j->empty()) {
										continue;
									}
									// Get the radius of the j-coloured sphere in the stack.
									auto rj_squared =
										(points_exact_q(all, (*verts_j)[0]) - centre).squaredNorm();
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
									auto&& cofacets = simplex->cofacets();
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
	}

	return delX;
}

// Create the chromatic Delaunay-Čech complex.
auto delaunay_cech(const MatrixXd& points, const vector<Colour>& colours) -> Filtration {
	Filtration delX(delaunay(points, colours));
	auto       points_exact   = points.template cast<Mpzf>();
	auto       points_exact_q = points.template cast<Quotient<Mpzf>>();
	const auto one_by_four    = Quotient<Mpzf>(1, 4);  // Used for the edge lengths.
	auto       to_double      = TypeConverter<cmb::SolutionExactType, double>{};
	if (delX.dimension() >= 1) {
		for (auto&& p = delX.dimension(); p > 1; p--) {
			for (auto&& [_, simplex]: delX.simplices()[p]) {
				auto&& verts = simplex->vertex_labels();
				auto&& [centre, sqRadius, success] =
					miniball<SolutionPrecision::EXACT>(points_exact(all, verts));
				assert(success && "Miniball failed.");
				simplex->value() = sqrt(to_double(sqRadius));
			}
		}
		// Fast version for dimension 1.
		for (auto&& [idx, edge]: delX.simplices()[1]) {
			auto&& verts = edge->vertex_labels();
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

// Create the chromatic Delaunay--Čech complex.
auto delaunay_cech_parallel(
	const MatrixXd&       points,
	const vector<Colour>& colours,
	const int             max_num_threads
) -> Filtration {
	Filtration delX(delaunay<CGAL::Parallel_tag>(points, colours));
	auto       points_exact   = points.template cast<Mpzf>();
	auto       points_exact_q = points.template cast<Quotient<Mpzf>>();
	const auto one_by_four    = Quotient<Mpzf>(1, 4);  // Used for the edge lengths.
	if (delX.dimension() >= 1) {
		task_arena arena(
			constraints{
				automatic,                                                       // numa_id
				max_num_threads == 0 ? task_arena::automatic : max_num_threads,  // max_concurrency
				automatic,                                                       // core_type
				1  // max_threads_per_core
			}
		);

		enumerable_thread_specific<TypeConverter<cmb::SolutionExactType, double>>
			thread_local_to_double;

		vector<Filtration::Simplex*> simplices;
		simplices.reserve(delX.size());
		for (auto&& p = delX.dimension(); p >= 1; p--) {
			for (auto&& [_, simplex_ptr]: delX.simplices()[p]) {
				simplices.push_back(simplex_ptr);
			}
		}

		arena.execute([&] {
			parallel_for(
				blocked_range<size_t>(0, simplices.size(), 500),
				[&](const blocked_range<size_t>& r) {
					auto&& to_double = thread_local_to_double.local();
					for (size_t idx = r.begin(); idx < r.end(); idx++) {
						auto&& simplex = simplices[idx];
						auto&& verts   = simplex->vertex_labels();
						if (simplex->dimension() > 1) {
							auto&& [centre, sqRadius, success] =
								miniball<SolutionPrecision::EXACT>(points_exact(all, verts).eval());
							simplex->value() = sqrt(to_double(sqRadius));
						} else {
							simplex->value() = sqrt(to_double(
								(points_exact_q.col(verts[0]).eval() -
						         points_exact_q.col(verts[1]).eval())
									.squaredNorm() *
								one_by_four
							));
						}
					}
				}
			);
		});
	}

	return delX;
}

}  // namespace chalc
