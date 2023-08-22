/*
    This file is part of Chalc.

    Chalc: Chromatic Alpha Complexes.
    Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha Complexes”. 
    Online preprint available at http://arxiv.org/abs/2212.03128. 
    Accessed: 2023-02-28 22:07:23 UTC. 
    DOI: 10.48550/arXiv.2212.03128.

    Project homepage:    http://github.com/abhinavnatarajan/Chalc

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
/*
DO NOT MOVE THE NEXT INCLUDE.
ConstrainedMiniball.h includes (by proxy) mpreal.h
mpreal.h defines a macro MPFR_USE_NO_MACRO before including 
mpfr.h, without which there are compilation issues. Therefore 
mpreal.h must be included BEFORE any CGAL headers, which also
include mpfr.h.
*/
#include "ConstrainedMiniball/ConstrainedMiniball.h"

#include <CGAL/Epick_d.h>
#include <CGAL/Dimension.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Triangulation.h>
#include <algorithm>
#include <iostream>

namespace chalc {
    using Kernel_d = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
    using Triangulation_data_structure = CGAL::Triangulation_data_structure< Kernel_d::Dimension,
        CGAL::Triangulation_vertex<Kernel_d, size_t>,
        CGAL::Triangulation_full_cell<Kernel_d>>;
    using DelaunayTriangulation = CGAL::Delaunay_triangulation<Kernel_d, Triangulation_data_structure>;
    using Point_d = Kernel_d::Point_d;
    
    template <class Real_t>
    using RealVector = Eigen::Matrix<Real_t, Eigen::Dynamic, 1>;
    
    template <class Real_t>
    using RealMatrix = Eigen::Matrix<Real_t, Eigen::Dynamic, Eigen::Dynamic>;

    // Convert a collection of coordinate vectors to a vector of CGAL points
    vector<Point_d> coordvecs_to_points(const RealMatrix<double>& x_arr) {
        vector<Point_d> points(x_arr.cols());
        for (size_t i = 0; i < x_arr.cols(); i++) {
            points[i] = Point_d(x_arr.col(i).begin(), x_arr.col(i).end());
        }
        return points;
    }

    template <typename T>
    vector<T> reorder(const vector<T>& v, const vector<size_t>& idx) {
        vector<T> result(v.size());
        for (auto i = 0; i < v.size(); i++) {
            result[i] = v[idx[i]];
        }
        return result;
    }

    template <typename T>
    tuple<vector<T>, vector<size_t>> sort_with_indices(const vector<T>& v, bool (*compare)(const T& a, const T& b)) {
        vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(),
            [&v, &compare](size_t i1, size_t i2) { return compare(v[i1], v[i2]); });
        return tuple{ reorder(v, idx), idx };
    }

    template <typename T>
    int compare(const void* a, const void* b)
    {
        const auto& arg1 = *(static_cast<const T*>(a));
        const auto& arg2 = *(static_cast<const T*>(b));
        return arg1 < arg2 ? -1
            : arg1 > arg2 ? +1
            : 0;
    }

    // make a vector contiguous and start at zero
    // returns new vector and number of distinct elements
    tuple<vector<size_t>, size_t> canonicalise(const vector<size_t>& vec) {
        vector<size_t> new_vec(vec.size());
        map<size_t, size_t> m;
        for (auto [c, i] = tuple{vec.begin(), 0}; c != vec.end(); c++) {
            if (!m.contains(*c)) {
                m[*c] = i++;
            }
        }
        for (size_t i = 0; i < vec.size(); i++) {
            new_vec[i] = m[vec[i]];
        }
        return tuple{new_vec, m.size()};
    }


    // Stratify a coloured point set
    // Points are provided as columns of a matrix
    // Colours are provided as a vector
    RealMatrix<double> stratify(const RealMatrix<double>& points, const vector<size_t>& colours, size_t num_colours) {
        size_t dim = points.rows();
        vector<size_t> new_colours;
        // make sure colours are contiguous and start at zero
        if (num_colours == -1) {
            std::tie(new_colours, num_colours) = canonicalise(colours);
        }
        else {
            new_colours = colours;
        }
        RealMatrix<double> result(dim + num_colours - 1, points.cols());
        result.topRows(dim) = points;
        //vector<size_t> new_colours(colours);
        if (num_colours != 1) {
            result.bottomRows(num_colours - 1).setZero();
            for (auto i = 0; i < result.cols(); i++) {
                if (new_colours[i] != 0) {
                    result(dim - 1 + new_colours[i], i) = 1.0;
                }
            }
        }
        return result;
    }

    // Create a Delaunay triangulation from a collection of coordinate vectors
    FilteredComplex delaunay_complex(const RealMatrix<double>& X) {
        auto dim = X.rows();
        FilteredComplex result(X.cols(), dim);
        if (X.cols() != 0) {
            auto points = coordvecs_to_points(X);
            auto [sorted_points, sorted_indices] = sort_with_indices<Point_d>(points, [](const Point_d& a, const Point_d& b) -> bool { return a < b; });
            auto delX = DelaunayTriangulation(dim);
            delX.insert(points.begin(), points.end());
            // iterate over all finite vertices and associate them with their original label
            for (auto vert_it = delX.finite_vertices_begin(); vert_it != delX.finite_vertices_end(); vert_it++) {
                // find the index of its associated point
                auto point = vert_it->point();
                const Point_d* p = static_cast<Point_d*>(
                    std::bsearch(&point, sorted_points.data(), sorted_points.size(), sizeof(Point_d), compare<Point_d>)
                    );
                vert_it->data() = sorted_indices[p - sorted_points.data()];
            }
            // iterate over the top dimensional cells and add them to the filtration
            vector<size_t> max_cell_vertex_labels(dim + 1);
            for (auto cell_it = delX.finite_full_cells_begin();
                cell_it != delX.finite_full_cells_end();
                cell_it++) {
                // iterate over the vertices of the cell and get their labels
                for (auto [label_it, vert_it] = tuple{ max_cell_vertex_labels.begin(), cell_it->vertices_begin() };
                    vert_it != cell_it->vertices_end();
                    label_it++, vert_it++) {
                    *label_it = (*vert_it)->data();
                }
                result.add_simplex(max_cell_vertex_labels, 0.0);
            }
        }
        return result;
    }

    // Create the weak chromatic alpha complex
    FilteredComplex weak_chromatic_alpha_complex(const RealMatrix<double>& points, const vector<size_t>& colours) {
        auto [new_colours, num_colours] = canonicalise(colours);
        RealMatrix<double> stratified_points = stratify(points, new_colours, num_colours);
        auto delX = delaunay_complex(stratified_points);
        // modify the colours of the vertices
        for (auto& [idx, vert] : delX.get_simplices()[0]) {
            vert->colours.set(new_colours[idx]);
        }
        delX.propagate_colours();
        // modify the filtration values
        if (delX.max_dim >= 1) { 
            for (auto& [idx, edge] : delX.get_simplices()[1]) {
                auto verts = edge->get_vertex_labels();
                edge->value = (points.col(verts[0]) - points.col(verts[1])).norm();
            }
            delX.propagate_filt_values(1, true);
        }
        return delX;
    }

    template <class Derived>
    tuple<RealMatrix<typename Derived::Scalar>, RealVector<typename Derived::Scalar>> equidistant_subspace(const Eigen::MatrixBase<Derived>& X) {
        int n = X.cols();
        typedef Derived::Scalar Real_t;
        RealMatrix<Real_t> E(n-1, X.rows());
        RealVector<Real_t> b(n-1);
        if (n > 1) {
            b = (X.rightCols(n-1).colwise().squaredNorm().array() - X.col(0).squaredNorm()).transpose();
            E = (X.rightCols(n-1).colwise() - X.col(0)).transpose();
        }
        return tuple{E, b};
    }

    // Create the chromatic alpha complex
    FilteredComplex chromatic_alpha_complex(const RealMatrix<double>& points, const vector<size_t>& colours) {
        typedef double Real_t;
        auto tol = Eigen::NumTraits<double>::dummy_precision();
        auto [new_colours, num_colours] = canonicalise(colours);
        RealMatrix<double> stratified_points = stratify(points, new_colours, num_colours);
        auto delX = delaunay_complex(stratified_points);
        // modify the colours of the vertices
        for (auto& [idx, vert] : delX.get_simplices()[0]) {
            vert->colours.set(new_colours[idx]);
        }
        delX.propagate_colours();
        // sort the points by colour
        map<size_t, vector<size_t>> verts_by_colour;
        for (size_t i = 0; i < points.cols(); i++) {
            verts_by_colour[new_colours[i]].push_back(i);
        }
        // modify the filtration values
        bool numerical_issues = false;
        if (delX.max_dim >= 1) {
            for (int p = delX.max_dim; p >= 1; p--) {
                for (auto& [idx, simplex] : delX.get_simplices()[p]) {
                    auto verts = simplex->get_vertex_labels();
                    map<size_t, vector<size_t>> verts_by_colour_in_simplex;
                    RealMatrix<Real_t> E(0, points.rows());
                    RealVector<Real_t> b;
                    for (auto v : verts) {
                        verts_by_colour_in_simplex[new_colours[v]].push_back(v);
                    }
                    for (auto& [j, verts_j] : verts_by_colour_in_simplex) {
                        size_t num_new_rows = verts_j.size()-1;
                        E.conservativeResize(E.rows() + num_new_rows, Eigen::NoChange);
                        b.conservativeResize(b.rows() + num_new_rows);
                        RealMatrix<Real_t>::BlockXpr E_new_rows = E.bottomRows(num_new_rows);
                        Eigen::VectorBlock<RealVector<Real_t>> b_new_rows = b(Eigen::lastN(num_new_rows));
                        std::tie(E_new_rows, b_new_rows) = equidistant_subspace(points(Eigen::all, verts_j));
                    }
                    auto [centre, sqRadius, success] = cmb::constrained_miniball<double>(points.rows(), points(Eigen::all, verts), E, b);
                    numerical_issues |= !success;
                    #ifndef NDEBUG
                        if (!success) {
                            std::cerr << "Error : dim = " << p << ", simplexID = " << idx << ", squared radius : " << sqRadius << std::endl;
                        }
                    #endif
                    // check if the stack is empty
                    bool stack_is_empty = true;
                    for (auto& [j, verts_j] : verts_by_colour_in_simplex) {
                        // radius of sphere of colour j
                        double rj = (points(Eigen::all, verts_j).colwise() - centre).colwise().squaredNorm().minCoeff();
                        double distance_to_nearest_point_of_colour_j = (points(Eigen::all, verts_by_colour[j]).colwise() - centre).colwise().squaredNorm().minCoeff();
                        stack_is_empty &= distance_to_nearest_point_of_colour_j - rj >= -tol;
                    }
                    if (stack_is_empty) {
                        simplex->value = sqrt(sqRadius);
                    }
                }
                delX.propagate_filt_values(p, false);
            }            
        }
        if (numerical_issues) {
            std::cerr << "Warning: encountered numerical problems. Filtration values may be inaccurate." << std::endl;
        }
        return delX;
    }

    // Create the chromatic Del-Cech complex
    FilteredComplex chromatic_delcech_complex(const RealMatrix<double>& points, const vector<size_t>& colours) {
        auto [new_colours, num_colours] = canonicalise(colours);
        RealMatrix<double> stratified_points = stratify(points, colours, num_colours);
        auto delX = delaunay_complex(stratified_points);
        // modify the colours of the vertices
        for (auto& [idx, vert] : delX.get_simplices()[0]) {
            vert->colours.set(new_colours[idx]);
        }
        delX.propagate_colours();
        // modify the filtration values
        bool numerical_issues = false;
        if (delX.max_dim >= 1) {
            for (int p = delX.max_dim; p >= 1; p--) {
                for (auto& [idx, simplex] : delX.get_simplices()[p]) {
                    auto verts = simplex->get_vertex_labels();
                    auto [centre, sqRadius, success] = cmb::miniball<double>(points.rows(), points(Eigen::all, verts));
                    numerical_issues |= !success;
                    simplex->value = sqrt(sqRadius);
                }
            }            
        }
        if (numerical_issues) {
            std::cerr << "Warning: encountered numerical problems. Filtration values may be inaccurate." << std::endl;
        }
        return delX;
    }
}