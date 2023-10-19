/*
    This file is part of Chalc.

    Chalc: Chromatic Alpha Complexes.
    Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha Complexes”.
    Online preprint available at http://arxiv.org/abs/2212.03128.
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
/*
DO NOT MOVE THE NEXT INCLUDE.
mpreal must be included before CGAL due to
preprocessor macro clashes
*/
#include "../../mpreal/mpreal.h"
#include <unsupported/Eigen/MPRealSupport>
#include "../../ConstrainedMiniball/ConstrainedMiniball.h"
#include <CGAL/Epick_d.h>
#include <CGAL/Dimension.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Triangulation.h>
#include <stdexcept>
#include <algorithm>

namespace
{
    using namespace chalc::stl;
    using 
        chalc::index_t,
        chalc::MAX_NUM_COLOURS,
        Eigen::lastN,
        Eigen::all,
        std::min,
        std::runtime_error,
        std::bad_weak_ptr,
        std::sort,
        std::unique;
    using Kernel_d = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
    using Triangulation_data_structure = CGAL::Triangulation_data_structure<Kernel_d::Dimension,
                                                                            CGAL::Triangulation_vertex<Kernel_d, size_t>,
                                                                            CGAL::Triangulation_full_cell<Kernel_d>>;
    using DelaunayTriangulation = CGAL::Delaunay_triangulation<Kernel_d, Triangulation_data_structure>;
    using Point_d = Kernel_d::Point_d;

    template <class Real_t>
    using RealVector = Eigen::Matrix<Real_t, Eigen::Dynamic, 1>;

    template <class Real_t>
    using RealMatrix = Eigen::Matrix<Real_t, Eigen::Dynamic, Eigen::Dynamic>;

    template <class Derived>
    concept MatrixXpr = requires { typename Eigen::MatrixBase<Derived>; };

    template <class Derived, class Real_t>
    concept RealMatrixXpr = MatrixXpr<Derived> && std::same_as<typename Derived::Scalar, Real_t>;

    // Convert a collection of coordinate vectors to a vector of CGAL points
    vector<Point_d> coordvecs_to_points(const RealMatrix<double> &x_arr)
    {
        vector<Point_d> points(x_arr.cols());
        for (index_t i = 0; i < x_arr.cols(); i++)
        {
            points[i] = Point_d(x_arr.col(i).begin(), x_arr.col(i).end());
        }
        return points;
    }

    template <typename T>
    vector<T> reorder(const vector<T> &v, const vector<index_t> &idx)
    {
        vector<T> result(v.size());
        for (auto i = 0; i < v.size(); i++)
        {
            result[i] = v[idx[i]];
        }
        return result;
    }

    template <typename T>
    tuple<vector<T>, vector<index_t>> sort_with_indices(const vector<T> &v, bool (*compare)(const T &a, const T &b))
    {
        vector<index_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(),
                         [&v, &compare](index_t i1, index_t i2)
                         { return compare(v[i1], v[i2]); });
        return tuple{reorder(v, idx), idx};
    }

    template <typename T>
    int compare(const void *a, const void *b)
    {
        const auto &arg1 = *(static_cast<const T *>(a));
        const auto &arg2 = *(static_cast<const T *>(b));
        return arg1 < arg2   ? -1
               : arg1 > arg2 ? +1
                             : 0;
    }

    // make a vector contiguous and start at zero
    // returns new vector and number of distinct elements
    tuple<vector<index_t>, index_t> canonicalise(const vector<index_t> &vec)
    {
        vector<index_t> new_vec(vec.size());
        map<index_t, index_t> m;
        for (auto [c, i] = tuple{vec.begin(), 0}; c != vec.end(); c++)
        {
            if (!m.contains(*c))
            {
                m[*c] = i++;
            }
        }
        for (index_t i = 0; i < vec.size(); i++)
        {
            new_vec[i] = m[vec[i]];
        }
        return tuple{new_vec, m.size()};
    }

    /*
    Return the affine subspace of n points equidistant from a given
    set of points in R^d. The points are input as columns of a (d, n)
    matrix or matrix expression, and the subspace is returned as
    a tuple (E, b), where
        E is an (n - 1, d) matrix
        b is a vector in R^(n-1)
    such that
        E * x = b
    defines the desired subspace.

    If (x_1, ..., x_n) are the points, then for any point x in the
    desired subspace we have for all 1 < i <= n

        |x_1 - x|^2 = |x_i - x|^2
    --> |x_1|^2 + |x|^2 - 2<x_1, x> = |x_i|^2 + |x|^2 - 2<x_i, x>
    --> 2<x_i - x_1, x> = |x_i|^2 - |x_1|^2

    The last line is a linear constraint, so we add a row
        x_i - x_1
    to the matrix E, and we add a corresponding entry
        0.5 * (|x_i|^2 - |x_1|^2)
    to the vector b.
    */
    template <typename Real_t, RealMatrixXpr<Real_t> Derived>
    tuple<
        RealMatrix<Real_t>,
        RealVector<Real_t>>
    equidistant_subspace(const Eigen::MatrixBase<Derived> &X)
    {
        int n = X.cols(); // number of points
        int d = X.rows(); // dimension of ambient space
        RealMatrix<Real_t> E(n - 1, d);
        RealVector<Real_t> b(n - 1);
        if (n > 1)
        {
            b = 0.5 * (X.rightCols(n - 1).colwise().squaredNorm().array() - X.col(0).squaredNorm()).transpose();
            E = (X.rightCols(n - 1).colwise() - X.col(0)).transpose();
        }
        return tuple{E.eval(), b.eval()};
    }

    /*
    Stratify a coloured point set.
    Points are provided as columns of a matrix or matrix expression.
    Colours are provided as a vector.
    */
    RealMatrix<double> stratify(const RealMatrix<double> &points, const vector<index_t> &colours)
    {
        index_t dim = points.rows();
        // Make sure colours are contiguous and start at zero
        auto [new_colours, num_colours] = canonicalise(colours);
        RealMatrix<double> result(dim + num_colours - 1, points.cols());
        result.topRows(dim) = points;
        if (num_colours != 1)
        {
            result.bottomRows(num_colours - 1).setZero();
            for (auto i = 0; i < result.cols(); i++)
            {
                if (new_colours[i] != 0)
                {
                    result(dim - 1 + new_colours[i], i) = 1.0;
                }
            }
        }
        return result;
    }

    bool bad_colour(const index_t& colour)
    {
        return (colour >= MAX_NUM_COLOURS || colour < 0);
    }
}

namespace chalc::chromatic
{
    // Create a Delaunay triangulation from a collection of coordinate vectors
    FilteredComplex delaunay(const RealMatrix<double> &X)
    {
        auto dim = X.rows();
        FilteredComplex result(X.cols(), dim);
        if (X.cols() != 0)
        {
            auto points = coordvecs_to_points(X);
            auto [sorted_points, sorted_indices] =
                sort_with_indices<Point_d>(
                    points, [](const Point_d &a, const Point_d &b) -> bool
                    { return a < b; });
            auto delX = DelaunayTriangulation(dim);
            delX.insert(points.begin(), points.end());
            // iterate over all finite vertices and associate them with their original label
            for (auto vert_it = delX.finite_vertices_begin(); vert_it != delX.finite_vertices_end(); vert_it++)
            {
                // find the index of its associated point
                auto point = vert_it->point();
                const Point_d *p = static_cast<Point_d *>(
                    std::bsearch(&point, sorted_points.data(), sorted_points.size(), sizeof(Point_d), compare<Point_d>));
                vert_it->data() = sorted_indices[p - sorted_points.data()];
            }
            // iterate over the top dimensional cells and add them to the filtration
            vector<index_t> max_cell_vertex_labels(dim + 1);
            for (auto cell_it = delX.finite_full_cells_begin();
                 cell_it != delX.finite_full_cells_end();
                 cell_it++)
            {
                // iterate over the vertices of the cell and get their labels
                for (auto [label_it, vert_it] = tuple{max_cell_vertex_labels.begin(), cell_it->vertices_begin()};
                     vert_it != cell_it->vertices_end();
                     label_it++, vert_it++)
                {
                    *label_it = (*vert_it)->data();
                }
                result.add_simplex(max_cell_vertex_labels, 0.0);
            }
        }
        return result;
    }

    // Create the chromatic Del-VR complex
    FilteredComplex delrips(
        const RealMatrix<double> &points,
        const vector<index_t> &colours)
    {
        // input checks
        if (std::any_of(colours.begin(), colours.end(), bad_colour))
        {
            throw std::domain_error("All colours must be between 0 and " + std::to_string(MAX_NUM_COLOURS-1));
        }
        if (colours.size() != points.cols())
        {
            throw std::domain_error("len(colours) must equal number of points.");
        }

        // Start
        // Stratify the points by colour
        RealMatrix<double> stratified_points = stratify(points, colours);
        auto delX = delaunay(stratified_points);
        // modify the colours of the vertices
        for (auto &[idx, vert] : delX.get_simplices()[0])
        {
            vert->set_colour(colours[idx]);
        }
        delX.propagate_colours();
        // modify the filtration values
        if (delX.dimension() >= 1)
        {
            for (auto &[idx, edge] : delX.get_simplices()[1])
            {
                auto verts = edge->get_vertex_labels();
                edge->value = (points.col(verts[0]) - points.col(verts[1])).norm();
            }
            delX.propagate_filt_values(1, true);
        }
        return delX;
    }

    // Create the chromatic alpha complex
    FilteredComplex alpha(
        const RealMatrix<double> &points,
        const vector<index_t> &colours,
        std::ostream &ostream)
    {
        // input checks
        if (std::any_of(colours.begin(), colours.end(), bad_colour))
        {
            throw std::domain_error("All colours must be between 0 and " + std::to_string(MAX_NUM_COLOURS-1));
        }
        if (colours.size() != points.cols())
        {
            throw std::domain_error("len(colours) must equal number of points.");
        }

        // Start
        // Stratify the points by colour
        RealMatrix<double> stratified_points = stratify(points, colours);
        // Get the delaunay triangulation
        auto delX = delaunay(stratified_points);
        // Modify the colours of the vertices and simplices
        for (auto &[idx, vert] : delX.get_simplices()[0])
        {
            vert->set_colour(colours[idx]);
        }
        delX.propagate_colours();

        // Partition the vertices by colour
        // We will need this later to check if stacks are empty
        map<index_t, vector<index_t>> verts_by_colour;
        for (index_t i = 0; i < points.cols(); i++)
        {
            verts_by_colour[colours[i]].push_back(i);
        }

        if (delX.dimension() >= 1)
        {
            // Now comes the hard bit
            // Modify the filtration values
            // We will use an arbitrary precision scalar type
            // for the geometric computation
            using mpfr::mpreal;
            mpreal::set_default_prec(mpfr::digits2bits(25));
            RealMatrix<mpreal> X = points.template cast<mpreal>();

            bool numerical_instability = false; // flag to check numerical instability

            // Start at the maximum dimension
            for (int p = delX.dimension(); p >= 1; p--)
            {
                // Iterate over p-simplices
                for (auto &[idx, simplex] : delX.get_simplices()[p])
                {
                    auto verts = simplex->get_vertex_labels();

                    // Partition the vertices of this simplex by colour
                    map<index_t, vector<index_t>> verts_by_colour_in_simplex;
                    for (auto v : verts)
                    {
                        verts_by_colour_in_simplex[colours[v]].push_back(v);
                    }

                    /*
                    For each colour j in the simplex, find the affine subspace
                    of points equidistant to the vertices of colour j in the simplex.
                    Suppose this subspace is E_j and is defined by a matrix equation:
                        E_j * x = b_j
                    We append E_j and b_j to the bottom of a matrix E and a vector b.
                    Then the solution set E of the matrix equation
                        E * x = b
                    is the intersection over all j of the affine subspaces E_j
                    */
                    RealMatrix<mpreal> E(0, X.rows());
                    RealVector<mpreal> b;
                    for (auto &[j, verts_j] : verts_by_colour_in_simplex)
                    {
                        index_t num_new_rows = verts_j.size() - 1;
                        E.conservativeResize(E.rows() + num_new_rows, Eigen::NoChange);
                        b.conservativeResize(b.rows() + num_new_rows);
                        RealMatrix<mpreal>::BlockXpr E_new_rows = E.bottomRows(num_new_rows);
                        Eigen::VectorBlock<RealVector<mpreal>> b_new_rows = b(lastN(num_new_rows));
                        tie(E_new_rows, b_new_rows) =
                            equidistant_subspace<mpreal>(X(all, verts_j));
                    }

                    /*
                    Get the smallest bounding ball of the points
                    in the simplex, with the added constraint that
                    the centre x of the ball must satisfy the equation
                        E * x = b
                    i.e., it lies in the affine subspace E
                    */
                    auto [centre, sqRadius, success] =
                        cmb::constrained_miniball<mpreal>(
                            X(all, verts), E, b);
                    // Check if there were any numerical issues
                    numerical_instability |= !success;

                    // Check if the stack is empty
                    bool stack_is_empty = true;
                    for (auto &[j, verts_j] : verts_by_colour_in_simplex)
                    {
                        /*
                        Compute the radius r_j of the sphere of colour j.
                        Take all points of colour j in the simplex
                        and find the minimum distance of any of those
                        points from the centre of the bounding ball.
                        Theoretically they should all be equal, but
                        we have to account for floating point errors.
                        We choose the radius this way to ensure that
                        if any j-coloured point from the dataset lies
                        in the sphere of radius r_j, then it is closer
                        to the centre than every  j-coloured point in the
                        simplex. In this case there is a high chance that
                        the stack is not empty.
                        */
                        mpreal rj_squared =
                            (X(all, verts_j).colwise() - centre)
                                .colwise()
                                .squaredNorm()
                                .minCoeff();
                        mpreal squared_dist_to_nearest_pt_of_colour_j =
                            (X(all, verts_by_colour[j]).colwise() - centre).colwise().squaredNorm().minCoeff();

                        stack_is_empty &=
                            (squared_dist_to_nearest_pt_of_colour_j >= rj_squared);
                    }
                    // If the stack is empty, assign the filtration value
                    if (stack_is_empty)
                    {
                        simplex->value = static_cast<double>(sqrt(sqRadius));
                        // We need to take into account floating point issues
                        // so we make sure that we satisfy the filtration property
                        try
                        {
                            for (auto &cofacet : simplex->get_cofacets())
                            {
                                simplex->value = min(simplex->value, cofacet.lock()->value);
                            }
                        }
                        catch (bad_weak_ptr &e)
                        {
                            throw runtime_error("Tried to dereference expired cofacet handle.");
                        }
                    }
                }
                // Propagate filtration values down
                delX.propagate_filt_values(p, false);
            }
            if (numerical_instability)
            {
                ostream << "Warning: encountered numerical problems. \
                            Filtration values may be inaccurate."
                        << std::endl;
            }
        }
        return delX;
    }

    // Create the chromatic Del-Cech complex
    FilteredComplex delcech(const RealMatrix<double> &points, const vector<index_t> &colours, std::ostream &ostream)
    {
        // input checks
        if (std::any_of(colours.begin(), colours.end(), bad_colour))
        {
            throw std::domain_error("All colours must be between 0 and " + std::to_string(MAX_NUM_COLOURS-1));
        }
        if (colours.size() != points.cols())
        {
            throw std::domain_error("len(colours) must equal number of points.");
        }

        // Start
        // Stratify the points by colour
        RealMatrix<double> stratified_points = stratify(points, colours);
        auto delX = delaunay(stratified_points);
        // modify the colours of the vertices
        for (auto &[idx, vert] : delX.get_simplices()[0])
        {
            vert->set_colour(colours[idx]);
        }
        delX.propagate_colours();
        // modify the filtration values
        bool numerical_instability = false;
        if (delX.dimension() >= 1)
        {
            for (int p = delX.dimension(); p >= 1; p--)
            {
                for (auto &[idx, simplex] : delX.get_simplices()[p])
                {
                    auto verts = simplex->get_vertex_labels();
                    auto [centre, sqRadius, success] =
                        cmb::miniball<double>(points(all, verts));
                    numerical_instability |= !success;
                    simplex->value = sqrt(sqRadius);
                    auto cofacets = simplex->get_cofacets();
                    if (cofacets.size() == 0)
                    {
                        continue;
                    }
                    try
                    {
                        for (auto &cofacet : simplex->get_cofacets())
                        {
                            simplex->value = min(simplex->value, cofacet.lock()->value);
                        }
                    }
                    catch (bad_weak_ptr &e)
                    {
                        throw runtime_error("Tried to dereference expired cofacet handle.");
                    }
                }
                // Propagate filtration values down
                delX.propagate_filt_values(p, false);
            }
        }
        if (numerical_instability)
        {
            ostream << "Warning: encountered numerical problems. \
                        Filtration values may be inaccurate."
                    << std::endl;
        }
        return delX;
    }
}
