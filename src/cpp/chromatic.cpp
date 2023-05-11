#include <CGAL/Epick_d.h>
#include <CGAL/Dimension.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Triangulation.h>
#include <algorithm>

#include "chromatic.h"

namespace chalc {
    using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
    using Triangulation_data_structure = CGAL::Triangulation_data_structure< Kernel::Dimension,
        CGAL::Triangulation_vertex<Kernel, size_t>,
        CGAL::Triangulation_full_cell<Kernel>>;
    using DelaunayTriangulation = CGAL::Delaunay_triangulation<Kernel, Triangulation_data_structure>;
    using Point = Kernel::Point_d;
    using Colouring = vector<size_t>;

    // Convert a collection of coordinate vectors to a vector of CGAL points
    vector<Point> coordvecs_to_points(const MatrixXd& x_arr) {
        vector<Point> points(x_arr.cols());
        for (size_t i = 0; i < x_arr.cols(); i++) {
            points[i] = Point(x_arr.col(i).begin(), x_arr.col(i).end());
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
    tuple<Colouring, size_t> canonicalise(const Colouring& colours) {
        Colouring new_colours(colours.size());
        map<size_t, size_t> colour_map;
        for (auto [c, i] = tuple{colours.begin(), 0}; c != colours.end(); c++) {
            if (!colour_map.contains(*c)) {
                colour_map[*c] = i++;
            }
        }
        for (size_t i = 0; i < colours.size(); i++) {
            new_colours[i] = colour_map[colours[i]];
        }
        return tuple{new_colours, colour_map.size()};
    }


    // Stratify a coloured point set
    // Points are provided as columns of a matrix
    // Colours are provided as a vector
    MatrixXd stratify(const MatrixXd& points, const Colouring& colours, size_t num_colours) {
        size_t dim = points.rows();
        Colouring new_colours(colours.size());
        // make sure colours are contiguous and start at zero
        if (num_colours == -1) {
            auto tmp = canonicalise(colours);
            new_colours = std::get<0>(tmp);
            num_colours = std::get<1>(tmp);
        }
        else {
            new_colours = colours;
        }
        MatrixXd result(dim + num_colours - 1, points.cols());
        result.topRows(dim) = points;
        //Colouring new_colours(colours);
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
    FilteredComplex delaunay_complex(const MatrixXd& X) {
        auto dim = X.rows();
        FilteredComplex result(X.cols(), dim);
        if (X.cols() != 0) {
            auto points = coordvecs_to_points(X);
            auto [sorted_points, sorted_indices] = sort_with_indices<Point>(points, [](const Point& a, const Point& b) -> bool { return a < b; });
            auto delX = DelaunayTriangulation(dim);
            delX.insert(points.begin(), points.end());
            // iterate over all finite vertices and associate them with their original label
            for (auto vert_it = delX.finite_vertices_begin(); vert_it != delX.finite_vertices_end(); vert_it++) {
                // find the index of its associated point
                auto point = vert_it->point();
                const Point* p = static_cast<Point*>(
                    std::bsearch(&point, sorted_points.data(), sorted_points.size(), sizeof(Point), compare<Point>)
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
    FilteredComplex weak_chromatic_alpha_complex(const MatrixXd& points, const Colouring& colours) {
        auto [new_colours, num_colours] = canonicalise(colours);
        MatrixXd stratified_points = stratify(points, colours, num_colours);
        auto delX = delaunay_complex(stratified_points);
        // modify the colours of the vertices
        for (auto& [idx, vert] : delX.get_simplices()[0]) {
            vert->colours = new_colours[idx];
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

    // Create the chromatic alpha complex
    // FilteredComplex chromatic_alpha_complex(const Matrix<double, 2, Dynamic>& points, const Colouring& colours) {
    //     [stratified_points, new_colours, num_colours] = stratify(points, colours);
    //     delX = delaunay_complex(stratified_points);
    //     CGAL::Bbox_2 plane();
    //     // modify the filtration values
    //     for (auto p = delX.max_dim; p > 0; p--) {
    //         CGAL:: = plane;
    //         vector<size_t> verts(p+1);
    //         for (auto& [ind, val] : delX.simplices[p]) {
    //             verts = delX.get_simplex_vertices(p, ind);
    //             auto E = CGAL::
    //             for(j = 0; j < num_colours; j++) { 
    //                 // find all points of the simplex with colour j
    //                 vector<int> vj;
    //                 for (auto& v : verts) {
    //                     if (new_colours[v] == j) {
    //                         vj.push_back(v);
    //                     }
    //                 }
    //                 // compute E_j
    //                 // ...
    //                 // intersect E_j with E - must get a point
    //             }
    //         }

    //     }
    // }

}