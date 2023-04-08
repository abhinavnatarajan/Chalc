#include <pybind11/pybind11.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Dimension.h>
#include <CGAL/Delaunay_triangulation.h>

#include <vector>

namespace py = pybind11;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef Kernel::Point_d Point;
typedef std::vector<double> CoordinateVector;
typedef CGAL::Delaunay_triangulation<Kernel> DelComplex;

template <typename T>
concept Iterable = requires(T x) {
    x.begin();
    x.end();
};

template <typename T, typename U>
concept IterableOfType = Iterable<T> && std::same_as<typename T::value_type, U>;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

// Convert a vector of coordinates to a CGAL point
Point coordvec_to_point(const CoordinateVector &x) {
    return Point(x.begin(), x.end());
}

// Convert a collection of coordinate vectors to a vector of CGAL points
template <typename Container_t>
requires IterableOfType<Container_t, CoordinateVector>
std::vector<Point> coordvecs_to_points(const Container_t& x_arr) {
    std::vector<Point> points;
    for (auto x = x_arr.begin(); x != x_arr.end(); x++) {
        points.push_back(coords_to_point(*x));
    }
    return points;
}

// Convert a CGAL point to a coordinate vector
std::vector<double> point_to_coordvec(const Point &p) {
    CoordinateVector x;
    for(auto coord = p.cartesian_begin(); coord != p.cartesian_end(); coord++) {
        x.push_back(*coord);
    }
    return x;
}

// Convert a collection of CGAL points to a vector of coordinate vectors
template <typename Container_t>
requires IterableOfType<Container_t, Point>
std::vector<CoordinateVector> points_to_coordvecs(const Container_t& points) {
    std::vector<CoordinateVector> x_arr;
    for (auto p = points.begin(); p != points.end(); p++) {
        x_arr.push_back(point_to_coords(*p));
    }
    return x_arr;
}

// Create a Delaunay triangulation from a vector of coordinate vectors
// template <typename Container_t>
// requires IterableOfType<Container_t, CoordinateVector>
// DelComplex del_complex(const Container_t& x_arr) {
//     if (x_arr.size() == 0) {
//         return 
//     }
//     for (auto it = x_arr.begin(); it != x_arr.end(); it++) {
//         if it->size
//     }
//     auto points = coordvecs_to_points(x_arr);
    
// }

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Chalc
        -----------------------
        .. currentmodule:: chalc
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}