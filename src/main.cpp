#include <pybind11/pybind11.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Dimension.h>
#include <CGAL/Delaunay_triangulation.h>

#include <vector>

namespace py = pybind11;

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
typedef K::Point_d Point;
typedef std::vector<double> CartesianCoordinates;

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

// Convert a vector of 
Point coords_to_point(const CartesianCoordinates &x) {
    return Point(x.begin(), x.end());
}

template <typename Container_t>
requires IterableOfType<Container_t, CartesianCoordinates>
std::vector<Point> coords_to_point(const Container_t& x_arr) {
    std::vector<Point> points;
    for (auto x = x_arr.begin(); x != x_arr.end(); x++) {
        points.push_back(coords_to_point(*x));
    }
    return points;
}

std::vector<double> point_to_coords(const Point &p) {
    CartesianCoordinates x;
    for(auto coord = p.cartesian_begin(); coord != p.cartesian_end(); coord++) {
        x.push_back(*coord);
    }
    return x;
}

template <typename Container_t>
requires IterableOfType<Container_t, Point>
std::vector<CartesianCoordinates> point_to_coords(const Container_t& points) {
    std::vector<CartesianCoordinates> x_arr;
    for (auto p = points.begin(); p != points.end(); p++) {
        x_arr.push_back(point_to_coords(*p));
    }
    return x_arr;
}


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