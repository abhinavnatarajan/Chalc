#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/attr.h>
#include <string>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#include "filtration.h"
#include "chromatic.h"

namespace py = pybind11;
using namespace chalc;

PYBIND11_MODULE(core, m) {
    m.doc() = 
        R"docstring(
            Module containing the core functionality of this package. 
        )docstring";

    py::class_<FilteredComplex::Simplex, shared_ptr<FilteredComplex::Simplex>> simplex_class(m, "SimplexObject");
        simplex_class.doc() = 
            R"docstring(
                Class representing a simplex in a filtered simplicial complex.
            )docstring";
        simplex_class.def_readonly("dimension", 
            &FilteredComplex::Simplex::dim, 
            R"docstring(
                Dimension of the simplex.
            )docstring")
        .def_readonly("label", 
            &FilteredComplex::Simplex::label,
            R"docstring(
                Label of the simplex in its parent filtered complex.
            )docstring")
        .def_property_readonly("vertices", 
            static_cast<vector<size_t>(FilteredComplex::Simplex::*)() const>(&FilteredComplex::Simplex::get_vertex_labels),
            R"docstring(
                List of (sorted, ascending) vertex labels of the simplex.
            )docstring")
        .def_readwrite("filtration_value", 
            &FilteredComplex::Simplex::value,
            R"docstring(
                Filtration value of the simplex. If you change the value, 
                you are responsible for checking whether the parent filtered simplicial complex
                satisfies the filtration property.
            )docstring")
        .def_readwrite("colours", 
            &FilteredComplex::Simplex::colours,
            R"docstring(
                Integer whose binary representation is the bitmask of colours
                in the simplex. For example, this value is a power of two for
                monochromatic simplices. If you change this value for a vertex, 
                you should call FilteredComplex.propagate_colours() afterwards 
                to ensure that colours of higher dimensional simplices are consistent
                with the colours of their vertices.
            )docstring")
        .def_property_readonly("facets", 
            &FilteredComplex::Simplex::get_facets, 
            R"docstring(
                Returns a read-only list of handles to the facets of the simplex.
            )docstring")
        .def("__repr__",
            [](const shared_ptr<FilteredComplex::Simplex> s_ptr) {
                return "<" + std::to_string(s_ptr->dim) + "-simplex>";
            });

    py::class_<FilteredComplex> filtered_complex_class(m, "FilteredComplex");
        filtered_complex_class.doc() =
            R"docstring(
                Class representing a filtered simplicial complex.
            )docstring";
        filtered_complex_class.def(py::init<const size_t, const size_t>(),
            R"docstring(
                Construct a discrete filtered simplicial complex with default filtration time of 0.
                
                Parameters
                ----------
                n : 
                Number of vertices. Cannot be changed after initialisation.
                k : 
                Maximum dimension of a simplex that the complex can have.
                This parameter is required for memory efficiency, and cannot
                be changed after initialisation.
            )docstring",
            py::arg("n"), py::arg("k"))
        .def(py::init<const FilteredComplex&, const size_t>(),
           R"docstring(
                Returns a handle to the k-skeleton of an existing filtered simplicial complex.
                
                Parameters
                ----------
                other : 
                Filtered simplicial complex whose k-skeleton is required.
                k : 
                Integer determing the dimension of the skeleton.
            )docstring",
            py::arg("other"), py::arg("k"))
        .def("add_simplex", &FilteredComplex::add_simplex,
            R"docstring(
                Add a simplex to a filtered simplicial complex.
                
                Parameters
                ----------
                vertices : 
                List of vertex labels corresponding to existing vertices in the complex
                filt_value : 
                Filtration value to associate to the new simplex. 
                Faces of the added simplex that are already present 
                in the simplicial complex will have their filtration 
                values reduced if necessary.                 
            )docstring",
            py::arg("vertices"), py::arg("filt_value"))
        .def_property_readonly("num_simplices", &FilteredComplex::size,
            R"docstring(
                The total number of simplices in the complex.
            )docstring")
        .def("count_simplices_in_dim", &FilteredComplex::size_in_dim,
            R"docstring(
                The number of k-simplices in the complex.
            )docstring",
            py::arg("k"))
        .def_property_readonly("dimension", &FilteredComplex::dimension,
            R"docstring(
                Current maximum dimension of a maximal simplex in the complex.
            )docstring")
        .def_readonly("max_dimension", &FilteredComplex::max_dim,
            R"docstring(
                Maximum dimension of simplex that this complex can store.
                Set during initialisation. 
            )docstring")
        .def_readonly("num_vertices", &FilteredComplex::N,
            R"docstring(
                Number of vertices in the simplicial complex.
                Set during initialisation.
            )docstring")
        .def("propagate_filt_values", &FilteredComplex::propagate_filt_values,
            R"docstring(
                Propagate filtration values upwards or downwards to ensure that
                every simplex appears after its faces.
                For example, setting the filtration values in dimension 1
                and propagating upwards is akin to the Rips filtration.
                
                Parameters
                ----------
                start_dim : 
                Dimension from which to start propagating (exclusive). 
                upwards : 
                If true then values are propagated upwards, downwards otherwise.
            )docstring",
            py::arg("start_dim"), py::arg("upwards"))
        .def("has_simplex", static_cast<bool (FilteredComplex::*)(vector<size_t>&) const>(&FilteredComplex::has_simplex),
            R"docstring(
                Check for membership of a simplex in the complex.
                
                Parameters
                ----------
                vertices : 
                Vertex labels of the simplex to check for. 
            )docstring",
            py::arg("vertices"))
        .def_property_readonly("simplices", &FilteredComplex::get_simplices,
            R"docstring(
                A list of dictionaries, where list[k] contains
                the k-simplices of the complex. 

                The key of a k-simplex in the [k]th dictionary is 
                the lexicographic index of that simplex with respect 
                to its vertex labels sorted in ascending order, 
                counting all possible sorted subsequences of (0, ..., N-1)
                of length k.
            )docstring")
        .def("get_label_from_vertex_labels", &FilteredComplex::get_label_from_vertex_labels,
            R"docstring(
                Returns the dictionary key of a simplex with respect to 
                its vertex labels sorted in ascending order, counting all 
                possible sorted subsequences of (0, ..., N-1) of length k, 
                where N is the number of vertices in the complex. 
                The simplex need not be present in the simplicial complex.
                
                Parameters
                ----------
                vertices : 
                List of vertex labels of the simplex.
            )docstring",
            py::arg("vertices"))
        .def("propagate_colours", &FilteredComplex::propagate_colours,
            R"docstring(
                Function to make sure that simplex colours are consistent with the colours 
                of their vertices. You should call this whenever you change the colour of 
                any vertex. 
            )docstring")
        .def("flat_representation", &FilteredComplex::flat_representation,
            R"docstring(
                Returns a serialised representation of the simplicial complex.
                Each list element is a tuple of the following signature:
                (f : list[int], idx: int, v : float, c: int)
                This corresponds to a simplex whose facets are the elements of the list 
                at indices f, whose original label in the FilteredComplex is idx, whose 
                filtration value is v, and whose colour bitmask is c. 
                Simplices appear in ascending order of dimension, then filtration value. 
            )docstring")
        .def("__repr__",
            [](const FilteredComplex& K) {
                return "<" + std::to_string(K.dimension()) + "-dimensional simplicial complex with " + std::to_string(K.N) + " vertices>";
            });

    m.def("clique_complex", &FilteredComplex::clique_complex, 
        R"docstring(
            Returns the k-skeleton of the complete simplicial complex
            on n vertices, with filtration values initialised to zero.
        )docstring",
        py::arg("n"), py::arg("k"));
    m.def("standard_simplex", &standard_simplex, 
        R"docstring(
            Returns the simplicial complex corresponding to 
            the standard abstract n-simplex.
        )docstring",
        py::arg("n"));
    m.def("delaunay_complex", &delaunay_complex,
        R"docstring(
            Returns the Delaunay triangulation of a point cloud in Euclidean space.

            Parameters
            ----------
            x :
            A numpy matrix whose columns are points in the point cloud.
        )docstring",
        py::arg("x"));
    m.def("stratify", 
        [](const MatrixXd& M, const Colouring& c) -> MatrixXd {
            return stratify(M, c); // default third argument
        },
        R"docstring(
            Returns the stratification of a point cloud by a sequence of colours.
            If s is the number of colours and d is the ambient dimension of the point cloud,
            then the points are translated along the vertices of the s-simplex whose
            vertices are the origin and the last s-1 basis elements of 
            R^(d + s)-dimensional Euclidean space. 

            Parameters
            ----------
            x :
            A numpy matrix whose columns are points in the point cloud.
            colours :
            A list of integers.
        )docstring",
        py::arg("x"), py::arg("colours"));
    m.def("weak_chromatic_alpha_complex", &weak_chromatic_alpha_complex,
        R"docstring(
            Returns the weak chromatic alpha complex of a coloured point cloud.
            This is a filtration of the Delaunay complex of the point cloud
            after stratification by the colours, where the filtration time
            of a simplex is the filtration time of that simplex in the 
            Vietoris-Rips filtration of the non-stratified point cloud. 

            Parameters
            ----------
            x :
            A numpy matrix whose columns are points in the point cloud.
            colours :
            A list of integers representing the colours of vertices.
            Note that the actual colours of simplices in the output filtration
            may not correspond to the input colours unless the set of values in
            colours is contiguous and colours[0] = 0.
        )docstring",
        py::arg("x"), py::arg("colours"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}