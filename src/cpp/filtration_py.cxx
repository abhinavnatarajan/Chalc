#include "filtration.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

PYBIND11_MODULE(filtration, m)
{
    using namespace chalc;
    using namespace chalc::common;
    namespace py = pybind11;
    m.doc() =
        R"docstring(
            Module containing utilities to store and manipulate abstract filtered simplicial complexes. 
        )docstring";

    py::class_<FilteredComplex::Simplex, shared_ptr<FilteredComplex::Simplex>> simplex_class(m, "SimplexObject");
    simplex_class.doc() =   
        R"docstring(
            Class representing a simplex in a filtered simplicial complex.
        )docstring";
    simplex_class.def_readonly("dimension", &FilteredComplex::Simplex::dim,
        R"docstring(
            Dimension of the simplex.
        )docstring")
        .def_readonly("label", &FilteredComplex::Simplex::label,
            R"docstring(
                Label of the simplex in its parent filtered complex.
                A :math:`k`-simplex :math:`\sigma` is labelled 
                by the lexicographic index of :math:`\sigma` with respect 
                to its vertex labels sorted in ascending order, 
                counting all possible sorted subsequences of :math:`(0, ..., N-1)`
                of length :math:`k`.
            )docstring")
        .def_property_readonly("vertices", static_cast<vector<index_t> (FilteredComplex::Simplex::*)() const>(&FilteredComplex::Simplex::get_vertex_labels),
            R"docstring(
                List of (sorted, ascending) vertex labels of the simplex.
            )docstring")
        .def_readwrite("filtration_value", &FilteredComplex::Simplex::value,
            R"docstring(
                Filtration value of the simplex. If you change the value, 
                you are responsible for checking whether the parent complex
                satisfies the filtration property.
            )docstring")
        .def_property_readonly(
            "colours",
            [](const shared_ptr<FilteredComplex::Simplex> &s_ptr)
            {
                return s_ptr->colours.to_ulong();
            },
            R"docstring(
                Bitmask of colours in the simplex, where the rightmost
                bit represents the smallest colour index. 
            )docstring")
        .def(
            "set_colour",
            [](const shared_ptr<FilteredComplex::Simplex> &s_ptr, unsigned long c)
            {
                if (s_ptr->dim == 0)
                {
                    if (c <= MAX_NUM_COLOURS)
                    {
                        s_ptr->colours.reset().set(c);
                    }
                    else
                    {
                        throw std::domain_error("Colour index too large.");
                    }
                }
                else
                {
                    throw std::domain_error("Can't change colour unless simplex is a vertex.");
                }
            },
            R"docstring(
                Change the colour of a vertex.
                
                Raises
                ------
                ValueError
                    If the simplex is not a vertex or if the colour is too large.
                
                Tip
                ---
                It is recommended to call the member function 
                `propagate_colours` of the parent simplicial complex 
                after changing the colour of a vertex.
            )docstring",
            py::arg("colour"))
        .def_property_readonly("facets",
                               &FilteredComplex::Simplex::get_facets,
            R"docstring(
                Returns a read-only list of handles to the facets of the simplex.
            )docstring")
        .def("__repr__",
             [](const shared_ptr<FilteredComplex::Simplex> &s_ptr)
             {
                 return "<" + std::to_string(s_ptr->dim) + "-simplex>";
             });

    py::class_<FilteredComplex> filtered_complex_class(m, "FilteredComplex");
    filtered_complex_class.doc() =
        R"docstring(
            Class representing a filtered simplicial complex.
        )docstring";
    filtered_complex_class.def(py::init<const index_t, const index_t>(),
        R"docstring(
            Construct a discrete filtered simplicial complex 
            with default filtration time of 0.
            
            Parameters
            ----------
            n : int
                Number of vertices. Cannot be changed after initialisation.
            k : int
                Maximum dimension of a simplex that the complex can have.
                This parameter is required for memory efficiency, and cannot
                be changed after initialisation.
        )docstring",
        py::arg("n"), py::arg("k"))
    .def(py::init<const FilteredComplex &, const index_t>(),
        R"docstring(
            Returns a handle to the :math:`k`-skeleton 
            of an existing (filtered) simplicial complex.
            
            Parameters
            ----------
            other : FilteredComplex
                Filtered simplicial complex.
            k : int
                Skeleton dimension.
        )docstring",
        py::arg("other"), py::arg("k"))
        .def("add_simplex", &FilteredComplex::add_simplex,
            R"docstring(
                Add a simplex to a filtered simplicial complex.
                
                Parameters
                ----------
                vertices : list[int]
                    List of vertex labels corresponding to existing vertices 
                    in the complex.
                filt_value : float
                    Filtration value to associate to the new simplex. 
                
                Note
                ----
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
                Count the number of simplices of a given dimension.
            )docstring",
            py::arg("dimension"))
        .def_property_readonly("dimension", &FilteredComplex::dimension,
            R"docstring(
                Current maximum dimension of a maximal simplex in the complex.
            )docstring")
        .def_property_readonly("max_dimension", &FilteredComplex::max_dim,
            R"docstring(
                Maximum dimension of simplex that this complex can store.
                Set during initialisation. 
            )docstring")
        .def_property_readonly("num_vertices", &FilteredComplex::N,
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
                start_dim : int
                    Dimension from which to start propagating (exclusive). 
                upwards : bool
                    If true then values are propagated upwards, downwards otherwise.
            )docstring",
            py::arg("start_dim"), py::arg("upwards"))
        .def("has_simplex", static_cast<bool (FilteredComplex::*)(vector<index_t> &) const>(&FilteredComplex::has_simplex),
            R"docstring(
                Check for membership of a simplex in the complex.
                
                Parameters
                ----------
                vertices : list[int]
                    Vertex labels of the simplex to check for. 
            )docstring",
            py::arg("vertices"))
        .def_property_readonly("simplices", &FilteredComplex::get_simplices,
            R"docstring(
                A list such that ``simplices[k]`` contains
                the dictionary of :math:`k`-simplices of the complex. 
                The key of a :math:`k`-simplex in ``simplices[k]`` is 
                the lexicographic index of that simplex with respect 
                to its vertex labels sorted in ascending order, 
                counting all possible sorted subsequences of :math:`(0, ..., N-1)`
                of length :math:`k`.
            )docstring")
        .def("get_label_from_vertex_labels", &FilteredComplex::get_label_from_vertex_labels,
             R"docstring(
                Returns the dictionary key of a simplex with respect to 
                its vertex labels sorted in ascending order, counting all 
                possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`, 
                where :math:`N` is the number of vertices in the complex. 
                The simplex need not be present in the simplicial complex.
                
                Parameters
                ----------
                vertices : list[int]
                    List of vertex labels of the simplex.
            )docstring",
            py::arg("vertices"))
        .def("propagate_colours", &FilteredComplex::propagate_colours,
             R"docstring(
                Function to make sure that simplex colours are consistent 
                with the colours of their vertices. You should call this 
                whenever you change the colour of any vertices. 
            )docstring")
        .def("serialised", &FilteredComplex::serialised,
             R"docstring(
                Serialised representation of the simplicial complex in a format
                suitable for persistent homology computations.

                Returns
                -------
                x : list[tuple[list[int], int, float, int]] 
                    List of simplices in the simplicial complex ordered by dimension 
                    followed by filtration time. Each simplex :math:`\sigma` is 
                    represented by a tuple containing the following items.

                    1.  A list containing the indices in ``x`` of the facets of :math:`\sigma`.
                    2.  The lexicographic key of :math:`\sigma` in the simplicial complex.
                    3.  The filtration time of :math:`\sigma`.
                    4.  The colour bitmask of :math:`\sigma`.
            )docstring")
        .def("__repr__",
            [](const FilteredComplex &K)
            {
                return "<" + std::to_string(K.dimension()) + "-dimensional simplicial complex with " + std::to_string(K.N) + " vertices>";
            });
    m.def("clique_complex", &FilteredComplex::clique_complex,
        R"docstring(
            Returns the :math:`k`-skeleton of the complete simplicial complex
            on :math:`n` vertices, with filtration values initialised to zero.
        )docstring",
        py::arg("n"), py::arg("k"));
    m.def("standard_simplex", &standard_simplex,
        R"docstring(
            Returns the simplicial complex corresponding to 
            the standard abstract :math:`n`-simplex.
        )docstring",
        py::arg("n"));
}