#include <chalc/filtration/filtration.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(filtration, m) {
	using namespace chalc;
	using namespace chalc::stl;
	namespace py = pybind11;
	m.doc() =
		R"docstring(
Module containing utilities to store and manipulate
abstract filtered simplicial complexes.
		)docstring";

	// forward declare classes
	py::class_<FilteredComplex> filtered_complex(m, "FilteredComplex");
	py::class_<FilteredComplex::Simplex, shared_ptr<FilteredComplex::Simplex>> simplex(m,
	                                                                                   "Simplex");

	/* FilteredComplex Interface */
	filtered_complex.doc() =
		R"docstring(
Class representing a filtered simplicial complex.
)docstring";
	filtered_complex
		.def(py::init<const index_t, const index_t>(),
	         R"docstring(
Construct a discrete filtered simplicial complex with default filtration time of 0.

Args:
	n : Number of vertices. Cannot be changed after initialisation.
	k : Maximum dimension of a simplex that the complex can have. This parameter is required for memory efficiency, and cannot be changed after initialisation.
)docstring",
	         py::arg("n"),
	         py::arg("k"))
		.def("add_simplex",
	         &FilteredComplex::add_simplex,
	         R"docstring(
Add a simplex to a filtered simplicial complex.

Args:
	vertices : List of vertex labels corresponding to existing vertices in the complex.
	filt_value : Filtration value to associate to the new simplex.

Note:
	Faces of the added simplex that are already present in the simplicial complex will have their filtration values reduced if necessary.
)docstring",
	         py::arg("vertices"),
	         py::arg("filt_value"))
		.def_property_readonly("num_simplices",
	                           &FilteredComplex::size,
	                           R"docstring(
The total number of simplices in the complex.
)docstring")
		.def_property_readonly("dimension",
	                           &FilteredComplex::dimension,
	                           R"docstring(
Current maximum dimension of a maximal simplex in the complex.
)docstring")
		.def_property_readonly("max_filtration_time",
	                           &FilteredComplex::max_filt_value,
	                           R"docstring(
Current maximum dimension of a maximal simplex in the complex.
)docstring")
		.def_readonly("max_dimension",
	                  &FilteredComplex::max_dim,
	                  R"docstring(
Maximum dimension of simplex that this complex can store. \
Set during initialisation.
)docstring")
		.def_readonly("num_vertices",
	                  &FilteredComplex::N,
	                  R"docstring(
Number of vertices in the simplicial complex. \
Set during initialisation.
)docstring")
		.def("propagate_filt_values",
	         &FilteredComplex::propagate_filt_values,
	         R"docstring(
Propagate filtration values upwards or downwards to ensure that every simplex appears after its faces. For example, setting the filtration values in dimension 1 and propagating upwards is akin to the Rips filtration.

Args:
	start_dim : Dimension from which to start propagating (exclusive).
	upwards : If true then values are propagated upwards, downwards otherwise. Defaults to true.
			 )docstring",
	         py::arg("start_dim"),
	         py::arg("upwards") = true)
		.def("has_simplex",
	         static_cast<bool (FilteredComplex::*)(vector<index_t>& v) const>(
				 &FilteredComplex::has_simplex),
	         R"docstring(
				Check for membership of a simplex in the complex.

				Args:
					vertices : Vertex labels of the simplex to check for.
			)docstring",
	         py::arg("vertices"))
		.def_property_readonly("simplices",
	                           &FilteredComplex::get_simplices,
	                           R"docstring(
				A list such that ``simplices[k]`` is a dictionary of handles to the :math:`k`-simplices in the complex. The key of a :math:`k`-simplex in ``simplices[k]`` is the lexicographic index of that simplex with respect to its vertex labels sorted in ascending order, counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.
			)docstring")
		.def("get_label_from_vertex_labels",
	         &FilteredComplex::get_label_from_vertex_labels,
	         R"docstring(
				Returns the dictionary key of a simplex with respect to its vertex labels sorted in ascending order, counting all
				possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`, where :math:`N` is the number of vertices in the complex. The simplex need not be present in the simplicial complex.

				Args:
					vertices : List of vertex labels of the simplex.
			)docstring",
	         py::arg("vertices"))
		.def("propagate_colours",
	         &FilteredComplex::propagate_colours,
	         R"docstring(
				Function to make sure that simplex colours are consistent with the colours of their vertices. You should call this whenever you change the colour of any vertices.
			)docstring")
		.def("serialised",
	         &FilteredComplex::serialised,
	         R"docstring(
				Serialised representation of the simplicial complex in a format suitable for persistent homology computations.

				:return:
					A list `x` of simplices in the simplicial complex ordered by dimension followed by filtration time. Each simplex :math:`\sigma` is represented by a tuple containing the following items.

					1.  A list containing the indices in `x` of the facets of :math:`\sigma`, sorted in ascending order.
					2.  The lexicographic key of :math:`\sigma` in the simplicial complex.
					3.  The filtration time of :math:`\sigma`.
					4.  The colour bitmask of :math:`\sigma`.
			)docstring")
		.def("is_filtration",
	         &FilteredComplex::is_filtration,
	         R"docstring(
			 Returns true if the filtration property is satisfied; that is, if each simplex has a filtration value at least as large as each of its faces.
			 )docstring")
		.def("__repr__", [](const FilteredComplex& K) {
			return "<" + std::to_string(K.dimension()) + "-dimensional simplicial complex with " +
		           std::to_string(K.N) + " vertices>";
		});

	/* Simplex Interface */
	simplex.doc() =
		R"docstring(
			Class representing a simplex in a filtered simplicial complex.
		)docstring";
	simplex
		.def_readonly("dimension",
	                  &FilteredComplex::Simplex::dim,
	                  R"docstring(
			Dimension of the simplex.
		)docstring")
		.def_readonly("label",
	                  &FilteredComplex::Simplex::label,
	                  R"docstring(
				Label of the simplex in its parent filtered complex. A :math:`k`-simplex :math:`\sigma` is labelled by the lexicographic index of :math:`\sigma` with respect to its vertex labels sorted in ascending order, counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.
			)docstring")
		.def_property_readonly("vertices",
	                           static_cast<vector<index_t> (FilteredComplex::Simplex::*)() const>(
								   &FilteredComplex::Simplex::get_vertex_labels),
	                           R"docstring(
				List of (sorted, ascending) vertex labels of the simplex.
			)docstring")
		.def_readwrite("filtration_value",
	                   &FilteredComplex::Simplex::value,
	                   R"docstring(
				Filtration value of the simplex. If you modify this value, you should call :meth:`propagate_filt_values() <chalc.filtration.FilteredComplex.propagate_filt_values>` from the parent complex to ensure that filtration times remain monotonic.
			)docstring")
		.def_property_readonly("colours",
	                           &FilteredComplex::Simplex::get_colours_as_int,
	                           R"docstring(
				Bitmask of colours in the simplex, where the rightmost bit represents the smallest colour index. More precisely, :math:`\mathrm{bitmask} = \sum_{c \in \mathrm{colours}(\sigma)} 2^c`. For example, a simplex having vertices of colours 0 and 1 has a colour bitmask of 3.
			)docstring")
		.def(
			"set_colour",
			[](const shared_ptr<FilteredComplex::Simplex>& s_ptr, index_t c) {
				if (s_ptr->dim == 0) {
					if (c < MAX_NUM_COLOURS) {
						s_ptr->set_colour(c);
					} else {
						throw std::domain_error("Colour index too large.");
					}
				} else {
					throw std::domain_error("Can't change colour unless simplex is a vertex.");
				}
			},
			R"docstring(
				Change the colour of a vertex.

				Raises:
					ValueError: If the simplex is not a vertex or if `colour >=` :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>`.

				Tip:
					It is recommended to call the member function :meth:`propagate_colours() <chalc.filtration.FilteredComplex.propagate_colours>` from the parent simplicial complex after changing the colour of a vertex.
			)docstring",
			py::arg("colour"))
		.def_property_readonly("facets",
	                           &FilteredComplex::Simplex::get_facets,
	                           R"docstring(
				Read-only list of handles to the facets of the simplex.
			)docstring")
		.def("__repr__", [](const shared_ptr<FilteredComplex::Simplex>& s_ptr) {
			return "<" + std::to_string(s_ptr->dim) + "-simplex>";
		});

	m.def("clique_complex",
	      &FilteredComplex::clique_complex,
	      R"docstring(
			Returns the :math:`k`-skeleton of the complete simplicial complex on :math:`n` vertices, with filtration values initialised to zero and all vertices coloured with the colour 0.
		)docstring",
	      py::arg("n"),
	      py::arg("k"));
	m.def("standard_simplex",
	      &standard_simplex,
	      R"docstring(
			Returns the filtered simplicial complex corresponding to the standard abstract :math:`n`-simplex, with filtration values initialised to zero and all vertices coloured with the colour 0.
		)docstring",
	      py::arg("n"));
}
