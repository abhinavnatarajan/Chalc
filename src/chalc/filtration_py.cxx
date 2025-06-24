#include <chalc/filtration/filtration.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(filtration, m) {  // NOLINT
	using chalc::FilteredComplex;
	using chalc::index_t;
	using chalc::MAX_NUM_COLOURS;
	using chalc::standard_simplex;
	using std::domain_error;
	using std::shared_ptr;
	using std::to_string;
	using std::vector;
	namespace py = pybind11;

	m.doc() =
		"Module containing utilities to store and manipulate abstract filtered simplicial complexes.";

	// forward declare classes
	py::class_<FilteredComplex> filtered_complex(m, "FilteredComplex");
	py::class_<FilteredComplex::Simplex, shared_ptr<FilteredComplex::Simplex>> simplex(
		m,
		"Simplex"
	);

	/* FilteredComplex Interface */
	filtered_complex.doc() = "Class representing a filtered simplicial complex.";
	filtered_complex
		.def(
			py::init<const index_t, const index_t>(),
			R"docstring(Construct a discrete filtered simplicial complex with default filtration time of 0.

Args:
	n : Number of vertices. Cannot be changed after initialisation.
	k : Maximum dimension of a simplex that the complex can have.
		This parameter is required for memory efficiency,
		and cannot be changed after initialisation.

)docstring",
			py::arg("n"),
			py::arg("k")
		)
		.def(
			"add_simplex",
			&FilteredComplex::add_simplex,
			R"docstring(Add a simplex to a filtered simplicial complex.

Args:
	vertices : List of vertex labels corresponding to existing vertices in the complex.
	filt_value : Filtration value to associate to the new simplex.

Note:
	Faces of the added simplex that are already present
	in the simplicial complex will have their filtration values reduced if necessary.

)docstring",
			py::arg("vertices"),
			py::arg("filt_value")
		)
		.def_property_readonly(
			"num_simplices",
			&FilteredComplex::size,
			"The total number of simplices in the complex."
		)
		.def_property_readonly(
			"dimension",
			&FilteredComplex::dimension,
			"Current maximum dimension of a maximal simplex in the complex."
		)
		.def_property_readonly(
			"max_filtration_time",
			&FilteredComplex::max_filt_value,
			"Current maximum filtration value in the complex."
		)
		.def_property_readonly(
			"max_dimension",
			&FilteredComplex::max_dimension,
			R"docstring(Maximum dimension of simplex that this complex can store.

Set during initialisation.

)docstring"
		)
		.def_property_readonly(
			"num_vertices",
			&FilteredComplex::num_vertices,
			R"docstring(Number of vertices in the simplicial complex.

Set during initialisation.

)docstring"
		)
		.def(
			"propagate_filt_values",
			&FilteredComplex::propagate_filt_values,
			R"docstring(Propagate filtration values upwards or downwards.

If propagating upwards, the filtration value of each simplex
will be set to the maximum filtration value of any of its faces.
If propagating downwards, the filtration value of each simplex
will be set to the minimum filtration value of any of its cofaces.
For example, setting the filtration values in dimension 1
to edge lengths and propagating upwards gives Rips-type filtration.

Args:
	start_dim : Dimension from which to start propagating (exclusive).
	upwards : If true then values are propagated upwards, downwards otherwise. Defaults to true.

)docstring",
			py::arg("start_dim"),
			py::arg("upwards") = true
		)
		.def(
			"has_simplex",
			static_cast<bool (FilteredComplex::*)(vector<index_t>& v) const>(
				&FilteredComplex::has_simplex
			),
			R"docstring(Check for membership of a simplex in the complex.

Args:
	vertices : Vertex labels of the simplex to check for.

)docstring",
			py::arg("vertices")
		)
		.def_property_readonly(
			"simplices",
			&FilteredComplex::get_simplices,
			R"docstring(A list such that ``simplices[k]`` is a dictionary of handles
to the :math:`k`-simplices in the complex.

The key of a :math:`k`-simplex in ``simplices[k]`` is the lexicographic index
of that simplex with respect to its vertex labels sorted in ascending order,
counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.

)docstring"
		)
		.def(
			"get_label_from_vertex_labels",
			&FilteredComplex::get_label_from_vertex_labels,
			R"docstring(Get the dictionary key of a simplex.

The key of a simplex is its lexicographic index with respect to
its vertex labels sorted in ascending order, counting all
possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`,
where :math:`N` is the number of vertices in the complex.
The simplex need not be present in the simplicial complex.

Args:
	vertices : List of vertex labels of the simplex.

)docstring",
			py::arg("vertices")
		)
		.def(
			"propagate_colours",
			&FilteredComplex::propagate_colours,
			R"docstring(Ensure that simplex colours are consistent
with the colours of their vertices.

You should call this whenever you change the colour of any vertices.

)docstring"
		)
		.def(
			"serialised",
			&FilteredComplex::serialised,
			R"docstring(Compute the boundary matrix of the simplicial complex.

:return:
	A list `x` of simplices in the simplicial complex ordered by
	filtration time, dimension, and label (in that order).
	Each simplex :math:`\sigma` is represented by a tuple containing the following items.

	1.  A list containing the indices in `x` of the facets of :math:`\sigma`, sorted in ascending order.
	2.  The lexicographic key of :math:`\sigma` in the simplicial complex.
	3.  The filtration time of :math:`\sigma`.
	4.  The set of colours of the vertices of :math:`\sigma`.

)docstring"
		)
		.def(
			"is_filtration",
			&FilteredComplex::is_filtration,
			R"docstring(Check if the filtration property is satisfied.

Returns true if each simplex has a filtration value at least as large as each of its faces.

)docstring"
		)
		.def("__repr__", [](const FilteredComplex& K) {
			return "<" + to_string(K.dimension()) + "-dimensional simplicial complex with " +
		           to_string(K.num_vertices()) + " vertices>";
		});

	/* Simplex Interface */
	simplex.doc() = "Class representing a simplex in a filtered simplicial complex.";
	simplex
		.def_property_readonly(
			"dimension",
			&FilteredComplex::Simplex::get_dim,
			R"docstring(
			Dimension of the simplex.
		)docstring"
		)
		.def_property_readonly(
			"label",
			&FilteredComplex::Simplex::get_label,
			R"docstring(Label of the simplex in its parent filtered complex.

A :math:`k`-simplex :math:`\sigma` is labelled by the lexicographic index of :math:`\sigma`
with respect to its vertex labels sorted in ascending order,
counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.

)docstring"
		)
		.def_property_readonly(
			"vertices",
			static_cast<vector<index_t> (FilteredComplex::Simplex::*)() const>(
				&FilteredComplex::Simplex::get_vertex_labels
			),
			"List of (sorted, ascending) vertex labels of the simplex."
		)
		.def_property(
			"filtration_value",
			&FilteredComplex::Simplex::get_value,
			&FilteredComplex::Simplex::set_value,
			R"docstring(Filtration value of the simplex.

If you modify this value, you should call
:meth:`propagate_filt_values() <chalc.filtration.FilteredComplex.propagate_filt_values>`
from the parent complex to ensure that filtration times remain monotonic.

)docstring"
		)
		.def_property_readonly(
			"colours",
			&FilteredComplex::Simplex::get_colours_as_vec,
			"List of colours of the vertices of the simplex."
		)
		.def(
			"set_colour",
			[](const shared_ptr<FilteredComplex::Simplex>& s_ptr, index_t c) {
				if (s_ptr->get_dim() == 0) {
					if (c < MAX_NUM_COLOURS) {
						s_ptr->set_colour(c);
					} else {
						throw domain_error("Colour index too large.");
					}
				} else {
					throw domain_error("Can't change colour unless simplex is a vertex.");
				}
			},
			R"docstring(Change the colour of a vertex.

Raises:
	ValueError: If the simplex is not a vertex or if
	`colour >=` :attr:`MaxColoursChromatic <chalc.chromatic.MaxColoursChromatic>`.

Tip:
	It is recommended to call the member function
	:meth:`propagate_colours() <chalc.filtration.FilteredComplex.propagate_colours>`
	from the parent simplicial complex after changing the colour of a vertex.

)docstring",
			py::arg("colour")
		)
		.def_property_readonly(
			"facets",
			&FilteredComplex::Simplex::get_facets,
			"Read-only list of handles to the facets of the simplex."
		)
		.def("__repr__", [](const shared_ptr<FilteredComplex::Simplex>& s_ptr) {
			return "<" + to_string(s_ptr->get_dim()) + "-simplex>";
		});

	m.def(
		"complete_complex",
		&FilteredComplex::complete_complex,
		R"docstring(Compute the :math:`k`-skeleton of the complete simplicial complex on :math:`n` vertices.

Filtration values are initialised to zero and all vertices coloured with the colour 0.

)docstring",
		py::arg("n"),
		py::arg("k")
	);
	m.def(
		"standard_simplex",
		&standard_simplex,
		R"docstring(Compute the filtered simplicial complex corresponding to the standard abstract :math:`n`-simplex.

Filtration values are initialised to zero and all vertices coloured with the colour 0.

)docstring",
		py::arg("n")
	);
}
