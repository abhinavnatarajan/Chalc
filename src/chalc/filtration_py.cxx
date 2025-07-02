#include <chalc/filtration/filtration.h>
#include <iterator>
#include <memory>
#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <type_traits>

namespace {
using chalc::Filtration;
using std::shared_ptr;

// A robust, flattening iterator for the vector of maps of simplices.
// It iterates over the Simplex smart pointers.
class SimplexIterator {
  private:
	using Vec = std::remove_reference_t<
		std::invoke_result_t<decltype(&Filtration::get_simplices), Filtration>>;
	using Map        = Vec::value_type;
	using VecConstIt = Vec::const_iterator;
	using MapConstIt = Map::const_iterator;

	VecConstIt vec_it, vec_end;
	MapConstIt map_it;

	// Advances the iterator to the next valid element, skipping empty maps.
	void advance_to_valid() {
		while (vec_it != vec_end && map_it == vec_it->end()) {
			++vec_it;
			if (vec_it != vec_end) {
				map_it = vec_it->begin();
			}
		}
	}

  public:
	using iterator_category = std::input_iterator_tag;
	using value_type        = shared_ptr<chalc::Filtration::Simplex>;
	using difference_type   = std::ptrdiff_t;
	using pointer           = const value_type*;
	using reference         = const value_type&;

	// Constructor for begin and other valid iterators.
	explicit SimplexIterator(const VecConstIt& vec_begin, const VecConstIt& vec_end) :
		vec_it(vec_begin),
		vec_end(vec_end) {
		if (vec_it != vec_end) {
			map_it = vec_it->begin();
		}
		advance_to_valid();
	}

	auto operator++() -> SimplexIterator& {
		++map_it;
		advance_to_valid();
		return *this;
	}

	auto operator*() const -> reference {
		return map_it->second;
	}

	auto operator->() const -> pointer {
		return &map_it->second;
	}

	auto operator==(const SimplexIterator& other) const -> bool {
		// If one is past-the-end then they are equal only if both are.
		if (vec_it == vec_end) {
			return (other.vec_it == other.vec_end && vec_it == other.vec_it);
		}
		return (vec_it == other.vec_it && vec_end == other.vec_end && map_it == other.map_it);
	}

	auto operator!=(const SimplexIterator& other) const -> bool {
		return !(*this == other);
	}
};
}  // namespace

PYBIND11_MODULE(filtration, m) {  // NOLINT
	using chalc::index_t;
	using chalc::MAX_NUM_COLOURS;
	using chalc::standard_simplex;
	using std::domain_error;
	using std::to_string;
	using std::vector;
	namespace py = pybind11;

	m.doc() =
		"Module containing utilities to store and manipulate abstract filtered simplicial complexes.";

	// Forward declare our classes.
	py::class_<Filtration> filtered_complex(m, "Filtration");
	py::class_<Filtration::Simplex, shared_ptr<Filtration::Simplex>> simplex(m, "Simplex");

	// Subclass Filtration from Sized, Iterable, and Container.
	auto Sized     = py::module_::import("collections.abc").attr("Sized");
	auto Iterable  = py::module_::import("collections.abc").attr("Iterable");
	auto Container = py::module_::import("collections.abc").attr("Container");
	Sized.attr("register")(filtered_complex);
	Iterable.attr("register")(filtered_complex);
	Container.attr("register")(filtered_complex);

	/* Filtration Interface */
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
			"skeleton",
			&Filtration::skeleton,
			R"docstring(Get a copy of the k-skeleton of the filtration.

Args:
	k: Dimension of the skeleton to return.

)docstring",
			py::arg("k")
		)
		.def(
			"add_simplex",
			&Filtration::add_simplex,
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
		// Requirements from Sized
		.def("__len__", &Filtration::size, "The total number of simplices in the complex.")
		// Requirements from Container
		.def(
			"__contains__",
			static_cast<bool (Filtration::*)(vector<index_t>& v) const>(&Filtration::has_simplex),
			R"docstring(Check for membership of a simplex in the complex.

Args:
	vertices : Vertex labels of the simplex to check for.

)docstring",
			py::arg("vertices")
		)
		// Requirements from Iterable
		.def(
			"__iter__",
			[](const Filtration& self) {
				const auto& simplices = self.get_simplices();
				return py::make_iterator(
					SimplexIterator(simplices.begin(), simplices.end()),
					SimplexIterator(simplices.end(), simplices.end())
				);
			},
			R"docstring(Iterate over the simplices in the complex, ordered by dimension.

There are no guarantees on the relative order of simplices with the same dimension.
)docstring",
			py::keep_alive<0, 1>()  // Keep the Filtration alive while iterating
		)
		.def_property_readonly(
			"dimension",
			&Filtration::dimension,
			"Current maximum dimension of a maximal simplex in the complex."
		)
		.def_property_readonly(
			"max_dimension",
			&Filtration::max_dimension,
			R"docstring(Maximum dimension of simplex that this complex can store.

Set during initialisation.

)docstring"
		)
		.def_property_readonly(
			"num_vertices",
			&Filtration::num_vertices,
			R"docstring(Number of vertices in the simplicial complex.

Set during initialisation.

)docstring"
		)
		.def(
			"propagate_filt_values",
			&Filtration::propagate_filt_values,
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
		.def_property_readonly(
			"simplices",
			&Filtration::get_simplices,
			R"docstring(A list such that ``simplices[k]`` is a dictionary of handles
to the :math:`k`-simplices in the complex.

The key of a :math:`k`-simplex in ``simplices[k]`` is the lexicographic index
of that simplex with respect to its vertex labels sorted in ascending order,
counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.

)docstring"
		)
		.def(
			"get_label_from_vertex_labels",
			&Filtration::get_label_from_vertex_labels,
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
			&Filtration::propagate_colours,
			R"docstring(Ensure that simplex colours are consistent
with the colours of their vertices.

You should call this whenever you change the colour of any vertices.

)docstring"
		)
		.def(
			"boundary_matrix",
			&Filtration::boundary_matrix,
			R"docstring(Compute the boundary matrix of the filtration.

Args:
	max_dimension: The maximum dimension of simplices to be considered.

Returns:
	A list `x` of simplices in the simplicial complex ordered by
	filtration time, dimension, and label (in that order).
	Each simplex :math:`\sigma` has dimension at most ``max_dimension`` and is
	represented by a tuple containing the following items.

	1.  A list containing the indices in `x` of the facets of :math:`\sigma`, sorted in ascending order.
	2.  The lexicographic key of :math:`\sigma` in the simplicial complex.
	3.  The filtration time of :math:`\sigma`.
	4.  The set of colours of the vertices of :math:`\sigma`.

)docstring",
			py::arg("max_dimension") = -1
		)
		.def(
			"is_filtration",
			&Filtration::is_filtration,
			R"docstring(Check if the filtration property is satisfied.

Returns true if each simplex has a filtration value at least as large as each of its faces.

)docstring"
		)
		.def("__repr__", [](const Filtration& self) {
			return "<" + to_string(self.dimension()) + "-dimensional simplicial complex with " +
		           to_string(self.num_vertices()) + " vertices>";
		});

	/* Simplex Interface */
	simplex.doc() = "Class representing a simplex in a filtered simplicial complex.";
	simplex
		.def_property_readonly(
			"dimension",
			&Filtration::Simplex::get_dim,
			R"docstring(
			Dimension of the simplex.
		)docstring"
		)
		.def_property_readonly(
			"label",
			&Filtration::Simplex::get_label,
			R"docstring(Label of the simplex in its parent filtered complex.

A :math:`k`-simplex :math:`\sigma` is labelled by the lexicographic index of :math:`\sigma`
with respect to its vertex labels sorted in ascending order,
counting all possible sorted subsequences of :math:`(0, ..., N-1)` of length :math:`k`.

)docstring"
		)
		.def_property_readonly(
			"vertices",
			static_cast<vector<index_t> (Filtration::Simplex::*)() const>(
				&Filtration::Simplex::get_vertex_labels
			),
			"List of (sorted, ascending) vertex labels of the simplex."
		)
		.def_property(
			"filtration_value",
			&Filtration::Simplex::get_value,
			&Filtration::Simplex::set_value,
			R"docstring(Filtration value of the simplex.

If you modify this value, you should call
:meth:`propagate_filt_values() <chalc.filtration.Filtration.propagate_filt_values>`
from the parent complex to ensure that filtration times remain monotonic.

)docstring"
		)
		.def_property_readonly(
			"colours",
			&Filtration::Simplex::get_colours_as_vec,
			"Set of colours of the vertices of the simplex."
		)
		.def(
			"set_colour",
			[](const shared_ptr<Filtration::Simplex>& s_ptr, index_t c) {
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
	:meth:`propagate_colours() <chalc.filtration.Filtration.propagate_colours>`
	from the parent simplicial complex after changing the colour of a vertex.

)docstring",
			py::arg("colour")
		)
		.def_property_readonly(
			"facets",
			[](const shared_ptr<Filtration::Simplex>& self) {
				auto                                         facets = self->get_facets();
				std::vector<shared_ptr<Filtration::Simplex>> result;
				result.reserve(facets.size());
				for (const auto& f: facets) {
					result.push_back(f->get_handle());
				}
				return result;
			},
			"Read-only list of handles to the facets of the simplex."
		)
		.def("__repr__", [](const shared_ptr<Filtration::Simplex>& s_ptr) {
			return "<" + to_string(s_ptr->get_dim()) + "-simplex>";
		});

	m.def(
		"complete_complex",
		&Filtration::complete_complex,
		R"docstring(Compute the :math:`k`-skeleton of the complete simplicial complex on :math:`n` vertices.

Filtration values are initialised to zero and all vertices coloured with the colour 0.

Raises:
	RuntimeError:
		If ``n<= 0`` or ``k >= n`` or ``k < 0``.

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
