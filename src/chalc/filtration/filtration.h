#pragma once
#ifndef FILTRATION_H
	#define FILTRATION_H

	#include <bitset>
	#include <cstdint>
	#include <iterator>
	#include <memory>
	#include <type_traits>
	#include <unordered_map>
	#include <vector>

namespace chalc {
// The type for filtration values.
using Value = double;
// The type for simplex labels.
using Label = uint64_t;
// The maximum number of colours that can be represented.
// 16 colours means that we can represent a colour label
// in a uint16_t.
constexpr uint8_t MAX_NUM_COLOURS = 16;
// The type for colour bitmasks.
using Colours = std::bitset<MAX_NUM_COLOURS>;
// The type for colour labels.
using Colour = uint16_t;
// Type for dimension, number of vertices, number of simplices,
// and most other counts that should not become too big.
using Index = int64_t;

namespace detail {
class BinomialCoeffTable {
	std::vector<std::vector<Label>> B;

  public:
	// Constructs a binomial coefficient table that will hold all values of
	// i_C_j for i = 0, ..., n and j = 0, ..., min(k, floor(i/2)) for i <= n.
	BinomialCoeffTable(Index n, Index k);

	auto operator()(Index n, Index k) const -> Label {
		return B[n][std::min(k, n - k)];
	}
};
}  // namespace detail

struct Filtration {
	struct Simplex;
	using SimplexHandle = std::unique_ptr<Simplex>;
	class SimplicesIterator;


  private:
	/* PRIVATE MEMBERS OF Filtration */

	detail::BinomialCoeffTable            binomial_;   // binomial coefficients
	std::vector<SimplexHandle> simplices_;  // container for the simplices
	std::vector<std::unordered_map<Label, Simplex*>>
		simplices_map_;      // std::vector whose kth element is a table of k-simplex handles,
	                         // labelled by their lexicographic index
	Index cur_dim_;        // current maximum dimension of a maximal simplex
	Index num_vertices_;   // number of vertices, labelled from 0 to n-1
	Index max_dimension_;  // maximum dimension of any simplex in the filtration

  public:
	// Constructor - create a vertex set.
	Filtration(const Index num_vertices, const Index max_dimension);
	// vector<unique_ptr> has a copy constructor that will fail to compile
	// but is nevertheless defined. Hence the copy constructor of Filtration
	// is not deleted (i.e., std::is_copy_constructible<Filtration>::value == true),
	// but this is bad because Pybind11 will try to make a copy constructor for our class.

	Filtration(const Filtration& other)                   = delete;
	auto operator=(const Filtration& other) -> Filtration = delete;
	Filtration(Filtration&& other)                        = default;
	auto operator=(Filtration&& other) -> Filtration&     = default;
	~Filtration()                                         = default;

	// Factory method - get k-skeleton
	[[nodiscard]]
	auto skeleton(const Index k) const -> Filtration;  // exported

	// Get label of a simplex from the labels of its vertices.
	[[nodiscard]]
	auto get_label_from_vertex_labels(const std::vector<Index>& verts) const
		-> Label;  // exported

	// Check if the filtration has a specific simplex.
	[[nodiscard]]
	auto has_simplex(const Index dim, const Label label) const -> bool;

	// Check if the filtration has a specific simplex.
	[[nodiscard]]
	auto has_simplex(std::vector<Index>& verts) const -> bool;  // exported

	// Add a new simplex to the filtration by its vertices.
	// Returns true if the simplex was added, false if it already exists.
	auto add_simplex(const std::vector<Index>& verts, const Value filt_value)
		-> bool;  // exported

	// Propagate filtration values from start_dim.
	void propagate_filt_values(const Index start_dim, const bool upwards);  // exported

	// Returns a handle to the simplices.
	[[nodiscard]]
	auto simplices() const noexcept  // exported
		-> const std::vector<std::unordered_map<Label, Simplex*>>& {
		return simplices_map_;
	}

	// Returns an iterator over the simplices in the filtration.
	[[nodiscard]]
	auto simplices_begin() -> SimplicesIterator;

	// Returns a past-the-end iterator over the simplices in the filtration.
	[[nodiscard]]
	auto simplices_end() -> SimplicesIterator;

	// Returns the number of simplices in the filtration.
	[[nodiscard]]
	auto size() const noexcept -> size_t {  // exported
		return simplices_.size();
	}

	// Returns the current dimension.
	[[nodiscard]]
	auto dimension() const noexcept -> Index {  // exported
		return cur_dim_;
	}

	// Returns the maximum dimension of simplex that can be stored in this filtration.
	[[nodiscard]]
	auto max_dimension() const noexcept -> Index {  // exported
		return max_dimension_;
	}

	// Returns the number of vertices of the filtration.
	[[nodiscard]]
	auto num_vertices() const noexcept -> Index {  // exported
		return num_vertices_;
	}

	// Returns the boundary matrix of the filtration.
	[[nodiscard]]
	auto boundary_matrix(Index max_dimension = -1) const -> std::vector<
		std::tuple<std::vector<Index>, Label, Value, std::vector<Colour>>>;  // exported

	// Returns the k-skeleton of this filtration.
	// Returns the k-skeleton of the complete simplicial complex on n vertices.
	[[nodiscard]]
	static auto complete_complex(const Index n, const Index k) -> Filtration;  // exported

	// Bitwise OR accumulates colours upwards from vertices.
	void propagate_colours() noexcept;  // exported

	// Check if filtration values are monotonic.
	[[nodiscard]]
	auto is_filtration() const noexcept -> bool;  // exported

  private:
	// Checks that verts is a non-empty subset of {0, ..., n_vertices-1} and
	// verts.size() <= max_dim + 1.
	// Returns a sorted copy of verts.
	[[nodiscard]]
	auto validated_vertex_sequence(const std::vector<Index>& verts) const -> std::vector<Index>;

	// Get label of a simplex (possibly not in the filtration) from the labels of
	// its vertices.
	// Assumes that verts is valid.
	[[nodiscard]]
	auto lex_label(const std::vector<Index>& verts) const -> Label;

	// Check if the filtration has a specific simplex.
	// Assumes that dim is valid.
	[[nodiscard]]
	auto has_simplex_unchecked(const Index dim, const Label label) const -> bool;

	// Check if the filtration has a specific simplex.
	// Assumes that verts is valid.
	[[nodiscard]]
	auto has_simplex_unchecked(const std::vector<Index>& verts) const noexcept -> bool;

	// Add a simplex to the filtration with the specified vertices and filtration value.
	auto add_simplex_unchecked(const std::vector<Index>& verts, const Value filt_value)
		-> Simplex*;

	// Min-fold filtration values downwards from start_dim.
	// If start_dim > cur_dim, does nothing.
	void propagate_filt_values_up(const Index start_dim) noexcept;

	// Max-fold filtration values upwards from start_dim.
	// If start_dim == 0, does nothing.
	void propagate_filt_values_down(const Index start_dim);
};

struct Filtration::Simplex {
	friend class Filtration;  // allow Filtration to access private members
	/* PUBLIC MEMBERS OF Simplex */

	static constexpr Value DEFAULT_FILT_VALUE = 0.0;

	// Delete the default constructor
	Simplex() = delete;

	// Return the sorted vertex labels of the simplex.
	// Assumes that the simplex has valid faces.
	[[nodiscard]]
	auto vertex_labels() const -> std::vector<Index>;  // exported

	// Return the indices of the facets of the simplex, [i]th element of the
	// result is the [i]th facet.
	[[nodiscard]]
	auto facet_labels() const -> std::vector<Label>;

	// Return a const reference to the facets.
	[[nodiscard]]
	auto facets() const noexcept -> const std::vector<Simplex*>& {  // exported
		return facets_;
	}

	// Return a const reference to the cofacets.
	[[nodiscard]]
	auto cofacets() const noexcept -> const std::vector<Simplex*>& {
		return cofacets_;
	}

	// Set a monochromatic colour for the simplex.
	void set_colour(Colour c) {  // exported
		colours_.reset().set(c);
	}

	// Colours of the simplex as a vector of colour labels.
	[[nodiscard]]
	auto colours() const -> std::vector<Colour>;  // exported

	// Dimension of this simplex.
	[[nodiscard]]
	auto dimension() const noexcept -> Index {  // exported
		return dim_;
	}

	// Get label of the simplex.
	[[nodiscard]]
	auto label() const noexcept -> Label {  // exported
		return label_;
	}

	// Reference to the filtration value of the simplex.
	[[nodiscard]]
	auto value() noexcept -> Value& {
		return filt_value_;
	}

	// Set the filtration value of the simplex.
	auto set_value(Value v) noexcept -> void {  // exported
		filt_value_ = v;
	}

	// Remove all colours from the simplex.
	void make_colourless() noexcept {
		colours_.reset();
	}

	Simplex(const Simplex& other)                   = delete;
	auto operator=(const Simplex& other) -> Simplex = delete;
	Simplex(Simplex&& other)                        = default;
	auto operator=(Simplex&& other) -> Simplex&     = default;
	~Simplex()                                      = default;

	/* PRIVATE MEMBERS OF Simplex */
  private:
	Label               label_;       // label for the simplex
	Index               max_vertex_;  // largest vertex label
	Index               dim_;         // number of vertices minus 1
	Value               filt_value_;  // filtration value
	std::vector<Simplex*> facets_;      // pointers to the [i]th facets of the simplex
	std::vector<Simplex*> cofacets_;    // pointers to the cofacets of the simplex
	Colours             colours_;     // bitmask representing the colours of its vertices

	// Constructor for internal use only.
	Simplex(Label label, Index max_vertex, Value value, const std::vector<Simplex*>& facets);

	// Factory method for internal use.
	static auto make_simplex(
		Label label,
		Index max_vertex,
		Value value                       = DEFAULT_FILT_VALUE,
		const std::vector<Simplex*>& facets = std::vector<Simplex*>{}
	) -> SimplexHandle;

	// Writes the sorted vertex labels of the simplex into a buffer
	// Assumes that the simplex has valid faces
	template <typename OutputIterator> void vertex_labels_(OutputIterator&& buf) const;

	// Set the colours of the simplex from a bitmask.
	void set_colours_bitmask(Colours c) noexcept {
		colours_.reset();
		add_colours_bitmask(c);
	}

	// Add multiple colours to the simplex from a bitmask.
	void add_colours_bitmask(Colours c) noexcept {
		colours_ |= c;
	}

	// Get the colours of the simplex as a bitmask.
	[[nodiscard]]
	auto colours_bitmask() const noexcept -> const Colours& {
		int i = 5;
		return colours_;
	};
};

// A robust, flattening iterator for the vector of maps of simplices.
// It iterates over the Simplex smart pointers.
class Filtration::SimplicesIterator {
  private:
	using Vec =
		std::remove_reference_t<std::invoke_result_t<decltype(&Filtration::simplices), Filtration>>;
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
	using value_type        = Filtration::Simplex;
	using difference_type   = std::ptrdiff_t;
	using pointer           = const value_type*;
	using reference         = const value_type&;

	// Constructor for begin and other valid iterators.
	explicit SimplicesIterator(const VecConstIt& vec_begin, const VecConstIt& vec_end) :
		vec_it(vec_begin),
		vec_end(vec_end) {
		if (vec_it != vec_end) {
			map_it = vec_it->begin();
		}
		advance_to_valid();
	}

	auto operator++() -> SimplicesIterator& {
		++map_it;
		advance_to_valid();
		return *this;
	}

	auto operator*() const -> reference {
		return *map_it->second;
	}

	auto operator->() const -> pointer {
		return map_it->second;
	}

	auto operator==(const SimplicesIterator& other) const -> bool {
		// If one is past-the-end then they are equal only if both are.
		if (vec_it == vec_end) {
			return (other.vec_it == other.vec_end && vec_it == other.vec_it);
		}
		return (vec_it == other.vec_it && vec_end == other.vec_end && map_it == other.map_it);
	}

	auto operator!=(const SimplicesIterator& other) const -> bool {
		return !(*this == other);
	}
};

// The simplicial complex associated to the standard n-simplex.
auto standard_simplex(const Index n) -> Filtration;  // exported

}  // namespace chalc

#endif
