/*
This file is part of Chalc.

Chalc: Chromatic Alpha Complexes.
Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha
Complexes”. Online preprint available at http://arxiv.org/abs/2212.03128.
Accessed: 2023-02-28 22:07:23 UTC.
DOI: 10.48550/arXiv.2212.03128.

Project homepage:    http://github.com/abhinavnatarajan/Chalc

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

#pragma once
#ifndef FILTRATION_H
	#define FILTRATION_H

	#include <bitset>
	#include <cstdint>
	#include <memory>
	#include <unordered_map>
	#include <vector>

namespace chalc {
// The type for filtration values.
using value_t = double;
// The type for simplex labels.
using label_t = uint64_t;
// The maximum number of colours that can be represented.
// 16 colours means that we can represent a colour label
// in a uint16_t.
constexpr uint8_t MAX_NUM_COLOURS = 16;
// The type for colour bitmasks.
using colours_t = std::bitset<MAX_NUM_COLOURS>;
// The type for colour labels.
using colour_t = uint16_t;
// Type for dimension, number of vertices, number of simplices,
// and most other counts that should not become too big.
using index_t = int64_t;

namespace detail {
class BinomialCoeffTable {
	std::vector<std::vector<label_t>> B;

  public:
	// Constructs a binomial coefficient table that will hold all values of
	// i_C_j for i = 0, ..., n and j = 0, ..., min(k, floor(i/2)) for i <= n.
	BinomialCoeffTable(index_t n, index_t k);

	auto operator()(index_t n, index_t k) const -> label_t {
		return B[n][std::min(k, n - k)];
	}
};
}  // namespace detail

struct Filtration {
	struct Simplex;

	// Constructor - create a vertex set.
	Filtration(const index_t num_vertices, const index_t max_dimension);

	// Factory method - get k-skeleton
	[[nodiscard]]
	auto skeleton(const index_t k) const -> Filtration;  // exported

	// Get label of a simplex from the labels of its vertices.
	[[nodiscard]]
	auto get_label_from_vertex_labels(const std::vector<index_t>& verts) const
		-> label_t;  // exported

	// Check if the filtration has a specific simplex.
	[[nodiscard]]
	auto has_simplex(const index_t dim, const label_t label) const -> bool;

	// Check if the filtration has a specific simplex.
	[[nodiscard]]
	auto has_simplex(std::vector<index_t>& verts) const -> bool;  // exported

	// Add a new simplex to the filtration by its vertices.
	// Returns true if the simplex was added, false if it already exists.
	auto add_simplex(const std::vector<index_t>& verts, const value_t filt_value)
		-> bool;  // exported

	// Propagate filtration values from start_dim.
	void propagate_filt_values(const index_t start_dim, const bool upwards);  // exported

	// Returns a handle to the simplices.
	[[nodiscard]]
	auto get_simplices() const noexcept  // exported
		-> const std::vector<std::unordered_map<label_t, std::shared_ptr<Simplex>>>& {
		return simplices;
	}

	// Returns the number of simplices in the filtration.
	[[nodiscard]]
	auto size() const noexcept -> size_t {  // exported
		return num_simplices;
	}

	// Returns the current dimension.
	[[nodiscard]]
	auto dimension() const noexcept -> index_t {  // exported
		return cur_dim;
	}

	// Returns the maximum dimension of simplex that can be stored in this filtration.
	[[nodiscard]]
	auto max_dimension() const noexcept -> index_t {  // exported
		return max_dim;
	}

	// Returns the number of vertices of the filtration.
	[[nodiscard]]
	auto num_vertices() const noexcept -> index_t {  // exported
		return n_vertices;
	}

	// Returns the boundary matrix of the filtration.
	[[nodiscard]]
	auto boundary_matrix(index_t max_dimension = -1) const -> std::vector<
		std::tuple<std::vector<index_t>, label_t, value_t, std::vector<colour_t>>>;  // exported

	// Returns the k-skeleton of this filtration.
	// Returns the k-skeleton of the complete simplicial complex on n vertices.
	[[nodiscard]]
	static auto complete_complex(const index_t n, const index_t k) -> Filtration;  // exported

	// Bitwise OR accumulates colours upwards from vertices.
	void propagate_colours() noexcept;  // exported

	// Check if filtration values are monotonic.
	[[nodiscard]]
	auto is_filtration() const noexcept -> bool;  // exported

  private:
	/* PRIVATE MEMBERS OF Filtration */

	detail::BinomialCoeffTable binomial;  // binomial coefficients
	std::vector<std::unordered_map<label_t, std::shared_ptr<Simplex>>>
		simplices;          // std::vector whose kth element is a table of k-simplices,
	                        // labelled by their lexicographic index
	size_t  num_simplices;  // total number of simplices
	index_t cur_dim;        // current maximum dimension of a maximal simplex
	index_t n_vertices;     // number of vertices, labelled from 0 to n-1
	index_t max_dim;        // maximum dimension of any simplex in the filtration

	/* PRIVATE METHODS OF Filtration */

	// Checks that verts is a non-empty subset of {0, ..., n_vertices-1} and
	// verts.size() <= max_dim + 1.
	// Returns a sorted copy of verts.
	[[nodiscard]]
	auto validated_vertex_sequence(const std::vector<index_t>& verts) const -> std::vector<index_t>;

	// Get label of a simplex (possibly not in the filtration) from the labels of
	// its vertices.
	// Assumes that verts is valid.
	[[nodiscard]]
	auto lex_label(const std::vector<index_t>& verts) const -> label_t;

	// Check if the filtration has a specific simplex.
	// Assumes that dim is valid.
	[[nodiscard]]
	auto has_simplex_unchecked(const index_t dim, const label_t label) const -> bool;

	// Check if the filtration has a specific simplex.
	// Assumes that verts is valid.
	[[nodiscard]]
	auto has_simplex_unchecked(const std::vector<index_t>& verts) const noexcept -> bool;

	// Add a simplex to the filtration with the specified vertices and filtration value.
	auto add_simplex_unchecked(const std::vector<index_t>& verts, const value_t filt_value)
		-> std::shared_ptr<Simplex>;

	// Min-fold filtration values downwards from start_dim.
	// If start_dim > cur_dim, does nothing.
	void propagate_filt_values_up(const index_t start_dim) noexcept;

	// Max-fold filtration values upwards from start_dim.
	// If start_dim == 0, does nothing.
	void propagate_filt_values_down(const index_t start_dim);
};

struct Filtration::Simplex : public std::enable_shared_from_this<Filtration::Simplex> {
	friend class Filtration;  // allow Filtration to access private members
	/* PUBLIC MEMBERS OF Simplex */

	static constexpr value_t DEFAULT_FILT_VALUE = 0.0;

	// Delete the default constructor
	Simplex() = delete;

	// Get a handle to this simplex
	[[nodiscard]]
	auto get_handle() -> std::shared_ptr<Simplex>;

	// Return the sorted vertex labels of the simplex.
	// Assumes that the simplex has valid faces.
	[[nodiscard]]
	auto get_vertex_labels() const -> std::vector<index_t>;  // exported

	// Return the indices of the facets of the simplex, [i]th element of the
	// result is the [i]th facet.
	[[nodiscard]]
	auto get_facet_labels() const -> std::vector<label_t>;

	// Return a const reference to the facets.
	[[nodiscard]]
	auto get_facets() const noexcept -> const std::vector<Simplex*>& {  // exported
		return facets;
	}

	// Return a const reference to the cofacets.
	[[nodiscard]]
	auto get_cofacets() const noexcept -> const std::vector<Simplex*>& {
		return cofacets;
	}

	// Set a monochromatic colour for the simplex.
	void set_colour(colour_t c) {  // exported
		colours.reset().set(c);
	}

	// Colours of the simplex as a vector of colour labels.
	[[nodiscard]]
	auto get_colours_as_vec() const -> std::vector<colour_t>;  // exported

	// Dimension of this simplex.
	[[nodiscard]]
	auto get_dim() const noexcept -> index_t {  // exported
		return m_dim;
	}

	// Get label of the simplex.
	[[nodiscard]]
	auto get_label() const noexcept -> label_t {  // exported
		return m_label;
	}

	// Reference to the filtration value of the simplex.
	[[nodiscard]]
	auto value() noexcept -> value_t& {
		return m_filt_value;
	}

	// Filtration value of the simplex.
	[[nodiscard]]
	auto get_value() const noexcept -> value_t {  // exported
		return m_filt_value;
	}

	// Set the filtration value of the simplex.
	auto set_value(value_t v) noexcept -> void {  // exported
		m_filt_value = v;
	}

	// Add a colour to the simplex from a colour label.
	void add_colour(colour_t c) {
		colours.set(c);
	}

	// Remove all colours from the simplex.
	void make_colourless() noexcept {
		colours.reset();
	}

	/* PRIVATE MEMBERS OF Simplex */
  private:
	label_t               m_label;       // label for the simplex
	index_t               m_max_vertex;  // largest vertex label
	index_t               m_dim;         // number of vertices minus 1
	value_t               m_filt_value;  // filtration value
	std::vector<Simplex*> facets;        // pointers to the [i]th facets of the simplex
	std::vector<Simplex*> cofacets;      // pointers to the cofacets of the simplex
	colours_t             colours;       // bitmask representing the colours of its vertices

	// Constructor for internal use only.
	Simplex(label_t label, index_t max_vertex, value_t value, const std::vector<Simplex*>& facets);

	// Factory method for internal use.
	static auto make_simplex(
		label_t label,
		index_t max_vertex,
		value_t value                       = DEFAULT_FILT_VALUE,
		const std::vector<Simplex*>& facets = std::vector<Simplex*>{}
	) -> std::shared_ptr<Simplex>;

	// Writes the sorted vertex labels of the simplex into a buffer
	// Assumes that the simplex has valid faces
	template <typename OutputIterator> void _get_vertex_labels(OutputIterator&& buf) const;

	// Set the colours of the simplex from a bitmask.
	void set_colours_bitmask(colours_t c) noexcept {
		colours.reset();
		add_colours_bitmask(c);
	}

	// Add multiple colours to the simplex from a bitmask.
	void add_colours_bitmask(colours_t c) noexcept {
		colours |= c;
	}

	// Get the colours of the simplex as a bitmask.
	auto get_colours_bitmask() const noexcept -> const colours_t& {
		return colours;
	};
};

// The simplicial complex associated to the standard n-simplex.
auto standard_simplex(const index_t n) -> Filtration;  // exported

}  // namespace chalc

#endif
