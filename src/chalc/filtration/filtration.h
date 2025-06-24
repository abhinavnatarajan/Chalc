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
	#include <chalc/common.h>
	#include <map>
	#include <memory>
	#include <vector>

namespace chalc {
// The maximum number of colours that can be represented.
constexpr index_t MAX_NUM_COLOURS = 16;
using colours_t                   = std::bitset<MAX_NUM_COLOURS>;

class BinomialCoeffTable;

struct FilteredComplex {
	/* PUBLIC CLASSES OF FilteredComplex */

	struct Simplex;

	/* PUBLIC MEMBERS OF FilteredComplex */

	/* PUBLIC METHODS OF FilteredComplex */

	// Constructors
	// 1. Create a vertex set.
	FilteredComplex(const index_t num_vertices, const index_t max_dimension);

	// Get label of a simplex from the labels of its vertices
	auto get_label_from_vertex_labels(std::vector<index_t>& verts) const -> index_t;

	// Check if the complex has a specific simplex
	[[nodiscard]]
	auto has_simplex(const index_t dim, const index_t label) const -> bool;

	// Check if the complex has a specific simplex
	auto has_simplex(std::vector<index_t>& verts) const -> bool;

	// Add a new simplex to the complex by its vertices
	auto add_simplex(std::vector<index_t>& verts, const value_t filt_value) -> bool;

	// Propagate filtration values from start_dim
	void propagate_filt_values(const index_t start_dim, const bool upwards);

	// Returns a handle to the simplices
	[[nodiscard]]
	auto get_simplices() const noexcept
		-> const std::vector<std::map<index_t, std::shared_ptr<Simplex>>>&;

	// Returns the number of simplices in the complex
	[[nodiscard]]
	auto size() const noexcept -> index_t;

	// Returns the number of simplices of a specific dimension
	[[nodiscard]]
	auto size_in_dim(const index_t dim) const -> index_t;

	// Returns the current dimension
	[[nodiscard]]
	auto dimension() const noexcept -> index_t;

	// Returns the maximum dimension that can be stored
	[[nodiscard]]
	auto max_dimension() const noexcept -> index_t;

	// Returns the number of vertices of the complex
	[[nodiscard]]
	auto num_vertices() const noexcept -> index_t;

	// Returns the current maximum filtration value
	[[nodiscard]]
	auto max_filt_value() const noexcept -> value_t;

	// Returns a flat vectorised representation of the complex
	[[nodiscard]]
	auto serialised() const
		-> std::vector<std::tuple<std::vector<index_t>, index_t, value_t, std::vector<index_t>>>;

	// Returns the k-skeleton of the complete simplicial complex on n vertices.
	static auto complete_complex(const index_t n, const index_t k) -> FilteredComplex;

	// bitwise OR accumulates colours upwards from vertices
	void propagate_colours();

	// check if filtration property is satisfied
	[[nodiscard]]
	auto is_filtration() const -> bool;

  private:
	/* PRIVATE MEMBERS OF FilteredComplex */

	std::shared_ptr<const BinomialCoeffTable> binomial;  // binomial coefficients
	std::vector<std::map<index_t, std::shared_ptr<Simplex>>>
		simplices;               // std::vector whose kth element is a table of k-simplices,
	                             // labelled by their lexicographic index
	index_t num_simplices;       // total number of simplices
	index_t cur_dim;             // current maximum dimension of a maximal simplex
	value_t cur_max_filt_value;  // current maximum filtration value
	index_t n_vertices;          // number of vertices, labelled from 0 to n-1
	index_t max_dim;             // maximum dimension of any simplex in the complex

	/* PRIVATE METHODS OF FilteredComplex */

	// Checks that verts is a non-empty subsequence of (0, ..., n_vertices-1) and
	// verts.size() <= max_dim + 1.
	// WARNING: modifies verts.
	void validate_vertex_sequence(std::vector<index_t>& verts) const;

	// Get label of a simplex (possibly not in the complex) from the labels of
	// its vertices.
	// Assumes that verts is valid.
	[[nodiscard]]
	auto _get_label_from_vertex_labels(const std::vector<index_t>& verts) const -> index_t;

	// Check if the complex has a specific simplex.
	// Assumes that dim is valid.
	[[nodiscard]]
	auto _has_simplex(const index_t dim, const index_t label) const -> bool;

	// Check if the complex has a specific simplex.
	// Assumes that verts is valid.
	[[nodiscard]]
	auto _has_simplex(const std::vector<index_t>& verts) const -> bool;

	// Add a simplex to the complex with the specified vertices and filtration value.
	auto _add_simplex(const std::vector<index_t>& verts, const value_t filt_value)
		-> std::shared_ptr<Simplex>;

	// Min-fold filtration values downwards from start_dim.
	// Assumes that start_dim is valid.
	void propagate_filt_values_up(const index_t start_dim);

	// Max-fold filtration values upwards from start_dim.
	// Assumes that start_dim is valid.
	void propagate_filt_values_down(const index_t start_dim);
};

struct FilteredComplex::Simplex : public std::enable_shared_from_this<FilteredComplex::Simplex> {
	// friend class FilteredComplex;  // allow FilteredComplex to access private members
	/* PUBLIC MEMBERS OF Simplex */

	const index_t            label;       // label for the simplex
	const index_t            max_vertex;  // largest vertex label
	const index_t            dim;         // number of vertices minus 1
	value_t                  value;       // filtration value
	static constexpr value_t DEFAULT_FILT_VALUE = 0.0;

	/* PUBLIC METHODS OF Simplex */
	// Delete the default constructor
	Simplex() = delete;
	// Factory method - only way to create new simplex
	static auto _make_simplex(
		index_t label,
		index_t max_vertex,
		value_t value = DEFAULT_FILT_VALUE,
		const std::vector<std::shared_ptr<Simplex>>& facets =
			std::vector<std::shared_ptr<Simplex>>{}
	) -> std::shared_ptr<Simplex>;

	// Get a handle to this simplex
	auto get_handle() -> std::shared_ptr<Simplex>;

	// Return the sorted vertex labels of the simplex
	// Assumes that the simplex has valid faces
	auto get_vertex_labels() const -> std::vector<index_t>;

	// Return the indices of the facets of the simplex, [i]th element of the
	// result is the [i]th facet
	auto get_facet_labels() const -> std::vector<index_t>;

	// Return a const reference to the facets
	auto get_facets() const -> const std::vector<std::shared_ptr<Simplex>>&;

	// Return a const reference to the cofacets
	auto get_cofacets() const -> const std::vector<std::weak_ptr<Simplex>>&;

	void set_colour(index_t c);  // not inline since we export this

	auto get_colours_as_vec() -> std::vector<index_t>;

	inline void _set_colours(colours_t c);

	inline void add_colour(index_t c);

	inline void _add_colours(colours_t c);

	auto _get_colours() -> colours_t;

	inline void make_colourless();

	/* PRIVATE MEMBERS OF Simplex */
  private:
	const std::vector<std::shared_ptr<Simplex>>
										facets;    // pointers to the [i]th facets of the simplex
	std::vector<std::weak_ptr<Simplex>> cofacets;  // pointers to the cofacets of the simplex
	colours_t colours;  // bitmask representing the colours of its vertices
	// Constructor
	Simplex(
		index_t                                      label,
		index_t                                      max_vertex,
		value_t                                      value,
		const std::vector<std::shared_ptr<Simplex>>& facets
	);

	// Writes the sorted vertex labels of the simplex into a buffer
	// Assumes that the simplex has valid faces
	template <typename OutputIterator> void _get_vertex_labels(OutputIterator&& buf) const;
};

// The simplicial complex associated to the standard n-simplex.
auto standard_simplex(const index_t n) -> FilteredComplex;

}  // namespace chalc

#endif
