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

namespace chalc {
constexpr index_t                    MAX_NUM_COLOURS = 64;
typedef std::bitset<MAX_NUM_COLOURS> colours_t;

class BinomialCoeffTable;

struct FilteredComplex {
	/* PUBLIC CLASSES OF FilteredComplex */

	struct Simplex;

	/* PUBLIC MEMBERS OF FilteredComplex */

	index_t N;        // number of vertices, labelled from 0 to n-1
	index_t max_dim;  // maximum dimension of any simplex in the complex

	/* PUBLIC METHODS OF FilteredComplex */

	// Constructors
	// 1. Create a vertex set.
	FilteredComplex(const index_t num_vertices, const index_t max_dimension);

	// Get label of a simplex from the labels of its vertices
	index_t get_label_from_vertex_labels(std::vector<index_t>& verts) const;

	// Check if the complex has a specific simplex
	bool has_simplex(const index_t dim, const index_t label) const;

	// Check if the complex has a specific simplex
	bool has_simplex(std::vector<index_t>& verts) const;

	// Add a new simplex to the complex by its vertices
	bool add_simplex(std::vector<index_t>& verts, const value_t filt_value);

	// Propagate filtration values from start_dim
	void propagate_filt_values(const index_t start_dim, const bool upwards);

	// Returns a handle to the simplices
	const std::vector<std::map<index_t, std::shared_ptr<Simplex>>>& get_simplices() const noexcept;

	// Returns the number of simplices in the complex
	index_t size() const noexcept;

	// Returns the number of simplices of a specific dimension
	index_t size_in_dim(const index_t dim) const;

	// Returns the current dimension
	index_t dimension() const noexcept;

	// Returns the current maximum filtration value
	value_t max_filt_value() const noexcept;

	// Returns a flat vectorised representation of the complex
	std::vector<std::tuple<std::vector<index_t>, index_t, value_t, unsigned long long int>>
	serialised() const;

	// Returns the k-skeleton of the clique complex on n vertices
	static FilteredComplex clique_complex(const index_t n, const index_t k);

	// bitwise OR accumulates colours upwards from vertices
	void propagate_colours();

	// check if filtration property is satisfied
	bool is_filtration() const;

  private:
	/* PRIVATE MEMBERS OF FilteredComplex */

	std::shared_ptr<const BinomialCoeffTable> binomial;  // binomial coefficients
	std::vector<std::map<index_t, std::shared_ptr<Simplex>>>
		simplices;               // std::vector whose kth element is a table of k-simplices,
	                             // labelled by their lexicographic index
	index_t num_simplices;       // total number of simplices
	index_t cur_dim;             // current maximum dimension of a maximal simplex
	value_t cur_max_filt_value;  // current maximum filtration value

	/* PRIVATE METHODS OF FilteredComplex */

	// Checks that verts is a non-empty subsequence of (0, ..., N-1) and
	// verts.size() <= max_dim + 1 does not leave verts unmodified
	void check_vertex_sequence_is_valid(std::vector<index_t>& verts) const;

	// Check that 0 <= dim <= max_dim
	void check_dimension_is_valid(const index_t dim) const;

	// Get label of a simplex (possibly not in the complex) from the labels of
	// its vertices Assumes that verts is valid
	index_t _get_label_from_vertex_labels(const std::vector<index_t>& verts) const;

	// Check if the complex has a specific simplex
	// Assumes that dim is valid
	bool _has_simplex(const index_t dim, const index_t label) const;

	// Check if the complex has a specific simplex
	// Assumes that verts is valid
	bool _has_simplex(const std::vector<index_t>& verts) const;

	std::shared_ptr<Simplex> _add_simplex(const std::vector<index_t>& verts,
	                                      const value_t               filt_value);

	// min accumulate filtration values downwards from start_dim
	// Assumes that start_dim is valid
	void propagate_filt_values_up(const index_t start_dim);

	// max accumulate filtration values upwards from start_dim
	// Assumes that start_dim is valid
	void propagate_filt_values_down(const index_t start_dim);
};

struct FilteredComplex::Simplex : public std::enable_shared_from_this<FilteredComplex::Simplex> {

	/* PUBLIC MEMBERS OF Simplex */

	const index_t            label;       // label for the simplex
	const index_t            max_vertex;  // largest vertex label
	const index_t            dim;         // number of vertices minus 1
	value_t                  value;       // filtration value
	static constexpr value_t DEFAULT_FILT_VALUE = 0.0;

	/* PUBLIC METHODS OF Simplex */
	// Factory method - only way to create new simplex
	static std::shared_ptr<Simplex>
	make_Simplex(index_t label,
	             index_t max_vertex,
	             value_t value = DEFAULT_FILT_VALUE,
	             const std::vector<std::shared_ptr<Simplex>>& facets =
	                 std::vector<std::shared_ptr<Simplex>>{});

	// Get a handle to this simplex
	std::shared_ptr<Simplex> get_handle();

	// Return the sorted vertex labels of the simplex
	// Assumes that the simplex has valid faces
	std::vector<index_t> get_vertex_labels() const;

	// Writes the sorted vertex labels of the simplex into a buffer
	// Assumes that the simplex has valid faces
	template <typename OutputIterator> void get_vertex_labels(OutputIterator&& buf) const;

	// Return the indices of the facets of the simplex, [i]th element of the
	// result is the [i]th facet
	std::vector<index_t> get_facet_labels() const;

	// Return a const reference to the facets
	const std::vector<std::shared_ptr<Simplex>>& get_facets() const;

	// Return a const reference to the cofacets
	const std::vector<std::weak_ptr<Simplex>>& get_cofacets() const;

	void set_colour(index_t c);  // not inline since we export this

	unsigned long long int get_colours_as_int();

	inline void set_colours(colours_t c);

	inline void add_colour(index_t c);

	inline void add_colours(colours_t c);

	colours_t get_colours();

	inline void make_colourless();

	/* PRIVATE MEMBERS OF Simplex */
  private:
	const std::vector<std::shared_ptr<Simplex>>
										facets;    // pointers to the [i]th facets of the simplex
	std::vector<std::weak_ptr<Simplex>> cofacets;  // pointers to the cofacets of the simplex
	colours_t colours;  // bitmask representing the colours of its vertices
	// Constructor
	Simplex(index_t                                      label,
	        index_t                                      max_vertex,
	        value_t                                      value,
	        const std::vector<std::shared_ptr<Simplex>>& facets);
	// Delete the default constructor
	Simplex() = delete;
};

// The simplicial complex associated to the standard n-simplex.
FilteredComplex standard_simplex(const index_t n);

}  // namespace chalc

#endif
