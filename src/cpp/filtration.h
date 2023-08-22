/*
    This file is part of Chalc.

    Chalc: Chromatic Alpha Complexes.
    Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha Complexes”. 
    Online preprint available at http://arxiv.org/abs/2212.03128. 
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

#include <vector>
#include <map>
#include <tuple>
#include <memory>
#include <bitset>

namespace chalc {
    using std::vector, std::map, std::shared_ptr, std::tuple, std::bitset;
    typedef double value_t;
    constexpr size_t MAX_NUM_COLOURS = 4;
    typedef bitset<MAX_NUM_COLOURS> colours_t;

    class BinomialCoeffTable {
        vector<vector<size_t>> B;
    public:
        BinomialCoeffTable(size_t n, size_t k);
        size_t operator()(size_t n, size_t k) const;
    };

    class FilteredComplex {

    public:

        /* PUBLIC CLASSES OF FilteredComplex */

        class Simplex;

        /* PUBLIC MEMBERS OF FilteredComplex */

        const size_t N; // number of vertices, labelled from 0 to n-1
        const size_t max_dim; // maximum dimension of any simplex in the complex

        /* PUBLIC METHODS OF FilteredComplex */

        // Constructors
        // 1. Create a vertex set.
        FilteredComplex(const size_t num_vertices, const size_t max_dimension);
        // 2. k-skeleton of another complex
        FilteredComplex(const FilteredComplex& other, const size_t k);

        // Get label of a simplex from the labels of its vertices
        size_t get_label_from_vertex_labels(vector<size_t>& verts) const;

        // Check if the complex has a specific simplex
        bool has_simplex(const size_t dim, const size_t label) const;

        // Check if the complex has a specific simplex
        bool has_simplex(vector<size_t>& verts) const;

        // Add a new simplex to the complex by its vertices
        bool add_simplex(vector<size_t>& verts, const value_t filt_value);

        // Propagate filtration values from start_dim
        void propagate_filt_values(const size_t start_dim, const bool upwards);

        // Returns a handle to the simplices
        const vector<map<size_t, shared_ptr<Simplex>>>& get_simplices() const noexcept;

        // Returns the number of simplices in the complex
        size_t size() const noexcept;

        // Returns the number of simplices of a specific dimension
        size_t size_in_dim(const size_t dim) const;

        // Returns the current dimension
        size_t dimension() const noexcept;

        // Returns the current maximum filtration value
        value_t max_filt_value() const noexcept;

        // Returns a flat vectorised representation of the complex
        vector<tuple<vector<size_t>, size_t, value_t, unsigned long>> flat_representation() const;

        // Returns the k-skeleton of the clique complex on n vertices
        static FilteredComplex clique_complex(const size_t n, const size_t k);

        // bitwise OR accumulates colours upwards from vertices
        void propagate_colours() const;

    private:

        /* PRIVATE MEMBERS OF FilteredComplex */

        const BinomialCoeffTable binomial; // binomial coefficients
        vector<map<size_t, shared_ptr<Simplex>>> simplices; // vector whose kth element is a table of k-simplices, labelled by their lexicographic index
        size_t num_simplices; // total number of simplices
        size_t cur_dim; // current maximum dimension of a maximal simplex
        value_t cur_max_filt_value; // current maximum filtration value

        /* PRIVATE METHODS OF FilteredComplex */

        // Checks that verts is a non-empty subsequence of (0, ..., N-1) and verts.size() <= max_dim + 1
        // does not leave verts unmodified
        void check_vertex_sequence_is_valid(vector<size_t>& verts) const;

        // Check that 0 <= dim <= max_dim
        void check_dimension_is_valid(const size_t dim) const;

        // Get label of a simplex (possibly not in the complex) from the labels of its vertices
        // Assumes that verts is valid
        size_t _get_label_from_vertex_labels(const vector<size_t>& verts) const;

        // Check if the complex has a specific simplex
        // Assumes that dim is valid
        bool _has_simplex(const size_t dim, const size_t label) const;

        // Check if the complex has a specific simplex
        // Assumes that verts is valid
        bool _has_simplex(const vector<size_t>& verts) const;

        shared_ptr<Simplex> _add_simplex(const vector<size_t>& verts, const value_t filt_value);

        // min accumulate filtration values downwards from start_dim
        // Assumes that start_dim is valid
        void propagate_filt_values_up(const size_t start_dim) const;

        // max accumulate filtration values upwards from start_dim
        // Assumes that start_dim is valid
        void propagate_filt_values_down(const size_t start_dim) const;
    };

    class FilteredComplex::Simplex {

        /* PRIVATE MEMBERS OF Simplex */

        const vector<shared_ptr<Simplex>> facets; // pointers to the [i]th facets of the simplex
        // vector<std::weak_ptr<Simplex>> cofacets; // pointers to the [i]th cofacets of the simplex

    public:

        /* PUBLIC MEMBERS OF Simplex */

        value_t value; // filtration value
        const size_t label; // label for the simplex
        const size_t max_vertex; // largest vertex label
        const size_t dim; // number of vertices - 1
        colours_t colours; // bitmask representing the colours of its vertices
        static constexpr value_t DEFAULT_FILT_VALUE = 0.0;

        /* PUBLIC METHODS OF Simplex */

        // Constructor
        Simplex(size_t label, size_t max_vertex, size_t dim = 0, 
        colours_t colours = colours_t(), value_t value = DEFAULT_FILT_VALUE, 
        const vector<shared_ptr<Simplex>>& facets = vector<shared_ptr<Simplex>>{});

        // Return the sorted vertex labels of the simplex
        // Assumes that the simplex has valid faces
        vector<size_t> get_vertex_labels() const;

        // Writes the sorted vertex labels of the simplex into a buffer
        // Assumes that the simplex has valid faces
        template <typename OutputIterator>
        void get_vertex_labels(OutputIterator&& buf) const;

        // Return the indices of the facets of the simplex, [i]th element of the result is the [i]th facet
        vector<size_t> get_facet_labels() const;

        // Return a const reference to the facets
        const vector<shared_ptr<Simplex>>& get_facets() const;

    };

    // The simplicial complex associated to the standard n-simplex.
    FilteredComplex standard_simplex(const size_t n);

}
#endif