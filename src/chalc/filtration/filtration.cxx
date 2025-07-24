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
#include <algorithm>
#include <cassert>
#include <chalc/filtration/filtration.h>
#include <cmath>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
using std::domain_error;
using std::invalid_argument;
using std::max;
using std::min;
using std::numeric_limits;
using std::ranges::adjacent_find;
using std::ranges::fill;
using std::ranges::prev_permutation;
using std::ranges::sort;
using std::ranges::stable_sort;
using std::runtime_error;
using std::to_string;
using std::tuple;
using std::unordered_map;
using std::vector;
}  // namespace

namespace chalc {

// Constructs a binomial coefficient table that will hold all values of
// i_C_j for i = 0, ..., n and j = 0, ..., min(k, floor(i/2)) for i <= n.
detail::BinomialCoeffTable::BinomialCoeffTable(Index n, Index k) :
	B(n + 1) {
	if (n < k || n < 0 || k < 0) {
		throw invalid_argument(
			"Binomial coefficient table cannot be constructed with n < k or n < 0 or k < 0."
		);
	}
	B[0].resize(1, 0);
	Index j_max = 0;
	for (Index i = 1; i <= n; ++i) {
		j_max = min(i / 2, k);
		B[i].resize(j_max + 1, 0);
		B[i][0] = 1;
		for (Index j = 1; j <= j_max; ++j) {
			// Use recurrence relation to fill the entries
			// SAFETY:
			// If i == 1 then j_max == 0 and hence
			// this inner loop is not executed.
			// Therefore we can assume that inside this loop
			// 2 <= i <= n and i - 1 >= i/2.
			// Therefore 1 <= j <= j_max <= min(k, i/2), so that
			// i - 1 - j >= i/2 - j_max >= 0.
			B[i][j] = B[i - 1][min(j - 1, i - j)] + B[i - 1][min(j, i - 1 - j)];
		}
		// Bounds checking: only check at the largest entry i_C_floor(i/2)
		// We check by using the recurrence relation and the datatype max.
		if (B[i - 1][min(j_max - 1, i - j_max)] >
		    numeric_limits<Label>::max() - B[i - 1][min(j_max, i - 1 - j_max)]) {
			throw runtime_error("Simplex index is too large.");
		}
	}
}

auto Filtration::Simplex::make_simplex(
	Label                   label,
	Index                   max_vertex,
	Value                   value,
	const vector<Simplex*>& facets
) -> SimplexHandle {
	auto simplex = SimplexHandle(new Simplex{label, max_vertex, value, facets});
	for (auto&& f: simplex->facets_) {
		f->cofacets_.push_back(simplex.get());
		simplex->add_colours_bitmask(f->colours_);
	}
	return simplex;
}

Filtration::Simplex::Simplex(
	Label                   label,
	Index                   max_vertex,
	Value                   value,
	const vector<Simplex*>& facets
) :
	label_(label),
	max_vertex_(max_vertex),
	// Narrowing cast here is not a problem in 64-bit systems
    // unless we are using a ridiculously large number of vertices,
    // in which case we have other problems to worry about.
	filt_value_(value),
	facets_(facets) {
	assert(facets.size() < std::numeric_limits<index_t>::max);
	dim_ = (facets.size() == 0 ? 0 : static_cast<Index>(facets.size()) - 1);
}

auto Filtration::Simplex::vertex_labels() const -> vector<Index> {
	vector<Index> result(dim_ + 1);
	vertex_labels_(result.begin());
	return result;
}

template <typename OutputIterator>
void Filtration::Simplex::vertex_labels_(OutputIterator&& buf) const {
	if (dim_ > 0) {
		// vertices of a simplex are the vertices of its last face along with
		// its last vertex
		facets_.back()->vertex_labels_(std::forward<OutputIterator>(buf));
		buf++;
	}
	*buf = max_vertex_;
}

auto Filtration::Simplex::facet_labels() const -> vector<Label> {
	vector<Label> result;
	if (dim_ != 0) {
		result.reserve(facets_.size());
		for (auto& f: facets_) {
			result.push_back(f->label_);
		}
	}
	return result;
}

auto Filtration::Simplex::colours() const -> vector<Colour> {
	vector<Colour> result;
	result.reserve(colours_.count());
	for (Colour i = 0; i < MAX_NUM_COLOURS; i++) {
		if (colours_[i]) {
			result.push_back(i);
		}
	}
	return result;
}

Filtration::Filtration(const Index num_vertices, const Index max_dimension) :
	binomial_(
		num_vertices,
		max_dimension + 1
	),  // we want nCk for all 0 <= n <= num_vertices and 0 <= k <=
        // max_num_verts_in_a_simplex = max_dim + 1
	simplices_map_(max_dimension + 1),
	num_vertices_(num_vertices),
	max_dimension_(max_dimension),
	cur_dim_(0) {
	if (max_dimension_ >= num_vertices_) {
		throw invalid_argument("Dimension must be less than number of points.");
	}
	simplices_.reserve(num_vertices);
	for (Index i = 0; i < num_vertices_; i++) {
		simplices_.push_back(Simplex::make_simplex(i, i, 0.0));
		simplices_map_[0][i] = simplices_.back().get();
		// simplices are initialised with colours unset
		// so we set them here
		simplices_map_[0][i]->set_colour(0);
	}
}

// Factory method - get k-skeleton
[[nodiscard]]
auto Filtration::skeleton(const Index k) const -> Filtration {
	if (k < 0) {
		throw invalid_argument("k must be non-negative.");
	}
	Filtration ret(num_vertices_, k);
	for (Index i = min(cur_dim_, k); i >= 1; i--) {
		for (auto&& [key, simplex]: simplices_map_[i]) {
			if (simplex->cofacets().empty()) {
				ret.add_simplex_unchecked(simplex->vertex_labels(), simplex->get_value());
			}
		}
	}
	return ret;
}

auto Filtration::validated_vertex_sequence(const vector<Index>& verts) const -> vector<Index> {
	if (verts.size() == 0) {
		throw invalid_argument("Vertex sequence cannot be empty.");
	}
	if (verts.size() - 1 > max_dimension_) {
		throw invalid_argument("Vertex sequence is too long.");
	}
	vector<Index> verts_sorted(verts);
	sort(verts_sorted);
	if (verts_sorted.front() < 0 || verts_sorted.back() >= num_vertices_ ||
	    adjacent_find(verts_sorted) != verts_sorted.cend()) {
		throw invalid_argument(
			"Vertex sequence must be a subset of {0, ..., num_verts-1} without repetitions."
		);
	};
	return verts_sorted;
}

// Needs verts to be validated.
auto Filtration::lex_label(const vector<Index>& verts) const -> Label {
	auto num_verts = static_cast<Index>(verts.size());
	assert(num_verts != 0);  // DEBUG
	if (num_verts == 1) {
		return verts[0];     // if we have a vertex, return its label
	}
	Label label = 0;
	for (Index i = 0, v = 0; i < num_verts; v = verts[i++] + 1) {
		for (Index j = v; j < verts[i]; j++) {
			label += binomial_(num_vertices_ - (j + 1), num_verts - (i + 1));
		}
	}
	return label;
}

auto Filtration::get_label_from_vertex_labels(const vector<Index>& verts) const -> Label {
	auto verts_validated = validated_vertex_sequence(verts);
	return lex_label(verts_validated);
}

auto Filtration::has_simplex_unchecked(const Index dim, const Label label) const -> bool {
	if (simplices_map_[dim].find(label) == simplices_map_[dim].end()) {
		return false;
	} else {
		return true;
	}
}

auto Filtration::has_simplex(const Index dim, const Label label) const -> bool {
	if (dim > max_dimension_ || dim < 0) {
		throw invalid_argument("Invalid dimension.");
	}
	if (label >= binomial_(num_vertices_, dim)) {
		// The check above is valid because we have
		// binomial coefficients for all 0 <= n <=
		// n_vertices and 0 <= k <= max_dim + 1.
		throw invalid_argument("Invalid label.");
	}
	return has_simplex_unchecked(dim, label);
}

auto Filtration::has_simplex_unchecked(const vector<Index>& verts) const noexcept -> bool {
	// Narrowing cast here is not a problem since this is
	// called only from has_simplex, which performs validation.
	auto num_verts = static_cast<Index>(verts.size());
	assert(num_verts != 0);
	auto dim   = num_verts - 1;  // we assume that verts is valid
	auto label = lex_label(verts);
	return (has_simplex_unchecked(dim, label));
}

auto Filtration::has_simplex(vector<Index>& verts) const -> bool {
	auto verts_validated = validated_vertex_sequence(verts);
	return (has_simplex_unchecked(verts_validated));
}

auto Filtration::add_simplex_unchecked(const vector<Index>& verts, const Value filt_value)
	-> Filtration::Simplex* {  // Narrowing cast here is not a problem since this is
	// called only from add_simplex, which performs validation.
	auto num_verts = static_cast<Index>(verts.size());
	assert(num_verts != 0);
	Simplex* new_simplex    = nullptr;
	auto     dim            = num_verts - 1;  // safe because we assume we have validated verts
	auto     label          = lex_label(verts);
	auto     search_simplex = simplices_map_[dim].find(label);
	if (search_simplex != simplices_map_[dim].end()) {
		// the simplex already exists
		new_simplex              = search_simplex->second;
		new_simplex->filt_value_ = min(max(filt_value, 0.0), new_simplex->filt_value_);
	} else {
		// simplex does not exist so we need to add it
		// first we recursively add faces of the simplex
		vector<Simplex*> facets(dim + 1);
		vector<Index>    facet_verts;
		facet_verts.reserve(dim);
		for (auto i = 0; i <= dim; i++) {
			facet_verts.insert(facet_verts.cend(), verts.cbegin(), verts.cbegin() + i);
			facet_verts.insert(
				facet_verts.cend(),
				verts.cbegin() + i + 1,
				verts.cend()
			);  // copy all except the i-th vertex
			facets[i] = add_simplex_unchecked(facet_verts, filt_value);
			facet_verts.clear();
		}
		auto max_vertex = verts.back();
		simplices_.push_back(Simplex::make_simplex(label, max_vertex, filt_value, facets));
		new_simplex                = simplices_.back().get();
		simplices_map_[dim][label] = new_simplex;
		cur_dim_                   = max(cur_dim_,
                       static_cast<Index>(verts.size() - 1));  // need verts to be valid
	}
	return new_simplex;
}

auto Filtration::add_simplex(
	const vector<Index>& verts,
	Value                filt_value = Filtration::Simplex::DEFAULT_FILT_VALUE
) -> bool {
	auto verts_validated = validated_vertex_sequence(verts);
	if (has_simplex_unchecked(verts_validated)) {
		return false;
	} else {
		add_simplex_unchecked(verts_validated, filt_value);
		return true;
	}
}

void Filtration::propagate_filt_values_up(const Index start_dim) noexcept {
	auto p = max(start_dim + static_cast<Index>(1), static_cast<Index>(1));
	while (p <= cur_dim_) {
		Value max_facet_filt_value = NAN;
		// iterate over the p-simplices
		for (auto& [label, simplex]: simplices_map_[p]) {
			max_facet_filt_value = simplex->facets()[0]->filt_value_;
			// iterate over faces of simplex and get the maximum filtration
			// value
			for (const auto& facet: simplex->facets()) {
				max_facet_filt_value = max(max_facet_filt_value, facet->filt_value_);
			}
			simplex->filt_value_ = max_facet_filt_value;
		}
		p++;
	}
}

void Filtration::propagate_filt_values_down(const Index start_dim) {
	if (start_dim <= 0) {
		return;
	}
	for (Index p = min(start_dim, cur_dim_) - 1; p >= 1; p--) {
		// iterate over the p-simplices
		for (const auto& [label, simplex]: simplices_map_[p]) {
			// iterate over cofacets of simplex and
			// modify the filtration value of simplex if needed
			auto& cofacets = simplex->cofacets();
			if (cofacets.size() == 0) {
				continue;
			}
			simplex->filt_value_ = cofacets[0]->filt_value_;
			for (auto&& cofacet: cofacets) {
				simplex->filt_value_ = min(simplex->filt_value_, cofacet->filt_value_);
			}
		}
	}
}

void Filtration::propagate_colours() noexcept {
	for (Index d = 1; d <= cur_dim_; d++) {
		for (auto& [idx, s]: simplices_map_[d]) {
			s->make_colourless();
			for (auto& f: s->facets()) {
				s->add_colours_bitmask(f->colours_bitmask());
			}
		}
	}
}

void Filtration::propagate_filt_values(const Index start_dim, const bool up) {
	if (up) {
		propagate_filt_values_up(start_dim);
	} else {
		propagate_filt_values_down(start_dim);
	}
}

auto Filtration::simplices_begin() -> Filtration::SimplicesIterator {
	return SimplicesIterator(simplices_map_.begin(), simplices_map_.end());
}

auto Filtration::simplices_end() -> Filtration::SimplicesIterator {
	return SimplicesIterator(simplices_map_.end(), simplices_map_.end());
}

auto Filtration::boundary_matrix(Index max_dimension) const
	-> vector<tuple<vector<Index>, Label, Value, vector<Colour>>> {
	if (max_dimension < -1) {
		throw domain_error("max_dimension must be >= -1, received " + to_string(max_dimension));
	}

	// Sort the simplices by their filtration value, dimension, and label.
	size_t length_return = 0;
	switch (max_dimension) {
		case -1 :
			max_dimension = cur_dim_;
			length_return = simplices_.size();
			break;
		case 0 :
			length_return = num_vertices_;
			break;
		default :
			max_dimension = min(max_dimension, cur_dim_);
			for (Index p = 0; p <= max_dimension; p++) {
				length_return += simplices_map_[p].size();
			}
	}
	vector<const Simplex*> sort_by_val;
	sort_by_val.reserve(length_return);
	for (Index p = 0; p <= max_dimension; p++) {
		for (auto&& s: simplices_map_[p]) {
			sort_by_val.push_back(s.second);
		}
	}
	stable_sort(sort_by_val, [](const Simplex* s1, const Simplex* s2) {
		return (
			s1->filt_value_ < s2->filt_value_ ||
			(s1->filt_value_ == s2->filt_value_ && s1->dim_ < s2->dim_) ||
			(s1->filt_value_ == s2->filt_value_ && s1->dim_ == s2->dim_ && s1->label_ < s2->label_)
		);
	});

	// Create mapping from labels to sorted index.
	vector<tuple<vector<Index>, Label, Value, vector<Colour>>> result(length_return);
	vector<unordered_map<Label, Index>>                        indices(max_dimension + 1);
	for (auto&& [i, s] = tuple{0L, sort_by_val.begin()}; s != sort_by_val.end(); s++, i++) {
		auto simplex = *s;
		auto dim     = simplex->dim_;
		auto label   = simplex->label_;
		// Add the label to the indices map for this dimension
		indices[dim][label] = i;
		// Add this simplex to the result
		auto          value        = simplex->filt_value_;
		auto          colours      = simplex->colours();
		auto          facet_labels = simplex->facet_labels();
		vector<Index> facet_idx(facet_labels.size());
		for (auto&& [idx_it, label_it] = tuple{facet_idx.begin(), facet_labels.cbegin()};
		     idx_it != facet_idx.cend();
		     idx_it++, label_it++) {
			*idx_it = indices[dim - 1][*label_it];
		}
		sort(facet_idx.begin(), facet_idx.end());
		result[i] = tuple{facet_idx, label, value, colours};
	}
	return result;
}

auto Filtration::complete_complex(const Index n, const Index k) -> Filtration {
	if (k >= n || n <= 0 || k < 0) {
		throw invalid_argument("k must satisfy 0 <= k < n");
	}
	Filtration   result(n, k);
	vector<bool> verts_mask(n);
	auto         end = verts_mask.begin();
	for (Index i = 0; i < k + 1; i++) {
		end++;
	}
	fill(verts_mask.begin(), end, true);
	vector<Index> verts(k + 1);
	bool          perm_found = true;
	while (perm_found) {
		for (auto&& [i, verts_it, mask_it] = tuple{0L, verts.begin(), verts_mask.begin()};
		     mask_it != verts_mask.end();
		     mask_it++, i++) {
			if (*mask_it) {
				*verts_it = i;
				verts_it++;
			}
		}
		result.add_simplex_unchecked(verts, 0.0);
		perm_found = prev_permutation(verts_mask).found;
	}
	return result;
}

auto Filtration::is_filtration() const noexcept -> bool {
	for (auto&& simplices_i: simplices_map_) {
		for (auto& s: simplices_i) {
			for (auto& f: s.second->facets()) {
				if (f->filt_value_ > s.second->filt_value_) {
					return false;
				}
			}
		}
	}
	return true;
}

auto standard_simplex(const Index n) -> Filtration {
	return chalc::Filtration::complete_complex(n + 1, n);
}
}  // namespace chalc
