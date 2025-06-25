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
#include <utility>

namespace {
using std::bad_weak_ptr;
using std::invalid_argument;
using std::make_shared;
using std::map;
using std::max;
using std::min;
using std::numeric_limits;
using std::ranges::adjacent_find;
using std::ranges::any_of;
using std::ranges::fill;
using std::ranges::prev_permutation;
using std::ranges::sort;
using std::ranges::stable_sort;
using std::runtime_error;
using std::shared_ptr;
using std::tuple;
using std::vector;
}  // namespace

namespace chalc {

class BinomialCoeffTable {
	vector<vector<label_t>> B;

  public:
	// Constructs a binomial coefficient table that will hold all values of
	// i_C_j for i = 0, ..., n and j = 0, ..., min(k, floor(i/2)) for i <= n.
	BinomialCoeffTable(index_t n, index_t k) :
		B(n + 1) {
		if (n < k || n < 0 || k < 0) {
			throw invalid_argument(
				"Binomial coefficient table cannot be constructed with n < k or n < 0 or k < 0."
			);
		}
		B[0].resize(1, 0);
		index_t j_max = 0;
		for (index_t i = 1; i <= n; ++i) {
			j_max = min(i / 2, k);
			B[i].resize(j_max + 1, 0);
			B[i][0] = 1;
			for (index_t j = 1; j <= j_max; ++j) {
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
			    numeric_limits<label_t>::max() - B[i - 1][min(j_max, i - 1 - j_max)]) {
				throw runtime_error("Simplex index is too large.");
			}
		}
	}

	auto operator()(index_t n, index_t k) const -> label_t {
		return B.at(n).at(min(k, n - k));
	}
};

auto Filtration::Simplex::_make_simplex(
	label_t                            label,
	index_t                            max_vertex,
	value_t                            value,
	const vector<shared_ptr<Simplex>>& facets
) -> shared_ptr<Filtration::Simplex> {
	auto self = shared_ptr<Simplex>(new Simplex(label, max_vertex, value, facets));
	for (auto& f: self->facets) {
		f->cofacets.push_back(self->get_handle());
		self->_add_colours(f->colours);
	}
	return self;
}

Filtration::Simplex::Simplex(
	label_t                            label,
	index_t                            max_vertex,
	value_t                            value,
	const vector<shared_ptr<Simplex>>& facets
) :
	m_label(label),
	m_max_vertex(max_vertex),
	// Narrowing cast here is not a problem in 64-bit systems
    // unless we are using a ridiculously large number of vertices,
    // in which case we have other problems to worry about.
	m_dim(facets.size() == 0 ? 0 : facets.size() - 1),  // NOLINT
	m_filt_value(value),
	facets(facets) {}

auto Filtration::Simplex::get_handle() -> shared_ptr<Filtration::Simplex> {
	return shared_from_this();
}

auto Filtration::Simplex::get_vertex_labels() const -> vector<index_t> {
	vector<index_t> result(m_dim + 1);
	_get_vertex_labels(result.begin());
	return result;
}

template <typename OutputIterator>
void Filtration::Simplex::_get_vertex_labels(OutputIterator&& buf) const {
	if (m_dim > 0) {
		// vertices of a simplex are the vertices of its last face along with
		// its last vertex
		facets.back()->_get_vertex_labels(std::forward<OutputIterator>(buf));
		buf++;
	}
	*buf = m_max_vertex;
}

auto Filtration::Simplex::get_facet_labels() const -> vector<label_t> {
	vector<label_t> result;
	if (m_dim != 0) {
		result.reserve(facets.size());
		for (auto& f: facets) {
			result.push_back(f->m_label);
		}
	}
	return result;
}

auto Filtration::Simplex::get_colours_as_vec() const -> vector<colour_t> {
	vector<colour_t> result;
	result.reserve(colours.count());
	for (colour_t i = 0; i < MAX_NUM_COLOURS; i++) {
		if (colours[i]) {
			result.push_back(i);
		}
	}
	return result;
}

Filtration::Filtration(const index_t num_vertices, const index_t max_dimension) :
	binomial(
		make_shared<BinomialCoeffTable>(
			num_vertices,
			max_dimension + 1
		)
	),  // we want nCk for all 0 <= n <= num_vertices and 0 <= k <=
        // max_num_verts_in_a_simplex = max_dim + 1
	simplices(max_dimension + 1),
	n_vertices(num_vertices),
	max_dim(max_dimension),
	num_simplices(num_vertices),
	cur_max_filt_value(0),
	cur_dim(0) {
	if (max_dim >= n_vertices) {
		throw invalid_argument("Dimension must be less than number of points.");
	}
	for (index_t i = 0; i < n_vertices; i++) {
		simplices[0][i] = Simplex::_make_simplex(i, i, 0.0);
		// simplices are initialised with colours unset
		// so we set them here
		simplices[0][i]->set_colour(0);
	}
}

void Filtration::validate_vertex_sequence(vector<index_t>& verts) const {
	if (verts.size() == 0) {
		throw invalid_argument("Vertex sequence cannot be empty.");
	}
	if (verts.size() - 1 > max_dim) {
		throw invalid_argument("Vertex sequence is too long.");
	}
	sort(verts);
	if (!(verts.back() < n_vertices && adjacent_find(verts) == verts.cend())) {
		throw invalid_argument("Vertex sequence cannot have repetitions.");
	};
	if (any_of(verts, [this](index_t v) {
			return v < 0 || v >= n_vertices;
		})) {
		throw invalid_argument("Vertex sequence must be a subsequence of (0, ..., n_vertices-1).");
	}
}

auto Filtration::_get_label_from_vertex_labels(const vector<index_t>& verts) const -> label_t {
	auto num_verts = static_cast<index_t>(verts.size());
	assert(num_verts != 0);  // DEBUG
	if (num_verts == 1) {
		return verts[0];     // if we have a vertex, return its label
	}
	label_t label = 0;
	for (index_t i = 0, v = 0; i < num_verts; v = verts[i++] + 1) {
		for (index_t j = v; j < verts[i]; j++) {
			label += (*binomial)(n_vertices - (j + 1), num_verts - (i + 1));
		}
	}
	return label;
}

auto Filtration::get_label_from_vertex_labels(vector<index_t>& verts) const -> label_t {
	validate_vertex_sequence(verts);
	return _get_label_from_vertex_labels(verts);
}

auto Filtration::_has_simplex(const index_t dim, const label_t label) const -> bool {
	if (simplices[dim].find(label) == simplices[dim].end()) {
		return false;
	} else {
		return true;
	}
}

auto Filtration::has_simplex(const index_t dim, const label_t label) const -> bool {
	if (dim > max_dim || dim < 0) {
		throw invalid_argument("Invalid dimension.");
	}
	if (label >= (*binomial)(n_vertices, dim)) {
		// The check above is valid because we have
		// binomial coefficients for all 0 <= n <=
		// n_vertices and 0 <= k <= max_dim + 1.
		throw invalid_argument("Invalid label.");
	}
	return _has_simplex(dim, label);
}

auto Filtration::_has_simplex(const vector<index_t>& verts) const noexcept -> bool {
	// Narrowing cast here is not a problem since this is
	// called only from has_simplex, which performs validation.
	index_t num_verts = verts.size();  // NOLINT
	assert(num_verts != 0);
	auto dim   = num_verts - 1;        // we assume that verts is valid
	auto label = _get_label_from_vertex_labels(verts);
	return (_has_simplex(dim, label));
}

auto Filtration::has_simplex(vector<index_t>& verts) const -> bool {
	validate_vertex_sequence(verts);
	return (_has_simplex(verts));
}

auto Filtration::_add_simplex(const vector<index_t>& verts, const value_t filt_value)
	-> shared_ptr<Filtration::Simplex> {
	// Narrowing cast here is not a problem since this is
	// called only from add_simplex, which performs validation.
	index_t num_verts = verts.size();  // NOLINT
	assert(num_verts != 0);
	shared_ptr<Simplex> new_simplex;
	auto                dim            = num_verts - 1;
	auto                label          = _get_label_from_vertex_labels(verts);
	auto                search_simplex = simplices[dim].find(label);
	if (search_simplex != simplices[dim].end()) {
		// the simplex already exists
		new_simplex               = search_simplex->second;
		new_simplex->m_filt_value = min(max(filt_value, 0.0), new_simplex->m_filt_value);
	} else {
		// simplex does not exist so we need to add it
		// first we recursively add faces of the simplex
		vector<shared_ptr<Simplex>> facets(dim + 1);
		vector<index_t>             facet_verts;
		facet_verts.reserve(dim);
		for (auto i = 0; i <= dim; i++) {
			facet_verts.insert(facet_verts.cend(), verts.cbegin(), verts.cbegin() + i);
			facet_verts.insert(
				facet_verts.cend(),
				verts.cbegin() + i + 1,
				verts.cend()
			);  // copy all except the i-th vertex
			facets[i] = _add_simplex(facet_verts, filt_value);
			facet_verts.clear();
		}
		auto max_vertex       = verts.back();
		new_simplex           = Simplex::_make_simplex(label, max_vertex, filt_value, facets);
		simplices[dim][label] = new_simplex;
		num_simplices++;
	}
	return new_simplex;
}

auto Filtration::add_simplex(
	vector<index_t>& verts,
	value_t          filt_value = Filtration::Simplex::DEFAULT_FILT_VALUE
) -> bool {
	validate_vertex_sequence(verts);
	if (_has_simplex(verts)) {
		return false;
	} else {
		_add_simplex(verts, filt_value);
		cur_dim = max(cur_dim, static_cast<index_t>(verts.size() - 1));  // need verts to be valid
		cur_max_filt_value = max(cur_max_filt_value, filt_value);
		return true;
	}
}

void Filtration::propagate_filt_values_up(const index_t start_dim) noexcept {
	auto p = max(start_dim + static_cast<index_t>(1), static_cast<index_t>(1));
	while (p <= cur_dim) {
		value_t max_facet_filt_value = NAN;
		// iterate over the p-simplices
		for (auto& [label, simplex]: simplices[p]) {
			max_facet_filt_value = simplex->get_facets()[0]->m_filt_value;
			// iterate over faces of simplex and get the maximum filtration
			// value
			for (const auto& facet: simplex->get_facets()) {
				max_facet_filt_value = max(max_facet_filt_value, facet->m_filt_value);
			}
			simplex->m_filt_value = max_facet_filt_value;
		}
		p++;
	}
}

void Filtration::propagate_filt_values_down(const index_t start_dim) {
	if (start_dim <= 0) {
		return;
	}
	for (index_t p = min(start_dim, cur_dim) - 1; p >= 1; p--) {
		// iterate over the p-simplices
		for (const auto& [label, simplex]: simplices[p]) {
			// iterate over cofacets of simplex and
			// modify the filtration value of simplex if needed
			auto& cofacets = simplex->get_cofacets();
			if (cofacets.size() == 0) {
				continue;
			}
			try {
				simplex->m_filt_value = cofacets[0].lock()->m_filt_value;
				for (auto&& cofacet: cofacets) {
					simplex->m_filt_value =
						min(simplex->m_filt_value, cofacet.lock()->m_filt_value);
				}
			} catch (const bad_weak_ptr& e) {
				throw runtime_error("Tried to dereference expired cofacet handle.");
			}
		}
	}
}

void Filtration::propagate_colours() noexcept {
	for (index_t d = 1; d <= cur_dim; d++) {
		for (auto& [idx, s]: simplices[d]) {
			s->make_colourless();
			for (auto& f: s->get_facets()) {
				s->_add_colours(f->_get_colours());
			}
		}
	}
}

void Filtration::propagate_filt_values(const index_t start_dim, const bool up) {
	if (up) {
		propagate_filt_values_up(start_dim);
	} else {
		propagate_filt_values_down(start_dim);
	}
}

auto Filtration::serialised() const
	-> vector<tuple<vector<index_t>, label_t, value_t, vector<colour_t>>> {
	// Sort the simplices by their filtration value, dimension, and label.
	vector<shared_ptr<Simplex>> sort_by_val;
	sort_by_val.reserve(num_simplices);
	for (auto&& simplices_d: simplices) {
		for (auto&& s: simplices_d) {
			sort_by_val.push_back(s.second);
		}
	}
	stable_sort(sort_by_val, [](const shared_ptr<Simplex>& s1, const shared_ptr<Simplex>& s2) {
		return (
			s1->m_filt_value < s2->m_filt_value ||
			(s1->m_filt_value == s2->m_filt_value && s1->m_dim < s2->m_dim) ||
			(s1->m_filt_value == s2->m_filt_value && s1->m_dim == s2->m_dim &&
		     s1->m_label < s2->m_label)
		);
	});
	// Create mapping from labels to sorted index.
	vector<tuple<vector<index_t>, label_t, value_t, vector<colour_t>>> result(num_simplices);
	vector<map<label_t, index_t>>                                      indices(cur_dim + 1);
	for (auto&& [i, s] = tuple{0L, sort_by_val.begin()}; s != sort_by_val.end(); s++, i++) {
		auto simplex = *s;
		auto dim     = simplex->m_dim;
		// Add the label to the indices map for this dimension
		indices[simplex->m_dim][simplex->m_label] = i;
		// Add this simplex to the result
		auto            value        = simplex->m_filt_value;
		auto            label        = simplex->m_label;
		auto            dimension    = simplex->m_dim;
		auto            colours      = simplex->get_colours_as_vec();
		auto            facet_labels = simplex->get_facet_labels();
		vector<index_t> facet_idx(facet_labels.size());
		for (auto&& [idx_it, label_it] = tuple{facet_idx.begin(), facet_labels.cbegin()};
		     idx_it != facet_idx.cend();
		     idx_it++, label_it++) {
			*idx_it = indices[dimension - 1][*label_it];
		}
		sort(facet_idx.begin(), facet_idx.end());
		result[i] = tuple{facet_idx, label, value, colours};
	}
	return result;
}

auto Filtration::complete_complex(const index_t n, const index_t k) -> Filtration {
	if (k >= n || n <= 0 || k < 0) {
		throw invalid_argument("k must satisfy 0 <= k < n");
	}
	Filtration   result(n, k);
	vector<bool> verts_mask(n);
	auto         end = verts_mask.begin();
	for (index_t i = 0; i < k + 1; i++) {
		end++;
	}
	fill(verts_mask.begin(), end, true);
	vector<index_t> verts(k + 1);
	bool            perm_found = true;
	while (perm_found) {
		for (auto&& [i, verts_it, mask_it] = tuple{0L, verts.begin(), verts_mask.begin()};
		     mask_it != verts_mask.end();
		     mask_it++, i++) {
			if (*mask_it) {
				*verts_it = i;
				verts_it++;
			}
		}
		result._add_simplex(verts, 0.0);
		perm_found = prev_permutation(verts_mask).found;
	}
	result.cur_dim = k;
	return result;
}

auto Filtration::is_filtration() const noexcept -> bool {
	for (auto&& simplices_i: simplices) {
		for (auto& s: simplices_i) {
			for (auto& f: s.second->get_facets()) {
				if (f->m_filt_value > s.second->m_filt_value) {
					return false;
				}
			}
		}
	}
	return true;
}

auto standard_simplex(const index_t n) -> Filtration {
	return chalc::Filtration::complete_complex(n + 1, n);
}
}  // namespace chalc
