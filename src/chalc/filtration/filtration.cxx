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
#include <limits>
#include <ranges>
#include <stdexcept>

namespace {
using namespace chalc::stl;
using std::sort, std::stable_sort, std::min, std::max, std::fill, std::adjacent_find,
	std::prev_permutation, std::invalid_argument, std::runtime_error, std::bad_weak_ptr;
}  // namespace

namespace chalc {
class BinomialCoeffTable {
	vector<vector<index_t>> B;

  public:
	// Constructs a binomial coefficient table that will hold all values of
	// i_C_j for i = 0, ..., n and j = 0, ..., min(k, floor(i/2)) for i <= n.
	BinomialCoeffTable(index_t n, index_t k) :
		B(n + 1) {
		assert(n >= k);
		B[0].resize(1, 0);
		index_t j_max;
		for (index_t i = 1; i <= n; ++i) {
			j_max = min(i >> 1, k);
			B[i].resize(j_max + 1, 0);
			B[i][0] = 1;
			for (index_t j = 1; j <= j_max; ++j) {
				// Use recurrence relation to fill the entries
				B[i][j] = B[i - 1][min(j - 1, i - j)] + B[i - 1][min(j, i - 1 - j)];
			}
			// Bounds checking: only check at the largest entry i_C_floor(i/2)
			// We check by using the recurrence relation and the datatype max.
			if (B[i - 1][min(j_max - 1, i - j_max)] >
			    std::numeric_limits<index_t>::max() - B[i - 1][min(j_max, i - 1 - j_max)]) {
				throw runtime_error("Simplex index is too large.");
			}
		}
	}

	index_t operator()(index_t n, index_t k) const {
		return B.at(n).at(min(k, n - k));
	}
};

shared_ptr<FilteredComplex::Simplex> FilteredComplex::Simplex::_make_simplex(
	index_t                            label,
	index_t                            max_vertex,
	value_t                            value,
	const vector<shared_ptr<Simplex>>& facets
) {
	auto self = shared_ptr<Simplex>(new Simplex(label, max_vertex, value, facets));
	for (auto& f: self->facets) {
		f->cofacets.push_back(self->get_handle());
		self->_add_colours(f->colours);
	}
	return self;
}

FilteredComplex::Simplex::Simplex(
	index_t                            label,
	index_t                            max_vertex,
	value_t                            value,
	const vector<shared_ptr<Simplex>>& facets
) :
	label(label),
	max_vertex(max_vertex),
	dim(facets.size() == 0 ? 0 : facets.size() - 1),
	value(value),
	facets(facets) {}

shared_ptr<FilteredComplex::Simplex> FilteredComplex::Simplex::get_handle() {
	return shared_from_this();
}

vector<index_t> FilteredComplex::Simplex::get_vertex_labels() const {
	vector<index_t> result(dim + 1);
	_get_vertex_labels(result.begin());
	return result;
}

template <typename OutputIterator>
void FilteredComplex::Simplex::_get_vertex_labels(OutputIterator&& buf) const {
	if (dim > 0) {
		// vertices of a simplex are the vertices of its last face along with
		// its last vertex
		facets.back()->_get_vertex_labels(buf);
		buf++;
	}
	*buf = max_vertex;
}

vector<index_t> FilteredComplex::Simplex::get_facet_labels() const {
	vector<index_t> result;
	if (dim != 0) {
		result.reserve(facets.size());
		for (auto& f: facets) {
			result.push_back(f->label);
		}
	}
	return result;
}

const vector<shared_ptr<FilteredComplex::Simplex>>& FilteredComplex::Simplex::get_facets() const {
	return facets;
}

const vector<weak_ptr<FilteredComplex::Simplex>>& FilteredComplex::Simplex::get_cofacets() const {
	return cofacets;
}

inline void FilteredComplex::Simplex::add_colour(index_t c) {
	colours.set(c);
}

inline void FilteredComplex::Simplex::_add_colours(colours_t c) {
	colours |= c;
}

inline void FilteredComplex::Simplex::_set_colours(colours_t c) {
	colours.reset();
	_add_colours(c);
}

inline void FilteredComplex::Simplex::make_colourless() {
	colours.reset();
}

void FilteredComplex::Simplex::set_colour(index_t c) {
	colours.reset().set(c);
}

colours_t FilteredComplex::Simplex::_get_colours() {
	return colours;
}

vector<index_t> FilteredComplex::Simplex::get_colours_as_vec() {
	vector<index_t> result;
	result.reserve(colours.count());
	for (size_t i = 0; i < colours.size(); i++) {
		if (colours[i]) {
			result.push_back(static_cast<index_t>(i));
		}
	}
	return result;
}

FilteredComplex::FilteredComplex(const index_t num_vertices, const index_t max_dimension) :
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

void FilteredComplex::validate_vertex_sequence(vector<index_t>& verts) const {
	if (verts.size() == 0) {
		throw invalid_argument("Vertex sequence cannot be empty.");
	}
	if (static_cast<index_t>(verts.size() - 1) > max_dim) {
		throw invalid_argument("Vertex sequence is too long.");
	}
	sort(verts.begin(), verts.end());
	if (!(verts.back() < n_vertices &&
	      adjacent_find(verts.cbegin(), verts.cend()) == verts.cend())) {
		throw invalid_argument("Vertex sequence cannot have repetitions.");
	};
}

index_t FilteredComplex::_get_label_from_vertex_labels(const vector<index_t>& verts) const {
	index_t num_verts = static_cast<index_t>(verts.size());
	assert(num_verts != 0);
	if (num_verts == 1) {
		return verts[0];  // if we have a vertex, return its label
	}
	index_t label = 0;
	for (index_t i = 0, v = 0; i < static_cast<index_t>(num_verts); v = verts[i++] + 1) {
		for (index_t j = v; j < verts[i]; j++) {
			label += (*binomial)(n_vertices - (j + 1), num_verts - (i + 1));
		}
	}
	return label;
}

index_t FilteredComplex::get_label_from_vertex_labels(vector<index_t>& verts) const {
	validate_vertex_sequence(verts);
	return _get_label_from_vertex_labels(verts);
}

bool FilteredComplex::_has_simplex(const index_t dim, const index_t label) const {
	if (simplices[dim].find(label) == simplices[dim].end()) {
		return false;
	} else {
		return true;
	}
}

bool FilteredComplex::has_simplex(const index_t dim, const index_t label) const {
	if (dim > max_dim) {
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

bool FilteredComplex::_has_simplex(const vector<index_t>& verts) const {
	index_t num_verts = static_cast<index_t>(verts.size());
	assert(num_verts != 0);
	auto dim   = num_verts - 1;  // we assume that verts is valid
	auto label = _get_label_from_vertex_labels(verts);
	return (_has_simplex(dim, label));
}

bool FilteredComplex::has_simplex(vector<index_t>& verts) const {
	validate_vertex_sequence(verts);
	return (_has_simplex(verts));
}

shared_ptr<FilteredComplex::Simplex>
FilteredComplex::_add_simplex(const vector<index_t>& verts, const value_t filt_value) {
	index_t num_verts = static_cast<index_t>(verts.size());
	assert(num_verts != 0);
	shared_ptr<Simplex> new_simplex;
	auto                dim            = num_verts - 1;
	auto                label          = _get_label_from_vertex_labels(verts);
	auto                search_simplex = simplices[dim].find(label);
	if (search_simplex != simplices[dim].end()) {
		// the simplex already exists
		new_simplex        = search_simplex->second;
		new_simplex->value = min(max(filt_value, 0.0), new_simplex->value);
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
		auto max_vertex = verts.back();
		new_simplex     = Simplex::_make_simplex(label, max_vertex, filt_value, std::move(facets));
		simplices[dim][label] = new_simplex;
		num_simplices++;
	}
	return new_simplex;
}

bool FilteredComplex::add_simplex(
	vector<index_t>& verts,
	value_t          filt_value = FilteredComplex::Simplex::DEFAULT_FILT_VALUE
) {
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

void FilteredComplex::propagate_filt_values_up(const index_t start_dim) {
	auto p = start_dim + 1;
	while (p <= cur_dim) {
		value_t max_facet_filt_value;
		// iterate over the p-simplices
		for (auto& [label, simplex]: simplices[p]) {
			max_facet_filt_value = simplex->get_facets()[0]->value;
			// iterate over faces of simplex and get the maximum filtration
			// value
			for (const auto& facet: simplex->get_facets()) {
				max_facet_filt_value = max(max_facet_filt_value, facet->value);
			}
			simplex->value = max_facet_filt_value;
		}
		p++;
	}
}

void FilteredComplex::propagate_filt_values_down(const index_t start_dim) {
	if (start_dim == 0) {
		return;
	}
	int p = static_cast<int>(start_dim) - 1;
	for (index_t p = start_dim - 1; p >= 1; p--) {
		// iterate over the p-simplices
		for (const auto& [label, simplex]: simplices[p]) {
			// iterate over cofacets of simplex and
			// modify the filtration value of simplex if needed
			auto& cofacets = simplex->get_cofacets();
			if (cofacets.size() == 0) {
				continue;
			}
			try {
				simplex->value = cofacets[0].lock()->value;
				for (auto& cofacet: cofacets) {
					simplex->value = min(simplex->value, cofacet.lock()->value);
				}
			} catch (const bad_weak_ptr& e) {
				throw runtime_error("Tried to dereference expired cofacet handle.");
			}
		}
	}
}

void FilteredComplex::propagate_colours() {
	for (index_t d = 1; d <= cur_dim; d++) {
		for (auto& [idx, s]: simplices[d]) {
			s->make_colourless();
			for (auto& f: s->get_facets()) {
				s->_add_colours(f->_get_colours());
			}
		}
	}
}

void FilteredComplex::propagate_filt_values(const index_t start_dim, const bool up) {
	if (start_dim > cur_dim) {
		throw invalid_argument("Invalid starting dimension.");
	}
	if (up) {
		propagate_filt_values_up(start_dim);
	} else {
		propagate_filt_values_down(start_dim);
	}
}

const vector<map<index_t, shared_ptr<FilteredComplex::Simplex>>>&
FilteredComplex::get_simplices() const noexcept {
	return simplices;
}

index_t FilteredComplex::size() const noexcept {
	return num_simplices;
}

index_t FilteredComplex::dimension() const noexcept {
	return cur_dim;
}

index_t FilteredComplex::max_dimension() const noexcept {
	return max_dim;
}

index_t FilteredComplex::num_vertices() const noexcept {
	return n_vertices;
}

value_t FilteredComplex::max_filt_value() const noexcept {
	return cur_max_filt_value;
}

vector<tuple<vector<index_t>, index_t, value_t, vector<index_t>>>
FilteredComplex::serialised() const {
	// Sort the simplices by their filtration value, dimension, and label.
	vector<shared_ptr<Simplex>> sort_by_val;
	sort_by_val.reserve(num_simplices);
	for (index_t d = 0; d <= cur_dim; d++) {
		for (auto&& s: simplices[d]) {
			sort_by_val.push_back(s.second);
		}
	}
	stable_sort(
		sort_by_val.begin(),
		sort_by_val.end(),
		[](const shared_ptr<Simplex>& s1, const shared_ptr<Simplex>& s2) {
			return (
				s1->value < s2->value || (s1->value == s2->value && s1->dim < s2->dim) ||
				(s1->value == s2->value && s1->dim == s2->dim && s1->label < s2->label)
			);
		}
	);
	// Create mapping from labels to sorted index.
	vector<tuple<vector<index_t>, index_t, value_t, vector<index_t>>> result(num_simplices);
	vector<map<index_t, index_t>>                                     indices(cur_dim + 1);
	for (auto&& [i, s]: std::views::enumerate(sort_by_val)) {
		auto dim = s->dim;
		// Add the label to the indices map for this dimension
		indices[s->dim][s->label] = i;
		// Add this simplex to the result
		auto value     = s->value;
		auto label     = s->label;
		auto dimension = s->dim;
		auto colours   = s->get_colours_as_vec();
		auto faces     = s->get_facet_labels();
		for (auto&& f: faces) {
			f = indices[dimension - 1][f];
		}
		sort(faces.begin(), faces.end());
		result[i] = tuple{faces, label, value, colours};
	}
	return result;
}

FilteredComplex FilteredComplex::complete_complex(const index_t n, const index_t k) {
	if (k >= n) {
		throw invalid_argument("k must satisfy 0 <= k < n");
	}
	FilteredComplex result(n, k);
	vector<bool>    v(n);
	fill(v.begin(), v.begin() + k + 1, true);
	vector<index_t> verts(k + 1);
	do {
		for (index_t i = 0, j = 0; i < n; ++i) {
			if (v[i]) {
				verts[j++] = i;
			}
		}
		result._add_simplex(verts, 0.0);
	} while (prev_permutation(v.begin(), v.end()));
	result.cur_dim = k;
	return result;
}

bool FilteredComplex::is_filtration() const {
	for (int i = 1; i <= cur_dim; i++) {
		for (auto& s: simplices[i]) {
			for (auto& f: s.second->get_facets()) {
				if (f->value > s.second->value) {
					return false;
				}
			}
		}
	}
	return true;
}

FilteredComplex standard_simplex(const index_t n) {
	return chalc::FilteredComplex::complete_complex(n + 1, n);
}
}  // namespace chalc
