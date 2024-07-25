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
#include <chalc/filtration/filtration.h>
#include <stdexcept>
#ifndef NDEBUG
	#define NDEBUG
#endif
#include <algorithm>
#include <cassert>

namespace {
using namespace chalc::stl;
using std::sort, std::stable_sort, std::min, std::max, std::fill, std::adjacent_find,
	std::prev_permutation, std::invalid_argument, std::runtime_error, std::bad_weak_ptr;
}  // namespace

namespace chalc {
class BinomialCoeffTable {
	vector<vector<index_t>> B;

  public:
	BinomialCoeffTable(index_t n, index_t k) :
		B(n + 1) {
		for (index_t i = 0; i <= n; ++i) {
			B[i].resize(i + 1, 0);
			B[i][0] = 1;
			for (index_t j = 1; j < min(i, k + 1); ++j) {
				B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
			}
			if (i <= k) {
				B[i][i] = 1;
			}
			assert(("Simplex index is too large.", B[i][min(i >> 1, k)] >= 0));
		}
	}

	index_t operator()(index_t n, index_t k) const {
		assert(n < B.size() && k < B[n].size() && n >= k - 1);
		return B[n][k];
	}
};

shared_ptr<FilteredComplex::Simplex>
FilteredComplex::Simplex::make_Simplex(index_t                            label,
                                       index_t                            max_vertex,
                                       value_t                            value,
                                       const vector<shared_ptr<Simplex>>& facets) {
	auto self = shared_ptr<Simplex>(new Simplex(label, max_vertex, value, facets));
	for (auto& f: self->facets) {
		f->cofacets.push_back(self->get_handle());
		self->add_colours(f->colours);
	}
	return self;
}

FilteredComplex::Simplex::Simplex(index_t                            label,
                                  index_t                            max_vertex,
                                  value_t                            value,
                                  const vector<shared_ptr<Simplex>>& facets) :
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
	get_vertex_labels(result.begin());
	return result;
}

template <typename OutputIterator>
void FilteredComplex::Simplex::get_vertex_labels(OutputIterator&& buf) const {
	if (dim > 0) {
		// vertices of a simplex are the vertices of its last face along with
		// its last vertex
		facets.back()->get_vertex_labels(buf);
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

inline void FilteredComplex::Simplex::add_colours(colours_t c) {
	colours |= c;
}

inline void FilteredComplex::Simplex::set_colours(colours_t c) {
	colours.reset();
	add_colours(c);
}

inline void FilteredComplex::Simplex::make_colourless() {
	colours.reset();
}

void FilteredComplex::Simplex::set_colour(index_t c) {
	colours.reset().set(c);
}

colours_t FilteredComplex::Simplex::get_colours() {
	return colours;
}

unsigned long long int FilteredComplex::Simplex::get_colours_as_int() {
	return colours.to_ullong();
}

FilteredComplex::FilteredComplex(const index_t num_vertices, const index_t max_dimension) :
	binomial(make_shared<BinomialCoeffTable>(
		num_vertices,
		max_dimension + 1)),  // we want nCk for all 0 <= n <= num_vertices and 0 <= k <=
                              // max_num_verts_in_a_simplex = max_dim + 1
	simplices(max_dimension + 1),
	N(num_vertices),
	max_dim(max_dimension),
	num_simplices(num_vertices),
	cur_dim(0) {
	if (max_dim < 0) {
		throw invalid_argument("Dimension cannot be negative.");
	}
	if (N <= 0) {
		throw invalid_argument("Number of points must be positive.");
	}
	if (max_dim >= N) {
		throw invalid_argument("Dimension must be less than number of points.");
	}
	for (index_t i = 0; i < N; i++) {
		simplices[0][i] = Simplex::make_Simplex(i, i, 0.0);
		// simplices are initialised with colours unset
		// so we set them here
		simplices[0][i]->set_colour(0);
	}
}

void FilteredComplex::check_vertex_sequence_is_valid(vector<index_t>& verts) const {
	check_dimension_is_valid(verts.size() - 1);
	sort(verts.begin(), verts.end());
	if (!(verts.back() < N && verts.front() >= 0 &&
	      adjacent_find(verts.begin(), verts.end()) == verts.end())) {
		throw invalid_argument("Invalid vertex sequence.");
	};
}

void FilteredComplex::check_dimension_is_valid(const index_t dim) const {
	if (dim > max_dim || dim < 0) {
		throw invalid_argument("Invalid dimension.");
	}
}

index_t FilteredComplex::_get_label_from_vertex_labels(const vector<index_t>& verts) const {
	index_t label = 0;
	auto    dim   = verts.size() - 1;
	for (index_t i = 0, prev_vert = -1; i < verts.size(); i++) {
		for (index_t j = prev_vert + 1; j < verts[i]; j++) {
			label += (*binomial)(N - j - 1, dim - i);
		}
		prev_vert = verts[i];
	}
	return label;
}

index_t FilteredComplex::get_label_from_vertex_labels(vector<index_t>& verts) const {
	check_vertex_sequence_is_valid(verts);
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
	check_dimension_is_valid(dim);
	if (label < 0 || label >= (*binomial)(N, dim)) {
		throw invalid_argument("Invalid label.");
	}
	return _has_simplex(dim, label);
}

bool FilteredComplex::_has_simplex(const vector<index_t>& verts) const {
	auto dim   = verts.size() - 1;
	auto label = _get_label_from_vertex_labels(verts);
	return (_has_simplex(dim, label));
}

bool FilteredComplex::has_simplex(vector<index_t>& verts) const {
	check_vertex_sequence_is_valid(verts);
	return (_has_simplex(verts));
}

shared_ptr<FilteredComplex::Simplex> FilteredComplex::_add_simplex(const vector<index_t>& verts,
                                                                   const value_t filt_value) {
	shared_ptr<Simplex> new_simplex;
	auto                dim            = verts.size() - 1;
	auto                label          = _get_label_from_vertex_labels(verts);
	auto                search_simplex = simplices[dim].find(label);
	if (search_simplex != simplices[dim].end()) {
		// the simplex already exists
		new_simplex        = search_simplex->second;
		new_simplex->value = min(max(filt_value, static_cast<value_t>(0)), new_simplex->value);
	} else {
		// simplex does not exist so we need to add it
		// first we recursively add faces of the simplex
		vector<shared_ptr<Simplex>> facets(dim + 1);
		vector<index_t>             facet_verts(dim);
		for (auto i = 0; i <= dim; i++) {
			for (auto j = 0, k = 0; j < verts.size(); j++) {
				if (j == i) {
					continue;
				}
				facet_verts[k] = verts[j];
				k++;
			}
			facets[i] = _add_simplex(facet_verts, filt_value);
		}
		auto max_vertex = verts.back();
		new_simplex     = Simplex::make_Simplex(label, max_vertex, filt_value, std::move(facets));
		simplices[dim][label] = new_simplex;
		num_simplices++;
	}
	return new_simplex;
}

bool FilteredComplex::add_simplex(
	vector<index_t>& verts,
	value_t          filt_value = FilteredComplex::Simplex::DEFAULT_FILT_VALUE) {
	check_vertex_sequence_is_valid(verts);
	if (_has_simplex(verts)) {
		return false;
	} else {
		_add_simplex(verts, filt_value);
		cur_dim            = max(static_cast<size_t>(cur_dim), verts.size() - 1);
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
	int p = static_cast<int>(start_dim) - 1;
	while (p >= 1) {
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
		p--;
	}
}

void FilteredComplex::propagate_colours() {
	for (index_t d = 1; d <= cur_dim; d++) {
		for (auto& [idx, s]: simplices[d]) {
			s->make_colourless();
			for (auto& f: s->get_facets()) {
				s->add_colours(f->get_colours());
			}
		}
	}
}

void FilteredComplex::propagate_filt_values(const index_t start_dim, const bool up) {
	if (start_dim > cur_dim || start_dim < 0) {
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

value_t FilteredComplex::max_filt_value() const noexcept {
	return cur_max_filt_value;
}

vector<tuple<vector<index_t>, index_t, value_t, unsigned long long int>>
FilteredComplex::serialised() const {
	vector<tuple<vector<index_t>, index_t, value_t, unsigned long long int>> result(num_simplices);
	vector<map<index_t, index_t>>                                            indices(cur_dim + 1);
	for (index_t d = 0, i = 0; d <= cur_dim; d++) {
		// sort the d-dimensional simplices by filtration value
		vector<shared_ptr<Simplex>> sort_by_val;
		sort_by_val.reserve(simplices[d].size());
		for (auto& s: simplices[d]) {
			sort_by_val.push_back(s.second);
		}
		stable_sort(sort_by_val.begin(),
		            sort_by_val.end(),
		            [](const shared_ptr<Simplex>& s1, const shared_ptr<Simplex>& s2) {
						return (s1->value < s2->value);
					});
		// iterate over the sorted simplices
		for (auto& simplex: sort_by_val) {
			// replace the labels of the faces with their corresponding indices
			// in our result faces is empty if simplex is a vertex
			vector<index_t> faces = simplex->get_facet_labels();
			for (auto& f: faces) {
				f = indices[d - 1][f];
			}
			indices[d][simplex->label] = i;
			sort(faces.begin(), faces.end());
			result[i++] =
				tuple{faces, simplex->label, simplex->value, simplex->get_colours_as_int()};
		}
	}
	return result;
}

FilteredComplex FilteredComplex::clique_complex(const index_t n, const index_t k) {
	if (n <= 0) {
		throw invalid_argument("number of vertices must be >= 0.");
	}
	if (k < 0 || k >= n) {
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
	for (int i = cur_dim; i >= 0; i--) {
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
	return chalc::FilteredComplex::clique_complex(n + 1, n);
}
}  // namespace chalc
