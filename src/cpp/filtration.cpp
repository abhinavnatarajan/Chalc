#include "filtration.h"

#include <stdexcept>
#include <algorithm>
namespace chalc {
    constexpr value_t DEFAULT_FILT_VALUE = 0.0;

    void check_overflow(size_t i) {
        if (i < 0)
            throw std::overflow_error("simplex index is too large");
    }

    BinomialCoeffTable::BinomialCoeffTable(size_t n, size_t k) : B(n + 1) {
        for (size_t i = 0; i <= n; ++i) {
            B[i].resize(i + 1, 0);
            B[i][0] = 1;
            for (size_t j = 1; j < std::min(i, k + 1); ++j) {
                B[i][j] = B[i - 1][j - 1] + B[i - 1][j];
            }
            if (i <= k) { B[i][i] = 1; }
            check_overflow(B[i][std::min(i >> 1, k)]);
        }
    }

    size_t BinomialCoeffTable::operator()(size_t n, size_t k) const {
        //assert(n < B.size() && k < B[n].size() && n >= k - 1);
        return B[n][k];
    }

    FilteredComplex::Simplex::Simplex(size_t label = 0, size_t max_vertex = 0, size_t dim = 0, value_t value = DEFAULT_FILT_VALUE,
        const vector<shared_ptr<Simplex>>& facets = vector<shared_ptr<Simplex>>{}) :
        label(label), max_vertex(max_vertex), dim(dim), value(value), facets(facets) {}

    vector<size_t> FilteredComplex::Simplex::get_vertex_labels() const {
        vector<size_t> result(dim + 1);
        get_vertex_labels(result.begin());
        return result;
    }

    template <typename OutputIterator>
    void FilteredComplex::Simplex::get_vertex_labels(OutputIterator&& buf) const {
        if (dim != 0) {
            // vertices of a simplex are the vertices of its last face along with its last vertex
            facets.back()->get_vertex_labels(buf);
            buf++;
        }
        *buf = max_vertex;
    }

    vector<size_t> FilteredComplex::Simplex::get_facet_labels() const {
        vector<size_t> result;
        if (dim != 0) {
            result.reserve(facets.size());
            for (auto& f : facets) {
                result.push_back(f->label);
            }
        }
        return result;
    }

    const vector<shared_ptr<FilteredComplex::Simplex>>& FilteredComplex::Simplex::get_facets() const {
        return facets;
    }

    FilteredComplex::FilteredComplex(const size_t num_vertices, const size_t max_dimension) :
        binomial(num_vertices, max_dimension + 1), // we want nCk for all 0 <= n <= num_vertices and 0 <= k <= max_num_verts_in_a_simplex = max_dim + 1
        simplices(max_dimension + 1),
        N(num_vertices),
        max_dim(max_dimension),
        num_simplices(num_vertices),
        cur_dim(0) {
        if (max_dim < 0) {
            throw std::invalid_argument("Dimension cannot be negative.");
        }
        if (N <= 0) {
            throw std::invalid_argument("Number of points must be positive.");
        }
        if (max_dim >= N) {
            throw std::invalid_argument("Dimension must be less than number of points.");
        }
        for (size_t i = 0; i < N; i++) {
            simplices[0][i] = std::make_shared<Simplex>(i, i);
        }
    }

    FilteredComplex::FilteredComplex(const FilteredComplex& other, const size_t k) :
        binomial(other.size_in_dim(0), other.max_dim + 1),
        simplices(other.max_dim + 1),
        N(other.size_in_dim(0)),
        max_dim(other.max_dim),
        cur_dim(k) {
        if (k > max_dim) {
            throw std::invalid_argument("Dimension of subcomplex must be less than max dimension of original complex.");
        }
        if (k < 0) {
            throw std::invalid_argument("Dimension must be non-negative.");
        }
        num_simplices = 0;
        for (size_t d = 0; d <= k; d++) {
            auto& other_dsimplices = other.get_simplices()[d];
            for (auto simplex = other_dsimplices.cbegin(); simplex != other_dsimplices.cend(); simplex++) {
                simplices[d][simplex->first] = simplex->second;
            }
            num_simplices += simplices[d].size();
        }
    }

    void FilteredComplex::check_vertex_sequence_is_valid(vector<size_t>& verts) const {
        check_dimension_is_valid(verts.size() - 1);
        std::sort(verts.begin(), verts.end());
        if (!(verts.back() < N &&
            verts.front() >= 0 &&
            std::adjacent_find(verts.begin(), verts.end()) == verts.end())) {
            throw std::invalid_argument("Invalid vertex sequence.");
        };
    }

    void FilteredComplex::check_dimension_is_valid(const size_t dim) const {
        if (dim > max_dim || dim < 0) {
            throw std::invalid_argument("Invalid dimension.");
        }
    }

    size_t FilteredComplex::_get_label_from_vertex_labels(const vector<size_t>& verts) const {
        size_t label = 0;
        auto dim = verts.size() - 1;
        for (size_t i = 0, prev_vert = -1; i < verts.size(); i++) {
            for (size_t j = prev_vert + 1; j < verts[i]; j++) {
                label += binomial(N - j - 1, dim - i);
            }
            prev_vert = verts[i];
        }
        return label;
    }

    size_t FilteredComplex::get_label_from_vertex_labels(vector<size_t>& verts) const {
        check_vertex_sequence_is_valid(verts);
        return _get_label_from_vertex_labels(verts);
    }

    bool FilteredComplex::_has_simplex(const size_t dim, const size_t label) const {
        if (simplices[dim].find(label) == simplices[dim].end()) { return false; }
        else { return true; }
    }

    bool FilteredComplex::has_simplex(const size_t dim, const size_t label) const {
        check_dimension_is_valid(dim);
        if (label < 0 || label >= binomial(N, dim)) {
            throw std::invalid_argument("Invalid label.");
        }
        return _has_simplex(dim, label);
    }


    bool FilteredComplex::_has_simplex(const vector<size_t>& verts) const {
        auto dim = verts.size() - 1;
        auto label = _get_label_from_vertex_labels(verts);
        return (_has_simplex(dim, label));
    }

    bool FilteredComplex::has_simplex(vector<size_t>& verts) const {
        check_vertex_sequence_is_valid(verts);
        return (_has_simplex(verts));
    }


    shared_ptr<FilteredComplex::Simplex> FilteredComplex::_add_simplex(const vector<size_t>& verts, const value_t filt_value) {
        shared_ptr<Simplex> new_simplex;
        auto dim = verts.size() - 1;
        auto label = _get_label_from_vertex_labels(verts);
        auto search_simplex = simplices[dim].find(label);
        if (search_simplex != simplices[dim].end()) { // the simplex already exists
            new_simplex = search_simplex->second;
            new_simplex->value = std::min(filt_value, new_simplex->value);
        }
        else { // simplex does not exist so we need to add it
            // first we recursively add faces of the simplex
            vector<shared_ptr<Simplex>> facets(dim + 1);
            vector<size_t> facet_verts(dim);
            for (auto i = 0; i <= dim; i++) {
                for (auto j = 0, k = 0; j < verts.size(); j++) {
                    if (j == i) continue;
                    facet_verts[k] = verts[j];
                    k++;
                }
                facets[i] = _add_simplex(facet_verts, filt_value);
            }
            auto max_vertex = verts.back();
            new_simplex = std::make_shared<Simplex>(label, max_vertex, dim, filt_value, std::move(facets));
            simplices[dim][label] = new_simplex;
            num_simplices++;
        }
        return new_simplex;
    }

    bool FilteredComplex::add_simplex(vector<size_t>& verts, const value_t filt_value = DEFAULT_FILT_VALUE) {
        check_vertex_sequence_is_valid(verts);
        if (_has_simplex(verts)) {
            return false;
        }
        else {
            _add_simplex(verts, filt_value);
            cur_dim = std::max(cur_dim, verts.size() - 1);
            return true;
        }
    }

    void FilteredComplex::propagate_filt_values_up(const size_t start_dim) {
        auto p = start_dim + 1;
        while (p <= cur_dim) {
            value_t tmp;
            // iterate over the p-simplices
            for (auto& [label, simplex] : simplices[p]) {
                tmp = -1;
                // iterate over faces of simplex and get the maximum filtration value
                for (const auto& facet : simplex->get_facets()) {
                    tmp = std::max(tmp, facet->value);
                }
                simplex->value = tmp;
            }
            p++;
        }
    }

    void FilteredComplex::propagate_filt_values_down(const size_t start_dim) {
        auto p = start_dim;
        while (p >= 1) {
            // iterate over the p-simplices
            for (const auto& [label, simplex] : simplices[p - 1]) {
                // iterate over faces of simplex and modify the filtration value if needed
                for (auto& facet : simplex->get_facets()) {
                    facet->value = std::min(simplex->value, facet->value);
                }
            }
            p--;
        }
    }

    void FilteredComplex::propagate_filt_values(const size_t start_dim, const bool up) {
        if (start_dim > cur_dim || start_dim < 0) {
            throw std::invalid_argument("Invalid starting dimension.");
        }
        if (up) { propagate_filt_values_up(start_dim); }
        else { propagate_filt_values_down(start_dim); }
    }

    const vector<map<size_t, shared_ptr<FilteredComplex::Simplex>>>& FilteredComplex::get_simplices() const noexcept {
        return simplices;
    }

    size_t FilteredComplex::size() const noexcept { return num_simplices; }

    size_t FilteredComplex::size_in_dim(const size_t dim) const {
        check_dimension_is_valid(dim);
        return simplices[dim].size();
    }

    size_t FilteredComplex::dimension() const noexcept {
        return cur_dim;
    }

    vector<tuple<vector<size_t>, value_t>> FilteredComplex::flat_representation() const {
        vector<tuple<vector<size_t>, value_t>> result(num_simplices);
        vector<map<size_t, size_t>> indices(cur_dim + 1);
        // add all the vertices
        size_t i = 0;
        for (auto& [idx, vertex] : simplices[0]) {
            indices[0][idx] = i;
            result[i++] = tuple<vector<size_t>, value_t>{ vector<size_t>{}, vertex->value };
        }
        for (size_t dim = 1; dim <= cur_dim; dim++) {
            for (auto& [label, simplex] : simplices[dim]) {
                // replace the labels of the faces with their corresponding indices in our result
                vector<size_t> faces = simplex->get_facet_labels();
                for (auto& f : faces) { f = indices[dim - 1][f]; }
                indices[dim][label] = i;
                result[i++] = tuple{ faces, simplex->value };
            }
        }
        return result;
    }

    FilteredComplex FilteredComplex::clique_complex(const size_t n, const size_t k) {
        if (n <= 0) {
            throw std::invalid_argument("number of vertices must be >= 0.");
        }
        if (k < 0 || k >= n) {
            throw std::invalid_argument("k must satisfy 0 <= k < n");
        }
        FilteredComplex result(n, k);
        vector<bool> v(n);
        std::fill(v.begin(), v.begin() + k + 1, true);
        vector<size_t> verts(k + 1);
        do {
            for (size_t i = 0, j = 0; i < n; ++i) {
                if (v[i]) {
                    verts[j++] = i;
                }
            }
            result._add_simplex(verts, 0.0);
        } while (std::prev_permutation(v.begin(), v.end()));
        result.cur_dim = k;
        return result;
    }

    FilteredComplex standard_simplex(const size_t n) {
        return chalc::FilteredComplex::clique_complex(n + 1, n);
    }
}