#pragma once

#include <functional>
#include <limits>

#include "index_map.h"
#include "maybe.h"
#include "sequence.h"

using namespace std;

template <class data> struct vertexSubsetData {
    using S = tuple<uintE, data>;
    using D = tuple<bool, data>;

    // An empty vertex set.
    vertexSubsetData(size_t _n) : n(_n), m(0), d(NULL), s(NULL), isDense(0) {}

    // A vertexSubset from array of vertex indices.
    vertexSubsetData(long _n, long _m, S *indices)
        : n(_n), m(_m), s(indices), d(NULL), isDense(0) {}

    // A vertexSubset from boolean array giving number of true values.
    vertexSubsetData(long _n, long _m, D *_d)
        : n(_n), m(_m), s(NULL), d(_d), isDense(1) {}

    // A vertexSubset from boolean array giving number of true values. Calculate
    // number of nonzeros and store in m.
    vertexSubsetData(long _n, D *_d) : n(_n), s(NULL), d(_d), isDense(1) {
        auto d_map = make_in_imap<size_t>(
            n, [&](size_t i) { return (size_t)get<0>(_d[i]); });
        m = pbbs::reduce_add(d_map);
    }

    vertexSubsetData() : n(0), m(0), s(NULL), d(NULL), isDense(0) {}

    void del() {
        if (d != NULL)
            free(d);
        if (s != NULL)
            free(s);
    }

    // Sparse
    inline uintE &vtx(const uintE &i) const {return std::get<0>(s[i]); }
    inline data &vtxData(const uintE &i) const { return std::get<1>(s[i]); }
    inline tuple<uintE, data> vtxAndData(const uintE &i) const { return s[i]; }

    // Dense
    inline bool isIn(const uintE &v) const { return std::get<0>(d[v]); }
    inline data &ithData(const uintE &v) const { return std::get<1>(d[v]); }

    // Returns (uintE) -> Maybe<tuple<vertex, vertex-data>>.
    auto get_fn_repr() const {
        std::function<Maybe<tuple<uintE, data>>(const uintE &)> fn;
        if (isDense) {
            fn = [&](const uintE &v) -> Maybe<tuple<uintE, data>> {
                auto ret =
                    Maybe<tuple<uintE, data>>(make_tuple(v, std::get<1>(d[v])));
                ret.exists = std::get<0>(d[v]);
                return ret;
            };
        } else {
            fn = [&](const uintE &i) -> Maybe<tuple<uintE, data>> {
                return Maybe<tuple<uintE, data>>(s[i]);
            };
        }
        return fn;
    }

    long size() const { return m; }
    long numVertices() const { return n; }

    long numRows() const { return n; }
    long numNonzeros() const { return m; }

    bool isEmpty() const { return m == 0; }
    bool dense() const { return isDense; }

    void toSparse() {
        if (s == NULL && m > 0) {
            auto f = make_in_imap<D>(
                n, [&](size_t i) -> tuple<bool, data> { return d[i]; });
            auto out = pbbs::pack_index_and_data<uintE, data>(f, n);
            out.alloc = false;
            s = out.s;
            if (out.size() != m) {
                cout << "bad stored value of m" << endl;
                abort();
            }
        }
        isDense = false;
    }

    // Convert to dense but keep sparse representation if it exists.
    void toDense() {
        if (d == NULL) {
            d = newA(D, n);
            { parallel_for(long i = 0; i < n; i++) std::get<0>(d[i]) = false; }
            {
                parallel_for(long i = 0; i < m; i++) d[std::get<0>(s[i])] =
                    make_tuple(true, std::get<1>(s[i]));
            }
        }
        isDense = true;
    }

    S *s;
    D *d;
    size_t n, m;
    bool isDense;
};

// Specialized version where data = pbbs::empty.
template <> struct vertexSubsetData<pbbs::empty> {
    using S = uintE;

    // An empty vertex set.
    vertexSubsetData<pbbs::empty>(size_t _n)
        : n(_n), m(0), d(NULL), s(NULL), isDense(0) {}

    vertexSubsetData<pbbs::empty>(size_t _n, bool _isDense)
        : n(_n), m(0), d(NULL), s(NULL), isDense(_isDense) {}

    // A vertexSubset with a single vertex.
    vertexSubsetData<pbbs::empty>(long _n, uintE v)
        : n(_n), m(1), d(NULL), isDense(0) {
        s = newA(uintE, 1);
        s[0] = v;
    }

    // A vertexSubset from array of vertex indices.
    vertexSubsetData<pbbs::empty>(long _n, long _m, S *indices)
        : n(_n), m(_m), s(indices), d(NULL), isDense(0) {}

    // A vertexSubset from array of vertex indices.
    vertexSubsetData<pbbs::empty>(long _n, long _m,
                                  tuple<uintE, pbbs::empty> *indices)
        : n(_n), m(_m), s((uintE *)indices), d(NULL), isDense(0) {}

    // A vertexSubset from boolean array giving number of true values.
    vertexSubsetData<pbbs::empty>(long _n, long _m, bool *_d)
        : n(_n), m(_m), s(NULL), d(_d), isDense(1) {}

    // A vertexSubset from boolean array giving number of true values. Calculate
    // number of nonzeros and store in m.
    vertexSubsetData<pbbs::empty>(long _n, bool *_d)
        : n(_n), s(NULL), d(_d), isDense(1) {
        auto d_map = make_in_imap<size_t>(n, [&](size_t i) { return _d[i]; });
        auto f = [&](size_t i, size_t j) { return i + j; };
        m = pbbs::reduce(d_map, f);
    }

    // A vertexSubset from boolean array giving number of true values. Calculate
    // number of nonzeros and store in m.
    vertexSubsetData<pbbs::empty>(long _n, tuple<bool, pbbs::empty> *_d)
        : n(_n), s(NULL), d((bool *)_d), isDense(1) {
        auto d_map =
            make_in_imap<size_t>(n, [&](size_t i) { return get<0>(_d[i]); });
        auto f = [&](size_t i, size_t j) { return i + j; };
        m = pbbs::reduce(d_map, f);
    }

    void del() {
        if (d != NULL)
            free(d);
        if (s != NULL)
            free(s);
    }

    // Sparse
    inline uintE &vtx(const uintE &i) const {return s[i]; }
    inline pbbs::empty vtxData(const uintE &i) const { return pbbs::empty(); }
    inline tuple<uintE, pbbs::empty> vtxAndData(const uintE &i) const {
        return make_tuple(s[i], pbbs::empty());
    }

    // Dense
    inline bool isIn(const uintE &v) const { return d[v]; }
    inline pbbs::empty ithData(const uintE &v) const { return pbbs::empty(); }

    // Returns (uintE) -> Maybe<tuple<vertex, vertex-data>>.
    auto get_fn_repr() const {
        std::function<Maybe<tuple<uintE, pbbs::empty>>(const uintE &)> fn;
        if (isDense) {
            fn = [&](const uintE &v) -> Maybe<tuple<uintE, pbbs::empty>> {
                auto ret = Maybe<tuple<uintE, pbbs::empty>>(
                    make_tuple(v, pbbs::empty()));
                ret.exists = d[v];
                return ret;
            };
        } else {
            fn = [&](const uintE &i) -> Maybe<tuple<uintE, pbbs::empty>> {
                return Maybe<tuple<uintE, pbbs::empty>>(
                    make_tuple(s[i], pbbs::empty()));
            };
        }
        return fn;
    }

    long size() { return m; }
    long numVertices() { return n; }

    long numRows() { return n; }
    long numNonzeros() { return m; }

    bool isEmpty() { return m == 0; }
    bool dense() { return isDense; }

    void toSparse() {
        if (s == NULL && m > 0) {
            auto _d = d;
            auto f = [&](size_t i) { return _d[i]; };
            auto f_in = make_in_imap<bool>(n, f);
            auto out = pbbs::pack_index<uintE>(f_in);
            out.alloc = false;
            s = out.s;
            if (out.size() != m) {
                cout << "bad stored value of m" << endl;
                cout << "out.size = " << out.size() << " m = " << m
                     << " n = " << n << endl;
                abort();
            }
        }
        else{
            for (int i = 0; i < m; i++)
            if (s[i] > n){
                cout << s[i] << endl;
                break;
            }
        }
        isDense = false;
    }

    // Converts to dense but keeps sparse representation if it exists.
    void toDense() {
        if (d == NULL) {
            d = newA(bool, n);
            { parallel_for(long i = 0; i < n; i++) d[i] = 0; }
            { parallel_for(long i = 0; i < m; i++) d[s[i]] = 1; }
        }
        isDense = true;
    }

    void merge(vertexSubsetData<pbbs::empty> &other) {
        if (isDense || other.isDense) {
            if (!isDense)
            {
                toDense();
                free(s);
                s = NULL;
            }
            if (!other.isDense)
            {
                other.toDense();
                free(other.s);
                other.s = NULL;
            }
            // #pragma omp parallel for reduction(+:m)
            for (long i = 0; i < n; i++) {
                if (d[i] == 0 && other.d[i] == 1)
                {
                    d[i] = 1;
                    m++;
                }
            }
        } else {
            int new_m = m + other.m;
            S *new_s = new S[new_m];
            
            parallel_for (long i = 0; i < m; i++)
            {
                new_s[i] = s[i];
                // if (s[i] > n) cout << "error" << endl;
            }
            parallel_for (long i = m; i < new_m; i++)
            {
                new_s[i] = other.s[i - m];
                // if (other.s[i - m] > n) cout << "error" << endl;
            }
            free(s);
            sort(new_s, new_s + new_m);
            int pl = 1;
            for (int i = 1; i < new_m; i++)
            {
                if (new_s[i] != new_s[i-1])
                {
                    new_s[pl] = new_s[i];
                    pl++;
                }
            }
            s = new_s;
            m = pl;
            if (m > n/20)
            {
                toDense();
                free(s);
                s = NULL;
            }
        }
    }

    S *s;
    bool *d;
    size_t n, m;
    bool isDense;
};

using vertexSubset = vertexSubsetData<pbbs::empty>;