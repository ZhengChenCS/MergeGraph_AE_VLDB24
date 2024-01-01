#ifndef VERTEX_H
#define VERTEX_H
#include "vertexSubset.h"
using namespace std;

namespace decode_uncompressed {

// Used by edgeMapDense. Callers ensure cond(v_id). For each vertex, decode
// its in-edges, and check to see whether this neighbor is in the current
// frontier, calling update if it is. If processing the edges sequentially,
// break once !cond(v_id).
template <class vertex, class F, class G, class VS>
inline void decodeInNghBreakEarly(vertex* v,
                                  long v_id,
                                  VS& vertexSubset,
                                  F& f,
                                  G& g,
                                  bool parallel = 0) {
    uintE d = v->getInDegree();
    if (!parallel || d < 1000) {
        for (size_t j = 0; j < d; j++) {
            uintE ngh = v->getInNeighbor(j);
            if (vertexSubset.isIn(ngh)) {
#ifndef WEIGHTED
                auto m = f.update(ngh, v_id);
#else
                auto m = f.update(ngh, v_id, v->getInWeight(j));
#endif
                g(v_id, m);
            }
            if (!f.cond(v_id))
                break;
        }
    } else {
        parallel_for(size_t j = 0; j < d; j++) {
            uintE ngh = v->getInNeighbor(j);
            if (vertexSubset.isIn(ngh)) {
#ifndef WEIGHTED
                auto m = f.updateAtomic(ngh, v_id);
#else
                auto m = f.updateAtomic(ngh, v_id, v->getInWeight(j));
#endif
                g(v_id, m);
            }
        }
    }
}

// Used by edgeMapDenseForward. For each out-neighbor satisfying cond, call
// updateAtomic.
template <class V, class F, class G>
inline void decodeOutNgh(V* v, long i, F& f, G& g) {
    uintE d = v->getOutDegree();
    #ifndef OPTIMIZE_GET_NEIGHBOR
    granular_for(j, 0, d, (d > 1000), {
        uintE ngh = v->getOutNeighbor(j);
        if (f.cond(ngh)) {
#ifndef WEIGHTED
            auto m = f.updateAtomic(i, ngh);
#else
            auto m = f.updateAtomic(i,ngh,v->getOutWeight(j));
#endif
            g(ngh, m);
        }
    });
    #else
    uintT batch = 256;
    if (d > 4*batch)
    {
        parallel_for(int base = 0; base < d; base += batch)
        {
            uintT gid = 0;
            for (; gid < v->subgraphNum; gid++)
            {
                if (v->outOffsets[gid] <= base && v->outOffsets[gid + 1] > base)
                {
                    break;
                }
                else if (gid == (v->subgraphNum) - 1)
                {
                    break;
                }
            }
            uintT ed = min(base + batch, (uintT)d);
            for (int j = base; j < ed; j++)
            {
                uintE ngh = v->getOutNeighbor(j, gid);
                if (f.cond(ngh)) {
                    #ifndef WEIGHTED
                    auto m = f.updateAtomic(i,ngh);
                    #else
                    auto m = f.updateAtomic(i,ngh,v->getOutWeight(j, gid));
                    #endif
                    g(ngh, m);
                }
            }
        }
    }
    else
    {
        uintT gid = 0;
        for(int e = 0; e < d; e ++)
        {
            uintE ngh = v->getOutNeighbor(e, gid);
            if (f.cond(ngh)) {
                #ifndef WEIGHTED
                auto m = f.updateAtomic(i,ngh);
                #else
                auto m = f.updateAtomic(i,ngh,v->getOutWeight(e, gid));
                #endif
                g(ngh, m);
            }
        }
    }
    #endif
}

template <class V, class F, class G>
inline void decodeOutNgh_single(V* v, long i, F& f, G& g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
        uintE ngh = v->getOutNeighbor(j);
        if (f.cond(ngh)) {
#ifndef WEIGHTED
            auto m = f.updateAtomic(i, ngh);
#else
            auto m = f.updateAtomic(i,ngh,v->getOutWeight(j));
#endif
            g(ngh, m);
        }
    });
}

// Used by edgeMapSparse. For each out-neighbor satisfying cond, call
// updateAtomic.
template <class V, class F, class G>
inline void decodeOutNghSparse(V* v, long i, uintT o, F& f, G& g) {
    uintE d = v->getOutDegree();
    // #ifndef OPTIMIZE_GET_NEIGHBOR
    granular_for(j, 0, d, (d > 1000), {
        uintE ngh = v->getOutNeighbor(j);
        if (f.cond(ngh)) {
#ifndef WEIGHTED
            auto m = f.updateAtomic(i, ngh);
#else
            auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
            g(ngh, o + j, m);
        } else {
            g(ngh, o + j);
        }
    });
    // #else
    // uintT batch = 256;
    // if (d > 4*batch)
    // {
    //     parallel_for(int base = 0; base < d; base += batch)
    //     {
    //         uintT gid = 0;
    //         for (; gid < v->subgraphNum; gid++)
    //         {
    //             if (v->outOffsets[gid] <= base && v->outOffsets[gid + 1] > base)
    //             {
    //                 break;
    //             }
    //             else if (gid == (v->subgraphNum) - 1)
    //             {
    //                 break;
    //             }
    //         }
    //         uintT ed = min(base + batch, (uintT)d);
    //         for (int j = base; j < ed; j++)
    //         {
    //             uintE ngh = v->getOutNeighbor(j, gid);
    //             if (f.cond(ngh)) {
    //             #ifndef WEIGHTED
    //                 auto m = f.updateAtomic(i, ngh);
    //             #else
    //                 auto m = f.updateAtomic(i, ngh, v->getOutWeight(j, gid));
    //             #endif
    //                 g(ngh, o + j, m);
    //             } else {
    //                 g(ngh, o + j);
    //             }
    //         }
    //     }
    // }
    // else
    // {
    //     uintT gid = 0;
    //     for(int e = 0; e < d; e ++)
    //     {
    //         uintE ngh = v->getOutNeighbor(e, gid);
    //         if (f.cond(ngh)) {
    //         #ifndef WEIGHTED
    //             auto m = f.updateAtomic(i, ngh);
    //         #else
    //             auto m = f.updateAtomic(i, ngh, v->getOutWeight(e, gid));
    //         #endif
    //             g(ngh, o + e, m);
    //         } else {
    //             g(ngh, o + e);
    //         }
    //     }
    // }
    // #endif
}

// Used by edgeMapSparse_no_filter. Sequentially decode the out-neighbors,
// and compactly write all neighbors satisfying g().
template <class V, class F, class G>
inline size_t decodeOutNghSparseSeq(V* v, long i, uintT o, F& f, G& g) {
    uintE d = v->getOutDegree();
    size_t k = 0;
    for (size_t j = 0; j < d; j++) {
        uintE ngh = v->getOutNeighbor(j);
        if (f.cond(ngh)) {
#ifndef WEIGHTED
            auto m = f.updateAtomic(i, ngh);
#else
            auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
            bool wrote = g(ngh, o + k, m);
            if (wrote) {
                k++;
            }
        }
    }
    return k;
}

template <class V, class F, class G>
inline size_t decodeOutNghSparseSeq(V* v,
                                    long i,
                                    uintT o,
                                    F& f,
                                    G& g,
                                    uint64_t base) {
    uintE d = v->getOutDegree();
    size_t k = 0;
    for (size_t j = 0; j < d; j++) {
        uintE ngh = v->getOutNeighbor(j);
        if (f.cond(ngh)) {
#ifndef WEIGHTED
            auto m = f.updateAtomic(i, ngh);
#else
            auto m = f.updateAtomic(i, ngh, v->getOutWeight(j));
#endif
            bool wrote = g(ngh, o + k + base, m);
            if (wrote) {
                k++;
            }
        }
    }
    return k;
}

// Decode the out-neighbors of v, and return the number of neighbors
// that satisfy f.
template <class V, class F>
inline size_t countOutNgh(V* v, long vtx_id, F& f) {
    uintE d = v->getOutDegree();
    if (d < 2000) {
        size_t ct = 0;
        for (size_t i = 0; i < d; i++) {
            uintE ngh = v->getOutNeighbor(i);
#ifndef WEIGHTED
            if (f(vtx_id, ngh))
#else
            if (f(vtx_id, ngh, v->getOutWeight(i)))
#endif
                ct++;
        }
        return ct;
    } else {
        size_t b_size = 2000;
        size_t blocks = 1 + ((d - 1) / b_size);
        auto cts = array_imap<uintE>(blocks, [&](size_t i) { return 0; });
        parallel_for_1(size_t i = 0; i < blocks; i++) {
            size_t s = b_size * i;
            size_t e = std::min(s + b_size, (size_t)d);
            uintE ct = 0;
            for (size_t j = s; j < e; j++) {
                uintE ngh = v->getOutNeighbor(j);
#ifndef WEIGHTED
                if (f(vtx_id, ngh))
#else
                if (f(vtx_id, ngh, v->getOutNeighbor(j)))
#endif
                    ct++;
            }
            cts[i] = ct;
        }
        size_t count = 0;
        return pbbs::reduce_add(cts);
    }
}

// Decode the out-neighbors of v. Apply f(src, ngh) and store the result
// using g.
template <class V, class E, class F, class G>
inline void copyOutNgh(V* v, long src, uintT o, F& f, G& g) {
    uintE d = v->getOutDegree();
    granular_for(j, 0, d, (d > 1000), {
        uintE ngh = v->getOutNeighbor(j);
#ifdef WEIGHTED
        E val = f(src, ngh, v->getOutWeight(j));
#else
      E val = f(src, ngh);
#endif
        g(ngh, o + j, val);
    });
}

// TODO(laxmand): Add support for weighted graphs.
template <class V, class Pred>
inline size_t packOutNgh(V* v, long vtx_id, Pred& p, bool* bits, uintE* tmp) {
    uintE d = v->getOutDegree();
    if (d < 5000) {
        size_t k = 0;
        for (size_t i = 0; i < d; i++) {
            uintE ngh = v->getOutNeighbor(i);
            if (p(vtx_id, ngh)) {
                v->setOutNeighbor(k, ngh);
                k++;
            }
        }
        v->setOutDegree(k);
        return k;
    } else {
        parallel_for(size_t i = 0; i < d; i++) {
            uintE ngh = v->getOutNeighbor(i);
            tmp[i] = ngh;
            bits[i] = p(vtx_id, ngh);
        }
        size_t k = sequence::pack(tmp, v->getOutNeighbors(), bits, d);
        v->setOutDegree(k);
        return k;
    }
}

}  // namespace decode_uncompressed

struct symmetricVertex {
#ifndef WEIGHTED
    uintE* neighbors;
#else
    intE* neighbors;
#endif
    uintT degree;
    void del() { free(neighbors); }
#ifndef WEIGHTED
    symmetricVertex(uintE* n, uintT d)
#else
    symmetricVertex(intE* n, uintT d)
#endif
        : neighbors(n), degree(d) {
    }
#ifndef WEIGHTED
    uintE* getInNeighbors() { return neighbors; }
    const uintE* getInNeighbors() const { return neighbors; }
    uintE* getOutNeighbors() { return neighbors; }
    const uintE* getOutNeighbors() const { return neighbors; }
    uintE getInNeighbor(uintT j) const { return neighbors[j]; }
    uintE getOutNeighbor(uintT j) const { return neighbors[j]; }

    void setInNeighbor(uintT j, uintE ngh) { neighbors[j] = ngh; }
    void setOutNeighbor(uintT j, uintE ngh) { neighbors[j] = ngh; }
    void setInNeighbors(uintE* _i) { neighbors = _i; }
    void setOutNeighbors(uintE* _i) { neighbors = _i; }
#else
    // weights are stored in the entry after the neighbor ID
    // so size of neighbor list is twice the degree
    intE* getInNeighbors() { return neighbors; }
    const intE* getInNeighbors() const { return neighbors; }
    intE* getOutNeighbors() { return neighbors; }
    const intE* getOutNeighbors() const { return neighbors; }
    intE getInNeighbor(intT j) const { return neighbors[2 * j]; }
    intE getOutNeighbor(intT j) const { return neighbors[2 * j]; }
    intE getInWeight(intT j) const { return neighbors[2 * j + 1]; }
    intE getOutWeight(intT j) const { return neighbors[2 * j + 1]; }
    void setInNeighbor(uintT j, uintE ngh) { neighbors[2 * j] = ngh; }
    void setOutNeighbor(uintT j, uintE ngh) { neighbors[2 * j] = ngh; }
    void setInWeight(uintT j, intE wgh) { neighbors[2 * j + 1] = wgh; }
    void setOutWeight(uintT j, intE wgh) { neighbors[2 * j + 1] = wgh; }
    void setInNeighbors(intE* _i) { neighbors = _i; }
    void setOutNeighbors(intE* _i) { neighbors = _i; }
#endif

    uintT getInDegree() const { return degree; }
    uintT getOutDegree() const { return degree; }
    void setInDegree(uintT _d) { degree = _d; }
    void setOutDegree(uintT _d) { degree = _d; }
    void flipEdges() {}

    template <class VS, class F, class G>
    inline void decodeInNghBreakEarly(long v_id,
                                      VS& vertexSubset,
                                      F& f,
                                      G& g,
                                      bool parallel = 0) {
        decode_uncompressed::decodeInNghBreakEarly<symmetricVertex, F, G, VS>(
            this, v_id, vertexSubset, f, g, parallel);
    }

    template <class F, class G>
    inline void decodeOutNgh(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh<symmetricVertex, F, G>(this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNgh_single(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh_single<symmetricVertex, F, G>(this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNghSparse(long i, uintT o, F& f, G& g) {
        decode_uncompressed::decodeOutNghSparse<symmetricVertex, F>(this, i, o,
                                                                    f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i, uintT o, F& f, G& g) {
        return decode_uncompressed::decodeOutNghSparseSeq<symmetricVertex, F>(
            this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i,
                                        uintT o,
                                        F& f,
                                        G& g,
                                        uint64_t base) {
        return decode_uncompressed::decodeOutNghSparseSeq<symmetricVertex, F>(
            this, i, o, f, g, base);
    }

    template <class E, class F, class G>
    inline void copyOutNgh(long i, uintT o, F& f, G& g) {
        decode_uncompressed::copyOutNgh<symmetricVertex, E>(this, i, o, f, g);
    }

    template <class F>
    inline size_t countOutNgh(long i, F& f) {
        return decode_uncompressed::countOutNgh<symmetricVertex, F>(this, i, f);
    }

    template <class F>
    inline size_t packOutNgh(long i,
                             F& f,
                             bool* bits,
                             uintE* tmp1,
                             uintE* tmp2) {
        return decode_uncompressed::packOutNgh<symmetricVertex, F>(this, i, f,
                                                                   bits, tmp1);
    }
};

struct asymmetricVertex {
#ifndef WEIGHTED
    uintE *inNeighbors, *outNeighbors;
#else
    intE *inNeighbors, *outNeighbors;
#endif
    uintT outDegree;
    uintT inDegree;
    void del() {
        free(inNeighbors);
        free(outNeighbors);
    }
#ifndef WEIGHTED
    asymmetricVertex(uintE* iN, uintE* oN, uintT id, uintT od)
#else
    asymmetricVertex(intE* iN, intE* oN, uintT id, uintT od)
#endif
        : inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {
    }
#ifndef WEIGHTED
    uintE* getInNeighbors() { return inNeighbors; }
    const uintE* getInNeighbors() const { return inNeighbors; }
    uintE* getOutNeighbors() { return outNeighbors; }
    const uintE* getOutNeighbors() const { return outNeighbors; }
    uintE getInNeighbor(uintT j) const { return inNeighbors[j]; }
    uintE getOutNeighbor(uintT j) const { return outNeighbors[j]; }
    void setInNeighbor(uintT j, uintE ngh) { inNeighbors[j] = ngh; }
    void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[j] = ngh; }
    void setInNeighbors(uintE* _i) { inNeighbors = _i; }
    void setOutNeighbors(uintE* _i) { outNeighbors = _i; }
#else
    intE* getInNeighbors() { return inNeighbors; }
    const intE* getInNeighbors() const { return inNeighbors; }
    intE* getOutNeighbors() { return outNeighbors; }
    const intE* getOutNeighbors() const { return outNeighbors; }
    intE getInNeighbor(uintT j) const { return inNeighbors[2 * j]; }
    intE getOutNeighbor(uintT j) const { return outNeighbors[2 * j]; }
    intE getInWeight(uintT j) const { return inNeighbors[2 * j + 1]; }
    intE getOutWeight(uintT j) const { return outNeighbors[2 * j + 1]; }
    void setInNeighbor(uintT j, uintE ngh) { inNeighbors[2 * j] = ngh; }
    void setOutNeighbor(uintT j, uintE ngh) { outNeighbors[2 * j] = ngh; }
    void setInWeight(uintT j, uintE wgh) { inNeighbors[2 * j + 1] = wgh; }
    void setOutWeight(uintT j, uintE wgh) { outNeighbors[2 * j + 1] = wgh; }
    void setInNeighbors(intE* _i) { inNeighbors = _i; }
    void setOutNeighbors(intE* _i) { outNeighbors = _i; }
#endif

    uintT getInDegree() const { return inDegree; }
    uintT getOutDegree() const { return outDegree; }
    void setInDegree(uintT _d) { inDegree = _d; }
    void setOutDegree(uintT _d) { outDegree = _d; }
    void flipEdges() {
        swap(inNeighbors, outNeighbors);
        swap(inDegree, outDegree);
    }

    template <class VS, class F, class G>
    inline void decodeInNghBreakEarly(long v_id,
                                      VS& vertexSubset,
                                      F& f,
                                      G& g,
                                      bool parallel = 0) {
        decode_uncompressed::decodeInNghBreakEarly<asymmetricVertex, F, G, VS>(
            this, v_id, vertexSubset, f, g, parallel);
    }

    template <class F, class G>
    inline void decodeOutNgh(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh<asymmetricVertex, F, G>(this, i, f,
                                                                  g);
    }

    template <class F, class G>
    inline void decodeOutNgh_single(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh_single<asymmetricVertex, F, G>(this, i, f,
                                                                  g);
    }

    template <class F, class G>
    inline void decodeOutNghSparse(long i, uintT o, F& f, G& g) {
        decode_uncompressed::decodeOutNghSparse<asymmetricVertex, F>(this, i, o,
                                                                     f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i, uintT o, F& f, G& g) {
        return decode_uncompressed::decodeOutNghSparseSeq<asymmetricVertex, F>(
            this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i,
                                        uintT o,
                                        F& f,
                                        G& g,
                                        uint64_t base) {
        return decode_uncompressed::decodeOutNghSparseSeq<asymmetricVertex, F>(
            this, i, o, f, g, base);
    }

    template <class E, class F, class G>
    inline void copyOutNgh(long i, uintT o, F& f, G& g) {
        decode_uncompressed::copyOutNgh<asymmetricVertex, E>(this, i, o, f, g);
    }

    template <class F>
    inline size_t countOutNgh(long i, F& f) {
        return decode_uncompressed::countOutNgh<asymmetricVertex, F>(this, i,
                                                                     f);
    }

    template <class F>
    inline size_t packOutNgh(long i,
                             F& f,
                             bool* bits,
                             uintE* tmp1,
                             uintE* tmp2) {
        return decode_uncompressed::packOutNgh<asymmetricVertex, F>(this, i, f,
                                                                    bits, tmp1);
    }
};

struct symmetricSubgraphVertex {
    uintE subgraphNum;
#ifndef WEIGHTED
    uintE** neighbors;
#else
    intE** neighbors;
#endif
    uintT degree;
    uintT* degrees;  // length of vertex neighbor of each subgraph
    void del() {
        free(neighbors);
        free(degrees);
    }
#ifndef WEIGHTED
    symmetricSubgraphVertex(uintE sn, bool is, uintE** n, uintE* ns, uintT d, uintT* D)
#else
    symmetricSubgraphVertex(uintE sn, bool is, intE** n, intE* ns, uintT d, uintT* D)
#endif
        : subgraphNum(sn), /*isSmall(is),*/ neighbors(n), /*neighborsSmall(ns),*/ degree(d), degrees(D) {
    }
#ifndef WEIGHTED
    uintE* getInNeighbors() {
        uintE* res = newA(uintE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    const uintE* getInNeighbors() const {
        uintE* res = newA(uintE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    uintE* getOutNeighbors() {
        uintE* res = newA(uintE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    const uintE* getOutNeighbors() const {
        uintE* res = newA(uintE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j) const {
        uintE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                res = neighbors[i][j];
            } else {
                j -= degrees[i];
            }
        }
        return res;
    }
    uintE getOutNeighbor(uintT j) const {
        uintE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                res = neighbors[i][j];
            } else {
                j -= degrees[i];
            }
        }
        return res;
    }

    void setInNeighbor(uintT j, uintE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                neighbors[i][j] = ngh;
                return;
            } else {
                j -= degrees[i];
            }
        }
    }
    void setOutNeighbor(uintT j, uintE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                neighbors[i][j] = ngh;
                return;
            } else {
                j -= degrees[i];
            }
        }
    }
    void setInNeighbors(uintE** _i) { neighbors = _i; }
    void setOutNeighbors(uintE** _i) { neighbors = _i; }
#else
    // weights are stored in the entry after the neighbor ID
    // so size of neighbor list is twice the degree
    intE* getInNeighbors() {
        intE* res = newA(intE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    const intE* getInNeighbors() const {
        intE* res = newA(intE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    intE* getOutNeighbors() {
        intE* res = newA(intE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    const intE* getOutNeighbors() const {
        intE* res = newA(intE, degree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < degrees[i]; j++, k++) {
                res[k] = neighbors[i][j];
            }
        }
        return res;
    }
    intE getInNeighbor(intT j) const {
        intE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                res = neighbors[i][2 * j];
            } else {
                j -= degrees[i];
            }
        }
        return res;
    }
    intE getOutNeighbor(intT j) const {
        intE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                res = neighbors[i][2 * j];
            } else {
                j -= degrees[i];
            }
        }
        return res;
    }
    intE getInWeight(intT j) const {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                return neighbors[i][2 * j + 1];
            } else {
                j -= degrees[i];
            }
        }
    }
    intE getOutWeight(intT j) const {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                return neighbors[i][2 * j + 1];
            } else {
                j -= degrees[i];
            }
        }
    }
    void setInNeighbor(uintT j, uintE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                neighbors[i][2 * j] = ngh;
                return;
            } else {
                j -= degrees[i];
            }
        }
    }
    void setOutNeighbor(uintT j, uintE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                neighbors[i][2 * j] = ngh;
                return;
            } else {
                j -= degrees[i];
            }
        }
    }
    void setInWeight(uintT j, intE wgh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                neighbors[i][2 * j + 1] = wgh;
                return;
            } else {
                j -= degrees[i];
            }
        }
    }
    void setOutWeight(uintT j, intE wgh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < degrees[i]) {
                neighbors[i][2 * j + 1] = wgh;
                return;
            } else {
                j -= degrees[i];
            }
        }
    }
    void setInNeighbors(intE** _i) { neighbors = _i; }
    void setOutNeighbors(intE** _i) { neighbors = _i; }
#endif

    uintT getInDegree() const { return degree; }
    uintT getOutDegree() const { return degree; }
    uintT* getInDegrees() const { return degrees; }
    uintT* getOutDegrees() const { return degrees; }
    void setInDegree(uintT _d) { degree = _d; }
    void setOutDegree(uintT _d) { degree = _d; }
    void setInDegrees(uintT* _d) { degrees = _d; }
    void setOutDegrees(uintT* _d) { degrees = _d; }
    void flipEdges() {}

    template <class VS, class F, class G>
    inline void decodeInNghBreakEarly(long v_id,
                                      VS& vertexSubset,
                                      F& f,
                                      G& g,
                                      bool parallel = 0) {
        decode_uncompressed::decodeInNghBreakEarly<symmetricSubgraphVertex, F,
                                                   G, VS>(
            this, v_id, vertexSubset, f, g, parallel);
    }

    template <class F, class G>
    inline void decodeOutNgh(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh<symmetricSubgraphVertex, F, G>(
            this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNgh_single(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh_single<symmetricSubgraphVertex, F, G>(
            this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNghSparse(long i, uintT o, F& f, G& g) {
        decode_uncompressed::decodeOutNghSparse<symmetricSubgraphVertex, F>(
            this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i, uintT o, F& f, G& g) {
        return decode_uncompressed::decodeOutNghSparseSeq<
            symmetricSubgraphVertex, F>(this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i,
                                        uintT o,
                                        F& f,
                                        G& g,
                                        uint64_t base) {
        return decode_uncompressed::decodeOutNghSparseSeq<
            symmetricSubgraphVertex, F>(this, i, o, f, g, base);
    }

    template <class E, class F, class G>
    inline void copyOutNgh(long i, uintT o, F& f, G& g) {
        decode_uncompressed::copyOutNgh<symmetricSubgraphVertex, E>(this, i, o,
                                                                    f, g);
    }

    template <class F>
    inline size_t countOutNgh(long i, F& f) {
        return decode_uncompressed::countOutNgh<symmetricSubgraphVertex, F>(
            this, i, f);
    }

    template <class F>
    inline size_t packOutNgh(long i,
                             F& f,
                             bool* bits,
                             uintE* tmp1,
                             uintE* tmp2) {
        return decode_uncompressed::packOutNgh<symmetricSubgraphVertex, F>(
            this, i, f, bits, tmp1);
    }
};

struct asymmetricSubgraphVertex {
    uintE subgraphNum;
#ifndef WEIGHTED
    uintE **inNeighbors, **outNeighbors;
#else
    intE **inNeighbors, **outNeighbors;
#endif
    uintT outDegree;
    uintT inDegree;
    uintT* outDegrees;  // length of vertex outNeighbor of each subgraph
    uintT* inDegrees;   // length of vertex inNeighbor of each subgraph
    uintT *outOffsets, *inOffsets;

    void del() {
        free(inNeighbors);
        free(outNeighbors);
        free(inDegrees);
        free(outDegrees);
        free(outOffsets);
        free(inOffsets);
    }

#ifndef WEIGHTED
    asymmetricSubgraphVertex(uintE num,
                             uintE** iN,
                             uintE** oN,
                             uintT id,
                             uintT od,
                             uintT* iD,
                             uintT* oD,
                             uintT* iO,
                             uintT* oO)
#else
    asymmetricSubgraphVertex(uintE num,
                             intE** iN,
                             intE** oN,
                             uintT id,
                             uintT od,
                             uintT* iD,
                             uintT* oD,
                             uintT* iO,
                             uintT* oO)
#endif
        : subgraphNum(num),
          inNeighbors(iN),
          outNeighbors(oN),
          inDegree(id),
          outDegree(od),
          inDegrees(iD),
          outDegrees(oD),
          inOffsets(iO),
          outOffsets(oO)
          {}

#ifndef WEIGHTED
    uintE* getInNeighbors() {
        uintE* res = newA(uintE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    const uintE* getInNeighbors() const {
        uintE* res = newA(uintE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    uintE* getOutNeighbors() {
        uintE* res = newA(uintE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    const uintE* getOutNeighbors() const {
        uintE* res = newA(uintE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j) const {
        uintE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                res = inNeighbors[i][j];
                break;
            } else {
                j -= inDegrees[i];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j, uintT &gid) const {
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = inNeighbors[gid][j - inOffsets[gid]];
                break;
            }
        }
        return res;
    }
    uintE getOutNeighbor(uintT j) const {
        uintE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                res = outNeighbors[i][j];
                break;
            } else {
                j -= outDegrees[i];
            }
        }
        return res;
    }
    uintE getOutNeighbor(uintT j, uintT &gid) const {
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = outNeighbors[gid][j - outOffsets[gid]];
                break;
            }
        }
        return res;
    }
    void setInNeighbor(uintT j, uintE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                inNeighbors[i][j] = ngh;
                return;
            } else {
                j -= inDegrees[i];
            }
        }
    }
    void setInNeighbor(uintT j, uintE ngh, uintT &gid) {
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                inNeighbors[gid][j - inOffsets[gid]] = ngh;
                break;
            }
        }
    }
    void setOutNeighbor(uintT j, uintE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                outNeighbors[i][j] = ngh;
                return;
            } else {
                j -= outDegrees[i];
            }
        }
    }
    void setOutNeighbor(uintT j, uintE ngh, uintT &gid) {
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                outNeighbors[gid][j - outOffsets[gid]] = ngh;
                break;
            }
        }
    }
    void setInNeighbors(uintE** _i) { inNeighbors = _i; }
    void setOutNeighbors(uintE** _i) { outNeighbors = _i; }
#else
    intE* getInNeighbors() {
        intE* res = newA(intE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    const intE* getInNeighbors() const {
        intE* res = newA(intE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    intE* getOutNeighbors() {
        intE* res = newA(intE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    const intE* getOutNeighbors() const {
        intE* res = newA(intE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    intE getInNeighbor(uintT j) const {
        intE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                res = inNeighbors[i][2 * j];
                break;
            } else {
                j -= inDegrees[i];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j, uintT &gid) const {
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = inNeighbors[gid][(j - inOffsets[gid])*2];
                break;
            }
        }
        return res;
    }
    intE getOutNeighbor(uintT j) const {
        intE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                res = outNeighbors[i][2 * j];
                break;
            } else {
                j -= outDegrees[i];
            }
        }

        return res;
    }
    uintE getOutNeighbor(uintT j, uintT &gid) const {
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = outNeighbors[gid][(j - outOffsets[gid])*2];
                break;
            }
        }
        return res;
    }

    intE getInWeight(uintT j) const {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                return inNeighbors[i][2 * j + 1];
            } else {
                j -= inDegrees[i];
            }
        }
        return INT_E_MAX;
    }
    intE getInWeight(uintT j, uintT &gid) const {
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = inNeighbors[gid][(j - inOffsets[gid])*2 + 1];
                break;
            }
        }
        return res;
    }

    intE getOutWeight(uintT j) const {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                return outNeighbors[i][2 * j + 1];
            } else {
                j -= outDegrees[i];
            }
        }
        return INT_E_MAX;
    }

    intE getOutWeight(uintT j, uintT &gid) const {
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = outNeighbors[gid][(j - outOffsets[gid])*2 + 1];
                break;
            }
        }
        return res;
    }

    void setInNeighbor(uintT j, intE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                inNeighbors[i][2 * j] = ngh;
                return;
            } else {
                j -= inDegrees[i];
            }
        }
    }
    void setInNeighbor(uintT j, intE ngh, uintT &gid) const {
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                inNeighbors[gid][(j - inOffsets[gid])*2] = ngh;
                break;
            }
        }
    }
    void setOutNeighbor(uintT j, intE ngh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                outNeighbors[i][2 * j] = ngh;
                return;
            } else {
                j -= outDegrees[i];
            }
        }
    }

    void setOutNeighbor(uintT j, intE ngh, uintT &gid) const {
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                outNeighbors[gid][(j - outOffsets[gid])*2] = ngh;
                break;
            }
        }
    }

    void setInWeight(uintT j, uintE wgh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                inNeighbors[i][2 * j + 1] = wgh;
                return;
            } else {
                j -= inDegrees[i];
            }
        }
    }

    void setInWeight(uintT j, uintE wgh, uintT &gid) const {
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                inNeighbors[gid][(j - inOffsets[gid])*2 + 1] = wgh;
                break;
            }
        }
    }

    void setOutWeight(uintT j, uintE wgh) {
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                outNeighbors[i][2 * j + 1] = wgh;
                return;
            } else {
                j -= outDegrees[i];
            }
        }
    }

    void setOutWeight(uintT j, uintE wgh, uintT &gid) const {
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                outNeighbors[gid][(j - outOffsets[gid])*2 + 1] = wgh;
                break;
            }
        }
    }

    void setInNeighbors(intE** _i) { inNeighbors = _i; }
    void setOutNeighbors(intE** _i) { outNeighbors = _i; }
#endif

    uintT getInDegree() const { return inDegree; }
    uintT getOutDegree() const { return outDegree; }
    uintT* getInDegrees() const { return inDegrees; }
    uintT* getOutDegrees() const { return outDegrees; }
    void setInDegree(uintT _d) { inDegree = _d; }
    void setOutDegree(uintT _d) { outDegree = _d; }
    void setInDegrees(uintT* _d) { inDegrees = _d; }
    void setOutDegrees(uintT* _d) { outDegrees = _d; }

    void flipEdges() {
        swap(inNeighbors, outNeighbors);
        swap(inDegree, outDegree);
        swap(inDegrees, outDegrees);
        swap(inOffsets, outOffsets);
    }

    template <class VS, class F, class G>
    inline void decodeInNghBreakEarly(long v_id,
                                      VS& vertexSubset,
                                      F& f,
                                      G& g,
                                      bool parallel = 0) {
        decode_uncompressed::decodeInNghBreakEarly<asymmetricSubgraphVertex, F,
                                                   G, VS>(
            this, v_id, vertexSubset, f, g, parallel);
    }

    template <class F, class G>
    inline void decodeOutNgh(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh<asymmetricSubgraphVertex, F, G>(
            this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNgh_single(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh_single<symmetricSubgraphVertex, F, G>(
            this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNghSparse(long i, uintT o, F& f, G& g) {
        decode_uncompressed::decodeOutNghSparse<asymmetricSubgraphVertex, F>(
            this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i, uintT o, F& f, G& g) {
        return decode_uncompressed::decodeOutNghSparseSeq<
            asymmetricSubgraphVertex, F>(this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i,
                                        uintT o,
                                        F& f,
                                        G& g,
                                        uint64_t base) {
        return decode_uncompressed::decodeOutNghSparseSeq<
            asymmetricSubgraphVertex, F>(this, i, o, f, g, base);
    }

    template <class E, class F, class G>
    inline void copyOutNgh(long i, uintT o, F& f, G& g) {
        decode_uncompressed::copyOutNgh<asymmetricSubgraphVertex, E>(this, i, o,
                                                                     f, g);
    }

    template <class F>
    inline size_t countOutNgh(long i, F& f) {
        return decode_uncompressed::countOutNgh<asymmetricSubgraphVertex, F>(
            this, i, f);
    }

    template <class F>
    inline size_t packOutNgh(long i,
                             F& f,
                             bool* bits,
                             uintE* tmp1,
                             uintE* tmp2) {
        return decode_uncompressed::packOutNgh<asymmetricSubgraphVertex, F>(
            this, i, f, bits, tmp1);
    }
};

struct asymmetricSubgraphVertexWithSmall {
    uintE subgraphNum;
#ifndef WEIGHTED
    uintE **inNeighbors, **outNeighbors;
#else
    intE **inNeighbors, **outNeighbors;
#endif
    uintT outDegree;
    uintT inDegree;
    uintT* outDegrees;  // length of vertex outNeighbor of each subgraph
    uintT* inDegrees;   // length of vertex inNeighbor of each subgraph
    uintT *outOffsets, *inOffsets;
    bool isInSmall;
    bool isOutSmall;
#ifndef WEIGHTED
    uintE *inNeighborsSmall, *outNeighborsSmall;
#else
    intE *inNeighborsSmall, *outNeighborsSmall;
#endif

    void del() {
        free(inNeighbors);
        free(outNeighbors);
        free(inDegrees);
        free(outDegrees);
        free(outOffsets);
        free(inOffsets);
        if (isInSmall) {
            free(inNeighborsSmall);
        }
        if (isOutSmall) {
            free(outNeighborsSmall);
        }
    }

#ifndef WEIGHTED
    asymmetricSubgraphVertexWithSmall(uintE num,
                             uintE** iN,
                             uintE** oN,
                             uintT id,
                             uintT od,
                             uintT* iD,
                             uintT* oD,
                             uintT* iO,
                             uintT* oO,
                             bool iIS,
                             bool iOS,
                             uintE* iNS,
                             uintE* oNS)
#else
    asymmetricSubgraphVertexWithSmall(uintE num,
                             intE** iN,
                             intE** oN,
                             uintT id,
                             uintT od,
                             uintT* iD,
                             uintT* oD,
                             uintT* iO,
                             uintT* oO,
                             bool iIS,
                             bool iOS,
                             intE* iNS,
                             intE* oNS)
#endif
        : subgraphNum(num),
          inNeighbors(iN),
          outNeighbors(oN),
          inDegree(id),
          outDegree(od),
          inDegrees(iD),
          outDegrees(oD),
          inOffsets(iO),
          outOffsets(oO),
          isInSmall(iIS),
          isOutSmall(iOS),
          inNeighborsSmall(iNS),
          outNeighborsSmall(oNS)
          {}

#ifndef WEIGHTED
    uintE* getInNeighbors() {
        if(isInSmall) return inNeighborsSmall;
        uintE* res = newA(uintE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    const uintE* getInNeighbors() const {
        if(isInSmall) return inNeighborsSmall;
        uintE* res = newA(uintE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    uintE* getOutNeighbors() {
        if(isOutSmall) return outNeighborsSmall;
        uintE* res = newA(uintE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    const uintE* getOutNeighbors() const {
        if(isOutSmall) return outNeighborsSmall;
        uintE* res = newA(uintE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j) const {
        if(isInSmall) return inNeighborsSmall[j];
        uintE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                res = inNeighbors[i][j];
                break;
            } else {
                j -= inDegrees[i];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j, uintT &gid) const {
        if (isInSmall) return inNeighborsSmall[j];
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = inNeighbors[gid][j - inOffsets[gid]];
                break;
            }
        }
        return res;
    }
    uintE getOutNeighbor(uintT j) const {
        if (isOutSmall) return outNeighborsSmall[j];
        uintE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                res = outNeighbors[i][j];
                break;
            } else {
                j -= outDegrees[i];
            }
        }
        return res;
    }
    uintE getOutNeighbor(uintT j, uintT &gid) const {
        if (isOutSmall) return outNeighborsSmall[j];
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = outNeighbors[gid][j - outOffsets[gid]];
                break;
            }
        }
        return res;
    }
    void setInNeighbor(uintT j, uintE ngh) {
        if (isInSmall) {
            inNeighborsSmall[j] = ngh;
            return ;
        }
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                inNeighbors[i][j] = ngh;
                return;
            } else {
                j -= inDegrees[i];
            }
        }
    }
    void setInNeighbor(uintT j, uintE ngh, uintT &gid) {
        if (isInSmall) {
            inNeighborsSmall[j] = ngh;
            return ;
        }
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                inNeighbors[gid][j - inOffsets[gid]] = ngh;
                break;
            }
        }
    }
    void setOutNeighbor(uintT j, uintE ngh) {
        // naive implementation
        if (isOutSmall) {
            outNeighborsSmall[j] = ngh;
            return ;
        }
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                outNeighbors[i][j] = ngh;
                return;
            } else {
                j -= outDegrees[i];
            }
        }
    }
    void setOutNeighbor(uintT j, uintE ngh, uintT &gid) {
        if (isOutSmall) {
            outNeighborsSmall[j] = ngh;
            return ;
        }
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                outNeighbors[gid][j - outOffsets[gid]] = ngh;
                break;
            }
        }
    }
    void setInNeighbors(uintE** _i) { inNeighbors = _i; }
    void setOutNeighbors(uintE** _i) { outNeighbors = _i; }
#else
    intE* getInNeighbors() {
        if (isInSmall) return inNeighborsSmall;
        intE* res = newA(intE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    const intE* getInNeighbors() const {
        if (isInSmall) return inNeighborsSmall;
        intE* res = newA(intE, inDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < inDegrees[i]; j++, k++) {
                res[k] = inNeighbors[i][j];
            }
        }
        return res;
    }
    intE* getOutNeighbors() {
        if (isOutSmall) return outNeighborsSmall;
        intE* res = newA(intE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    const intE* getOutNeighbors() const {
        if (isOutSmall) return outNeighborsSmall;
        intE* res = newA(intE, outDegree);
        intE k = 0;
        for (int i = 0; i < subgraphNum; i++) {
            for (int j = 0; j < outDegrees[i]; j++, k++) {
                res[k] = outNeighbors[i][j];
            }
        }
        return res;
    }
    intE getInNeighbor(uintT j) const {
        if (isInSmall) return inNeighborsSmall[2*j];
        intE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                res = inNeighbors[i][2 * j];
                break;
            } else {
                j -= inDegrees[i];
            }
        }
        return res;
    }
    uintE getInNeighbor(uintT j, uintT &gid) const {
        if (isInSmall) return inNeighborsSmall[2*j];
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = inNeighbors[gid][(j - inOffsets[gid])*2];
                break;
            }
        }
        return res;
    }
    intE getOutNeighbor(uintT j) const {
        if (isOutSmall) return outNeighborsSmall[2*j];
        intE res;
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                res = outNeighbors[i][2 * j];
                break;
            } else {
                j -= outDegrees[i];
            }
        }

        return res;
    }
    uintE getOutNeighbor(uintT j, uintT &gid) const {
        if (isOutSmall) return outNeighborsSmall[2*j];
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = outNeighbors[gid][(j - outOffsets[gid])*2];
                break;
            }
        }
        return res;
    }

    intE getInWeight(uintT j) const {
        if (isInSmall) return inNeighborsSmall[2*j + 1];
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                return inNeighbors[i][2 * j + 1];
            } else {
                j -= inDegrees[i];
            }
        }
        return INT_E_MAX;
    }
    intE getInWeight(uintT j, uintT &gid) const {
        if (isInSmall) return inNeighborsSmall[2*j + 1];
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = inNeighbors[gid][(j - inOffsets[gid])*2 + 1];
                break;
            }
        }
        return res;
    }

    intE getOutWeight(uintT j) const {
        if (isOutSmall) return outNeighborsSmall[2*j + 1];
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                return outNeighbors[i][2 * j + 1];
            } else {
                j -= outDegrees[i];
            }
        }
        return INT_E_MAX;
    }

    intE getOutWeight(uintT j, uintT &gid) const {
        if (isOutSmall) return outNeighborsSmall[2*j + 1];
        uintE res;
        // naive implementation
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                res = outNeighbors[gid][(j - outOffsets[gid])*2 + 1];
                break;
            }
        }
        return res;
    }

    void setInNeighbor(uintT j, intE ngh) {
        if (isInSmall) {
            inNeighborsSmall[2*j] = ngh;
            return ;
        }
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                inNeighbors[i][2 * j] = ngh;
                return;
            } else {
                j -= inDegrees[i];
            }
        }
    }
    void setInNeighbor(uintT j, intE ngh, uintT &gid) const {
        if (isInSmall) {
            inNeighborsSmall[2*j] = ngh;
            return ;
        }
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                inNeighbors[gid][(j - inOffsets[gid])*2] = ngh;
                break;
            }
        }
    }
    void setOutNeighbor(uintT j, intE ngh) {
        if (isOutSmall) {
            outNeighborsSmall[2*j] = ngh;
            return ;
        }
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                outNeighbors[i][2 * j] = ngh;
                return;
            } else {
                j -= outDegrees[i];
            }
        }
    }

    void setOutNeighbor(uintT j, intE ngh, uintT &gid) const {
        if (isOutSmall) {
            outNeighborsSmall[2*j] = ngh;
            return ;
        }
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                outNeighbors[gid][(j - outOffsets[gid])*2] = ngh;
                break;
            }
        }
    }

    void setInWeight(uintT j, uintE wgh) {
        if (isInSmall) {
            inNeighborsSmall[2*j + 1] = wgh;
            return ;
        }
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < inDegrees[i]) {
                inNeighbors[i][2 * j + 1] = wgh;
                return;
            } else {
                j -= inDegrees[i];
            }
        }
    }

    void setInWeight(uintT j, uintE wgh, uintT &gid) const {
        if (isInSmall) {
            inNeighborsSmall[2*j + 1] = wgh;
            return ;
        }
        for (; gid < subgraphNum; gid++) {
            if ((j >= inOffsets[gid] && j < inOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                inNeighbors[gid][(j - inOffsets[gid])*2 + 1] = wgh;
                break;
            }
        }
    }

    void setOutWeight(uintT j, uintE wgh) {
        if (isOutSmall) {
            outNeighborsSmall[2*j + 1] = wgh;
            return ;
        }
        // naive implementation
        for (int i = 0; i < subgraphNum; i++) {
            if (j < outDegrees[i]) {
                outNeighbors[i][2 * j + 1] = wgh;
                return;
            } else {
                j -= outDegrees[i];
            }
        }
    }

    void setOutWeight(uintT j, uintE wgh, uintT &gid) const {
        if (isOutSmall) {
            outNeighborsSmall[2*j + 1] = wgh;
            return ;
        }
        for (; gid < subgraphNum; gid++) {
            if ((j >= outOffsets[gid] && j < outOffsets[gid + 1]) || (gid == subgraphNum - 1)) {
                outNeighbors[gid][(j - outOffsets[gid])*2 + 1] = wgh;
                break;
            }
        }
    }

    void setInNeighbors(intE** _i) { inNeighbors = _i; }
    void setOutNeighbors(intE** _i) { outNeighbors = _i; }
#endif

    uintT getInDegree() const { return inDegree; }
    uintT getOutDegree() const { return outDegree; }
    uintT* getInDegrees() const { return inDegrees; }
    uintT* getOutDegrees() const { return outDegrees; }
    void setInDegree(uintT _d) { inDegree = _d; }
    void setOutDegree(uintT _d) { outDegree = _d; }
    void setInDegrees(uintT* _d) { inDegrees = _d; }
    void setOutDegrees(uintT* _d) { outDegrees = _d; }

    void flipEdges() {
        swap(inNeighbors, outNeighbors);
        swap(inDegree, outDegree);
        swap(inDegrees, outDegrees);
        swap(inOffsets, outOffsets);
        swap(isInSmall, isOutSmall);
        swap(inNeighborsSmall, outNeighborsSmall);
    }

    template <class VS, class F, class G>
    inline void decodeInNghBreakEarly(long v_id,
                                      VS& vertexSubset,
                                      F& f,
                                      G& g,
                                      bool parallel = 0) {
        decode_uncompressed::decodeInNghBreakEarly<asymmetricSubgraphVertexWithSmall, F,
                                                   G, VS>(
            this, v_id, vertexSubset, f, g, parallel);
    }

    template <class F, class G>
    inline void decodeOutNgh(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh<asymmetricSubgraphVertexWithSmall, F, G>(
            this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNgh_single(long i, F& f, G& g) {
        decode_uncompressed::decodeOutNgh_single<symmetricSubgraphVertex, F, G>(
            this, i, f, g);
    }

    template <class F, class G>
    inline void decodeOutNghSparse(long i, uintT o, F& f, G& g) {
        decode_uncompressed::decodeOutNghSparse<asymmetricSubgraphVertexWithSmall, F>(
            this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i, uintT o, F& f, G& g) {
        return decode_uncompressed::decodeOutNghSparseSeq<
            asymmetricSubgraphVertexWithSmall, F>(this, i, o, f, g);
    }

    template <class F, class G>
    inline size_t decodeOutNghSparseSeq(long i,
                                        uintT o,
                                        F& f,
                                        G& g,
                                        uint64_t base) {
        return decode_uncompressed::decodeOutNghSparseSeq<
            asymmetricSubgraphVertexWithSmall, F>(this, i, o, f, g, base);
    }

    template <class E, class F, class G>
    inline void copyOutNgh(long i, uintT o, F& f, G& g) {
        decode_uncompressed::copyOutNgh<asymmetricSubgraphVertexWithSmall, E>(this, i, o,
                                                                     f, g);
    }

    template <class F>
    inline size_t countOutNgh(long i, F& f) {
        return decode_uncompressed::countOutNgh<asymmetricSubgraphVertexWithSmall, F>(
            this, i, f);
    }

    template <class F>
    inline size_t packOutNgh(long i,
                             F& f,
                             bool* bits,
                             uintE* tmp1,
                             uintE* tmp2) {
        return decode_uncompressed::packOutNgh<asymmetricSubgraphVertexWithSmall, F>(
            this, i, f, bits, tmp1);
    }
};

#endif
