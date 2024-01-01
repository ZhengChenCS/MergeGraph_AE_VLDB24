#ifndef EDGEMAP_ONCE
#define EDGEMAP_ONCE
#include "IO.h"
#include "compressedVertex.h"
#include "edgeMap_utils.h"
#include "gettime.h"
#include "graph.h"
#include "index_map.h"
#include "parallel.h"
#include "parseCommandLine.h"
#include "utils.h"
#include "vertex.h"
#include "vertexSubset.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <memory>
#include <vector>
// #define EDGEMAP_DEBUG

typedef uint32_t flags;

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDense_once(std::vector<graph<vertex>> &GA, VS &vertexSubset, F &f,
                                    const flags fl) {
    static int rounds = 0;
    rounds++;
    #ifdef EDGEMAP_DEBUG
    if(rounds <= 10)
        std::cout << "dense pull" << std::endl;
    #endif
    using D = tuple<bool, data>;
    long n = GA[0].n;
    long numGraphs = GA.size();
    std::vector<vertex*> G(numGraphs);
    for (int gid = 0; gid < numGraphs; gid++)
        G[gid] = GA[gid].V;
    if (should_output(fl)) {
        D *next = newA(D, n);
        auto g = get_emdense_gen<data>(next);
        parallel_for(long v = 0; v < n; v++) {
            std::get<0>(next[v]) = 0;
            if (f.cond(v)) {
                for (int gid = 0; gid < numGraphs; gid++)
                    G[gid][v].decodeInNghBreakEarly(v, vertexSubset, f, g,
                                            fl & dense_parallel);
            }
        }
        return vertexSubsetData<data>(n, next);
    } else {
        auto g = get_emdense_nooutput_gen<data>();
        parallel_for(long v = 0; v < n; v++) {
            if (f.cond(v)) {
                for (int gid = 0; gid < numGraphs; gid++)
                    G[gid][v].decodeInNghBreakEarly(v, vertexSubset, f, g,
                                           fl & dense_parallel);
            }
        }
        return vertexSubsetData<data>(n);
    }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDenseForward_once(std::vector<graph<vertex>> &GA, VS &vertexSubset,
                                           F &f, const flags fl) {
    using D = tuple<bool, data>;
    static int rounds = 0;
    rounds++;
    #ifdef EDGEMAP_DEBUG
    if(rounds <= 10)
        std::cout << "dense forward" << std::endl;
    #endif
    long n = GA[0].n;
    long numGraphs = GA.size();
    std::vector<vertex*> G(numGraphs);
    for(int gid= 0; gid < numGraphs; gid++)
        G[gid] = GA[gid].V;
    if (should_output(fl)) {
        D *next = newA(D, n);
        auto g = get_emdense_forward_gen<data>(next);
        parallel_for(long i = 0; i < n; i++) { std::get<0>(next[i]) = 0; }
        parallel_for(long i = 0; i < n; i++) {
            if (vertexSubset.isIn(i)) {
                for(int gid = 0; gid < numGraphs; gid++)
                {
                    G[gid][i].decodeOutNgh_single(i, f, g);
                }
            }
        }
        // parallel_for(int gid = 0; gid < numGraphs; gid++)
        // {
        //     parallel_for(long i = 0; i < n; i++)
        //     {
        //         if (vertexSubset.isIn(i)) {
        //             G[gid][i].decodeOutNgh_single(i, f, g);
        //         }
        //     }
        // }
        return vertexSubsetData<data>(n, next);
    } else {
        auto g = get_emdense_forward_nooutput_gen<data>();
        parallel_for(long i = 0; i < n; i++) {
            if (vertexSubset.isIn(i)) {
                for(int gid = 0; gid < numGraphs; gid++)
                {
                    G[gid][i].decodeOutNgh_single(i, f, g);
                }
            }
        }
        return vertexSubsetData<data>(n);
    }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data>
edgeMapSparse_once(std::vector<graph<vertex>> &GA, std::vector<vertex *> frontierVertices, VS &indices,
              std::vector<uintT*> degrees, uintT m, F &f, const flags fl) {
    using S = tuple<uintE, data>;
    static int rounds = 0;
    rounds++;
    #ifdef EDGEMAP_DEBUG
    if(rounds <= 10)
        std::cout << "sparse" << std::endl;
    #endif
    long n = indices.n;
    long numGraphs = GA.size();
    S * outEdges;
    long outEdgeCount = 0;
    std::vector<long> temp_outEdgeCount(numGraphs);
    std::vector<S*> temp_outEdges(numGraphs);

    if (should_output(fl)) {
        std::vector<uintT *> offsets = degrees;

        for(int i = 0; i < numGraphs; i++)
        {
            temp_outEdgeCount[i] = sequence::plusScan(offsets[i], offsets[i], m);
            temp_outEdges[i] = newA(S, temp_outEdgeCount[i]);
        }
        parallel_for(size_t i = 0; i < m; i++) {
            for(int gid = 0; gid < numGraphs; gid++)
            {
                auto g = get_emsparse_gen<data>(temp_outEdges[gid]);
                uintT v = indices.vtx(i), o = offsets[gid][i];
                vertex vert = frontierVertices[gid][i];
                vert.decodeOutNghSparse(v, o, f, g);
            }
        }
    } else {
        auto g = get_emsparse_nooutput_gen<data>();
        parallel_for(size_t i = 0; i < m; i++) {
            for(int gid = 0; gid < numGraphs; gid++)
            {
                uintT v = indices.vtx(i);
                vertex vert = frontierVertices[gid][i];
                vert.decodeOutNghSparse(v, 0, f, g);
            }
        }
    }

    if (should_output(fl)) {
        for (int i = 0; i < numGraphs; i++)
        {
            outEdgeCount += temp_outEdgeCount[i];
        }
        outEdges = newA(S, outEdgeCount);
        long base = 0;
        for (int gid = 0; gid < numGraphs; gid++)
        {
            parallel_for(int i = 0; i < temp_outEdgeCount[gid]; i++)
            {
                outEdges[base + i] = temp_outEdges[gid][i];
            }
            base += temp_outEdgeCount[gid];
        }
        for(int i = 0; i < numGraphs; i++)
        {
            free(temp_outEdges[i]);
        }
        S *nextIndices = newA(S, outEdgeCount);
        if (fl & remove_duplicates) {
            if (GA[0].flags == NULL) {
                GA[0].flags = newA(uintE, n);
                parallel_for(long i = 0; i < n; i++) {
                    GA[0].flags[i] = UINT_E_MAX;
                }
            }
            auto get_key = [&](size_t i) -> uintE & {
                return std::get<0>(outEdges[i]);
            };
            remDuplicates(get_key, GA[0].flags, outEdgeCount, n);
        }
        auto p = [](tuple<uintE, data> &v) {
            return std::get<0>(v) != UINT_E_MAX;
        };
        size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
        free(outEdges);
        return vertexSubsetData<data>(n, nextM, nextIndices);
    } else {
        return vertexSubsetData<data>(n);
    }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter_once(std::vector<graph<vertex>> &GA,
                                               std::vector<vertex *> frontierVertices,
                                               VS &indices, std::vector<uintT *> offsets,
                                               uintT m, F &f, const flags fl) {
    static int rounds = 0;
    rounds++;
    #ifdef EDGEMAP_DEBUG
    if(rounds <= 10)
        std::cout << "sparse no filter" << std::endl;
    #endif
    using S = tuple<uintE, data>;
    long n = indices.n;
    long outEdgeCount = 0;
    long numGraphs = GA.size();
    for(int i = 0; i < numGraphs; i++)
        outEdgeCount += sequence::plusScan(offsets[i], offsets[i], m);
    S *outEdges = newA(S, outEdgeCount);

    auto g = get_emsparse_no_filter_gen<data>(outEdges);

    // binary-search into scan to map workers->chunks
    size_t b_size = 10000;
    size_t n_blocks = nblocks(outEdgeCount, b_size);

    uintE *cts = newA(uintE, n_blocks + 1);
    size_t *block_offs = newA(size_t, n_blocks + 1);

    auto complete_off = newA(uintT, n);
    for (int gid = 0; gid < numGraphs; gid++)
    {
        parallel_for (int v = 0; v < n; v++)
        {
            complete_off[v] += offsets[gid][v];
        }
    }
    auto offsets_m =
        make_in_imap<uintT>(m, [&](size_t i) { return complete_off[i]; });
    auto lt = [](const uintT &l, const uintT &r) { return l < r; };
    parallel_for(size_t i = 0; i < n_blocks; i++) {
        size_t s_val = i * b_size;
        block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
    }
    block_offs[n_blocks] = m;
    parallel_for(size_t i = 0; i < n_blocks; i++) {
        if ((i == n_blocks - 1) || block_offs[i] != block_offs[i + 1]) {
            // start and end are offsets in [m]
            size_t start = block_offs[i];
            size_t end = block_offs[i + 1];
            uintT start_o = complete_off[start];
            uintT k = start_o;
            for (size_t j = start; j < end; j++) {
                uintE v = indices.vtx(j);
                size_t num_in = 0;
                uint64_t base = 0;
                for (int gid = 0; gid < numGraphs; gid++)
                {
                    num_in += frontierVertices[gid][j].decodeOutNghSparseSeq(v, k, f, g, base);
                    base += offsets[gid][j];
                }
                k += num_in;
            }
            cts[i] = (k - start_o);
        } else {
            cts[i] = 0;
        }
    }

    long outSize = sequence::plusScan(cts, cts, n_blocks);
    cts[n_blocks] = outSize;

    S *out = newA(S, outSize);

    parallel_for(size_t i = 0; i < n_blocks; i++) {
        if ((i == n_blocks - 1) || block_offs[i] != block_offs[i + 1]) {
            size_t start = block_offs[i];
            size_t start_o = complete_off[start];
            size_t out_off = cts[i];
            size_t block_size = cts[i + 1] - out_off;
            for (size_t j = 0; j < block_size; j++) {
                out[out_off + j] = outEdges[start_o + j];
            }
        }
    }
    free(outEdges);
    free(cts);
    free(block_offs);

    if (fl & remove_duplicates) {
        if (GA[0].flags == NULL) {
            GA[0].flags = newA(uintE, n);
            parallel_for(size_t i = 0; i < n; i++) { GA[0].flags[i] = UINT_E_MAX; }
        }
        auto get_key = [&](size_t i) -> uintE & { return std::get<0>(out[i]); };
        remDuplicates(get_key, GA[0].flags, outSize, n);
        S *nextIndices = newA(S, outSize);
        auto p = [](tuple<uintE, data> &v) {
            return std::get<0>(v) != UINT_E_MAX;
        };
        size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
        free(out);
        return vertexSubsetData<data>(n, nextM, nextIndices);
    }
    return vertexSubsetData<data>(n, outSize, out);
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapData_once(std::vector<graph<vertex>> &GA, VS &vs, F f,
                                   intT threshold = -1, const flags &fl = 0) {
    long numVertices = GA[0].n, numEdges = 0, m = vs.numNonzeros(), numGraphs = GA.size();
    for (int i = 0; i < numGraphs; i++)
        numEdges += GA[i].m;
    if (threshold == -1)
        threshold = numEdges / 20; // default threshold
    if (numVertices != vs.numRows()) {
        cout << "edgeMap: Sizes Don't match" << endl;
        abort();
    }
    if (m == 0)
        return vertexSubsetData<data>(numVertices);
    std::vector<uintT *>degrees(numGraphs);
    std::vector<vertex *> frontierVertices(numGraphs);
    static int rounds = 0;
    rounds++;
    uintT outDegrees = 0;
    if (threshold > 0) { // compute sum of out-degrees if threshold > 0
        vs.toSparse();
        for(int i = 0; i < numGraphs; i++)
        {
            degrees[i] = newA(uintT, m);
            frontierVertices[i] = newA(vertex, m);
        }
        {
            std::vector<vertex *> G(numGraphs);
            for (int i = 0; i < numGraphs; i++)
                G[i] = GA[i].V;
            parallel_for(size_t i = 0; i < m; i++) {
                uintE v_id = vs.vtx(i);
                // if (v_id > numVertices) cout << "error " << v_id << endl;
                for (int gid = 0; gid < numGraphs; gid++)
                {
                    vertex v = G[gid][v_id];
                    degrees[gid][i] = v.getOutDegree();
                    frontierVertices[gid][i] = v;
                }
            }
        }
        uintT* tt_degrees = newA(uintT, m);
        parallel_for(int i = 0; i < m; i++)
        {
            tt_degrees[i] = 0;
            for (int j = 0; j < numGraphs; j++)
                tt_degrees[i] += degrees[j][i];
        }
        outDegrees = sequence::plusReduce(tt_degrees, m);
        if (outDegrees == 0)
            return vertexSubsetData<data>(numVertices);
    }
    #ifdef EDGEMAP_DEBUG
    if(rounds <= 10)
    cout << "m=" << m << " outDegree=" << outDegrees << endl;
    #endif
    #ifndef ONLY_DENSE
    if (!(fl & no_dense) && m + outDegrees > threshold) {
        vs.toDense();
        return (fl & dense_forward)
                   ? edgeMapDenseForward_once<data, vertex, VS, F>(GA, vs, f, fl)
                   : edgeMapDense_once<data, vertex, VS, F>(GA, vs, f, fl);
    } else {
        auto vs_out =
            (should_output(fl) && fl & sparse_no_filter)
                ? // only call snof when we output
                edgeMapSparse_no_filter_once<data, vertex, VS, F>(
                    GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl)
                : edgeMapSparse_once<data, vertex, VS, F>(GA, frontierVertices, vs,
                                                     degrees, vs.numNonzeros(),
                                                     f, fl);
        return vs_out;
    }
    #else
    return edgeMapDenseForward_once<data, vertex, VS, F>(GA, vs, f, fl);
    #endif
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class vertex, class VS, class F>
vertexSubset edgeMap_once(std::vector<graph<vertex>> &GA, VS &vs, F f, intT threshold = -1,
                     const flags &fl = 0) {
    return edgeMapData_once<pbbs::empty>(GA, vs, f, threshold, fl);
}

#endif