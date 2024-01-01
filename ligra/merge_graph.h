#include <vector>
#include <chrono>
#include <iomanip>
#include "ligra.h"
#include "edgemap_once.h"

template <class Vertex>
graph<Vertex> test_merge_time(vector<graph<Vertex>> &GA)
{
    // auto all_start = std::chrono::high_resolution_clock::now();
    int numGraph = GA.size();
    long n = GA[0].n;
    long m = 0;
    for (int i = 0; i < GA.size(); i++) {
        m += GA[i].m;
    }
    Vertex *v = newA(Vertex, n);
    #ifndef WEIGHTED
    uintE *inNgh = newA(uintE, m);
    uintE *outNgh = newA(uintE, m);
    #else
    intE *inNgh = newA(intE, m*2);
    intE *outNgh = newA(intE, m*2);
    #endif
    uintE *in_offset = newA(uintE, n + 1);
    uintE *out_offset = newA(uintE, n + 1);
    in_offset[0] = 0;
    out_offset[0] = 0;
    for (int i = 0; i < n; i++)
    {
        uintT outDegree = 0;
        uintT inDegree = 0;
        for (int j = 0; j < numGraph; j++){
            outDegree += GA[j].V[i].getOutDegree();
            inDegree += GA[j].V[i].getInDegree();
        }
        in_offset[i + 1] = in_offset[i] + inDegree;
        out_offset[i + 1] = out_offset[i] + outDegree;
        v[i].setOutDegree(outDegree);
        v[i].setInDegree(inDegree);
    }
    parallel_for (int i = 0; i < n; i++){
        uintT l = 0;
        #ifndef WEIGHTED
        uintE *in = inNgh + in_offset[i], *out = outNgh + out_offset[i];
        #else
        intE *in = inNgh + in_offset[i]*2, *out = outNgh + out_offset[i]*2;
        #endif
        v[i].setInNeighbors(in);
        v[i].setOutNeighbors(out);
        for (int j = 0; j < numGraph; j++){
            uintE outDegree = GA[j].V[i].getOutDegree();
            #ifndef WEIGHTED
            memcpy(out + l, GA[j].V[i].getOutNeighbors(), outDegree * sizeof(uintE));
            l += outDegree;
            #else
            memcpy(out + l, GA[j].V[i].getOutNeighbors(), outDegree * 2 * sizeof(uintE));
            l += outDegree*2;
            #endif
            // for (int k = 0; k < GA[j].V[i].getOutDegree(); k++){
            //     out[l] = GA[j].V[i].getOutNeighbor(k);
            //     l++;
            // }
        }
        l = 0;
        for (int j = 0; j < numGraph; j++){
            uintE inDegree = GA[j].V[i].getInDegree();
            #ifndef WEIGHTED
            memcpy(in + l, GA[j].V[i].getInNeighbors(), inDegree * sizeof(uintE));
            l += inDegree;
            #else
            memcpy(in + l, GA[j].V[i].getInNeighbors(), inDegree * 2 * sizeof(uintE));
            l += inDegree*2;
            #endif
            // for (int k = 0; k < GA[j].V[i].getInDegree(); k++){
            //     inNgh[l] = GA[j].V[i].getInNeighbor(k);
            //     l++;
            // }
        }
        // v[i].setOutNeighbors(outNgh);
        // v[i].setInNeighbors(inNgh);
    }
    Uncompressed_Mem<Vertex> *mem;
    graph<Vertex> G(v, n, m, mem);
    // auto all_end = std::chrono::high_resolution_clock::now();
    // auto all_duration = std::chrono::duration<double>(all_end - all_start).count();
    // std::cout << std::fixed << std::setprecision(4);
    // std::cout << "merge graph time: " << all_duration << " seconds" << std::endl;
    return G;
}