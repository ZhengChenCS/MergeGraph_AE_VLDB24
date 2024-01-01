// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include <vector>
#include <chrono>
#include <iomanip>
// #define ONLY_DENSE
#define OPTIMIZE_GET_NEIGHBOR
#include "ligra.h"
#include "edgemap_once.h"
#include "merge.h"
#define INIT (uintE)0xffffff
#define InterSubgraph 1
#define EdgeLevel 2
#define MergeSmallDegree 3 
// #define MERGE_DEBUG

struct BFS_F {
    uintE* Results;
    BFS_F(uintE* _Results) : Results(_Results) {}
    inline bool update(uintE s, uintE d) {  // Update
        if (Results[d] == 0) {
            Results[d] = 1;
            return 1;
        } else
            return 0;
    }
    inline bool updateAtomic(uintE s, uintE d) {  // atomic version of Update
        return (CAS(&Results[d], (uintE)0, (uintE)1));
    }
    // cond function checks if vertex has been visited yet
    inline bool cond(uintE d) { return (cond_true(d)); }
};


struct LocalResult {
    uintE* value;
    LocalResult() {}
    LocalResult(long n, uintE* arr) {
        value = new uintE[n];
        for (int i = 0; i < n; i++)
            value[i] = *(arr + i);
    }
    // ~LocalResult(){
    // 	free(value);
    // }
};

void init(uintE* Results, vertexSubset& Frontier, vector<LocalResult> &LocalResults)
{
    long n = Frontier.n;
    uint64_t numGraph = LocalResults.size();
    parallel_for(long i = 0; i < n; i++) {
        Results[i] = 0;
        for (int j = 0; j < LocalResults.size(); j++) {
            if (LocalResults[j].value[i] < INIT)
            {
                Results[i] = 1;
                break;
            }
        }
    }
    bool *d = new bool[n];
    parallel_for (int i = 0; i < n; i++)
    {
        d[i] = false;
        if (Results[i] == 0) continue;
        for (int j = 0; j < numGraph; j++)
        {
            if (LocalResults[j].value[i] >= INIT)
            {
                d[i] = true;
                break;
            }
        }
    }
    Frontier = vertexSubset(n, d);
}

template <class origin_vertex, class merged_vertex>
void Compute(vector<graph<origin_vertex>>& G_vector, graph<merged_vertex> &G_merged, vector<LocalResult> &LocalResults, int parallel_method) 
{
    long n = G_vector[0].n;
    int numGraph = G_vector.size();
    // creates Results array, initialized to all -1, except for start
    uintE* Results = newA(uintE, n);
    uint64_t edge_num = 0;
    for (int i = 0; i < numGraph; i++)
        edge_num += G_vector[i].m;
    
    vertexSubset Frontier(n);  // creates initial frontier
    init(Results, Frontier, LocalResults);
    uint64_t round_count = 0;
    while (!Frontier.isEmpty() ||
           round_count == 0) {  // loop until frontier is empty
        round_count++;
        int active_cnt = 0, comp_cnt = 0;
        #ifdef MERGE_DEBUG
        if (round_count <= 10)
            cout << "round " << round_count << " : ";
        #endif
        vertexSubset output(n);
        #ifdef MERGE_DEBUG
        if (round_count <= 10)
            startTime();
        #endif
        if (parallel_method == InterSubgraph)
            output = edgeMap_once(G_vector, Frontier, BFS_F(Results));
        else if (parallel_method == EdgeLevel || parallel_method == MergeSmallDegree)
            output = edgeMap(G_merged, Frontier, BFS_F(Results));
        #ifdef MERGE_DEBUG
        if (round_count <= 10)
            nextTime("Running time");
        #endif
        uint64_t active_sum = 0;
        Frontier.del();
        Frontier = output;  // set new frontier
    }

    uint64_t sum = 0;    // use the sum of Results checking if our system gets right result.
    for (int i = 0; i < n; i++) {
        if (Results[i] == 1)
            sum ++;
    }
    cout << "Results sum = " << sum << endl;
    Frontier.del();
    free(Results);
}


LocalResult readResults(string path) {
    uint64_t* input = mmap_binary_read(path);
    long vertex_cnt = *reinterpret_cast<long*>(input);
    uintE* node_value = reinterpret_cast<uintE*>(input + 1);
    return LocalResult(vertex_cnt, node_value);
}

int parallel_main(int argc, char* argv[]) {
    commandLine P(argc, argv, " [-s] <inFile>");
    vector<string> iFile = P.getAllOptionValue("-g");
    int graph_num = iFile.size();
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool binary = P.getOptionValue("-b");
    bool mmap = P.getOptionValue("-m");
    // cout << "mmap = " << mmap << endl;
    long rounds = 3;
    {
        vector<graph<asymmetricVertex>> G_vector;
        vector<LocalResult> LocalResults;
        uint64_t edge_sum = 0;
        uint64_t vertex_num = 0;
        parallel_for (int i = 0; i < graph_num; i++) {
            string graph_name = iFile[i] + ".graph";
            string tree_name = iFile[i] + ".tree";
            graph<asymmetricVertex> subG = readGraph<asymmetricVertex>(
                (char*)graph_name.c_str(), compressed, symmetric, binary,
                mmap);  // asymmetric graph
            vertex_num = subG.n;
            #pragma omp critical
            {
                LocalResults.push_back(readResults(tree_name));
                G_vector.push_back(subG);
                edge_sum += subG.m;
            }
        }
        #ifdef MERGE_DEBUG
        cout << "total edge num = " << edge_sum << endl;
        #endif

        double duration = 0, merge_time = 0, compute_time = 0;
        double avg_degree = (double)edge_sum / vertex_num;
        int parallel_method;
        if (avg_degree >= 2 * graph_num)
            parallel_method = InterSubgraph;
        else if (avg_degree > graph_num)
            parallel_method = EdgeLevel;
        else
            parallel_method = MergeSmallDegree;
        for (int i = 0; i < rounds + 1; i++)
        {
            if (parallel_method == InterSubgraph)
            {
                auto start = std::chrono::high_resolution_clock::now();
                graph<asymmetricSubgraphVertex> G_merged = graph<asymmetricSubgraphVertex>(nullptr, vertex_num, edge_sum, nullptr);
                auto mid = std::chrono::high_resolution_clock::now();
                // return 0;
                Compute<asymmetricVertex,asymmetricSubgraphVertex>(G_vector, G_merged, LocalResults, parallel_method);
                auto end = std::chrono::high_resolution_clock::now();
                if (i > 0){
                    duration += std::chrono::duration<double>(end - start).count();
                    merge_time += std::chrono::duration<double>(mid - start).count();
                    compute_time += std::chrono::duration<double>(end - mid).count();
                }
            }
            if (parallel_method == EdgeLevel)
            {
                auto start = std::chrono::high_resolution_clock::now();
                graph<asymmetricSubgraphVertex> G_merged = mergeGraph<asymmetricVertex>(G_vector);
                auto mid = std::chrono::high_resolution_clock::now();
                // return 0;
                Compute<asymmetricVertex,asymmetricSubgraphVertex>(G_vector, G_merged, LocalResults, parallel_method);
                auto end = std::chrono::high_resolution_clock::now();
                if (i > 0){
                    duration += std::chrono::duration<double>(end - start).count();
                    merge_time += std::chrono::duration<double>(mid - start).count();
                    compute_time += std::chrono::duration<double>(end - mid).count();
                }
            }
            if (parallel_method == MergeSmallDegree)
            {
                auto start = std::chrono::high_resolution_clock::now();
                graph<asymmetricSubgraphVertexWithSmall> G_merged = mergeGraphWithSmall<asymmetricVertex>(G_vector);
                auto mid = std::chrono::high_resolution_clock::now();
                // return 0;
                Compute<asymmetricVertex,asymmetricSubgraphVertexWithSmall>(G_vector, G_merged, LocalResults, parallel_method);
                auto end = std::chrono::high_resolution_clock::now();
                if (i > 0){
                    duration += std::chrono::duration<double>(end - start).count();
                    merge_time += std::chrono::duration<double>(mid - start).count();
                    compute_time += std::chrono::duration<double>(end - mid).count();
                }
            }
        }
        std::cout << "merge time: " << merge_time / rounds << " seconds" << std::endl;
        std::cout << "compute time: " << compute_time / rounds << " seconds" << std::endl;
        std::cout << "execution time: " << duration / rounds << " seconds" << std::endl;

        for (int i = 0; i < graph_num; i++) {
            if (G_vector[i].transposed)
                G_vector[i].transpose();
            G_vector[i].del();
        }
    }
}
