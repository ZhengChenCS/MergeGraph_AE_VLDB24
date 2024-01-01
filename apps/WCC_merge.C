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
// #define WEIGHTED 1
#define OPTIMIZE_GET_NEIGHBOR
// #define ONLY_DENSE
#include "ligra.h"
#include "edgemap_once.h"
#define INIT (uintE)0xffffff
// #define MERGE_DEBUG

struct LocalResult {
    uintE* value;
    long n;
    LocalResult() {}
    LocalResult(long _n, uintE* arr) {
        value = arr;
        n = _n;
    }
    // ~LocalResult(){
    // 	free(value);
    // }
};

struct CC_F{
    int numGraph;
    long numVertices;
    uintE **subTreeRoot;
    uintE *father;
    uintT comp_cnt;
    CC_F(long n, int numG, vector<LocalResult> &LocalResults, uintE *_father){
        numVertices = n;
        numGraph = numG;
        subTreeRoot = newA(uintE*, numGraph);
        for (int i = 0; i < numGraph; i++)
        {
            subTreeRoot[i] = LocalResults[i].value;
        }
        father = _father;
        comp_cnt = 0;
    }
    uintE find_root(int vid)
    {
        // comp_cnt ++;
        if (father[vid] == vid) return vid;
        else{
            uintE root = find_root(father[vid]);
            if (root != father[vid]) CAS(&father[vid], father[vid], root);
            return root;
        }
    }
    bool operator () (uintE vid) {
        for (int i = 0; i < numGraph; i++)
        {
            uintE dst = subTreeRoot[i][vid];
            if (dst == vid) continue;
            uintE vid_root = find_root(vid);
            uintE dst_root = find_root(dst);
            if (vid_root != dst_root)
            {
                if (vid_root > dst_root)
                {
                    CAS(&father[vid_root], vid_root, dst_root);
                }
                else
                {
                    CAS(&father[dst_root], dst_root, vid_root);
                }
            }
        }
        return 1;
    }
};

void Compute(vector<LocalResult> &LocalResults, commandLine P) {
    long n = LocalResults[0].n;
    int numGraph = LocalResults.size();
    uintE* ID = newA(uintE, n);

    //init 
    {
        parallel_for(long i = 0; i < n; i++) {
            ID[i] = LocalResults[0].value[i];
            for (int j = 1; j < numGraph; j++) {
                if (LocalResults[j].value[i] < ID[i])
                {
                    ID[i] = LocalResults[j].value[i];
                }
            }
        }
    }
    uint64_t round_count = 0;
    uint64_t active_sum = 0;
    uintE comp_cnt = 0;
    auto F = CC_F(n, numGraph, LocalResults, ID);
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
        F(i);
    uint64_t sum = 0;
    for (int i = 0; i < n; i++) {
        if (ID[i] == i)
            sum ++;
    }
    cout << "CC num = " << sum << endl;
    free(ID);
}

vector<LocalResult> LocalResults;

LocalResult readTree(string path) {
    uint64_t* input = mmap_binary_read(path);
    long vertex_cnt = *reinterpret_cast<long*>(input);
    uintE* node_value = reinterpret_cast<uintE*>(input + 1);
    return LocalResult(vertex_cnt, node_value);
}

int parallel_main(int argc, char* argv[]) {
    commandLine P(argc, argv, " [-s] <inFile>");
    vector<string> iFile = P.getAllOptionValue("-g");
    int g_num = iFile.size();
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool binary = P.getOptionValue("-b");
    bool mmap = P.getOptionValue("-m");
    long rounds = 3;
    {
        vector<graph<asymmetricVertex>> GA;
        uint64_t edge_sum = 0;
        parallel_for (int i = 0; i < g_num; i++) {
            string tree_name = iFile[i] + ".tree";
            #pragma omp critical
            {
                LocalResults.push_back(readTree(tree_name));
            }
        }
        #ifdef MERGE_DEBUG
        cout << "total edge num = " << edge_sum << endl;
        #endif
        double duration = 0, merge_time = 0, compute_time = 0;
        for (int i = 0; i < rounds + 1; i++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            Compute(LocalResults, P);
            auto end = std::chrono::high_resolution_clock::now();
            if (i > 0)
            {
                compute_time += std::chrono::duration<double>(end - start).count();
            }
        }
        std::cout << "execution time: " << compute_time / rounds << " seconds" << std::endl;

        // for (int i = 0; i < g_num; i++) {
        //     if (GA[i].transposed)
        //         GA[i].transpose();
        //     GA[i].del();
        // }
    }
}
