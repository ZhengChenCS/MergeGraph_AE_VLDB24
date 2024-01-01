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
#include "ligra.h"

struct BFS_F {
    uintE* distance;
    BFS_F(uintE* _distance) : distance(_distance) {}
    inline bool update(uintE s, uintE d) {  // Update
        if (distance[d] == UINT_E_MAX) {
            distance[d] = distance[s] + 1;
            return 1;
        } else
            return 0;
    }
    inline bool updateAtomic(uintE s, uintE d) {  // atomic version of Update
        uintE new_dis = distance[s] + 1;
        return (CAS(&distance[d], UINT_E_MAX, new_dis));
    }
    // cond function checks if vertex has been visited yet
    inline bool cond(uintE d) { return (distance[d] == UINT_E_MAX); }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
    cout << "|V|=" << GA.n << " |E|=" << GA.m << endl;
    long start = P.getOptionLongValue("-r", 0);
    cout << "start=" << start << endl;
    FILE* output_file;
    bool output = false;
    if (P.getOption("-o"))
    {
        output = true;
        std::string out_path = P.getOptionValue("-o");
        if( !(output_file=fopen(out_path.c_str(), "w"))){
            fprintf(stderr, "Open out file failed.\n");
        }
    }

    long n = GA.n;
    // creates distance array, initialized to all -1, except for start
    uintE* distance = newA(uintE, n);
    parallel_for(long i = 0; i < n; i++) distance[i] = UINT_E_MAX;
    distance[start] = 0;
    uintE comp_cnt = 0;
    vertexSubset Frontier(n, start);  // creates initial frontier
    while (!Frontier.isEmpty()) {     // loop until frontier is empty
        // comp_cnt += Frontier.m;
        // if (Frontier.isDense)
        // {
        //     for (int i = 0; i < n; i++){
        //         if (Frontier.isIn(i))
        //         {
        //             comp_cnt += GA.V[i].outDegree;
        //         }
        //     }
        // }
        // else
        // {
        //     for (int i = 0; i < Frontier.m; i++)
        //     {
        //         uintE id = Frontier.vtx(i);
        //         comp_cnt += GA.V[id].outDegree;
        //     }
        // }
        vertexSubset output = edgeMap(GA, Frontier, BFS_F(distance));
        Frontier.del();
        Frontier = output;  // set new frontier
    }

    // cout << "comp cnt = " << comp_cnt << endl;
    uint sum = 0;
    for (int i = 0; i < GA.n; i++)
    {
        if (distance[i] < UINT_E_MAX)
            sum++;
    }
    cout << sum << endl;

    if (output == true)
    {
        fwrite(&GA.n, sizeof(GA.n), 1, output_file);
        fwrite(&distance[0], sizeof(uintE), GA.n, output_file);
        fclose(output_file);
    }

    Frontier.del();
    free(distance);
}

int parallel_main(int argc, char* argv[]) {
    commandLine P(argc, argv, " [-s] <inFile>");
    char* iFile = P.getArgument(0);
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool binary = P.getOptionValue("-b");
    bool mmap = P.getOptionValue("-m");
    // cout << "mmap = " << mmap << endl;
    long rounds = P.getOptionLongValue("-rounds", 0);
    {
        graph<asymmetricVertex> G = readGraph<asymmetricVertex>(
            iFile, compressed, symmetric, binary, mmap);  // asymmetric graph
        startTime();
        Compute(G, P);
        nextTime("Running time");
        if (G.transposed)
            G.transpose();
        G.del();
    }
}
