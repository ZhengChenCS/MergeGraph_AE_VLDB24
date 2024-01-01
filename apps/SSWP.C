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
#define WEIGHTED 1
#include "ligra.h"

struct SSWP_F {
    intE* WidestPath;
    int* Visited;
    SSWP_F(intE* _WidestPath, int* _Visited)
        : WidestPath(_WidestPath), Visited(_Visited) {}
    inline bool update(
        uintE s,
        uintE d,
        intE edgeLen) {  // Update WidestPath if found a shorter path
        intE newDist = min(WidestPath[s], edgeLen);
        if (WidestPath[d] < newDist) {
            WidestPath[d] = newDist;
            if (Visited[d] == 0) {
                Visited[d] = 1;
                return 1;
            }
        }
        return 0;
    }
    inline bool updateAtomic(uintE s, uintE d, intE edgeLen) {  // atomic Update
        intE newDist = min(WidestPath[s], edgeLen);
        return (writeMax(&WidestPath[d], newDist) &&
                CAS(&Visited[d], 0, 1));
    }
    inline bool cond(uintE d) { return cond_true(d); }
};

// reset visited vertices
struct SSWP_Vertex_F {
    int* Visited;
    SSWP_Vertex_F(int* _Visited) : Visited(_Visited) {}
    inline bool operator()(uintE i) {
        Visited[i] = 0;
        return 1;
    }
};

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
    long start = P.getOptionLongValue("-r", 0);
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
    // initialize WidestPath to "-infinity"
    intE* WidestPath = newA(intE, n);
    { parallel_for(long i = 0; i < n; i++) WidestPath[i] = -INT_MAX / 2; }
    WidestPath[start] = INT_MAX / 2;

    int* Visited = newA(int, n);
    { parallel_for(long i = 0; i < n; i++) Visited[i] = 0; }

    vertexSubset Frontier(n, start);  // initial frontier

    long round = 0;
    uintE comp_cnt = 0;
    while (!Frontier.isEmpty()) {
        vertexSubset output =
            edgeMap(GA, Frontier, SSWP_F(WidestPath, Visited), GA.m / 20,
                    dense_forward);
        vertexMap(output, SSWP_Vertex_F(Visited));
        Frontier.del();
        Frontier = output;
        round++;
    }
    cout << "comp cnt = " << comp_cnt << endl;

    if (output == true)
    {
        fwrite(&GA.n, sizeof(GA.n), 1, output_file);
        fwrite(&WidestPath[0], sizeof(intE), GA.n, output_file);
        fclose(output_file);
    }

    uintE sum = 0;
    for (int i = 0; i < n; i++)
    {
        if (WidestPath[i] > 0 && i != start)
            sum += WidestPath[i];
    }
    cout << "sum = " << sum << endl;

    Frontier.del();
    free(Visited);
    free(WidestPath);
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