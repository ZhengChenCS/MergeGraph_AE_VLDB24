#pragma once
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "graph.h"
#include "parallel.h"
using namespace std;
#define MERGE_THRESHOLD 8
#define InterSubgraph 1
#define EdgeLevel 2
#define MergeSmallDegree 3 

template <class oringin_vertex>
graph<asymmetricSubgraphVertex> mergeGraph(vector<graph<oringin_vertex>> &GA) {
    long graphNum = GA.size();
    long n = GA[0].n;

    long m = 0;
    for (int vid = 0; vid < graphNum; vid++) {
        m += GA[vid].m;
    }

    asymmetricSubgraphVertex* v = newA(asymmetricSubgraphVertex, n);
    uintT * inDegrees = newA(uintT, graphNum*n); 
    uintT * outDegrees = newA(uintT, graphNum*n);
    uintT * inOffsets = newA(uintT, graphNum*n);
    uintT * outOffsets = newA(uintT, graphNum*n);
#ifndef WEIGHTED
    uintE ** inNeighbors = newA(uintE*, graphNum*n);
    uintE ** outNeighbors = newA(uintE*, graphNum*n);
#else
    intE ** inNeighbors = newA(intE*, graphNum*n);
    intE ** outNeighbors = newA(intE*, graphNum*n);
#endif
    parallel_for(int vid = 0; vid < n; vid++) {
        uintT inDegree = 0, outDegree = 0;
        uintT base = vid*graphNum;
        for (int gid = 0; gid < graphNum; gid++) {
            inNeighbors[gid+base] = GA[gid].V[vid].inNeighbors;
            outNeighbors[gid+base] = GA[gid].V[vid].outNeighbors;
            inDegrees[gid+base] = GA[gid].V[vid].inDegree;
            outDegrees[gid+base] = GA[gid].V[vid].outDegree;
            inOffsets[gid+base] = inDegree;
            outOffsets[gid+base] = outDegree;
            inDegree += inDegrees[gid+base];
            outDegree += outDegrees[gid+base];
        }
        v[vid] = asymmetricSubgraphVertex(graphNum, inNeighbors+base, outNeighbors+base, inDegree,
                    outDegree, inDegrees+base, outDegrees+base, inOffsets+base, outOffsets+base);
    }
    MergeGraph_Mem<asymmetricSubgraphVertex>* mem = new MergeGraph_Mem<asymmetricSubgraphVertex>(v, n, m);
    return graph<asymmetricSubgraphVertex>(v, n, m, mem);
}

template <class oringin_vertex>
graph<asymmetricSubgraphVertexWithSmall> mergeGraphWithSmall(vector<graph<oringin_vertex>> &GA) {
    long graphNum = GA.size();
    long n = GA[0].n;

    long m = 0;
    for (int vid = 0; vid < graphNum; vid++) {
        m += GA[vid].m;
    }

    asymmetricSubgraphVertexWithSmall* v = newA(asymmetricSubgraphVertexWithSmall, n);
    uintT * inDegrees = newA(uintT, graphNum*n); 
    uintT * outDegrees = newA(uintT, graphNum*n);
    uintT * inOffsets = newA(uintT, graphNum*n);
    uintT * outOffsets = newA(uintT, graphNum*n);
#ifndef WEIGHTED
    uintE ** inNeighbors = newA(uintE*, graphNum*n);
    uintE ** outNeighbors = newA(uintE*, graphNum*n);
    uintE * inNeighborSmall, *outNeighborSmall;
#else
    intE ** inNeighbors = newA(intE*, graphNum*n);
    intE ** outNeighbors = newA(intE*, graphNum*n);
    intE * inNeighborSmall, *outNeighborSmall;
#endif
#ifndef WEIGHTED
    inNeighborSmall = newA(uintE, MERGE_THRESHOLD*n);
    outNeighborSmall = newA(uintE, MERGE_THRESHOLD*n);
#else
    inNeighborSmall = newA(intE, MERGE_THRESHOLD*n*2);
    outNeighborSmall = newA(intE, MERGE_THRESHOLD*n*2);
#endif
    parallel_for(int vid = 0; vid < n; vid++) {
        uintT inDegree = 0, outDegree = 0;
        uintT base = vid*graphNum;
        for (int gid = 0; gid < graphNum; gid++) {
            inNeighbors[gid+base] = GA[gid].V[vid].inNeighbors;
            outNeighbors[gid+base] = GA[gid].V[vid].outNeighbors;
            inDegrees[gid+base] = GA[gid].V[vid].inDegree;
            outDegrees[gid+base] = GA[gid].V[vid].outDegree;
            inOffsets[gid+base] = inDegree;
            outOffsets[gid+base] = outDegree;
            inDegree += inDegrees[gid+base];
            outDegree += outDegrees[gid+base];
        }
        bool inIsSmall = false, outIsSmall = false;
        if (inDegree <= MERGE_THRESHOLD)
        {
            inIsSmall = true;
        #ifndef WEIGHTED
            uintT offset = vid*MERGE_THRESHOLD;
        #else
            uintT offset = vid*MERGE_THRESHOLD*2;
        #endif
            for (int gid = 0; gid < graphNum; gid++)
            {
                for (int eid = 0; eid < GA[gid].V[vid].inDegree; eid++)
                {
                #ifndef WEIGHTED
                    inNeighborSmall[offset] = GA[gid].V[vid].inNeighbors[eid];
                    offset++;
                #else
                    inNeighborSmall[offset] = GA[gid].V[vid].inNeighbors[eid*2];
                    inNeighborSmall[offset + 1] = GA[gid].V[vid].inNeighbors[eid*2 + 1];
                    offset += 2;
                #endif
                }
            }
        }
        if (outDegree <= MERGE_THRESHOLD)
        {
            outIsSmall = true;
        #ifndef WEIGHTED
            uintT offset = vid*MERGE_THRESHOLD;
        #else
            uintT offset = vid*MERGE_THRESHOLD*2;
        #endif
            for (int gid = 0; gid < graphNum; gid++)
            {
                for (int eid = 0; eid < GA[gid].V[vid].outDegree; eid++)
                {
                #ifndef WEIGHTED
                    outNeighborSmall[offset] = GA[gid].V[vid].outNeighbors[eid];
                    offset++;
                #else
                    outNeighborSmall[offset] = GA[gid].V[vid].outNeighbors[eid*2];
                    outNeighborSmall[offset + 1] = GA[gid].V[vid].outNeighbors[eid*2 + 1];
                    offset += 2;
                #endif
                }
            }
    }
#ifndef WEIGHTED
    uintT offset = vid*MERGE_THRESHOLD;
#else
    uintT offset = vid*MERGE_THRESHOLD*2;
#endif
    v[vid] = asymmetricSubgraphVertexWithSmall(graphNum, inNeighbors+base, outNeighbors+base, inDegree, 
                outDegree, inDegrees+base, outDegrees+base, inOffsets+base, outOffsets+base,
                inIsSmall, outIsSmall, inNeighborSmall+offset, outNeighborSmall+offset);
    }
    MergeGraph_Mem<asymmetricSubgraphVertexWithSmall>* mem = new MergeGraph_Mem<asymmetricSubgraphVertexWithSmall>(v, n, m);
    return graph<asymmetricSubgraphVertexWithSmall>(v, n, m, mem);
}