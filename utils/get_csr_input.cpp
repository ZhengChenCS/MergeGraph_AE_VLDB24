#include "ligra.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

struct edge
{
    uint64_t src, dst, weight, ts;
    edge(){}
    edge(uint64_t _src, uint64_t _dst): src(_src), dst(_dst)
    {
        weight = (src + dst) % 128;
    }
    // edge(uint64_t _src, uint64_t _dst, uint64_t _weight): src(_src), dst(_dst), weight(_weight){}
    edge(uint64_t _src, uint64_t _dst, uint64_t _ts): src(_src), dst(_dst), ts(_ts){weight = (src + dst) % 128;}
    edge(uint64_t _src, uint64_t _dst, uint64_t _weight, uint64_t _ts): src(_src), dst(_dst), weight(_weight), ts(_ts){}
    void build_weight(uint64_t _src, uint64_t _dst, uint64_t _weight)
    {
        src = _src;
        dst = _dst;
        weight = _weight;
    }
};

bool cmp_dst(const edge &a, const edge &b)
{
    return a.dst < b.dst;
}

struct du_data
{
    uint64_t du, id;
};

bool cmp(du_data &a, du_data &b)
{
    return a.du > b.du;
}

int main(int argc, char **argv)
{
    commandLine P(argc, argv);
    std::string file_path;
    if(P.getOption("-f")){
        file_path = P.getOptionValue("-f");
    }else{
        fprintf(stderr, "No input data.\n");
        return 0;
    }
    std::string graph_type;
    if(P.getOption("-t")){
        graph_type = P.getOptionValue("-t");
    }else{
        fprintf(stderr, "No graph type.\n");
        return 0;
    }
    std::string output_path;
    if(P.getOption("-o")){
        output_path = P.getOptionValue("-o");
    }else{
        fprintf(stderr, "No output path.\n");
        return 0;
    }
    bool gene_weight;
    if(P.getOption("-g")){
        gene_weight = true;
    }else{
        gene_weight = false;
    }
    bool undirect = false;
    if (P.getOption("-undirect"))
    {
        undirect = true;
    }
    bool create_edge = false;
    if (P.getOption("-create"))
    {
        create_edge = true;
    }
    std::string weighted_file, unweighted_file;
    if (undirect)
    {
        weighted_file = output_path + "_undirect_weighted.txt", unweighted_file = output_path + "_undirect_unweighted.txt";
    }
    else
    {
        weighted_file = output_path + "_weighted.txt", unweighted_file = output_path + "_unweighted.txt";
    }
    FILE *out_weighted, *out_unweighted;
    if( !(out_weighted=fopen(weighted_file.c_str(), "w"))){
        fprintf(stderr, "Open weighted file failed.\n");
    } 
    if(gene_weight == true)
    {
        if( !(out_unweighted=fopen(unweighted_file.c_str(), "w"))){
            fprintf(stderr, "Open unweighted file failed.\n");
        } 
    }
    // out = fopen(output_path.c_str(), "w");

    uint64_t* input;
    input = mmap_binary_read(file_path);
    uint64_t raw_len = get_num_vertices(file_path);
    uint64_t num_vertices = 0;
    uint64_t num_edges = 0;
    std::vector<std::vector<edge>> adj_graph;
    std::vector<uint64_t> vlist;
    std::vector<uint64_t> elist;
    std::vector<uint64_t> weight;
    std::vector<uint64_t> timestamp;
    uint64_t *du;

    if(graph_type == "plain") // no weight, no timestamp
    {
        num_edges = raw_len / 2;
        // for(uint64_t eid = 0; eid < num_edges; ++eid){
        //     num_vertices = std::max(num_vertices, input[eid*2]);
        //     num_vertices = std::max(num_vertices, input[eid*2+1]);
        // }
        num_vertices = *input; input++;
        printf("|V|=%lu,|E|=%lu\n",num_vertices, num_edges);
        if (create_edge)
            num_vertices++;
        adj_graph.resize(num_vertices);
        du = new uint64_t[num_vertices];
        for(uint64_t eid = 0; eid < raw_len; eid += 2)
        {
            uint64_t src = input[eid];
            uint64_t dst = input[eid+1];
            adj_graph[src].push_back(edge(src, dst));
            du[src]++;
            if (undirect)
            {
                // du[dst]++;
                adj_graph[dst].push_back(edge(dst, src));
            }
        }
        if (create_edge)
        {
            du_data *d = new du_data[num_vertices];
            for (int i = 0; i < num_vertices; i++)
            {
                d[i].du = du[i];
                d[i].id = i;
            }
            sort(d, d + num_vertices, cmp);
            uint64_t src = num_vertices - 1;
            for (int j = 0; j < 10; j++)
            {
                uint64_t dst = d[j].id;
                uint64_t weight = (src + dst) % 128;
                adj_graph[src].push_back(edge(src, dst));
                printf("create edge %lu -> %lu\n", src, dst);
                if (undirect)
                {
                    adj_graph[dst].push_back(edge(dst, src));
                    printf("create edge %lu -> %lu\n", dst, src);
                }
            }
        }
        uint64_t offset = 0;
        vlist.push_back(0);
        for(uint64_t vid = 0; vid < num_vertices; ++vid){
            for(uint64_t eid = 0; eid < adj_graph[vid].size(); ++eid){
                elist.push_back(adj_graph[vid][eid].dst);
                weight.push_back(adj_graph[vid][eid].weight);
            }
            offset += adj_graph[vid].size();
            vlist.push_back(offset);
        }
        num_edges = elist.size();
        {
            fprintf(out_unweighted, "AdjacencyGraph\n");
            fprintf(out_unweighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_unweighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_unweighted, "%lu\n", elist[eid]);
            }
        }
        if (gene_weight == true)
        {
            fprintf(out_weighted, "WeightedAdjacencyGraph\n");
            fprintf(out_weighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_weighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", elist[eid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", weight[eid]);
            }
        }
    }
    else if(graph_type == "only_weight")
    {
        num_edges = raw_len / 3;
        // for(uint64_t eid = 0; eid < num_edges; ++eid){
        //     num_vertices = std::max(num_vertices, input[eid*3]);
        //     num_vertices = std::max(num_vertices, input[eid*3+1]);
        // }
        num_vertices = *input; input++;
        printf("|V|=%lu,|E|=%lu\n",num_vertices, num_edges);
        if (create_edge)
            num_vertices++;
        du_data *d = new du_data[num_vertices];
        for (int i = 0; i < num_vertices; i++)
        {
            d[i].du = 0;
            d[i].id = i;
        }
        adj_graph.resize(num_vertices);
        for(uint64_t eid = 0; eid < num_edges; ++eid){
            edge e;
            uint64_t src = input[eid*3];
            uint64_t dst = input[eid*3+1];
            uint64_t _weight = input[eid*3+2];
            e.build_weight(src, dst, _weight);
            d[src].du++;
            adj_graph[src].push_back(e);
            if (undirect)
            {
                // d[dst].du++;
                e.build_weight(dst, src, _weight);
                adj_graph[dst].push_back(e);
            }
        }
        if (create_edge)
        {
            sort(d, d + num_vertices, cmp);
            uint64_t src = num_vertices - 1;
            for (int j = 0; j < 10; j++)
            {
                uint64_t dst = d[j].id;
                uint64_t weight = (src + dst) % 128;
                edge e;
                e.build_weight(src, dst, weight);
                adj_graph[src].push_back(e);
                printf("create edge %lu -> %lu\n", src, dst);
                if (undirect)
                {
                    e.build_weight(dst, src, weight);
                    adj_graph[dst].push_back(e);
                    printf("create edge %lu -> %lu\n", dst, src);
                }
            }
        }
        uint64_t offset = 0;
        vlist.push_back(0);
        for(uint64_t vid = 0; vid < num_vertices; ++vid){
            for(uint64_t eid = 0; eid < adj_graph[vid].size(); ++eid){
                elist.push_back(adj_graph[vid][eid].dst);
                weight.push_back(adj_graph[vid][eid].weight);
                timestamp.push_back(adj_graph[vid][eid].ts);
            }
            offset += adj_graph[vid].size();
            vlist.push_back(offset);
        }
        num_edges = elist.size();
        {
            fprintf(out_unweighted, "AdjacencyGraph\n");
            fprintf(out_unweighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_unweighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_unweighted, "%lu\n", elist[eid]);
            }
        }
        {
            fprintf(out_weighted, "WeightedAdjacencyGraph\n");
            fprintf(out_weighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_weighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", elist[eid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", weight[eid]);
            }
        }
    }
    else if(graph_type == "only_timestamp")
    {
        num_edges = raw_len / 3;
        num_vertices = *input; input++;
        printf("|V|=%lu,|E|=%lu\n",num_vertices, num_edges);
        if (create_edge)
            num_vertices++;
        // for(uint64_t eid = 0; eid < num_edges; ++eid){
        //     num_vertices = std::max(num_vertices, input[eid*3]);
        //     num_vertices = std::max(num_vertices, input[eid*3+1]);
        // }
        // printf("|V|=%lu,|E|=%lu\n", num_vertices, num_edges);
        adj_graph.resize(num_vertices);
        du_data *d = new du_data[num_vertices];
        for (int i = 0; i < num_vertices; i++)
        {
            d[i].du = 0;
            d[i].id = i;
        }
        uint64_t edge_num = 0;
        for(uint64_t eid = 0; eid < num_edges; ++eid){
            uint64_t src = input[eid*3];
            uint64_t dst = input[eid*3+1];
            uint64_t _timestamp = input[eid*3+2];
            if (src > num_vertices || dst > num_vertices)
            {
                // printf("src = %ld, dst = %ld\n", src, dst);
                edge_num++;
                continue;
            }
            adj_graph[src].push_back(edge(src, dst, _timestamp));
            d[src].du++;
            if (undirect)
            {
                // d[dst].du++;
                adj_graph[dst].push_back(edge(dst, src, _timestamp));
            }
        }
        cout << edge_num << endl;
        if (create_edge)
        {
            sort(d, d + num_vertices, cmp);
            uint64_t src = num_vertices - 1;
            for (int j = 0; j < 10; j++)
            {
                uint64_t dst = d[j].id;
                uint64_t weight = (src + dst) % 128;
                adj_graph[src].push_back(edge(src, dst));
                printf("create edge %lu -> %lu\n", src, dst);
                if (undirect)
                {
                    adj_graph[dst].push_back(edge(dst, src));
                    printf("create edge %lu -> %lu\n", dst, src);
                }
            }
        }
        uint64_t offset = 0;
        vlist.push_back(0);
        for(uint64_t vid = 0; vid < num_vertices; ++vid){
            for(uint64_t eid = 0; eid < adj_graph[vid].size(); ++eid){
                elist.push_back(adj_graph[vid][eid].dst);
                weight.push_back(adj_graph[vid][eid].weight);
                timestamp.push_back(adj_graph[vid][eid].ts);
            }
            offset += adj_graph[vid].size();
            vlist.push_back(offset);
        }
        num_edges = elist.size();
        {
            fprintf(out_unweighted, "AdjacencyGraph\n");
            fprintf(out_unweighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_unweighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_unweighted, "%lu\n", elist[eid]);
            }
        }
        if (gene_weight == true)
        {
            fprintf(out_weighted, "WeightedAdjacencyGraph\n");
            fprintf(out_weighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_weighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", elist[eid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", weight[eid]);
            }
        }
    }
    else if(graph_type == "weight_and_timestamp")
    {
        num_edges = raw_len / 4;
        // for(uint64_t eid = 0; eid < num_edges; ++eid){
        //     num_vertices = std::max(num_vertices, input[eid*4]);
        //     num_vertices = std::max(num_vertices, input[eid*4+1]);
        // }
        num_vertices = *input; input++;
        printf("|V|=%lu,|E|=%lu\n",num_vertices, num_edges);
        if (create_edge)
            num_vertices++;
        du_data *d = new du_data[num_vertices];
        for (int i = 0; i < num_vertices; i++)
        {
            d[i].du = 0;
            d[i].id = i;
        }
        adj_graph.resize(num_vertices);
        for(uint64_t eid = 0; eid < num_edges; ++eid){
            uint64_t src = input[eid*4];
            uint64_t dst = input[eid*4+1];
            uint64_t _weight = input[eid*4+2];
            uint64_t _timestamp = input[eid*4+3];
            d[src].du++;
            adj_graph[src].push_back(edge(src, dst, _weight, _timestamp));
            if (undirect)
            {
                // d[dst].du++;
                adj_graph[dst].push_back(edge(dst, src, _weight, _timestamp));
            }
        }
        if (create_edge)
        {
            sort(d, d + num_vertices, cmp);
            uint64_t src = num_vertices - 1;
            for (int j = 0; j < 10; j++)
            {
                uint64_t dst = d[j].id;
                uint64_t weight = (src + dst) % 128;
                adj_graph[src].push_back(edge(src, dst));
                printf("create edge %lu -> %lu\n", src, dst);
                if (undirect)
                {
                    adj_graph[dst].push_back(edge(dst, src));
                    printf("create edge %lu -> %lu\n", dst, src);
                }
            }
        }
        uint64_t offset = 0;
        vlist.push_back(0);
        for(uint64_t vid = 0; vid < num_vertices; ++vid){
            for(uint64_t eid = 0; eid < adj_graph[vid].size(); ++eid){
                elist.push_back(adj_graph[vid][eid].dst);
                weight.push_back(adj_graph[vid][eid].weight);
            }
            offset += adj_graph[vid].size();
            vlist.push_back(offset);
        }
        num_edges = elist.size();
        {
            fprintf(out_unweighted, "AdjacencyGraph\n");
            fprintf(out_unweighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_unweighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_unweighted, "%lu\n", elist[eid]);
            }
        }
        {
            fprintf(out_weighted, "WeightedAdjacencyGraph\n");
            fprintf(out_weighted, "%lu\n%lu\n", num_vertices, num_edges);
            for (int vid = 0; vid < num_vertices; vid++)
            {
                fprintf(out_weighted, "%lu\n", vlist[vid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", elist[eid]);
            }
            for (int eid = 0; eid < num_edges; eid++)
            {
                fprintf(out_weighted, "%lu\n", weight[eid]);
            }
        }
    }
    fclose(out_weighted);
    fclose(out_unweighted);
    return 0;
}