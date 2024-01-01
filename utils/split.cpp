#include "ligra.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

struct edge{
    uint64_t src;
    uint64_t dst;
    uint64_t ts;
    uint64_t weight;
    edge(uint64_t _src, uint64_t _dst, uint64_t _ts) : src(_src), dst(_dst), ts(_ts){}
    edge(uint64_t _src, uint64_t _dst, uint64_t _weight, uint64_t _ts) : src(_src), dst(_dst), weight(_weight), ts(_ts){}
};

struct du_data
{
    uint64_t du, id;
};

bool cmp(du_data &a, du_data &b)
{
    return a.du > b.du;
}

int RandomInt(int low, int high) {
    // 创建随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    // 定义随机数范围
    // 创建随机数分布范围
    std::uniform_int_distribution<> distr(low, high);
    // 生成随机数
    int num = distr(gen);
    return num;
}

double RandomDouble(double low, double high) {
    std::random_device rd;              // 用于获取随机数种子
    std::mt19937 gen(rd());            // 使用种子初始化Mersenne Twister生成器
    std::uniform_real_distribution<> dis(low, high); // 定义分布范围
    return dis(gen); // 生成并返回随机数
}

int main(int argc, char **argv){
    std::string file_path;
    int num_graph;
    commandLine P(argc, argv);
    if(P.getOption("-f")){
        file_path = P.getOptionValue("-f");
    }else{
        fprintf(stderr, "No input data.\n");
    }
    if(P.getOption("-n")){
        num_graph = P.getOptionIntValue("-n");
    }else{
        fprintf(stderr, "No num graph.\n");
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
    bool create_edge = false;
    if(P.getOption("-create"))
    {
        create_edge = true;
        fprintf(stderr, "create edge\n");
    }

    uint64_t *input;
    input = mmap_binary_read(file_path);

    uint64_t raw_len = get_num_vertices(file_path);
    uint64_t num_edges = 0;
    std::vector<struct edge> graph;
    std::vector<uint64_t> weight;
    
    uint64_t num_vertices = *input; input++;
    uint64_t *du = new uint64_t[num_vertices];

    if(graph_type == "has_weight"){
        num_edges = raw_len / 4;
        for(uint64_t i = 0; i < num_edges; ++i)
        {
            uint64_t src = input[i*4];
            uint64_t dst = input[i*4+1];
            uint64_t we = input[i*4+2];
            uint64_t ts = input[i*4+3];
            // num_vertices = std::max(num_vertices, std::max(src, dst));
            graph.emplace_back(edge(src, dst, we, ts));
            du[src]++;
        }
    }else{
        num_edges = raw_len / 3;
        for(uint64_t i = 0; i < num_edges; ++i)
        {
            uint64_t src = input[i*3];
            uint64_t dst = input[i*3+1];
            uint64_t ts = input[i*3+2];
            // num_vertices = std::max(num_vertices, std::max(src, dst));
            graph.emplace_back(edge(src, dst, ts));
            du[src]++;
        }
    }
    sort(graph.begin(), graph.end(), [](const edge &a, const edge &b){
        return a.ts < b.ts;
    }); 
    // std::cout << graph[0].ts << " " << graph[num_edges - 1].ts << endl;
    // return 0;
    du_data *d = new du_data[num_vertices];
    for (int i = 0; i < num_vertices; i++)
    {
        d[i].du = du[i];
        d[i].id = i;
    }
    sort(d, d + num_vertices, cmp);

    std::cout << "|V| = " << num_vertices 
              << " |E| = " << num_edges << " number of subgraph = " << num_graph << std::endl;

    bool avg = true;
    // if(P.getOption("-avg"))
    // {
    //     avg = true;
    //     fprintf(stderr, "avg split\n");
    // }

    uint64_t num_edges_per_graph = (num_edges + num_graph - 1) / num_graph;
    
    uint64_t offset = 0;
    if (create_edge)
        num_vertices++;
    if (avg == true)
    {
        for(uint64_t i = 0; i < num_graph; i++){
            FILE* fp;
            std::string output_sub_path = output_path + "/" + std::to_string(i) + ".bin";
            fp = fopen(output_sub_path.c_str(), "w");
            uint64_t start = i * num_edges_per_graph;
            uint64_t end = (i+1) * num_edges_per_graph > num_edges ? num_edges : (i+1) * num_edges_per_graph;

            fwrite(&num_vertices, sizeof(uint64_t), 1, fp);
            for(uint64_t j = start; j < end; ++j){
                if(graph_type == "has_weight"){
                    fwrite(&graph[j].src, sizeof(uint64_t), 1, fp);
                    fwrite(&graph[j].dst, sizeof(uint64_t), 1, fp);
                    fwrite(&graph[j].weight, sizeof(uint64_t), 1, fp);
                    fwrite(&graph[j].ts, sizeof(uint64_t), 1, fp);
                }else{
                    fwrite(&graph[j].src, sizeof(uint64_t), 1, fp);
                    fwrite(&graph[j].dst, sizeof(uint64_t), 1, fp);
                    fwrite(&graph[j].ts, sizeof(uint64_t), 1, fp);
                    // std::cout << graph[j].src << " " << graph[j].dst << std::endl;
                }
            }
            if (create_edge)
            {
                uint64_t src = num_vertices - 1;
                for (int j = 0; j < 10; j++)
                {
                    uint64_t dst = d[j].id;
                    if (i == 0 && num_graph == 2)
                    {
                        printf("create edge %lu -> %lu\n", src, dst);
                    }
                    if(graph_type == "has_weight"){
                        uint64_t weight = (src + dst) % 128;
                        fwrite(&src, sizeof(uint64_t), 1, fp);
                        fwrite(&dst, sizeof(uint64_t), 1, fp);
                        fwrite(&weight, sizeof(uint64_t), 1, fp);
                        fwrite(&graph[j].ts, sizeof(uint64_t), 1, fp);
                    }else{
                        fwrite(&src, sizeof(uint64_t), 1, fp);
                        fwrite(&dst, sizeof(uint64_t), 1, fp);
                        fwrite(&graph[j].ts, sizeof(uint64_t), 1, fp);
                        // std::cout << graph[j].src << " " << graph[j].dst << std::endl;
                    }
                }
            }
            fclose(fp);
        }
    }
    else
    {
        vector<vector<edge>> subgraph(num_graph);
        int eid = 0;
        double largest_size = 0.2;
        for (eid = 0; eid < num_edges * largest_size; eid++)
        {
            subgraph[0].push_back(graph[eid]);
        }
        double rest = 1 - largest_size;
        for (int i = 1; i < num_graph - 1; i++)
        {
            double this_size = RandomDouble(rest/num_graph/3, min(largest_size, rest - 0.025));
            while(subgraph[i].size() <= this_size * num_edges)
                subgraph[i].push_back(graph[eid++]);
            rest -= this_size;
        }
        while(eid < num_edges)
            subgraph[num_graph - 1].push_back(graph[eid++]);
        for (int i = 0; i < num_graph; i++)
        {
            int size = subgraph[i].size();
            cout << (double)size / num_edges << endl;
            FILE* fp;
            std::string output_sub_path = output_path + "/" + std::to_string(i) + ".bin";
            fp = fopen(output_sub_path.c_str(), "w");
            fwrite(&num_vertices, sizeof(uint64_t), 1, fp);
            for(uint64_t j = 0; j < size; ++j){
                if(graph_type == "has_weight"){
                    fwrite(&subgraph[i][j].src, sizeof(uint64_t), 1, fp);
                    fwrite(&subgraph[i][j].dst, sizeof(uint64_t), 1, fp);
                    fwrite(&subgraph[i][j].weight, sizeof(uint64_t), 1, fp);
                    fwrite(&subgraph[i][j].ts, sizeof(uint64_t), 1, fp);
                }else{
                    fwrite(&subgraph[i][j].src, sizeof(uint64_t), 1, fp);
                    fwrite(&subgraph[i][j].dst, sizeof(uint64_t), 1, fp);
                    fwrite(&subgraph[i][j].ts, sizeof(uint64_t), 1, fp);
                    // std::cout << graph[j].src << " " << graph[j].dst << std::endl;
                }
            }
            if (create_edge)
            {
                uint64_t src = num_vertices - 1;
                for (int j = 0; j < 10; j++)
                {
                    uint64_t dst = d[j].id;
                    if (i == 0 && num_graph == 2)
                    {
                        printf("create edge %lu -> %lu\n", src, dst);
                    }
                    if(graph_type == "has_weight"){
                        uint64_t weight = (src + dst) % 128;
                        fwrite(&src, sizeof(uint64_t), 1, fp);
                        fwrite(&dst, sizeof(uint64_t), 1, fp);
                        fwrite(&weight, sizeof(uint64_t), 1, fp);
                        fwrite(&subgraph[i][j].ts, sizeof(uint64_t), 1, fp);
                    }else{
                        fwrite(&src, sizeof(uint64_t), 1, fp);
                        fwrite(&dst, sizeof(uint64_t), 1, fp);
                        fwrite(&subgraph[i][j].ts, sizeof(uint64_t), 1, fp);
                        // std::cout << graph[j].src << " " << graph[j].dst << std::endl;
                    }
                }
            }
            fclose(fp);
        }
        
    }
    return 0;
}