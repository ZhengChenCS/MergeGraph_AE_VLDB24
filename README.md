# MergeGraph_AE_VLDB24

## 1. Description

This is the code for the paper "Enabling Window-Based Monotonic Graph Analytics with Reusable Transitional Results for Pattern-Consistent Queries" submitted to VLDB 2024.

## 2. Getting Started

### 2.1. System Dependencies

- [CMake](https://gitlab.kitware.com/cmake/cmake)
- OpenMP and C++17

### 2.2. Compilation

```bash
git clone https://github.com/ZhengChenCS/MergeGraph_AE_VLDB24
cd MergeGraph_AE_VLDB24
mkdir -p build
cd build
cmake ..
make -j
```

### 2.3. Graph Input Format
The initial input graph format should be in the adjacency graph format. For example, the [SNAP](https://snap.stanford.edu/data) format with timestamp and weight(optional) are shown below.

SNAP format:
```
0 1 1217567877
0 2 1217573801
2 0 1217606247
2 1 1217617639
```

### 2.4. Preprocessing

The graph format used by MergeGraph is binary edge lists with 64-bit vertex IDs.
MergeGraph provides some tools that converts edge lists in SNAP format to binary format.

#### 2.4.1. Split the graph into multiple files

If input graph is weighted, use "-t has_weight".

```bash
./split -f input_graph_path -t <has_weight|no_weight> -n number_of_subgraph -o output_graph_path
```

#### 2.4.2. Convert the splited graph into binary format

```bash
./convert_to_binary -n vertex_num <input_graph_path >output_graph_path
```

### 2.4.3. Get transitional results

```bash
./bfs -r source -o output_result_path input_graph_path
./sssp -r source -o output_result_path input_graph_path
./sswp -r source -o output_result_path input_graph_path
./wcc -r source -o output_result_path input_graph_path
```

### 2.5. Running

These applications will process the entire graph. 

### Breadth-First Search
```bash
./bfs_merge -g <graph_1_path graph_2_path ... graph_n_path>
```

### Single Source Shortest Path
```bash
./sssp_merge -g <graph_1_path graph_2_path ... graph_n_path>
```

### Single Source Widest Path
```bash
./sswp_merge -g <graph_1_path graph_2_path ... graph_n_path>
```

### Weakly Connected Components
```bash
./wcc_merge -g <graph_1_path graph_2_path ... graph_n_path>
```

### Run applications with script

```bash
cd script
bash get_input.sh
cd preprocess
bash bfs_preprocess.sh
bash sssp_preprocess.sh
bash sswp_preprocess.sh
bash wcc_preprocess.sh
cd ../merge
bash bfs_merge.sh
bash sssp_merge.sh
bash sswp_merge.sh
bash wcc_merge.sh
```