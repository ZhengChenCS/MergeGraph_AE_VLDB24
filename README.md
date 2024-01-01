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

```bash
./split input_graph_path output_graph_path
```

#### 2.4.2. Convert the splited graph into binary format

```bash
./convert_to_binary input_graph_path output_graph_path
```

### 2.4.3. Get transitional results

```bash
./get_bfs_transitional_results input_graph_path output_graph_path root
./get_sssp_transitional_results input_graph_path output_graph_path root
./get_sswp_transitional_results input_graph_path output_graph_path root
./get_wcc_transitional_results input_graph_path output_graph_path
```

### 2.5. Running

These applications will process the entire graph. 

### Breadth-First Search
```bash
./bfs binary_graph_path output_results_path root
```

### Single Source Shortest Path
```bash
./sssp binary_graph_path output_results_path root
```

### Single Source Widest Path
```bash
./sswp binary_graph_path output_results_path root
```

### Weakly Connected Components
```bash
./wcc binary_graph_path output_results_path
```

### Run applications with script

```bash
cd script
bash split.sh
bash convert_to_binary.sh
bash get_transitional_results.sh
bash bfs.sh
bash sssp.sh
bash sswp.sh
bash wcc.sh
```