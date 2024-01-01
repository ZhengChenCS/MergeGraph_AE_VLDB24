cmake_minimum_required(VERSION 3.12)

project(MergeGraph)

set(CMAKE_CXX_STANDARD 17)

include_directories(ligra)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_CXX_FLAGS "-fopenmp -Ofast -march=native -DOPENMP -DEDGELONG")
# set(CMAKE_CXX_FLAGS "-fopenmp -march=native -DOPENMP -DEDGELONG")
# set(CMAKE_CXX_FLAGS "-Ofast -march=native -DEDGELONG")
# set(CMAKE_CXX_FLAGS "-march=native -DEDGELONG")
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_executable(split utils/split.cpp)
add_executable(convert2binary utils/convert2binary.cpp)
add_executable(get_csr_input utils/get_csr_input.cpp)

add_executable(bfs apps/BFS.C)
add_executable(wcc apps/WCC.C)
add_executable(sssp apps/SSSP.C)
add_executable(sswp apps/SSWP.C)

add_executable(bfs_merge apps/BFS_merge.C)
add_executable(sssp_merge apps/SSSP_merge.C)
add_executable(sswp_merge apps/SSWP_merge.C)
add_executable(wcc_merge apps/WCC_merge.C)

