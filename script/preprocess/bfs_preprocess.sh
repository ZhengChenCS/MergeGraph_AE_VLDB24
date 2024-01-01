alg=bfs
exec=../../bin/${alg}
csr_data=(2 4 8 16 32)

data_set=HepPh source=0

length=${#csr_data[@]}

data_dir=../../example_datasets/$data_set/subgraph
# data_dir=/home/SSD/dataset/$data_set/subgraph
for ((i=0; i<$length; i++));
do
    sub_dir=$data_dir/${csr_data[$i]}
    for ((j=0; j<${csr_data[$i]}; j++));
    do
        $exec -r $source -o $sub_dir/${alg}_input/$j.tree $sub_dir/${alg}_input/$j.graph
    done
done
