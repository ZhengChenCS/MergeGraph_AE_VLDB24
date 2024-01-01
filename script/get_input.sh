dataset=HepPh
file=.txt
nums=(2 4 8 16 32)
../bin/convert2binary -n 28093 <../example_datasets/$dataset/${dataset}.txt >../example_datasets/$dataset/${dataset}.bin
for num in ${nums[@]}
do
    inpath=../example_datasets/$dataset/${dataset}.bin
    outpath=../example_datasets/$dataset/subgraph/${num}
    if [ ! -d ../example_datasets/$dataset/subgraph ]; then
        mkdir -p ../example_datasets/$dataset/subgraph
    fi
    if [ ! -d ${outpath} ]; then
        mkdir -p ${outpath}
    fi
    if [ ! -d ${outpath}/wcc_input ]; then
        mkdir -p ${outpath}/wcc_input
    fi
    if [ ! -d ${outpath}/bfs_input ]; then
        mkdir -p ${outpath}/bfs_input
    fi
    if [ ! -d ${outpath}/sssp_input ]; then
        mkdir -p ${outpath}/sssp_input
    fi
    if [ ! -d ${outpath}/sswp_input ]; then
        mkdir -p ${outpath}/sswp_input
    fi
    ../bin/split -f ${inpath} -t has_weight -n ${num} -o ${outpath}
    subDir=$num
    for ((k=0; k<$num; k++));
    do
        coo_bin=${outpath}/${k}.bin
        ../bin/get_csr_input -f $coo_bin -t weight_and_timestamp -o ${outpath}/${k} -g -undirect
        ../bin/get_csr_input -f $coo_bin -t weight_and_timestamp -o ${outpath}/${k} -g
        cp ${outpath}/${k}_unweighted.txt ${outpath}/bfs_input/${k}.graph
        cp ${outpath}/${k}_weighted.txt ${outpath}/sssp_input/${k}.graph
        cp ${outpath}/${k}_weighted.txt ${outpath}/sswp_input/${k}.graph
        cp ${outpath}/${k}_undirect_unweighted.txt ${outpath}/wcc_input/${k}.graph
    done
done