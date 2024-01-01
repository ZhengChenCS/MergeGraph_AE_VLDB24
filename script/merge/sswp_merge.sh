alg=sswp
exec=../../bin/${alg}_merge
data_set=HepPh source=28093
graph_num=(2 4 8 16 32)
graph_num_len=${#graph_num[@]}
term_width=$(tput cols)

for ((j=0; j<graph_num_len; j++));
do
    padding=$(( (term_width - ${#graph_num[$j]}) / 2 ))
    printf "%*s%s\n" $padding '' "graph_num = ${graph_num[$j]}"
    data_dir=../../example_datasets/$data_set/subgraph/${graph_num[$j]}/${alg}_input
    if [ ${graph_num[$j]} -eq 2 ];then 
        padding=$(( (term_width - ${#exec}) / 2 ))
        printf "%*s%s\n" $padding '' "${exec}"
        ${exec} -r $source -g $data_dir/0 $data_dir/1
    elif [ ${graph_num[$j]} -eq 4 ];then
        padding=$(( (term_width - ${#exec}) / 2 ))
        printf "%*s%s\n" $padding '' "${exec}"
        ${exec} -r $source -g $data_dir/0 $data_dir/1 $data_dir/2 $data_dir/3
    elif [ ${graph_num[$j]} -eq 8 ];then
        padding=$(( (term_width - ${#exec}) / 2 ))
        printf "%*s%s\n" $padding '' "${exec}"
        ${exec} -r $source -g $data_dir/0 $data_dir/1 $data_dir/2 $data_dir/3 $data_dir/4 $data_dir/5 $data_dir/6 $data_dir/7
    elif [ ${graph_num[$j]} -eq 16 ];then
        padding=$(( (term_width - ${#exec}) / 2 ))
        printf "%*s%s\n" $padding '' "${exec}"
        ${exec} -r $source -g $data_dir/0 $data_dir/1 $data_dir/2 $data_dir/3 $data_dir/4 $data_dir/5 $data_dir/6 $data_dir/7 $data_dir/8 $data_dir/9 $data_dir/10 $data_dir/11 $data_dir/12 $data_dir/13 $data_dir/14 $data_dir/15
    elif [ ${graph_num[$j]} -eq 32 ];then
        padding=$(( (term_width - ${#exec}) / 2 ))
        printf "%*s%s\n" $padding '' "${exec}"
        ${exec} -r $source -g $data_dir/0 $data_dir/1 $data_dir/2 $data_dir/3 $data_dir/4 $data_dir/5 $data_dir/6 $data_dir/7 $data_dir/8 $data_dir/9 $data_dir/10 $data_dir/11 $data_dir/12 $data_dir/13 $data_dir/14 $data_dir/15 $data_dir/16 $data_dir/17 $data_dir/18 $data_dir/19 $data_dir/20 $data_dir/21 $data_dir/22 $data_dir/23 $data_dir/24 $data_dir/25 $data_dir/26 $data_dir/27 $data_dir/28 $data_dir/29 $data_dir/30 $data_dir/31
    fi
done