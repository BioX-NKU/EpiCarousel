#!/bin/bash
time=`echo Begin time: $(date +"%H:%M.%S")`

tmp_fifofile="$$.fifo"

mkfifo $tmp_fifofile

exec 6<>$tmp_fifofile
rm $tmp_fifofile
thread_num=${10}  # threads

#create key
for i in `seq 1 $thread_num`
do
    echo
done >&6 

# beginTime=`date +%s`
declare -i folds
folds=${1}
data_name=${2}
chunk_dir=${3}
chunk_preprocessed_dir=${4}
chunk_mc_dir=${5}
chunk_mc_binarized_dir=${6}
chunk_method=${7}
metacell_resolution=${8}
steps=${9}
if_bi=${11}
if_mc_bi=${12}
mc_mode=${13}
threshold=${14}
filter_rate=${15}
neighbors_method=${16}
decomposition=${17}
n_components=${18}
svd_solver=${19}
random_state=${20}
n_iter=${21}
kernel=${22}

for i in `seq 1 $folds`
do
read -u6
    {
    echo "Start $i!"
    python metacell_construction.py ${data_name} $i ${chunk_dir} ${chunk_preprocessed_dir} ${chunk_mc_dir} ${chunk_mc_binarized_dir} ${chunk_method} ${metacell_resolution} ${steps} ${if_bi} ${if_mc_bi} ${mc_mode} ${threshold} ${filter_rate} ${neighbors_method} ${decomposition} ${n_components} ${svd_solver} ${random_state} ${n_iter} ${kernel}
    echo "Finish $i!"
    echo >&6
    }&
done
wait
echo $time 
echo Finish time: $(date +"%H:%M.%S")

exec 6>$-
# endTime=`date +%s`
# echo "Total runtime:" $(($endTime-$beginTime)) "seconds"
