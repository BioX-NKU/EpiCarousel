#!/bin/bash
time=`echo Begin time: $(date +"%H:%M.%S")`

tmp_fifofile="$$.fifo"

mkfifo $tmp_fifofile

exec 6<>$tmp_fifofile
rm $tmp_fifofile
thread_num=${9}  # threads

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
metacell_resolution=${7}
steps=${8}
if_bi=${10}
if_mc_bi=${11}
mc_mode=${12}
threshold=${13}
filter_rate=${14}
neighbors_method=${15}
n_components=${16}
svd_solver=${17}


for i in `seq 1 $folds`
do
read -u6
    {
    echo "Start $i!"
    python identify_metacells.py ${data_name} $i ${chunk_dir} ${chunk_preprocessed_dir} ${chunk_mc_dir} ${chunk_mc_binarized_dir} ${metacell_resolution} ${steps} ${if_bi} ${if_mc_bi} ${mc_mode} ${threshold} ${filter_rate} ${neighbors_method} ${n_components} ${svd_solver}
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
