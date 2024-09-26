#!/bin/bash -l

cwd=$(echo $PWD)

config=$1 # file name of configuration file in ../data_configs/
test_mode=$2 # should processing be run in test mode 0 (No) 1 (Yes, first 10000 interactions only)

if [ ! -f "${cwd}/../data_configs/$config" ] ; then
    echo "Error: config file $config not found"
    exit 1
fi
source $cwd/../data_configs/$config
outpath=$cwd/../Data/example_outputs/

python $cwd/../neteval/process_data.py \
    -A "$nodeA_col" -B "$nodeB_col" \
    -i $input_id_type -t $target_id_type -N $name \
    --species "$species_col" --species_value "$species_values" \
    --score "$score_col" --header $header \
    --sep $separator --testMode $test_mode  --prefix "$node_prefix_separator" \
    -o $outpath $input_datafile



