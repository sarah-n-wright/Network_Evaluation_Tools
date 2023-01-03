#!/bin/bash -l
#SBATCH --job-name=pcnet_data
#SBATCH --output /cellar/users/snwright/Data/SlurmOut/pcnet_data_%A.out
#SBATCH --error /cellar/users/snwright/Data/SlurmOut/pcnet_data_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --time=2:00:00

config=$1
source /cellar/users/snwright/Git/Network_Evaluation_Tools/data_configs/$config
test_mode=$2

srun -l python -m memory_profiler process_data.py -A "$nodeA_col" -B "$nodeB_col" \
    -i $input_id_type -t $target_id_type -N $name \
    --species "$species_col" --species_value "$species_values" \
    --score "$score_col" --header $header \
    --sep $separator --testMode $test_mode  --prefix "$node_prefix_separator" \
    -o $outpath $input_datafile



