#!/bin/bash -l
#SBATCH --job-name=shuffle
#SBATCH --output /cellar/users/snwright/Data/SlurmOut/pcnet_shuffle_%A_%a.out
#SBATCH --error /cellar/users/snwright/Data/SlurmOut/pcnet_shuffle_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2:00:00
#SBATCH --array=1-16

# Set the directory containing the input files
# Replace "/path/to/input/files" with the actual path
input_dir=$1
testmode=$2

# Set the output directory
# Replace "/path/to/output/files" with the actual path
output_dir=$input_dir/shuffled_networks/

# Identify a list of files in the input directory with the suffix "net.txt"
file_list=$(find $input_dir -maxdepth 1 -name "*net.txt")

# Determine the number of files in the list
num_files=$(echo $file_list | wc -w)

# Calculate the number of files per task
files_per_task=$((num_files / 16 + 1))

# Calculate the starting and ending indices for the current task
start_index=$(((SLURM_ARRAY_TASK_ID - 1) * files_per_task + 1))
end_index=$((SLURM_ARRAY_TASK_ID * files_per_task))

# Extract the subset of files for the current task
task_files=$(echo $file_list | cut -d " " -f $start_index-$end_index)

# Loop through each file in the list and apply a function to it
for file in $task_files; do
    echo $file
    # Replace this with your desired function
    # The output file name will be the same as the input file name, but in the output directory
    srun -l python -m memory_profiler shuffle_networks.py --testMode $testmode -o $output_dir $file
    # Add your function here
done