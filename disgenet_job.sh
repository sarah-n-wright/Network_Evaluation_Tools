#!/bin/bash -l
#SBATCH --job-name=disgenet
#SBATCH --output /cellar/users/snwright/Data/SlurmOut/pcnet_disgen_%A.out
#SBATCH --error /cellar/users/snwright/Data/SlurmOut/pcnet_disgen_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=4:00:00

email=$1
password=$2

srun -l python prepare_evaluation_data.py --email $email --password $password




