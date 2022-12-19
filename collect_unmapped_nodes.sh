#!/bin/bash -l
#SBATCH --job-name=pcnet_unmapped
#SBATCH --output /cellar/users/snwright/Data/SlurmOut/pcnet_unmapped.out
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0:02:00
jobid=$2

config=$1
source /cellar/users/snwright/Git/Network_Evaluation_Tools/data_configs/$config

outfile=${outpath}${name}_unmapped_nodes.txt
awk -v out=$outfile '/UNMAPPED NODES/, /END UNMAPPED NODES/ {print $2 > out}' /cellar/users/snwright/Data/SlurmOut/pcnet_data_$2.out
if [ ! -f "$outfile" ]; then
    touch $outfile
fi

sed -i '/UNMAPPED/d' $outfile
sed -i '/END/d' $outfile
