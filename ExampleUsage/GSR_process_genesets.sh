#!/bin/bash -l
#SBATCH --job-name=pcnet_genesets
#SBATCH --output /cellar/users/snwright/Data/SlurmOut/pcnet_genesets_%A.out
#SBATCH --error /cellar/users/snwright/Data/SlurmOut/pcnet_genesets_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00

# 1 = geneset file
# 2 = prefix file
# 3 = output file
# 4 = id type
# 5 = min genes
# 6 = strict mode ('strict' or none)

datadir=/cellar/users/snwright/Data/Network_Analysis
setfile=$datadir/Reference_Data/$1 #.genesets
gitdir=/cellar/users/snwright/Git/Network_Evaluation_Tools

script=$gitdir/neteval/prepare_evaluation_data.py
netdir=$datadir/Processed_Data/v2_fixed/
prefix_file=$gitdir/Data/$2
out_file=$gitdir/Data/$3

min_genes=$5 #m
max_genes=500 #M
id_type=$4 #i
strict_mode=$6 #--strict

if [ "$id_type" == "Entrez" ]; then #no conversion needed
    if [ "$strict_mode" == "strict" ]; then
        srun -l python $script -s $setfile -i $id_type -n $prefix_file \
            -F -m $min_genes -M $max_genes -o $out_file -d $netdir --strict
    else
        srun -l python $script -s $setfile -i $id_type -n $prefix_file \
            -F -m $min_genes -M $max_genes -o $out_file -d $netdir
    fi
else
    if [ "$strict_mode" == "strict" ]; then
        srun -l python $script -s $setfile -i $id_type -n $prefix_file \
            -C -F -m $min_genes -M $max_genes -o $out_file -d $netdir --strict
    else
        srun -l python $script -s $setfile -i $id_type -n $prefix_file \
            -C -F -m $min_genes -M $max_genes -o $out_file -d $netdir
    fi
fi



