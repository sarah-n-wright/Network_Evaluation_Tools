#!/bin/bash -l
#SBATCH --job-name=pcnet_eval
#SBATCH --output /cellar/users/snwright/Data/SlurmOut/pcnet_eval_%A.out
#SBATCH --error /cellar/users/snwright/Data/SlurmOut/pcnet_eval_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-00:00:00
# network=$1
net_name=$1
# network=biogrid.4.4.217_net.txt
set_file=GWAS_Catalog_genesets_entrez.txt
set_name=gwas_orig
# set_file=DisGeNET_genesets_entrez.txt
#network=$2
network=${net_name}_net.txt
# net_name=biogird.4.4.217

datadir=/cellar/users/snwright/Data/Network_Analysis
gitdir=/cellar/users/snwright/Git/Network_Evaluation_Tools
network_path=$datadir/Processed_Data/v2_2022/$network
node_sets_file=$gitdir/Data/$set_file
actual_AURPRCs_save_path=$datadir/Evaluation/AUPRCs/$net_name.$set_name.auprcs.csv

srun -l python run_network_evaluation.py \
	--cores 4 \
	--null_network_outdir $datadir/Evaluation/Null_Networks/${net_name}_ \
	--null_AUPRCs_save_path $datadir/Evaluation/AUPRCs/Null_AUPRCs/$net_name.$set_name.null_auprcs.csv \
	--performance_save_path $datadir/Evaluation/Performance/$net_name.$set_name.performance.csv \
	--performance_gain_save_path $datadir/Evaluation/Performance_Gain/$net_name.$set_name.performance_gain.csv \
	$network_path $node_sets_file $actual_AURPRCs_save_path

# all other parameters default
# --verbose -> not sure how to use this one
# --net_file_delim "\t"
# --set_file_delim "\t"
# --sample_p None 
# --alpha None
# --sub_sample_iter 30
# --background network
# --null_iter 30

