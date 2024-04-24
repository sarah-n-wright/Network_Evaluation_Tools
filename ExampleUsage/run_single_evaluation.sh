#!/bin/bash
#SBATCH --job-name=eval_net    # Job name
#SBATCH --output=eval_net_%A_%a.out 
#SBATCH --error=eval_net_%A_%a.err    # Standard output and error log
#SBATCH --cpus-per-task=2    # Number of CPUs per task
#SBATCH --mem-per-cpu=48G    # Memory per CPU
#SBATCH --time=3-00:00:00    # Maximum execution time

cwd=$(echo $PWD)
# Load any required modules or software here
network_full_path=${cwd}/../Data/example_outputs/dip_net.txt
net_name=dip
shuffs=50
samples=50
alpha=0.64
sampp=0.5

job_id=$SLURM_JOB_ID
set_files=(combined_genesets.genesets disgen_befree_Dec22_final5.txt gwas_catalog_Dec22_final.txt gwas_cat_2023_07_27.genesets)
set_names=(disgen gwas)
min_genes_list=(20 20)

datadir=$cwd/../Data/
gitdir=$cwd/../

for i in {1..2}; do
	set_file=${set_files[$i]}
	set_name=${set_names[$i]}
	min_genes=${min_genes_list[$i]}
	node_sets_file=$gitdir/Data/$set_file
	actual_AURPRCs_save_path=$datadir/example_outputs/$net_name.$set_name.auprcs.csv
	echo '>>>'$i $set_file $set_name

	srun -l python $gitdir/neteval/run_network_evaluation.py \
		--cores 2 -i $shuffs -n $samples -a $alpha -p $sampp --min_genes $min_genes \
		--null_network_outdir $datadir/Evaluation/Null_Networks/ \
		--null_AUPRCs_save_path $datadir/Evaluation/AUPRCs/Null_AUPRCs/$net_name.$set_name.null_auprcs.csv \
		--performance_save_path $datadir/Evaluation/Performance/$net_name.$set_name.performance.csv \
		--performance_gain_save_path $datadir/Evaluation/Performance_Gain/$net_name.$set_name.performance_gain.csv \
		$network_full_path $node_sets_file $actual_AURPRCs_save_path
done
