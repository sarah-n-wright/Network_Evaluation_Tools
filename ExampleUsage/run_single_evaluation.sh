#!/bin/bash
#SBATCH --job-name=eval_net    # Job name
#SBATCH --output=eval_net_%A_%a.out 
#SBATCH --error=eval_net_%A_%a.err    # Standard output and error log
#SBATCH --cpus-per-task=2    # Number of CPUs per task
#SBATCH --mem-per-cpu=48G    # Memory per CPU
#SBATCH --time=3-00:00:00    # Maximum execution time

cwd=$(echo $PWD)
net_prefix=$1
network_full_path=${cwd}/../Data/example_outputs/${net_prefix}_net.txt
net_name=$net_prefix
shuffs=50
samples=50
alpha=0.64
sampp=0.5

job_id=$SLURM_JOB_ID
set_files=(experimental.genesets disgen.genesets gwas.genestes gwas_20230727.genesets)
set_names=(exp_genesets disgen gwas gwas_230727)
min_genes_list=(10 20 20 10)

datadir=$cwd/../Data/
gitdir=$cwd/../

for i in {0..3}; do
	set_file=${set_files[$i]}
	set_name=${set_names[$i]}
	min_genes=${min_genes_list[$i]}
	node_sets_file=$gitdir/Data/$set_file
	actual_AURPRCs_save_path=$datadir/example_outputs/AUPRCs/$net_name.$set_name.auprcs.csv
	echo '>>>'$i $set_file $set_name

	srun -l python $gitdir/neteval/run_network_evaluation.py \
		--cores 2 -i $shuffs -n $samples -a $alpha -p $sampp --min_genes $min_genes \
		--null_network_outdir $datadir/example_outputs/Null_Networks/ \
		--null_AUPRCs_save_path $datadir/example_outputs/AUPRCs/Null_AUPRCs/$net_name.$set_name.null_auprcs.csv \
		--performance_save_path $datadir/example_outputs/Performance/$net_name.$set_name.performance.csv \
		--performance_gain_save_path $datadir/example_outputs/Performance_Gain/$net_name.$set_name.performance_gain.csv \
		$network_full_path $node_sets_file $actual_AURPRCs_save_path
done
