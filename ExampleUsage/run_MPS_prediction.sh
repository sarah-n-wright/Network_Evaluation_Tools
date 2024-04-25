#!/bin/bash
#SBATCH --job-name=mps_pred
#SBATCH --output=mps_pred_%A_%a.out
#SBATCH --error=mps_pred_%A_%a.err
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-10

cwd=$(echo $PWD)

network_list=$1
prefixes=($(cat $network_list))
net_name=${prefixes[$SLURM_ARRAY_TASK_ID-1]}
run_what=$2
datadir=/path/to/interactome/file/directory/
outdir=/path/to/write/results/
parallel_path=/path/to/parallel/executable

process_fold() {
  i=$1
  echo "FOLD:" $i
  python $cwd/../neteval/edge_prediction.py \
	--networkprefix $net_name --runwhat $run_what --networksuffix .fold$i \
        --benchmarks corum panther --outdir $outdir --pred_method MPS --datadir $datadir
}

export -f process_fold
export net_name outdir run_what

$parallel_path -j 3 process_fold ::: {1..10}
