#!/bin/bash
#SBATCH --job-name=l3
#SBATCH --output=l3_pred_%A_%a.out
#SBATCH --error=l3_pred_%A_%a.err
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-48
#SBATCH --time=1-12:00:00

cwd=$(echo $PWD)

network_list=$1
run_what=$2 # must be Prediction, Evaluation or Both
prefixes=($(cat $network_list))
net_name=${prefixes[$SLURM_ARRAY_TASK_ID-1]}
job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
datadir=/path/to/interactome/file/directory/
execdir=/path/to/L3/executable
outdir=/path/to/write/results/
parallel_path=/path/to/parallel/executable

process_fold() {
  i=$1
  echo "FOLD:" $i
  python $cwd/../neteval/edge_prediction.py \
	--networkprefix $net_name --runwhat $run_what --networksuffix .fold$i \
	--benchmarks corum panther --datadir $datadir --execdir $execdir \
	--outdir $outdir
}

export -f process_fold
export net_name run_what

$parallel_path -j 3 process_fold ::: {1..10}



