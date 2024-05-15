#!/bin/bash
#SBATCH --job-name=l3
#SBATCH --output=l3_pred_%A_%a.out
#SBATCH --error=l3_pred_%A_%a.err
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-3
#SBATCH --time=1-12:00:00

cwd=$(echo $PWD)

network_list=$1
prefixes=($(cat $network_list))
net_name=${prefixes[$SLURM_ARRAY_TASK_ID-1]}
job_id=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
datadir=${cwd}/../Data/example_outputs/
# !! Update to path of your L3 executable
execdir=${cwd}/../../../Data/Network_Analysis/Edge_Prediction/kpisti-L3-ed6b18f/
outdir=${cwd}/../Data/example_outputs/L3/
# !! Update to path of your parallel
parallel_path=~/anaconda3/bin/parallel

process_fold() {
  i=$1
  echo "FOLD:" $i
  python $cwd/../neteval/edge_prediction.py \
  --networkprefix $net_name --runwhat Prediction --networksuffix .fold$i \
  --benchmarks corum panther --datadir $datadir --execdir $execdir \
  --outdir $outdir
  python $cwd/../neteval/edge_prediction.py \
  --networkprefix $net_name --runwhat Evaluation --networksuffix .fold$i \
  --benchmarks corum panther --datadir $datadir --execdir $execdir \
  --outdir $outdir
  python $cwd/../neteval/edge_prediction.py \
	--networkprefix $net_name --runwhat GoldStandard --networksuffix .fold$i \
	--benchmarks corum panther --datadir $datadir --execdir $execdir \
	--outdir $outdir
}

export -f process_fold
export net_name run_what cwd datadir execdir outdir

$parallel_path -j 3 process_fold ::: {1..10}



