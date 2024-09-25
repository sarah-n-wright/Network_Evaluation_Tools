#!/bin/bash

cwd=$(echo $PWD)

network_list=$1
execdir=$2 # path to the L3 executable
parallel_path=$3 # path to the parallel executable
prefixes=($(cat $network_list))

datadir=${cwd}/../Data/example_outputs/Interaction_Prediction/
outdir=${cwd}/../Data/example_outputs/Interaction_Prediction/L3/


for net_name in ${prefixes[@]}; do
	# Perform prediction for all 10 folds
	process_fold() {
	  i=$1
	  echo "FOLD:" $i
	  python $cwd/../neteval/edge_prediction.py \
	  --networkprefix $net_name --runwhat Predict --networksuffix .fold$i \
	  --datadir $datadir --execdir $execdir \
	  --outdir $outdir --pred_method L3
	}
	export -f process_fold
	export net_name cwd datadir execdir outdir
	$parallel_path -j 3 process_fold ::: {1..10}

	# Perform prediction for the full network
	 python $cwd/../neteval/edge_prediction.py \
	 --networkprefix $net_name --runwhat Predict \
	 --datadir $datadir --execdir $execdir \
	 --outdir $outdir --pred_method L3

done

