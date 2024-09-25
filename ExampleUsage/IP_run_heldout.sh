#!/bin/bash

cwd=$(echo $PWD)

network_list=$1
pred_method=$2
parallel_path=$3 # path to the parallel executable
prefixes=($(cat $network_list))

datadir=${cwd}/../Data/example_outputs/InteractionPrediction/
outdir=${cwd}/../Data/example_outputs/InteractionPrediciton/L3/


for net_name in ${prefixes[@]}; do
	# Perform held out evaluation for all 10 folds
	process_fold() {
	  i=$1
	  echo "FOLD:" $i
	  python $cwd/../neteval/edge_prediction.py \
	  --networkprefix $net_name --runwhat EvaluateHeldOut --networksuffix .fold$i \
	  --datadir $datadir \
	  --outdir $outdir --pred_method $pred_method
	}
	export -f process_fold
	export net_name cwd datadir pred_method outdir
	$parallel_path -j 3 process_fold ::: {1..10}

done

