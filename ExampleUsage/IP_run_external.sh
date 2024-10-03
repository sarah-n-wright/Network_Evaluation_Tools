#!/bin/bash

cwd=$(echo $PWD)

network_list=$1
pred_method=$2
prefixes=($(cat $network_list))

datadir=${cwd}/../Data/example_outputs/InteractionPrediction/
outdir=${cwd}/../Data/example_outputs/InteractionPrediction/L3/


for net_name in ${prefixes[@]}; do
	python $cwd/../neteval/edge_prediction.py \
	  --networkprefix $net_name --runwhat EvaluateExternal \
	  --datadir $datadir --benchmarks corum panther \
	  --outdir $outdir --pred_method $pred_method

done

