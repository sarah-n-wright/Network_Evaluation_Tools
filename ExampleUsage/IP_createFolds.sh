#!/bin/bash

cwd=$(echo $PWD)

network_list=$1
prefixes=($(cat $network_list))
datadir=${cwd}/../Data/example_outputs/
outdir=${cwd}/../Data/example_outputs/InteractionPrediction/L3/

for net_name in ${prefixes[@]}; do
	python $cwd/../neteval/edge_prediction.py \
	--networkprefix $net_name --runwhat Folds --datadir $datadir \
	--outdir $outdir
done

