#!/bin/bash

pref_file=$1
datadir=$2
outdir=$3
mps_dir=$4
parallel_path=$5

prefixes=($(cat $pref_file))

for net_name in ${prefixes[@]}; do
	top_dir=$outdir/$net_name/Topological
	mkdir -p $top_dir
	cd "${mps_dir}/topological_feature_extractor/bin"
	# Run for all cross-validation folds
	process_fold() {
		i=$1
		cd "${mps_dir}/topological_feature_extractor/bin"
    		java -Xms12g -Xmx12g algos.TopologicalFeaturesExtractor e $datadir/csvs ${net_name}.fold$i.csv $top_dir
	}
	export -f process_fold
	export datadir net_name mps_dir top_dir
	$parallel_path -j 3 process_fold ::: {1..10}
	# Run for the full network
    	java -Xms12g -Xmx12g algos.TopologicalFeaturesExtractor e $datadir/csvs ${net_name}.csv $top_dir
done
