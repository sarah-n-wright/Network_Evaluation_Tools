#!/bin/bash

prefix_file=$1
stat=$2
min_genes=$3
max_genes=$4
nfolds=5
datadir=$5
gofile=$6

prefixes=($(cat $prefix_file))

for pref in "${prefixes[@]}":
	edgefile=$datadir/${pref}_net.txt
	outpath=$datadir/${pref}_GO_min${min_genes}_max${max_genes}_cv${nfolds}.${stat}

	srun -l Rscript egad_GBA.R $gofile $edgefile $min_genes $max_genes $stat $nfolds $outpath

