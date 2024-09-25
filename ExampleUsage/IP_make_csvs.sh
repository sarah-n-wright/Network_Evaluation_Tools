#!/bin/bash

pref_file=$1
datadir=$2

prefixes=($(cat $pref_file))

csv_dir=$datadir/csvs/
mkdir -p $csv_dir

for net_name in ${prefixes[@]}; do
	for i in {1..10}; do
		net_path=$datadir/${net_name}_net.txt.fold$i
		awk -F'\t' -v OFS=',' '{print $1, $2}' $net_path > $csv_dir/$net_name.fold$i.csv
		sed -i '1isource,target' $csv_dir/$net_name.fold$i.csv
	done
	full_net=$datadir/${net_name}_net.txt
	awk -F'\t' -v OFS=',' 'NR == 1 {gsub(/Entrez_A/, "source"); gsub(/Entrez_B/, "target"); print $1, $2; next} {print $1, $2}' "$full_net" > "$csv_dir/$net_name.csv"

done
