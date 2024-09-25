#!/bin/bash

# Parse the inputs
network_pref=$1
maxres=$2
datadir=$3
hidef_path=$4

# create full file paths
network_file=$datadir/${net_name}_net.txt
outpref=$datadir/$network_pref.$maxres

# Reformat the network file
awk -F'\t' '{print $1 "\t" $2}' $network_file > $network_file.hidef
sed -i 's/Entrez_/Node/g' $network_file.hidef

# Run HiDeF
/usr/bin/time -v python $hidef_path/hidef_finder.py --g $network_file.hidef \
	--maxres $maxres --o $outpref

rm $network_file.hidef