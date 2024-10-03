#!/bin/bash

cwd=$(echo $PWD)
net_name=$1

shuffs=50
samples=50
alpha=0.64
sampp=0.3

set_files=(experimental.genesets disgen.genesets gwas.genesets gwas_20230727.genesets)
min_genes_list=(10 20 20 10)

gitdir=$cwd/../
network_full_path=$gitdir/Data/example_outputs/${net_name}_net.txt


for i in {0..3}; do
	set_file=${set_files[$i]}
	min_genes=${min_genes_list[$i]}
	node_sets_file=$gitdir/Data/$set_file
	echo 'ANALYZING ' $set_file
	python $gitdir/neteval/run_network_evaluation.py \
		--cores 2 -i $shuffs -n $samples -a $alpha -p $sampp --min_genes $min_genes \
		-o $gitdir/Data/example_outputs/GeneSetRecovery/ $network_full_path $gitdir/Data/$set_file
done
