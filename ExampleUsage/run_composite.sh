#!/bin/bash -l
#SBATCH --job-name=pcnet_comp
#SBATCH --output pcnet_comp_%A.out
#SBATCH --error pcnet_comp_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2:00:00

cwd=$(echo $PWD)

# file containing prefixes of networks to use in composite construction. If method=ranked, prefixes should be listed in rank order.
pref_file=$1
# prefix for naming the output composite networks
name=$2
# method of composite construction ('global' or 'ranked')
method=$3
netdir=/path/to/network/directory/
outdir=/directory/path/to/write/composites/

if [ $method == 'global' ]; then
srun -l python $cwd/../neteval/network_constructor.py -d $netdir \
	-o $outdir -m parsimonious \
	-n 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 --name $name --nodepref Entrez_ \
	$pref_file
fi
if [ $method == 'ranked' ]; then
srun -l python $cwd/../neteval/network_constructor.py -d $netdir \
	-o $outdir -m ordered \
	-n 2 --name $name --nodepref Entrez_ \
	$pref_file
fi

