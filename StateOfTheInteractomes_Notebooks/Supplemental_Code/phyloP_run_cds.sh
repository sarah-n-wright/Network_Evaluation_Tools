#!/bin/bash
chr=$1
datadir=$2
chr_files=($(cat perGene_CDS/${chr}_files.txt))

task_start=$((SLURM_ARRAY_TASK_ID * 100))
task_end=$((task_start + 100))

if [ $task_end -gt ${#chr_files[@]} ]; then
	task_end=${#chr_files[@]}
fi

#fn=$datadir/${chr_files[$SLURM_ARRAY_TASK_ID]}

for ((i=task_start; i<task_end; i++)); do
	fn=$datadir/${chr_files[i]}
	if ! [ -f ${fn}.summary.bed ]; then
		echo ${fn}
		bedmap --faster --echo --delim '\t' --mean --median --stdev --mad --bases-uniq ${fn} $chr.phyloP30way.wigFix.gz.bed > ${fn}.summary.bed;
	fi
done
