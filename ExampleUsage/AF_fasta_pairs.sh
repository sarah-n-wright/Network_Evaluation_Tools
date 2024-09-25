#!/bin/bash -l

# path to the file containing list of pairs (.csv)
pair_file=$1
# path to the file containing protein sequences from uniprot
ref_file=$2

if [ -f $pairfile ]; then

while IFS=',' read -r id1 id2; do
	sh get_fastas.sh "$id1" "$id2" $pair_file
done < $pairfile

fi
