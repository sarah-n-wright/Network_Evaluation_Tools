# NCBI gene ids of each protein
id1=$1
id2=$2
# file path to the file containing protein sequences from uniprot
ref_file=$3

# save results to example_outputs
outdir=../Data/example_outputs/AlphaFold
outfile=$outdir/${id1}_${id2}.fasta

echo '>'${id1}_${id2} > $outfile.temp

# search reference file for id1 and append sequence
grep -w $id1 $ref_file | awk -F'\t' '{print $3}' >> $outfile.temp
#separate the two sequences with ':'
echo ':' >> $outfile.temp
# search  reference file for id2 and append sequence
grep -w $id2 $ref_file | awk -F'\t' '{print $3}' >> $outfile.temp

# create the final fasta file
awk 'NR == 1 {print; next} {printf "%s", $0}' $outfile.temp > $outfile

echo '' >> $outfile
# remove temporary files
rm $outfile.temp
