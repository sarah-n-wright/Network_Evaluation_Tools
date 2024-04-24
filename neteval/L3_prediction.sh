#!/bin/sh
file=$1
name=$2
outpath=$3
execdir=$4

cd $execdir

# check the number of columns
num_col=$(awk -F'\t' '{print NF; exit}' $file)
# if no score column, add a score column
# if number of columns == 3
if [[ "$num_col" -eq 3 ]]; then
    awk -F"\t" '(NR>1){print $1 "\t" $2 "\t" 1.0}' $file > ./${name}_temp_input_file.txt
else
    awk -F"\t" '(NR>1){print $1 "\t" $2 "\t" $3}' $file > ./${name}_temp_input_file.txt
# move to the directory. 
fi

./L3.out ${name}_temp_input_file.txt

mv ${name}_temp_input_file.txt $file

mv L3_predictions_${name}_temp_input_file.txt.dat $outpath/L3_predictions_${name}.dat

