datadir=$1
binary_file=$2
cocomp_file=$3

sed 's/\t/,/g' $datadir/$binary_file | awk -F',' '{print $0 ",binary"}' | sed '1s/binary/source/g' > $datadir/binary_temp

sed 's/\t/,/g' $datadir/$cocomp_file | awk -F',' -v out=$datadir/cocomp_temp '(NR>1){print $0 ",cocomp" > out}' 

### add new line to end of all lines
sed -i -e '$a\' $datadir/binary_temp
sed -i -e '$a\' $datadir/cocomp_temp

cat $datadir/binary_temp $datadir/cocomp_temp > $datadir/$binary_file.COMBINED.$cocomp_file

sed -i 's/,/\t/g' $datadir/$binary_file.COMBINED.$cocomp_file