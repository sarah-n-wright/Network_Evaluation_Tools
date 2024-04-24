net_dir=$1
net_pref=$2
temp_dir=$net_dir/csvs/
mkdir -p $temp_dir

net_path=$net_dir/${net_pref}_net.txt
awk -F'\t' -v OFS=',' 'NR == 1 {gsub(/Entrez_A/, "source"); gsub(/Entrez_B/, "target"); print $1, $2; next} {print $1, $2}' "$net_path" > "$temp_dir/$net_pref.csv"

echo $temp_dir
