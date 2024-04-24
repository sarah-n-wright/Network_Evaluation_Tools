#!/bin/bash
#SBATCH --job-name=mps
#SBATCH --output=mps_top_%A_%a.out
#SBATCH --error=mps_top_%A_%a.err
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:00:00

pref_file=$1
prefixes=($(cat $pref_file))
type=$2
name=${prefixes[$SLURM_ARRAY_TASK_ID-1]}
echo $name

output_dir=/directory/to/write/results/
input_dir=/directory/with/network/files/
mps_dir='/path/to/Git/PPI-Prediction-Project/Prediction methods/MPS'
top_dir=$output_dir/$type/$name/Topological
mkdir -p $top_dir


if [ "$type" = "external" ]; then

#for i in {1..10}; do

  #source activate MPS2
  process_fold() {
    i=$1
    echo 'FOLD:' $i
    #sh /cellar/users/snwright/Git/Network_Evaluation_Tools/run_files/fold_to_csv.sh $input_dir $name $i
    cd "${mps_dir}/topological_feature_extractor/bin"
    java -Xms12g -Xmx12g algos.TopologicalFeaturesExtractor e $input_dir/csvs ${name}.fold$i.csv $top_dir

  }
  export -f process_fold
  export input_dir name mps_dir top_dir

  ~/anaconda3/bin/parallel process_fold ::: {1..4}

#done
fi

if [ "$type" = "internal" ]; then
	cd "${mps_dir}/topological_feature_extractor/bin"
	java -Xms12g -Xmx12g algos.TopologicalFeaturesExtractor i $input_dir/csvs ${name}.csv $top_dir


fi
