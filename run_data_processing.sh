config=$1
test=$2

jobid=$(sbatch --parsable process_data_job.sh $config $test)

sbatch  --dependency=afterany:$jobid collect_unmapped_nodes.sh $config $jobid
