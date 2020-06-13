# meta script to submit jobs to the queue

time_sleep=2

n_start=0
n_patterns=3
for ((i = n_start;  i<n_patterns; i++))
do
echo "submitting job $i"
qsub jobfile_epth_ob_prelearning_single_pattern.sh -v x=$i
sleep $time_sleep
done

