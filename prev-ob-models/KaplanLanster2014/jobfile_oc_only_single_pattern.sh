# The name of the script is neuron_job
#PBS -N concInv

# 4 hour wall-clock time will be given to this job
#PBS -l walltime=0:05:00
# delm PBS -l walltime=0:30:00
# a response curve sweep takes about 6 min

# Number of cores to be allocated is 24
#PBS -l mppwidth=48

#PBS -e error_file_oc_only_$x.e
#PBS -o output_file_oc_only_$x.o

# Change to the work directory
cd $PBS_O_WORKDIR

PN=$x
PARAM_FN=/cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/Cluster_PatternCompletionComplexPatternsPostLearning_nGlom40_nHC12_nMC30_vqOvrlp4_np10_OcOnly/Parameters/simulation_params.hoc
echo "Starting pattern $PN at `date`"
aprun -n 48 /cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/x86_64/special -mpi  -c "x=$PN" -c "strdef param_file" -c "sprint(param_file, \"$PARAM_FN\")" /cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/start_file_oc_only_recognition_task.hoc > delme_oc_only_$PN \
echo "Stopping pattern $PN at `date`"
