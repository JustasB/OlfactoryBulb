# The name of the script is neuron_job
#PBS -N epth_ob_prel

# 4 hour wall-clock time will be given to this job
#PBS -l walltime=0:20:00
# delm PBS -l walltime=0:30:00
# a response curve sweep takes about 6 min

# Number of cores to be allocated is 24
#PBS -l mppwidth=384

#PBS -e error_file.e
#PBS -o output_file.o

# Change to the work directory
cd $PBS_O_WORKDIR

#PN=0
#PN=$x
# Run the executable named myexe and write the output into my_output_file
# aprun -n 24 ./i686/special -mpi test0.hoc > my_output_file 2>&1
# aprun -n 24 /cfs/emil/pdc/bkaplan/neuron/nrnmpi/x86_64/bin/nrniv -mpi /cfs/emil/pdf/bkaplan/test/test0.hoc > my_output_file 2>&1

#echo "Starting pattern $PN at `date`"
#rm current_pattern_nr
#echo "strdef param_file" > bullshit2
#echo "sprint(param_file, \"/cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/Cluster_ResponseCurvesEpthOb/Parameters/simulation_params.hoc\")" > bullshit3
#cat bullshit2 bullshit3 > current_simulation_params
#aprun -n 24 /cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/x86_64/special -mpi  -c "x=$PN" /cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/start_file_ob_response_curve.hoc > delme_all3_$PN 2>&1
#echo "Stopping pattern $PN at `date`"

PN=$x
PARAM_FN=/cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/Cluster_PaperData_nGlom40_nHC25_nMC16_rORN100_ORnoise0.0/Parameters/simulation_params.hoc
echo "Starting pattern $PN at `date`"
aprun -n 384 /cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/x86_64/special -mpi  -c "x=$PN" -c "strdef param_file" -c "sprint(param_file, \"$PARAM_FN\")" /cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files/start_file_epth_ob_prelearning.hoc > delme_$PN \
echo "Stopping pattern $PN at `date`"
