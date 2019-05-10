import os
import sys
import time
import numpy
import simulation_parameters # defines simulation parameters

param_tool = simulation_parameters.parameter_storage()
# params is the dictionary with all parameters
params = param_tool.params
param_tool.write_parameters_to_file(params["info_file"]) # human readable
param_tool.write_parameters_to_file() # 
param_tool.hoc_export() # write the simulation parameters to a NEURON executable file
param_tool.print_cell_gids()

t_all_0 = time.time()
for pn in xrange(params['n_patterns']):
    os.system("rm %s/*" % (params["spiketimes_folder"]))
    os.system("rm %s/*" % (params["volt_folder"]))


#for pn in [0]:
for pn in xrange(params['n_patterns']):
    t1 = time.time()

#    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
#            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_oc_only_recognition_task.hoc "\
#            % (params['n_proc'], pn, params['hoc_file'])
#    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
#            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_oc_only_recognition_task.hoc > delme%d" \
#            % (params['n_proc'], pn, params['hoc_file'], pn)

    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_full_system_recognition_task.hoc" \
            % (params['n_proc'], pn, params['hoc_file'])

    os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
    print 'NEURON simulation starts for pattern %d ...' % pn
    os.system(neuron_command)

    t2 = time.time() - t1
    print "Simulating %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_cells_oc'] , params["t_sim"], t2, t2/60.)

    os.chdir('../') # this is important to avoide problems with tabchannel files and the functions defined therein

t_all_1 = time.time() - t_all_0
print "Simulating %d patterns %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_patterns'], params['n_cells_oc'], params["t_sim"], t_all_1, t_all_1/60.)
