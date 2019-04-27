import os
import sys
import time
import numpy
import simulation_parameters # defines simulation parameters
import CreateObConnections
import CreateOrnParameters

param_tool = simulation_parameters.parameter_storage()
# params is the dictionary with all parameters
params = param_tool.params
param_tool.write_parameters_to_file(params["info_file"]) # human readable
param_tool.write_parameters_to_file() # 
param_tool.hoc_export() # write the simulation parameters to a NEURON executable file
param_tool.print_cell_gids()

prepare = True
if prepare:
    t1 = time.time()
    ConnectionClass = CreateObConnections.CreateObConnections(params)
    ConnectionClass.connect_orn_mit()
    ConnectionClass.connect_orn_pg()
    ConnectionClass.connect_pg_mit_serial()
    ConnectionClass.connect_pg_mit_reciprocal()
    ConnectionClass.connect_mt_gran_local()
    ConnectionClass.connect_mt_gran_global()

    print "Creating ORN parameters...."
    OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
    ok = OrnParamClass.create_single_odorant_patterns()

    t2 = time.time()
    print "Time: %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)


t_all_0 = time.time()
#for pn in xrange(params['n_patterns']):
#    os.system("rm %s/*" % (params["spiketimes_folder"]))
#    os.system("rm %s/*" % (params["volt_folder"]))

#for pn in [0]:
for pn in xrange(params['n_patterns']):
    t1 = time.time()
    param_file_orns = params['orn_params_fn_base'] + '%d.dat' % pn
    assert os.path.exists(param_file_orns), 'File does not exist: %s\nPlease run prepare_epth_response_curve.py before!' % param_file_orns

    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_epth_ob_prelearning.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)

    os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
    print 'NEURON simulation starts for pattern %d ...' % pn
    os.system(neuron_command)

    t2 = time.time() - t1
    print "Simulating %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_cells_ob'] , params["t_sim"], t2, t2/60.)

    os.chdir('../') # this is important to avoide problems with tabchannel files and the functions defined therein

t_all_1 = time.time() - t_all_0
print "Simulating %d patterns %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_patterns'], params['n_cells_ob'], params["t_sim"], t_all_1, t_all_1/60.)
