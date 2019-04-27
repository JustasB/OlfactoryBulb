"""
Make sure that prepare_epth_response_curve.py was run successfully before
calling this script.
"""

import simulation_parameters
import os
import time
import MergeSpikefiles
import SetOfCurvesPlotter
import CreateOrnParameters

param_tool = simulation_parameters.parameter_storage()
params = param_tool.params
param_tool.hoc_export()


prepare = True
if prepare:
    OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
    OrnParamClass.create_params_for_response_curve()

t1 = time.time()
sim_cnt = 0
pn = sim_cnt

param_file_orns = params['orn_params_fn_base'] + '%d.dat' % pn
assert os.path.exists(param_file_orns), 'File does not exist: %s\nPlease run prepare_epth_response_curve.py before!' % param_file_orns

if params['with_artificial_orns']:
    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_artificial_orn_spikes.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)
else:
    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_epth_response_curve.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)

param_tool.print_cell_gids()

# clean up data folder
os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
os.system("rm %s/*" % (params["spiketimes_folder"]))
os.system("rm %s/*" % (params["volt_folder"]))
print 'NEURON simulation starts ...'
os.system(neuron_command)

t2 = time.time() - t1

print "Simulating %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_orn'], params["t_sim"], t2, t2/60.)

os.chdir('../')

# ------- A N A L Y S I S --------------------
Merger = MergeSpikefiles.MergeSpikefiles(params)
Merger.merge_epth_spiketimes_file(pattern=pn)
Merger.merge_epth_nspike_files(pattern=pn)

SOCP = SetOfCurvesPlotter.SetOfCurvesPlotter(params)
output_fn = params['figure_folder'] + '/hand_tuned_%d.png' % sim_cnt
#print 'Saving figure to:', output_fn
x_data, y_data = SOCP.plot_set_of_curves(pn=sim_cnt, output_fn=output_fn)

print 'Opening with ristretto: %s' % (output_fn)
os.system('ristretto %s' % output_fn)
