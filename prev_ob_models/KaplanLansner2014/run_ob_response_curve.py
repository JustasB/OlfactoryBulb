import simulation_parameters
import os
import time
import MergeSpikefiles
import SetOfCurvesPlotter

param_tool = simulation_parameters.parameter_storage()
params = param_tool.params
param_tool.hoc_export()


"""
Make sure that prepare_ob_response_curve.py was run successfully before
calling this script.
"""

prepare = True
if prepare:
    os.system('python prepare_ob_response_curve.py')

t1 = time.time()
pn = 0

param_file_orns = params['orn_params_fn_base'] + '%d.dat' % pn
assert os.path.exists(param_file_orns), 'File does not exist: %s\nPlease run prepare_ob_response_curve.py before!' % param_file_orns

if params['with_artificial_orns']:
    # NOT YET IMPLEMENTED
    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_artificial_orn_spikes.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)
else:
    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_ob_response_curve.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)


param_tool.print_cell_gids()

# change folder and remove old data files to avoid possible confusion with old data
os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
os.system("rm %s/*" % (params["spiketimes_folder"]))
os.system("rm %s/*" % (params["volt_folder"]))
print 'Running the NEURON simulation ...'
print neuron_command
os.system(neuron_command)
t2 = time.time() - t1
print "Simulating %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_orn'], params["t_sim"], t2, t2/60.)
os.chdir('../')

# ------- A N A L Y S I S --------------------
Merger = MergeSpikefiles.MergeSpikefiles(params)
Merger.merge_ob_spiketimes_file(pattern=pn)
Merger.merge_ob_nspike_files(pattern=pn)

sim_cnt = 1

SOCP = SetOfCurvesPlotter.SetOfCurvesPlotter(params)
output_fn = params['figure_folder'] + '/ob_response_curve_%d.png' % sim_cnt
SOCP.plot_set_of_curves(output_fn=output_fn, cell_type='mit')
print 'Opening with ristretto: %s' % (output_fn)
os.system('ristretto %s' % output_fn)
