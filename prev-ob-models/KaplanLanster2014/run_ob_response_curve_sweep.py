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
import numpy as np
import CreateObConnections
import CreateMitParameters

sim_cnt_offset = 0

n_runs = 1

orn_mit_change_ids = range(8)
orn_mit_change_factors = np.ones((n_runs, len(orn_mit_change_ids)))
#                                0    1    2    3    4    5   6   7                                     
#orn_mit_change_factors[0, :] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#orn_mit_change_factors[0, :] = [2.0, 2.0, 1.8, 1.8, 1.3, 1.1, 1.0, 1.0]
#orn_mit_change_factors[0, :] = [2.5, 2.5, 1.8, 1.8, 1.3, 1.1, 1.0, 1.0]


orn_mit_pg_mult = [2., 3., 5., 6.]
n_runs = len(orn_mit_pg_mult)

for sim_cnt in xrange(n_runs):
    sim_cnt += sim_cnt_offset
    param_tool = simulation_parameters.parameter_storage()
    param_too.params['w_orn_mit_mult'] = orn_mit_pg_mult[sim_cnt]
    param_too.params['w_orn_pg_mult'] = orn_mit_pg_mult[sim_cnt]
#    param_tool.params['orn_mit_change_ids'] = orn_mit_change_ids
#    param_tool.params['orn_mit_change_factors'] = orn_mit_change_factors[sim_cnt, :].tolist()
    param_tool.hoc_export()
    params = param_tool.params
    param_fn = params['params_fn_json'].rsplit('.')[0] + '_%d.json' % (sim_cnt)
    param_tool.write_parameters_to_file(param_fn)

    # prepare
    OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
    OrnParamClass.create_params_for_response_curve()
    ConnectionClass = CreateObConnections.CreateObConnections(params)
    ConnectionClass.connect_orn_mit()
    ConnectionClass.connect_orn_pg()
    ConnectionClass.connect_pg_mit_serial()
    ConnectionClass.connect_pg_mit_reciprocal()
    ConnectionClass.connect_mt_gran_local()
    ConnectionClass.connect_mt_gran_global()
    # this is required by SetOfCurvesPlotter to store the information which ORN is connected to which MIT
    MitParamClass = CreateMitParameters.CreateMitParameters(params)
    ok = MitParamClass.create_parameters()

    t1 = time.time()
    pn = 0
    param_file_orns = params['orn_params_fn_base'] + '%d.dat' % pn
    assert os.path.exists(param_file_orns), 'File does not exist: %s\nPlease run prepare_ob_response_curve.py before!' % param_file_orns

    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_ob_response_curve.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)

    param_tool.print_cell_gids()
    print 'Running sim %d with out of %d' % (sim_cnt, n_runs), orn_mit_change_factors[sim_cnt, :]

    # clean up data folder
    os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
    os.system("rm %s/*" % (params["spiketimes_folder"]))
    os.system("rm %s/*" % (params["volt_folder"]))
    print 'NEURON SIM:\n', neuron_command
    os.system(neuron_command)
    os.chdir('../')

    t2 = time.time() - t1

    print "Simulating %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_orn'], params["t_sim"], t2, t2/60.)


    # ------- A N A L Y S I S --------------------
    Merger = MergeSpikefiles.MergeSpikefiles(params)
    Merger.merge_ob_spiketimes_file(pattern=pn)
    Merger.merge_ob_nspike_files(pattern=pn)

    SOCP = SetOfCurvesPlotter.SetOfCurvesPlotter(params)
    output_fn = params['figure_folder'] + '/ob_response_curve_%d.png' % sim_cnt
    print 'Saving figure to:', output_fn
    x_data, y_data = SOCP.plot_set_of_curves(output_fn=output_fn, cell_type='mit')
    data_fn = params['other_folder'] + '/response_curve_x_%d.dat' % (sim_cnt)
    header = '# %s\n' % (str(orn_mit_change_factors[sim_cnt, :]))
    f = open(data_fn, 'w')
    f.write(header)
    f.flush()
    print 'Saving data to:', data_fn
    np.savetxt(f, x_data)

    data_fn = params['other_folder'] + '/response_curve_y_%d.dat' % (sim_cnt)
    header = '# %s\n' % (str(orn_mit_change_factors[sim_cnt, :]))
    f = open(data_fn, 'w')
    f.write(header)
    f.flush()
    print 'Saving data to:', data_fn
    np.savetxt(f, y_data)


#    print 'Opening with ristretto: %s' % (output_fn)
#    os.system('ristretto %s' % output_fn)
