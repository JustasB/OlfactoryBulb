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


param_sets = [\
                [5.0e-05, -8.0e+05, 1.0e-04], \
                [5.0e-05, -7.0e+05, 1.0e-04], \
                [5.0e-05, -6.0e+05, 1.0e-04], \
                [6.0e-05, -8.0e+05, 1.0e-04], \
                [6.0e-05, -7.0e+05, 1.0e-04], \
                [6.0e-05, -6.0e+05, 1.0e-04], \
                [7.0e-05, -8.0e+05, 1.0e-04], \
                [7.0e-05, -7.0e+05, 1.0e-04], \
                [7.0e-05, -6.0e+05, 1.0e-04],  \

                [5.0e-05, -8.0e+05, 5.0e-05], \
                [5.0e-05, -7.0e+05, 5.0e-05], \
                [5.0e-05, -6.0e+05, 5.0e-05], \
                [6.0e-05, -8.0e+05, 5.0e-05], \
                [6.0e-05, -7.0e+05, 5.0e-05], \
                [6.0e-05, -6.0e+05, 5.0e-05], \
                [7.0e-05, -8.0e+05, 5.0e-05], \
                [7.0e-05, -7.0e+05, 5.0e-05], \
                [7.0e-05, -6.0e+05, 5.0e-05],  \

                [5.0e-05, -8.0e+05, 1.0e-05], \
                [5.0e-05, -7.0e+05, 1.0e-05], \
                [5.0e-05, -6.0e+05, 1.0e-05],  \
                [6.0e-05, -8.0e+05, 1.0e-05], \
                [6.0e-05, -7.0e+05, 1.0e-05], \
                [6.0e-05, -6.0e+05, 1.0e-05], \
                [7.0e-05, -8.0e+05, 1.0e-05], \
                [7.0e-05, -7.0e+05, 1.0e-05], \
                [7.0e-05, -6.0e+05, 1.0e-05],  \

                [5.0e-05, -8.0e+05, 1.0e-06], \
                [5.0e-05, -7.0e+05, 1.0e-06], \
                [5.0e-05, -6.0e+05, 1.0e-06], \
                [6.0e-05, -8.0e+05, 1.0e-06], \
                [6.0e-05, -7.0e+05, 1.0e-06], \
                [6.0e-05, -6.0e+05, 1.0e-06], \
                [7.0e-05, -8.0e+05, 1.0e-06], \
                [7.0e-05, -7.0e+05, 1.0e-06], \
                [7.0e-05, -6.0e+05, 1.0e-06]  \
                ]

param_name = 'gleak_params'

sim_cnt_offset = 0
for sim_cnt, param_set in enumerate(param_sets):
    sim_cnt += sim_cnt_offset
    param_tool = simulation_parameters.parameter_storage()
    param_tool.update_values({param_name : param_set})
    param_tool.hoc_export()
    params = param_tool.params
    param_fn = params['params_fn_json'].rsplit('.')[0] + '_%d.json' % (sim_cnt)
    param_tool.write_parameters_to_file(param_fn)

    prepare = True
    if prepare:
        OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
        OrnParamClass.create_params_for_response_curve()

    t1 = time.time()
    pn = 0

    param_file_orns = params['orn_params_fn_base'] + '%d.dat' % pn
    assert os.path.exists(param_file_orns), 'File does not exist: %s\nPlease run prepare_epth_response_curve.py before!' % param_file_orns

    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_epth_response_curve.hoc > delme%d" \
            % (params['n_proc'], pn, params['hoc_file'], pn)

    param_tool.print_cell_gids()
    print 'Running with ', param_name, param_set

    # clean up data folder
    os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
    os.system("rm %s/*" % (params["spiketimes_folder"]))
    os.system("rm %s/*" % (params["volt_folder"]))
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
    SOCP.plot_set_of_curves(output_fn)

#    print 'Opening with ristretto: %s' % (output_fn)
#    os.system('ristretto %s' % output_fn)
