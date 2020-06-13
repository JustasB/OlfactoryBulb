info_txt = \
"""
    This script is basically identical to prepare_pattern_completion
    but the training is with simple patterns (less OR activated)
    and the testing is with complex patterns (a superposition of training patterns)
"""

import os
import sys
import numpy as np
import simulation_parameters # defines simulation parameters
import CreateOrnParameters
import CreateObConnections
import CreateObOcConnections as CC
import AnalyseObOutput
from run_oc_only import add_first_line


def prepare_epth_ob_conn(params):

    # EPTH -> OB connections are not affected by the pattern
    print "Creating connections: orn -> mit"
    ConnectionClass = CreateObConnections.CreateObConnections(params)
    ConnectionClass.connect_orn_mit()
    ConnectionClass.connect_orn_pg()
    ConnectionClass.connect_pg_mit_serial()
    ConnectionClass.connect_pg_mit_reciprocal()
    ConnectionClass.connect_mt_gran_local()
    ConnectionClass.connect_mt_gran_global()


if __name__ == '__main__':
    print info_txt
    # ------------ I N I T -----------------------------
    # The simulation_parameters module defines a class for simulation parameter storage
    param_tool = simulation_parameters.parameter_storage()
    # params is the dictionary with all parameters
    params = param_tool.params

    print 'New folder:', params['folder_name']
    ok = raw_input('\nContinue to create this folder structure? Parameters therein will be overwritten\n\ty / Y / blank = OK; anything else --> exit\n')
    if not ((ok == '') or (ok.capitalize() == 'Y')):
        print 'quit'
        exit(1)

    OrnParamClass = CreateOrnParameters.CreateOrnParameters(params) # patterns for ORN activation must be recreated to add noise
    if params['train_pattern_completion']:
        n_active_OR = int(round(params['n_or'] * params['frac_min_active_OR']))
        OrnParamClass.create_pattern_completion_training(n_active_OR)
        required_connection_fns = [params['conn_list_orn_mit'], \
                            params['conn_list_orn_pg'], \
                            params['conn_list_pg_mit_serial'], \
                            params['conn_list_pg_mit_reciprocal'], \
                            params['conn_list_mit_pg_reciprocal'], \
                            params['conn_list_mit_pg_serial'], \
                            params['conn_list_mit_gran_local'], \
                            params['conn_list_mit_gran_global'], \
                            params['conn_list_gran_mit_local'], \
                            params['conn_list_gran_mit_global']]

        # the EPTH, OB activation will be composite patterns
        cmd = 'cp %s/NumberOfSpikes/mit_* %s' % (os.path.abspath(training_folder), params['nspikes_folder'])
        print cmd
        os.system(cmd)

        # training of OB-OC, OC-OC, OC-Readout connections will be based on the simple patterns
        # copy the mitral cell files
        ObAnalyser = AnalyseObOutput.AnalyseObOutput(params)
        ObAnalyser.get_output_activity() # optional, checks if activity for some patterns is missing
        ObAnalyser.get_mit_response_normalized()
        CC.mds_vq_ob_output(params)
        CC.bcpnn_ob_oc(params)
        CC.bcpnn_oc_oc(params)
        readout_activation = np.eye(params['n_patterns'])
        CC.bcpnn_oc_readout(params, readout_activation)
        CC.create_pyr_parameters(params)
        CC.create_connections(params)

    elif (params['test_pattern_completion']):
        assert (len(sys.argv) > 1), 'Please give the filename for the \'simple\' activation patterns with which the system has been trained'
        training_folder = sys.argv[1] # contains the EPTH and OB activity of simple patterns
        print 'Training folder:', training_folder

        # if you want to test on complex patterns, which are a composite of 'simpler' patterns with which the system has
        # been trained, then set oc_only to False!
        if params['oc_only']:
            cmd = 'cp %s/Connections/* %s' % (training_folder, params['conn_folder'])
            print cmd
            os.system(cmd)
            cmd = 'cp %s/Parameters/* %s' % (training_folder, params['params_folder'])
            print cmd
            os.system(cmd)
            # run with the mitral cell activity that is stored in training_folder
            cmd = 'cp %s/Spiketimes/mit_spiketimes_merged_* %s' % (training_folder, params['spiketimes_folder'])
            print cmd
            os.system(cmd)
            for pn in xrange(params['n_patterns']):
                add_first_line(pn)
            required_connection_fns = [ params['conn_list_mit_pyr'], \
                                        params['conn_list_mit_rsnp'], \
                                        params['conn_list_layer23'], \
                                        params['conn_list_pyr_pyr'], \
                                        params['conn_list_pyr_rsnp'], \
                                        params['conn_list_pyr_readout']]
        else:
#            prepare_epth_ob_conn(params)
            # full system
            cmd = 'cp %s/Connections/* %s' % (training_folder, params['conn_folder'])
            print cmd
            os.system(cmd)
            cmd = 'cp %s/Parameters/* %s' % (training_folder, params['params_folder'])
            print cmd
            os.system(cmd)
            activation_matrix_fn = os.path.abspath(training_folder) + '/Parameters/activation_matrix.dat'
            assert (os.path.exists(activation_matrix_fn)), 'Required activation matrix file not found: %s' % activation_matrix_fn

            # create simpler patterns created out of the training training patterns
            OrnParamClass.create_pattern_completion_test(activation_matrix_fn)
            required_connection_fns = params['all_connection_fns']
            for pn in xrange(params['n_patterns']):
                fn = params['orn_params_fn_base'] + '%d.dat' % pn
                assert os.path.exists(fn)

        for fn in [params['pyr_params_file'], params['readout_params_file']]:
            assert os.path.exists(fn), 'Required dile does not exist: %s' % fn

        if not params['oc_only']:
            # remove the mitral cell spike files that were necessary for training the system
            # but need to be removed in order to not confuse with the test results for full_system run
            cmd = 'rm %s/mit_* ' % (params['nspikes_folder'])
            print cmd
            os.system(cmd)

    for fn in required_connection_fns:
        assert os.path.exists(fn), 'Required dile does not exist: %s' % fn

    param_tool.write_parameters_to_file(params["info_file"])
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # 

    print 'Los, Los, Los!\t', params['folder_name']
