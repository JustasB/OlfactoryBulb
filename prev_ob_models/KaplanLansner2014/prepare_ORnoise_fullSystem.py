info_txt = \
"""
    Before running this script, you should
    - have run a full system without noise
    i.e. 
    - you have the simulation results of a 'pre-learning' run 
    - and run CreateObOcConnections.py and verified that (at least some patterns get recognized)

    This script basically copies the relevant files from the noise-free folders
    to the newly created with-noise folders
"""

import os
import sys
import numpy as np
import simulation_parameters # defines simulation parameters
import CreateOrnParameters
# classes for setting up connectivity and the individual cell parameters



if __name__ == '__main__':
    print info_txt
    # ------------ I N I T -----------------------------
    # The simulation_parameters module defines a class for simulation parameter storage
    param_tool = simulation_parameters.parameter_storage()
    # params is the dictionary with all parameters
    params = param_tool.params

    assert (len(sys.argv) > 1), 'Please give the name of the noise-free system (including all parameters to run a full-system simulation including EPTH, OB, OC)'
    # you want to create a new folder structure to not overwrite the results obtained from 
    # presenting patterns to EPTH and OB --> 'preLearning' results
    # it's convenient to write results into a seperate folder structure as OB->OC connectivity
    # and pattern recognition might require some tuning (as every model does ;))
    print 'New folder:', params['folder_name']

    ok = raw_input('\nContinue to create this folder structure? Parameters therein will be overwritten\n\ty / Y / blank = OK; anything else --> exit\n')
    if not ((ok == '') or (ok.capitalize() == 'Y')):
        print 'quit'
        exit(1)

    default_folder = os.path.abspath(sys.argv[1]) + '/' # this should include all the data ready to run a full system
    print 'Taking data from:', default_folder 

    srcs = []
    srcs.append('%sParameters/' % default_folder)
    srcs.append('%sConnections/' % default_folder)
    for src in srcs:
        cmd = 'cp -r %s %s' % (src, params['folder_name'])
        print cmd
        os.system(cmd)

    param_tool.write_parameters_to_file(params["info_file"])
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # 

    default_activation_matrix_fn = default_folder + 'Parameters/activation_matrix.dat'
    OrnParamClass = CreateOrnParameters.CreateOrnParameters(params) # patterns for ORN activation must be recreated to add noise
    default_activation_matrix = '%sParameters/activation_matrix.dat' % default_folder
    ok = OrnParamClass.create_single_odorant_patterns(given_activation_matrix=default_activation_matrix)

    for fn in params['all_connection_fns']:
        assert os.path.exists(fn), 'Required dile does not exist: %s' % fn
    for fn in [params['pyr_params_file'], params['readout_params_file']]:
        assert os.path.exists(fn), 'Required dile does not exist: %s' % fn

    for pn in xrange(params['n_patterns']):
        fn = params['orn_params_fn_base'] + '%d.dat' % pn
        assert os.path.exists(fn)

    print 'Ready - Go!\t', params['folder_name']
