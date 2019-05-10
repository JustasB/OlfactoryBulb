info_txt = \
"""
    This script trains takes noisy patterns as input to train the full system.
    This script basically copies the relevant files from the noisy folders
    to the newly created with-noise folders.
"""

import os
import sys
import numpy as np
import simulation_parameters
import CreateObOcConnections as CC
import time
from prepare_training_noisy_patterns import create_normalized_mit_response


def add_first_line(pn, input_fn, output_fn):
    """"
    This function is required if OC only is run 
    Adds a nrow, ncol first line in the mit_spiketimes_merged_ files 
    """
    print 'Loading', input_fn
    d = np.loadtxt(input_fn)
    n_row, n_col = d.shape[0], d.shape[1]
    # read the current contents of the file
    f = open(input_fn)
    text = f.read()
    f.close()
    # open the file again for writing
    f = open(output_fn, 'w')
    f.write('%d %d\n' % (n_row, n_col))
    # write the original contents
    f.write(text)
    f.close()




if __name__ == '__main__':
    t1 = time.time()

    # get the mitral cell activation
    assert (len(sys.argv) > 1), 'Please give the name of the start system (system which has been run with different concentrations'
    input_folder = sys.argv[1]
    # e.g. Cluster_ORnoise0.0_nGlom40_nHC12_nMC30_vqOvrlp4_np50_fullSystem_ConcInv/

    param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params
    print 'New folder:', params['folder_name']
    ok = raw_input('\nContinue to create this folder structure? Parameters therein will be overwritten\n\ty / Y / blank = OK; anything else --> exit\n')
    if not ((ok == '') or (ok.capitalize() == 'Y')):
        print 'quit'
        exit(1)

    # copy the EPTH, OB parameter and connectivity files to the new folder
    srcs = []
    srcs.append('%sParameters/' % input_folder)
    srcs.append('%sConnections/' % input_folder)
    for src in srcs:
        cmd = 'cp -r %s %s' % (src, params['folder_name'])
        print cmd
        os.system(cmd)
    param_tool.write_parameters_to_file(params["info_file"])
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # 

    mit_response_normalized_fn = os.path.abspath(input_folder) + '/' + 'NumberOfSpikes/mit_response_normalized_np50'
    print 'From input folder: mit_response_normalized_fn', mit_response_normalized_fn
    if not os.path.exists(mit_response_normalized_fn):
        cmd = 'python MergeSpikefiles.py %s mit' % input_folder
        print cmd
        os.system(cmd)
        create_normalized_mit_response(params, [input_folder], n_patterns=params['n_patterns'], copy=True)
    CC.mds_vq_ob_output(params)

    # create the readout activation for the training
    params['n_patterns_test_conc_inv']
#    n_patterns_total = n_patterns * len(input_folders)
    readout_activation = np.zeros((params['n_readout'], params['n_patterns']))
    row = 0 
    for readout_idx in xrange(params['n_patterns_test_conc_inv']):
        for conc_ in xrange(params['n_conc_check']):
            readout_activation[row, readout_idx] = 1.
            row += 1

    CC.bcpnn_ob_oc(params)
    CC.bcpnn_oc_oc(params)
    CC.bcpnn_oc_readout(params, readout_activation)
    CC.create_pyr_parameters(params)
    CC.create_connections(params)

    if params['oc_only']:
        params['all_connection_fns'] = [ params['conn_list_mit_pyr'], \
                                         params['conn_list_mit_rsnp'], \
                                         params['conn_list_layer23'], \
                                         params['conn_list_pyr_pyr'], \
                                         params['conn_list_pyr_rsnp'], \
                                         params['conn_list_pyr_readout']]
    else:
        for pn in xrange(params['n_patterns']):
            fn = params['orn_params_fn_base'] + '%d.dat' % pn
            assert os.path.exists(fn), 'Required file does not exist: %s' % fn

    for fn in params['all_connection_fns']:
        assert os.path.exists(fn), 'Required file does not exist: %s' % fn
    for fn in [params['pyr_params_file'], params['readout_params_file']]:
        assert os.path.exists(fn), 'Required file does not exist: %s' % fn

    print 'Ready - Go!\t', params['folder_name']
    t2 = time.time()
    t_diff = t2 - t1
    print 'Time: %.1f sec %.2f min' % (t_diff, t_diff / 60.)
