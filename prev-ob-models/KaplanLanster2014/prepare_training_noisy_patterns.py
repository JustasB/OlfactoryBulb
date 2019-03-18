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
import scipy.stats


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


def reorder_patterns(params, n_patterns_per_set, input_folders, mit_response_fn):
    """
    input_folders -- string of folder names containing the noisy data in ASCENDING order
    This function takes the noisy activation matrices from input_folders
    and re-orders them so that the same pattern with different degrees of noise 
    appears subsequent.

    Also, the mit_response_fn is reordered in the same way.
    """
    
    n_sets = len(input_folders)
    n_patterns_total = n_sets * n_patterns_per_set
    n_readout = params['n_readout']

    mit_response_unsorted = np.loadtxt(mit_response_fn)
    mit_response_sorted = np.zeros((n_patterns_total, params['n_mit']))
    activation_sorted = np.zeros((n_patterns_total, params['n_or']))


    # get the activation matrix files 
    new_pns = []
    for set_idx, folder in enumerate(input_folders):
        if set_idx == 0:
            fn = os.path.abspath(folder) + '/Parameters/activation_matrix.dat'
        else:
            fn = os.path.abspath(folder) + '/Parameters/activation_matrix_with_noise.dat'
        d = np.loadtxt(fn)
        for pn in xrange(n_patterns_per_set):
            pn_offset = pn * (n_sets - 1)
            new_pn = pn + set_idx + pn_offset
            print 'set_idx %d pn %d ----> %d' % (set_idx, pn, new_pn)
            new_pns.append(new_pn)
            activation_sorted[new_pn, :] = d[pn, :] 
            mit_response_sorted[new_pn, :] = mit_response_unsorted[pn, :]


    assert(np.unique(new_pns).size == n_patterns_total), 'Re-ordering failed'
    output_fn = params['params_folder'] + '/activation_matrix_reordered.dat'
    print 'Saving reordered action to:', output_fn
    np.savetxt(output_fn, activation_sorted)

    mit_response_output_fn = params['nspikes_folder'] + '/mit_response_normalized_reordered.dat'
    print 'Saving reordered mit response to:', mit_response_output_fn
    np.savetxt(mit_response_output_fn, mit_response_sorted)

    print 'Calculating pearson correlation coefficients...'
    pearsoncorr_act_matrix = np.zeros((n_patterns_total, n_patterns_total))
    pearsoncorr_mit_response = np.zeros((n_patterns_total, n_patterns_total))
    for i_ in xrange(n_patterns_total):
        for j_ in xrange(n_patterns_total):
            pearsoncorr_act_matrix[i_, j_] = scipy.stats.pearsonr(activation_sorted[i_, :], activation_sorted[j_, :])[0]
            pearsoncorr_mit_response[i_, j_] = scipy.stats.pearsonr(mit_response_sorted[i_, :], mit_response_sorted[j_, :])[0]

    output_fn = params['other_folder'] + '/pearsoncorr_activation_matrix.dat'
    print 'Saving reordered action to:', output_fn
    np.savetxt(output_fn, pearsoncorr_act_matrix)
    output_fn = params['other_folder'] + '/pearsoncorr_mitral_reponse.dat'
    print 'Saving reordered action to:', output_fn
    np.savetxt(output_fn, pearsoncorr_mit_response)
    return mit_response_output_fn








def create_normalized_mit_response(params, input_folders, n_patterns=None, copy=True):
    if n_patterns == None:
        n_patterns = params['n_patterns']
    # copy the mit spike files with a new pattern number to the new folder
    n_patterns_total = n_patterns * len(input_folders)
    if copy:
        for i_, folder in enumerate(input_folders):
            for pn_old in xrange(n_patterns):
                new_pn = i_ * n_patterns + pn_old
                old_fn = folder + '/Spiketimes/mit_spiketimes_merged_%d.dat' % (pn_old)
                new_fn = params['folder_name'] + '/Spiketimes/mit_spiketimes_merged_%d.dat' % (new_pn)
                cmd = 'cp %s %s' % (old_fn, new_fn)
                print cmd
                os.system(cmd)
                add_first_line(new_pn, new_fn, params['folder_name'] + '/Spiketimes/mit_spiketimes_merged_for_neuron_%d.dat' % (new_pn))

                old_fn = folder + '/NumberOfSpikes/mit_nspikes_merged_%d.dat' % (pn_old)
                new_fn = params['folder_name'] + '/NumberOfSpikes/mit_nspikes_merged_%d.dat' % (new_pn)
                cmd = 'cp %s %s' % (old_fn, new_fn)
                print cmd
                os.system(cmd)

    # create a new mitral cell response space from all the patterns combined
    mit_spikes = np.zeros((n_patterns_total, params['n_mit'])) # copy the number of spikes for each pattern into an array
    silent_mit = []
    for pn in xrange(n_patterns_total):
        fn = params['folder_name'] + '/NumberOfSpikes/mit_nspikes_merged_%d.dat' % (pn)
        print 'Loading', fn
        single_pattern_response = np.loadtxt(fn)
        if single_pattern_response.size > 0:
            for row in xrange(single_pattern_response[:, 0].size):
                gid = single_pattern_response[row, 0]
                cell_index = gid - params['mit_offset']
                mit_spikes[pn, cell_index] = single_pattern_response[row, 1]

    # 2a) Normalization: the sum of spikes fired by each cell during all patterns is normalized to 1 -> each mitral cell has a normalized activty
    normalized_activity = np.zeros((n_patterns_total, params['n_mit'])) # normalized_activity[:, cell_id].sum() = 1 for all cell_id (or = 0 if no spikes fired)
    for cell in xrange(params['n_mit']):
        spike_sum = mit_spikes[:, cell].sum() # sum of spikes for one cell over all patterns
        if (spike_sum == 0):
            normalized_activity[:, cell] = np.zeros(n_patterns_total)
            silent_mit.append(cell)
        else:
            normalized_activity[:, cell] = mit_spikes[:, cell] / spike_sum

    # 2b) the sum of normalized_activity within one glomerulus must not be > 1
    for pattern in xrange(n_patterns_total):
        for glom in xrange(params['n_mit_y']):
            # calculate which cells belong to one glomerulus
            id1 = glom * params['n_mit_x']
            id2 = (glom + 1) * params['n_mit_x']  
            # normalization within one glomerulus
            activity_sum = normalized_activity[pattern, id1:id2].sum()
            if (activity_sum > 1):
                normalized_activity[pattern, id1:id2] /= activity_sum

    if (len(silent_mit) > 0):
        silent_cells  = ['%d ' % mit for mit in silent_mit]
        silent_file = open(params['silent_mit_fn'], 'w')
        for s in silent_cells:
            silent_file.write(s)
        silent_file.close()

    print 'Silent mitral cells:', silent_mit
    print "AnalyseObOutput output file:", params["mit_response_normalized"]
    np.savetxt(params["mit_response_normalized"], normalized_activity)
    return params['mit_response_normalized']



if __name__ == '__main__':
    t1 = time.time()
#    input_folders = ['PaperData/Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/', 
#            'Cluster_ORnoise0.05_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/', 
#            'Cluster_ORnoise0.10_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/',
#            'Cluster_ORnoise0.15_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/',
#            'Cluster_ORnoise0.2_nGlom40_nHC12_nMC30_vqOvrlp4_np50_fullSystem/',
#            'Cluster_ORnoise0.25_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/',
#            'Cluster_ORnoise0.30_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/']
    input_folders = ['PaperData/Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/', 
            'Cluster_ORnoise0.05_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/', 
            'Cluster_ORnoise0.10_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/']


    n_pattern_sets = len(input_folders)
    n_patterns = 50
    n_patterns_total = n_pattern_sets * n_patterns
    param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params
    params['n_readout'] = n_patterns
    params['n_patterns'] = n_patterns_total
    param_tool.set_filenames()
    print 'New folder:', params['folder_name']
    ok = raw_input('\nContinue to create this folder structure? Parameters therein will be overwritten\n\ty / Y / blank = OK; anything else --> exit\n')
    if not ((ok == '') or (ok.capitalize() == 'Y')):
        print 'quit'
        exit(1)
    param_tool.write_parameters_to_file(params["info_file"])
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # 

#     first merge MIT spike files in all input_folders
#    for folder in input_folders:
#        cmd = 'python MergeSpikefiles.py %s mit' % folder
#        print cmd
#        os.system(cmd)

    mit_response_fn = create_normalized_mit_response(params, input_folders, n_patterns=n_patterns, copy=True)
#    mit_response_fn = create_normalized_mit_response(params, input_folders, n_patterns=n_patterns, copy=False)

    mit_response_fn = reorder_patterns(params, n_patterns, input_folders, mit_response_fn)


    CC.mds_vq_ob_output(params)

#    n_patterns_total = n_patterns * len(input_folders)
#    readout_activation = np.zeros((n_patterns_total, n_patterns))
#    for i_ in xrange(len(input_folders)):
#        idx0 = i_ * n_patterns
#        idx1 = (i_ + 1) * n_patterns
#        readout_activation[idx0:idx1, :] = np.eye(n_patterns)

#    CC.bcpnn_ob_oc(params)
#    CC.bcpnn_oc_oc(params)
#    CC.bcpnn_oc_readout(params, readout_activation)
#    CC.create_pyr_parameters(params)
#    CC.create_connections(params)

#    if params['oc_only']:
#        params['all_connection_fns'] = [ params['conn_list_mit_pyr'], \
#                                         params['conn_list_mit_rsnp'], \
#                                         params['conn_list_layer23'], \
#                                         params['conn_list_pyr_pyr'], \
#                                         params['conn_list_pyr_rsnp'], \
#                                         params['conn_list_pyr_readout']]
#    else:
#        for pn in xrange(params['n_patterns']):
#            fn = params['orn_params_fn_base'] + '%d.dat' % pn
#            assert os.path.exists(fn), 'Required file does not exist: %s' % fn

#    for fn in params['all_connection_fns']:
#        assert os.path.exists(fn), 'Required file does not exist: %s' % fn
#    for fn in [params['pyr_params_file'], params['readout_params_file']]:
#        assert os.path.exists(fn), 'Required file does not exist: %s' % fn


#    print 'Ready - Go!\t', params['folder_name']
#    t2 = time.time()
#    t_diff = t2 - t1
#    print 'Time: %.1f sec %.2f min' % (t_diff, t_diff / 60.)
