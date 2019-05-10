info_txt = \
"""
    Before running this script, you should have:
    - the simulation results of a 'pre-learning' run 

    - chosen a new folder name and run this script 
    only creating the new folder structure and parameter files by running (see below)
        param_tool.write_parameters_to_file() # 
    - then, copy the nspike files of mitral cells into the new NumberOfSpikes folder:
     cp PRE_LEARNING/NumberOfSpikes/mit*merged* POST_LEARNING/NumberOfSpikes/

    Usage:
        python CreateObConnections.py new PRE_LEANING_FOLDER
"""

import os
import sys
import numpy as np
import simulation_parameters # defines simulation parameters
import prepare_epth_ob_prelearning
import MergeSpikefiles
import MDSVQ
import BCPNN
import AnalyseObOutput
import CreatePyrReadoutParameters
import CreateHyperAndMinicolumns
import TransformAbstractToDetailedConnectivity
# classes for setting up connectivity and the individual cell parameters


from run_oc_only import add_first_line


def choose_initial_centroids_for_vq_ob_oc(params, input_fn, n_hc, mit_cell_coordinates_in_mi_space):
    import random
    random.seed(1)
    points = np.loadtxt(input_fn)

    rnd_points = []
    if os.path.exists(params['silent_mit_fn']):
        silent_mitral_cells = set(np.loadtxt(params['silent_mit_fn']))
        all_mit = set(range(params['mit_offset'], params['mit_offset'] + params['n_mit']))
        non_silent_mitral_cells = all_mit.difference(silent_mitral_cells)
    else:
        non_silent_mitral_cells = range(params['mit_offset'], params['mit_offset'] + params['n_mit'])

    while len(rnd_points) != n_hc:
        candidate = random.choice(non_silent_mitral_cells)
        while candidate in rnd_points:
            candidate = random.choice(non_silent_mitral_cells)
        rnd_points.append(candidate - params['mit_offset'])
    
    print 'rnd_points:', rnd_points
    rnd_centroids = mit_cell_coordinates_in_mi_space[rnd_points, :]
    
    return rnd_centroids




def mds_vq_ob_output(params, mds_output_fn=None):

    mdsvq = MDSVQ.MDSVQ(params)

    ob_activity_fn = params['mit_mds_input_fn'] # choose one of the different MIT - output files (with different normalizations)
    assert (os.path.exists(ob_activity_fn)), 'ERROR: File does not exist %s\n\t Have you created / copied the (normalized) mit_response file properly?\n' % (ob_activity_fn)
    # 1) calculate mutual information (MI) between MIT cells and their distances
    # 2) run MDS in this MI-space (OB pattern reponse or mutual information (MI) space) and save mitral cell coordinates
    if mds_output_fn == None:
        mds_output_fn = params['mds_ob_oc_output_fn']
        cell_type = 'mit'
        mit_cell_coordinates_in_mi_space = mdsvq.mds(ob_activity_fn, mds_output_fn, thresh=1e-6, cell_type=cell_type)
    else:
        mit_cell_coordinates_in_mi_space = np.loadtxt(mds_output_fn)
    print 'MDS-VQ-OB output:', mds_output_fn

    # 3) run a VQ in the MI space, with n_clusters = n_hc
    # i.e. assign one or more hypercolumns to each mitral cell
    vq_output_fn = params['vq_ob_oc_output_fn']

#    guessed_centroids = choose_initial_centroids_for_vq_ob_oc(params, mds_output_fn, params['n_hc'], mit_cell_coordinates_in_mi_space)
#    print 'guessed_centroids', guessed_centroids
#    np.random.seed(123)
#    guessed_centroids = 2 * np.random.random((params['n_hc'], 3))
    guessed_centroids = params['n_hc']
    mdsvq.vq(mds_output_fn, vq_output_fn, guessed_centroids, overlap=params['vq_ob_oc_overlap'], remove_silent_cells_fn=params['silent_mit_fn'], show=False)
#    mdsvq.vq(mds_output_fn, vq_output_fn, params['n_hc'], overlap=params['vq_ob_oc_overlap'], remove_silent_cells_fn=params['silent_mit_fn'])

    # 4) For each hypercolumn, create a new space spanned by the mitral cells projecting to the hc
    #   Each mitral cell represents one dimension and each pattern represents one vector or point in that space.
    #   The value of one component in such a pattern-vector is equal to the normalized activation of the respective mitral cell.
    #   The n_patterns vectors are clustered by VQ among the minicolumns in the hypercolumn.
    binary_oc_activation_fn = params['binary_oc_activation_fn']
    mit_mc_kmeans_trial = 0
    mdsvq.create_mitral_response_space(vq_output_fn, ob_activity_fn, binary_oc_activation_fn, \
            remove_silent_cells_fn=os.path.exists(params['silent_mit_fn']), mit_mc_kmeans_trial=mit_mc_kmeans_trial) 
    del mdsvq
    return (ob_activity_fn, binary_oc_activation_fn, vq_output_fn)


def bcpnn_ob_oc(params):
    print 'BCPNN OB -> OC'

    n_mit = params['n_mit']

    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    bcpnn = BCPNN.BCPNN(1, n_mit, n_hc, n_mc, n_patterns, params) # 1 src HC, n_mit MC --> n_hc, n_mc

    ob_activity_fn = params['mit_mds_input_fn']
    activity_fn = params['oc_abstract_activity_fn']
    
    bias_fn = params['ob_oc_abstract_bias_fn']
    binary_oc_activation_fn = params['binary_oc_activation_fn']
    w_ij_mit_hc = params['vq_ob_oc_output_fn']
    weights_fn = params['ob_oc_abstract_weights_fn']
    print "BCPNN OB -> OC loading files:", ob_activity_fn


    bcpnn.load_input_activity(ob_activity_fn)
    bcpnn.load_output_activity(binary_oc_activation_fn)
    bcpnn.load_mc_hc_mask(w_ij_mit_hc, silent_units_fn=os.path.exists(params['silent_mit_fn']))
#    bcpnn.initialize()

    n_steps = params['n_bcpnn_steps']
    for i in xrange(n_steps):
    #    bcpnn.train_network()
        bcpnn.train_network()#activity_fn, weights_fn, bias_fn)
    #bcpnn.train_network(activity_fn, weights_fn, bias_fn)

#    if os.path.exists(params['silent_mit_fn']):
#        bcpnn.silence_mit(params['silent_mit_fn'])
#    print 'bcpnn_ob_oc: input = %s\tweights = %s\tbias = %s\ttest_output = %s' % (test_input, weights_fn, bias_fn, test_output)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)
    del bcpnn


def bcpnn_oc_oc(params):
    print 'BCPNN OC -> OC'

    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    n_readout = params['n_readout']
    bcpnn = BCPNN.BCPNN(n_hc, n_mc, n_hc, n_mc, n_patterns, params)

    # train with the recurrent connections with OC output activity after learning OB -> OC
    #oc_activity_fn = params['oc_abstract_activity_fn']
    #bcpnn.load_input_activity(oc_activity_fn)
    #bcpnn.load_output_activity(oc_activity_fn)
    #bcpnn.initialize()

    # train with binary oc activation derived from WTA after 2nd VQ
    oc_oc_training_fn = params['binary_oc_activation_fn']
    # train with the output activity when learning the ob-oc connections
#    oc_oc_training_fn = params['oc_abstract_activity_fn']

    print 'BCPNN OC <-> OC training with:', oc_oc_training_fn
    bcpnn.load_input_activity(oc_oc_training_fn)
    bcpnn.load_output_activity(oc_oc_training_fn)

#    bcpnn.initialize()

    activity_fn = params['oc_oc_abstract_activity_fn'] # for output 
    weights_fn = params['oc_oc_abstract_weights_fn']
    bias_fn = params['oc_oc_abstract_bias_fn']
    n_steps = params['n_bcpnn_steps']
    for i in xrange(n_steps):
    #    bcpnn.train_network()
        bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)
    del bcpnn



def bcpnn_oc_readout(params, readout_activation=None):
    print 'BCPNN OC -> READOUT'
    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    if readout_activation == None:
        # same number of patterns and readout cells
        n_readout = params['n_patterns']
        readout_activation = np.eye(n_patterns)
    else: # e.g. noisy patterns
        n_readout = params['n_readout']

    bcpnn = BCPNN.BCPNN(n_hc, n_mc, 1, n_readout, n_patterns, params)
    # take the output activity of OB-OC Bcpnn as activity
    oc_activity_fn = params['oc_abstract_activity_fn']
    # take the output activity of OC-OC Bcpnn as activity
#    oc_activity_fn = params['oc_oc_abstract_activity_fn']
    print 'Loading as input activity', oc_activity_fn
    bcpnn.load_input_activity(oc_activity_fn)

    bcpnn.load_output_activity(readout_activation)
#    bcpnn.initialize()

    activity_fn = params['readout_abstract_activity_fn']
    weights_fn = params['oc_readout_abstract_weights_fn']
    bias_fn = params['oc_readout_abstract_bias_fn']

    #n_steps = 1
    n_steps = params['n_bcpnn_steps']
    for i in xrange(n_steps):
    #    bcpnn.train_network()
        bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    #bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)

    #if params['multiple_concentrations_per_pattern']:
    #    n_patterns = 10
    # testing 
    #del bcpnn
    bcpnn = BCPNN.BCPNN(n_hc, n_mc, 1, n_readout, n_patterns, params)
    test_input = oc_activity_fn
    test_output = params['readout_abstract_activity_fn'].rsplit('.dat')[0] + '_test.dat'
#    print 'BCPNN.testing(input = %s \nweights = %s \nbias = %s \ntest_output = %s' % (test_input, weights_fn, bias_fn, test_output)
    bcpnn.testing(test_input, weights_fn, bias_fn, output_fn=test_output)



def create_pyr_parameters(params):

    PyrReadoutParamClass = CreatePyrReadoutParameters.CreatePyrReadoutParameters(params)
    PyrReadoutParamClass.write_pyr_parameters()
    PyrReadoutParamClass.write_readout_parameters()
    del PyrReadoutParamClass 


def create_connections(params):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    print "Creating OB -> OC connections ..."
    GetConn = TransformAbstractToDetailedConnectivity.GetConnections(params, comm, rank, debug=1)

    if (rank == 0):
        GetConn.get_mit_rsnp_connections(random_conn=params['ob_oc_random_conns'])
        GetConn.get_mit_pyr_connections(random_conn=params['ob_oc_random_conns'])
        GetConn.get_oc_oc_connections(random_conn=params['oc_oc_random_conns'])
        GetConn.get_pyr_readout_connections()
    comm.Barrier()

    print "Creating hyper and minicolumns ..."
    ConnCreator = CreateHyperAndMinicolumns.CreateHyperAndMinicolumns(params)#, offset=0)
    ConnCreator.create_connections()
    comm.Barrier()

    print "Folder name:", params['folder_name']


if __name__ == '__main__':
    print info_txt
    # ------------ I N I T -----------------------------
    # The simulation_parameters module defines a class for simulation parameter storage
    param_tool = simulation_parameters.parameter_storage()
    # params is the dictionary with all parameters
    params = param_tool.params

    ok = raw_input('\nContinue to create this folder structure? Parameters therein will be overwritten\n\ty / Y / blank = OK; anything else --> exit\n')
    if not ((ok == '') or (ok.capitalize() == 'Y')):
        print 'quit'
        exit(1)

    if len(sys.argv) > 1: 
        # you want to create a new folder structure to not overwrite the results obtained from 
        # presenting patterns to EPTH and OB --> 'preLearning' results
        # it's convenient to write results into a seperate folder structure as OB->OC connectivity
        # and pattern recognition might require some tuning (as every model does ;))
        if sys.argv[1] == 'new':
            print 'New folder:', params['folder_name']

        try: 
            preLearning_folder = os.path.abspath(sys.argv[2])

            # merge the mit spike files in the preLearning_folder
            if not(os.path.exists(os.path.abspath(preLearning_folder) + '/NumberOfSpikes/mit_nspikes_merged_0.dat')):
                cmd = 'python MergeSpikefiles.py %s mit' % (preLearning_folder)
                print cmd
                os.system(cmd)
            cmd = 'cp %s/NumberOfSpikes/mit_nspikes_merged_* %s' % (preLearning_folder, params['nspikes_folder'])
            print cmd
            os.system(cmd)
            cmd = 'cp %s/Parameters/* %s' % (preLearning_folder, params['params_folder'])
            print cmd
            os.system(cmd)
            cmd = 'cp %s/Connections/* %s' % (preLearning_folder, params['conn_folder'])
            print cmd
            os.system(cmd)
            if params['oc_only']:
                cmd = 'cp %s/Spiketimes/mit_spiketimes_merged_* %s' % (preLearning_folder, params['spiketimes_folder'])
                print cmd
                os.system(cmd)
                for pn in xrange(params['n_patterns']):
                    add_first_line(pn)

        except:
            print '\n\tCould not copy MT spike times to new folder. You need to copy them by hand from the preLearning folder.\n'
            print '\tHave you merged the mit spike files in the preLearning folder?\n'
            exit(1)
        
    param_tool.write_parameters_to_file(params["info_file"])
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # 

    if not params['oc_only']:
        prepare_epth_ob_prelearning.prepare_epth_ob(params)

#     ------------ MDS + VQ of OB output ---------------
#    ObAnalyser = AnalyseObOutput.AnalyseObOutput(params)
#    ObAnalyser.get_output_activity() # optional, checks if activity for some patterns is missing
#    ObAnalyser.get_mit_response_normalized()
#    ObAnalyser.rescale_activity()
#    ObAnalyser.rescale_activity_cellwise()
#    ObAnalyser.rescale_activity_patternwise()
#    ObAnalyser.rescale_activity_glom_patterns()
#    mds_vq_ob_output(params)
#    mds_vq_ob_output(params, mds_output_fn='Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/Other/mds_ob_oc_output.dat')

    bcpnn_ob_oc(params)
    bcpnn_oc_oc(params)
    bcpnn_oc_readout(params)

    create_pyr_parameters(params)
    create_connections(params)

    for fn in params['all_connection_fns']:
        assert os.path.exists(fn), 'Required file does not exist: %s' % fn
    for fn in [params['pyr_params_file'], params['readout_params_file']]:
        assert os.path.exists(fn), 'Required file does not exist: %s' % fn

    for pn in xrange(params['n_patterns']):
        fn = params['orn_params_fn_base'] + '%d.dat' % pn
        assert os.path.exists(fn)

    print 'Ready - Go!\t', params['folder_name']
