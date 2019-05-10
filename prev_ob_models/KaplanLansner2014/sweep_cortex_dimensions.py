info_txt = \
"""
    Before running this script, you should have:
    - the simulation results of a 'pre-learning' run 

    - chosen a new folder name and run this script 
    only creating the new folder structure and parameter files by running (see below)
        param_tool.write_parameters_to_file() # 
    - then, copy the nspike files of mitral cells into the new NumberOfSpikes folder:
     cp PRE_LEARNING/NumberOfSpikes/mit*merged* POST_LEARNING/NumberOfSpikes/
"""

import os
import sys
import time
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
import CreateObOcConnections as ObOc
# classes for setting up connectivity and the individual cell parameters
t_init = time.time()


def check_readout_test(fn):
    d = np.loadtxt(fn)
    n_correct = 0
    for row in xrange(d.shape[0]):
        winner = d[row, :].argmax()
        if row == winner:
            n_correct += 1
    return n_correct






if __name__ == '__main__':



    # all setups should use the same MI-coordinates for the MT cells
#    mds_output_fn = 'Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/Other/mds_ob_oc_output.dat'
#    mit_response_normalized_fn = 'Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/NumberOfSpikes/mit_response_normalized_np50'

#    mds_output_fn = 'Cluster_trainingWithNoisyPatterns_nGlom40_nHC12_nMC30_vqOvrlp4_np350_OcOnly/Other/mds_ob_oc_output.dat'
#    mds_output_fn = 'Cluster_TrainingWithNoisyPatterns_nGlom40_nHC60_nMC6_vqOvrlp0_np150_OcOnly/Other/mds_ob_oc_output.dat'
#    mit_response_normalized_fn = 'Cluster_TrainingWithNoisyPatterns_nGlom40_nHC60_nMC6_vqOvrlp0_np150_OcOnly/NumberOfSpikes/mit_response_normalized_np150'
#    mit_response_normalized_fn = 'Cluster_TrainingWithNoisyPatterns_nGlom40_nHC60_nMC6_vqOvrlp0_np150_OcOnly/NumberOfSpikes/mit_response_normalized_np150'

    mds_output_fn = 'Cluster_SweepCtxDim_NoisyPatterns_BinAct_nGlom40_nHC12_nMC30_vqOvrlp4_np150_OcOnly/Other/mds_ob_oc_output.dat'
    mit_response_normalized_fn = 'Cluster_SweepCtxDim_NoisyPatterns_BinAct_nGlom40_nHC12_nMC30_vqOvrlp4_np150_OcOnly/NumberOfSpikes/mit_response_normalized_np150'


    n_patterns = 50
    n_pattern_sets = 3
    n_patterns_total = n_patterns * n_pattern_sets
    n_hc_mc = [(6, 60), (8, 45), (12, 30), (18,20), (20, 18), (30, 12), (45, 8), (60, 6)]

#    n_hc_mc = [(6, 60)]
    n_trials = 5
#    n_hc_mc = [(3, 160), (4, 120), (5, 96), (6, 80), (8, 60), (10, 48), (12, 40), (16, 30), (20, 24), \
#            (24, 20), (30, 16), (40, 12), (48, 10), (60, 8)]

#    n_hc_mc.reverse()

    log_file = open('ctx_params_sweep_noisy_%dpatterns_%dMCs_binactivation.txt' % (n_patterns_total, n_hc_mc[0][0] * n_hc_mc[0][1]), 'w')
    log_txt = '# n_hc, n_mc, vq_overlap, correct_trials.mean(), correct_trials.mean() / n_patterns,  correct_trials.std() / n_patterns\tn_dim_mds\tn_trials\n'
    log_file.write(log_txt)
    for i_ in xrange(len(n_hc_mc)):
        n_hc, n_mc = n_hc_mc[i_][0], n_hc_mc[i_][1]
#        for vq_overlap in [0]:
        max_overlap = max(0, n_hc)
#        max_overlap = min(max_overlap, 20)
        for vq_overlap in xrange(0, max_overlap):
            print '\nVQ overlap %d / %d\n' % (vq_overlap, max_overlap)
            correct_trials = np.zeros(n_trials)
            for trial in xrange(n_trials):

                print 'n_hc, n_mc', n_hc, n_mc
                param_tool = simulation_parameters.parameter_storage()
                params = param_tool.params
                params['n_hc'] = n_hc
                params['n_mc'] = n_mc
                params['n_readout'] = n_patterns
                params['vq_ob_oc_overlap'] = vq_overlap
                param_tool.set_filenames()

                cmd = 'cp %s %s' % (mit_response_normalized_fn, params['nspikes_folder'])
                print 'Copying', cmd
                os.system(cmd)

                ObOc.mds_vq_ob_output(params, mds_output_fn)
#                ObOc.mds_vq_ob_output(params)
                ObOc.bcpnn_ob_oc(params)
                ObOc.bcpnn_oc_oc(params)

                readout_activation = np.zeros((n_patterns_total, n_patterns))
                for i_ in xrange(n_pattern_sets):
                    idx0 = i_ * n_patterns
                    idx1 = (i_ + 1) * n_patterns
                    readout_activation[idx0:idx1, :] = np.eye(n_patterns)
                ObOc.bcpnn_oc_readout(params, readout_activation)
                
                test_output = params['readout_abstract_activity_fn'].rsplit('.dat')[0] + '_test.dat'
    #            wta_output_fn = output_fn.rsplit('.dat')[0] + '_wta.dat'
                correct_trials[trial] = check_readout_test(test_output)
            log_txt = '%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\n' % (n_hc, n_mc, vq_overlap, correct_trials.mean(), \
                    correct_trials.mean() / float(params['n_patterns']) * 100., correct_trials.std() / float(params['n_patterns']) * 100., params['n_dim_mds'], n_trials)
            print 'log_txt:', log_txt
            log_file.write(log_txt)
            log_file.flush()

    log_file.close()


