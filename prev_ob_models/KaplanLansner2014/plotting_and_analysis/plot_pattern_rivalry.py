import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import simulation_parameters
import numpy as np
import pylab
import MergeSpikefiles
from FigureCreator import plot_params
import matplotlib.cm as cm
import json
from plot_pattern_completion import plot_sniff_input

class Plotter(object):

    def __init__(self, params, training_params, pn_0, pn_1, cell_type='pyr', nspike_thresh=0):
        self.params = params
        self.training_params = training_params
        self.pn_min = pn_0
        self.pn_max = pn_1
        self.cell_type = cell_type
        self.nspike_thresh = nspike_thresh

        self.pyr_pattern_map = {} # key: gid, value: list of patterns for which cell responded
        gid_0 = params['%s_offset' % self.cell_type]
        gid_1 = params['n_%s' % self.cell_type] + params['%s_offset' % self.cell_type]
        for gid in xrange(gid_0, gid_1):
#            self.pyr_pattern_map[gid] = []
            self.pyr_pattern_map[str(int(gid))] = []
        self.pyr_pattern_map_loaded = False


    def save_pyr_pattern_map(self, output_fn):
#        output_fn = self.params['other_folder'] + '/' + 'pyr_pattern_map.json'
        print 'Saving pyr pattern map to:', output_fn
        output_file = file(output_fn, 'w')
        d = json.dump(self.pyr_pattern_map, output_file, indent=0)


    def load_pyr_pattern_map(self, input_fn):
#        input_fn = self.params['other_folder'] + '/' + 'pyr_pattern_map.json'
        print 'Loading pyr pattern map from :', input_fn
        input_file = file(input_fn, 'r')
        self.pyr_pattern_map = json.load(input_file)
        self.pyr_pattern_map_loaded = True


    def create_pyr_pattern_map(self, params):
        for pn in xrange(self.pn_min, self.pn_max):
            print 'create_pyr_pattern_map pn:', pn
            fn = params['%s_spikes_merged_fn_base' % (self.cell_type)] + str(pn) + '.dat'
            gid_0 = params['%s_offset' % self.cell_type]
            gid_1 = params['n_%s' % self.cell_type] + params['%s_offset' % self.cell_type]
            d = np.loadtxt(fn)
            for gid in xrange(gid_0, gid_1):
                if d.shape == (2,):
                    nspikes = d[1]
                else:
                    idx = (d[:, 0] == gid).nonzero()[0]
                    nspikes = d[idx, 1]
                if  nspikes > self.nspike_thresh:
                    self.pyr_pattern_map[str(int(gid))].append(pn)
#                    self.pyr_pattern_map[gid].append(pn)
        self.pyr_pattern_map_loaded = True


        
    def load_training_spiketrains(self, pn, params):

        fn = params['%s_spiketimes_merged_fn_base' % (self.cell_type)] + str(pn) + '.dat'
        if (os.path.exists(fn) == False):
            Merger = MergeSpikefiles.MergeSpikefiles(params)
            Merger.merge_spiketimes_files(params['%s_spiketimes_fn_base' % (self.cell_type)], params['%s_spiketimes_merged_fn_base' % (self.cell_type)], pn)

        d = np.loadtxt(fn)
        gid_0 = params['%s_offset' % self.cell_type]
        gid_1 = params['n_%s' % self.cell_type] + params['%s_offset' % self.cell_type]
        for gid in xrange(gid_0, gid_1):
            st_idx = (d[:, 1] == gid).nonzero()[0]
            spiketimes = d[st_idx, 0]


#    def plot_raster(self, params, fn, ax, pn, title='', alpha=1.):
#    def plot_raster(self, params, fn, fig, pn_morph, pn_0, pn_1, title='', alpha=1., output_fig=''):
#    def plot_raster(self, params, fn, ax, pn_morph, pn_0, pn_1, title='', alpha=1.):
    def plot_raster(self, params, fn, pn_morph, pn_0, pn_1, title='', alpha=1., output_fig=''):
        """
        pn_morph -- the 
        pn_0 -- the 
        pn_1 -- the 
        """

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        if (os.path.exists(fn) == False):
            Merger = MergeSpikefiles.MergeSpikefiles(self.params)
            Merger.merge_spiketimes_files(self.params['%s_spiketimes_fn_base' % (self.cell_type)], self.params['%s_spiketimes_merged_fn_base' % (self.cell_type)], pn_morph)

        print 'Loading ', fn
        d = np.loadtxt(fn)
        if (d.size == 0):
            print 'ERROR file %s has 0 size\nIf there was a problem when merging them, delete the empty one and rerun' % (fn)
            return

        assert self.pyr_pattern_map_loaded, 'Please call either create_pyr_pattern_map or load_pyr_pattern_map before plotting'
        if d.size == 2:
            print 'Found only one spike in file %s \nIf there was a problem when merging them, delete the empty one and rerun' % (fn)
            return 
        color_map = {}
        color_0 = 'b' # color if cell has only been active in the first pattern
        color_1 = 'r' # color if cell has only been active in the second
        color_2 = 'k' # color if cell has only been active in both patterns
        color_3 = '#A6A6A6' # color if cell has been active in none of the two

        spiking_cells = np.unique(d[:, 1])
        for gid in spiking_cells:
            idx = (d[:, 1] == gid).nonzero()[0]
            st = d[idx, 0]
            y = gid * np.ones(st.size)

            kid = str(int(gid))
            active_patterns = self.pyr_pattern_map[kid]
            if ((pn_0 in active_patterns) and ((pn_1) not in active_patterns)):
                c = color_0
                alpha = .8
            elif ((pn_1 in active_patterns) and (pn_0 not in active_patterns)):
                c = color_1
                alpha = .8
            elif ((pn_0 in active_patterns) and (pn_1 in active_patterns)):
                c = color_2
                alpha = 1.0
            else:
#                print 'gid %d active during test pattern %d %d times and training patterns' % (gid, st.size, pn_0), self.pyr_pattern_map[kid]
                c = color_3
                alpha = .4
#            ax.plot(st, y, 'o', markersize=4, c=c, alpha=alpha)
            ax.plot(st, y, 'o', markersize=4, markeredgewidth=.0, c=c, alpha=alpha)

        ax.set_title(title)
        ax.set_xlim((0, self.params['t_sim']))

#        ax.set_ylim((0 + params['%s_offset' % cell_type] - 10, params['n_%s' % cell_type] + params['%s_offset' % cell_type] - 10))
        ax.set_ylim((0 + params['%s_offset' % cell_type], params['n_%s' % cell_type] + params['%s_offset' % cell_type]))
        ylabels = ax.get_yticklabels()
        yticks = ax.get_yticks()
        new_ylabels = []
        for i_, y in enumerate(yticks[1:]):
            new_ylabels.append('%d' % (y - self.params['%s_offset' % self.cell_type]))

        ax.set_xticks(range(0, 1800, 300))
            
        if len(new_ylabels) > 0:
            ax.set_yticklabels(new_ylabels)
        ax.set_xlabel('Time [ms]')
        if output_fig != '':
            print 'Saving:', output_fig
            fig.savefig(output_fig, dpi=300)



if __name__ == '__main__':

    info_txt = \
    """
    Usage:
        python plot_pattern_completion_rivalry.py [TRAINING_FOLDER] [PATTERN_NUMBER_MIN] [PATTERN_NUMBER_MAX]
    """

#    assert (len(sys.argv) > 1), 'ERROR: pattern numbers not given\n' + info_txt
#    pn_0 = int(sys.argv[1])

    training_folder = 'Cluster_OcOcLearning_nGlom40_nHC12_nMC30_vqOvrlp4_np50_OcOnly/'
    # for ORN you need to have the actual ORN spikes
#    training_folder = 'Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/'

#    plot_folder = 'Cluster_PatternRivalryTestPostLearningWithSniff2_ORnoise0.00_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/' #Spiketimes/pyr_spiketimes_merged_0.dat
    plot_folder = 'Cluster_PatternRivalryMorphingPLWithSniff_wRsnpPyr3.0e-03_ORnoise0.00_nGlom40_nHC12_nMC30_vqOvrlp4_np350_FullSystem/'

    params_fn = os.path.abspath(plot_folder) + '/Parameters/simulation_parameters.json'
    param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    params = param_tool.params
    training_params_fn = os.path.abspath(training_folder) + '/Parameters/simulation_parameters.json'
    training_param_tool = simulation_parameters.parameter_storage(params_fn=training_params_fn)
    training_params = training_param_tool.params

    cell_type = sys.argv[1]
#    cell_type = 'pyr'
#    cell_type = 'pyr'
#    n_morph_sets = 2
    morph_set_start = 15
    morph_set_stop = 16
    nspike_thresh = 1

    pn_max = training_params['n_patterns']
    P = Plotter(params, training_params, 0, params['n_patterns_test_rivalry'], cell_type=cell_type, nspike_thresh=nspike_thresh)
#    P.create_pyr_pattern_map(training_params)
#    P.save_pyr_pattern_map(training_params['other_folder'] + '/%s_pattern_map.json' % cell_type)
    P.load_pyr_pattern_map(training_params['other_folder'] + '/%s_pattern_map.json' % cell_type)

    plot_params['figure.subplot.left'] = .14
    plot_params['figure.subplot.right'] = .91
    plot_params['figure.subplot.top'] = .82
    plot_params['figure.subplot.bottom'] = .15
#    plot_params['figure.subplot.top'] = .94
    plot_params['xtick.labelsize'] = 24
    plot_params['ytick.labelsize'] = 24

    plot_params['axes.labelsize'] = 32
    plot_params['axes.titlesize'] = 32
#    plot_params['lines.markeredgewidth'] = .0
    pylab.rcParams.update(plot_params)

    pn_per_set = len(params['rivalry_morph_stages'])
    for morph_set in xrange(morph_set_start, morph_set_stop):

        for pn_in_set in xrange(0, pn_per_set):

            pn_morph = morph_set * pn_per_set + pn_in_set
#            fig = pylab.figure()
#            ax = fig.add_subplot(111)
            fn = params['%s_spiketimes_merged_fn_base' % cell_type] + str(pn_morph) + '.dat'
            pn_0 = morph_set
            pn_1 = (morph_set + 1) % params['n_patterns_test_rivalry']
            output_fig = params['figure_folder'] + '/rivalry_raster_%s_pn%d.png' % (cell_type, pn_morph)
            P.plot_raster(params, fn, pn_morph, pn_0, pn_1, title='%s spikes for\nmixture: $%.1f$ B+$%.1f$ R' % \
                    (cell_type.capitalize(), params['rivalry_morph_stages'][pn_in_set], params['rivalry_morph_stages'][-(pn_in_set+1)]), alpha=1., output_fig=output_fig)
#            fig = P.plot_raster(params, fn, fig, pn_morph, pn_0, pn_1, title='%s spikes for\nmixture: $%.1f$ B+$%.1f$ R' % \
#                    (cell_type.capitalize(), params['rivalry_morph_stages'][pn_in_set], params['rivalry_morph_stages'][-(pn_in_set+1)]), alpha=1., output_fig=output_fig)
#            fig = P.plot_raster(params, fn, ax, pn_morph, pn_0, pn_1, title='%s spikes for\nmixture: $%.1f$ B+$%.1f$ R' % \
#                    (cell_type.capitalize(), params['rivalry_morph_stages'][pn_in_set], params['rivalry_morph_stages'][-(pn_in_set+1)]), alpha=1.)
#                    (params['rivalry_morph_stages'][pn_in_set], pn_morph - 1, params['rivalry_morph_stages'][-pn_in_set], pn_morph), alpha=1.)
#            output_fig = params['figure_folder'] + '/rivalry_raster_%s_pn%d.pdf' % (cell_type, pn_morph)
#            print 'Saving:', output_fig
#            fig.savefig(output_fig, dpi=300)
#            del fig
#            del ax

    #        plot_sniff_input(params, ax)
    #        pylab.show()

    #        P.plot_raster(params, fn, ax, pn_morph, title='Pattern Rivalry (mixture of patterns %d + %d)' % (pn_morph - 1, pn_morph), alpha=1.)
