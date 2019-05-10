import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import pylab
import numpy as np
import sys
import os
import MergeSpikefiles 
import simulation_parameters
from FigureCreator import plot_params

class Plotter(object):
    def __init__(self, params):

        self.params = params
        self.Merger = MergeSpikefiles.MergeSpikefiles(params)

    def plot(self, pn_max, cell_type, mc_clustering=False):
        self.cell_type = cell_type
        self.pn_max = pn_max
        d = np.zeros((pn_max, self.params['n_%s' % cell_type]))
        for pn in xrange(pn_max):

            print "Celltype: %s pattern_nr: %d" % (cell_type, pn)
            fn = self.params['%s_spikes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
            if not os.path.exists(fn):
                print 'DEBUG Merging ...\n\n'
                self.Merger.merge_nspike_files(self.params['%s_spike_fn_base' % (cell_type)], self.params['%s_spikes_merged_fn_base' % (cell_type)], pn)
                self.Merger.merge_spiketimes_files(self.params['%s_spiketimes_fn_base' % (cell_type)], self.params['%s_spiketimes_merged_fn_base' % (cell_type)], pn)

            print "Loading data ", fn
            data = np.loadtxt(fn)
            idx = np.array(data[:, 0], dtype=np.int) - self.params['%s_offset' % cell_type]
#            print 'debug', idx.shape, d.shape, data.shape
            d[pn, idx] = data[:, 1]

        if mc_clustering:
            d = self.cluster_nspikes_by_mc(d)
            xlabel = 'Minicolumn'
            clabel = 'Average number of spikes in MC'
            fig_fn = self.params['figure_folder'] + '/' + '%s_activity_%dpatterns_mcclustered.png' % (self.cell_type, self.pn_max)
        else:
            xlabel = 'Cells'
            clabel = 'Number of spikes'
            fig_fn = self.params['figure_folder'] + '/' + '%s_activity_%dpatterns.png' % (self.cell_type, self.pn_max)

        self.get_units_above_thresh(d)

        fig = pylab.figure()
        ax = fig.add_subplot(111)
#        ax.set_title('Activity for %s cells over %d patterns' % (cell_type.capitalize(), pn_max))
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Pattern number')
        print "plotting ...."
        cax = ax.pcolormesh(d)#, cmap='binary')
        ax.set_xlim((0, d.shape[1]))
        ax.set_ylim((0, d.shape[0]))

        cbar = pylab.colorbar(cax)
        cbar.set_label(clabel)
        print 'Saving to:', fig_fn
        pylab.savefig(fig_fn, dpi=300)
        pylab.show()


    def cluster_nspikes_by_mc(self, all_nspikes):

        n_col = self.params['n_hc'] * self.params['n_mc']
        d_out = np.zeros((all_nspikes.shape[0], n_col))
        n_cells_per_mc = self.params['n_%s_per_mc' % self.cell_type]
        for pn in xrange(self.pn_max):
            for mc in xrange(n_col):
                idx0 = mc * n_cells_per_mc
                idx1 = (mc + 1) * n_cells_per_mc
                d_out[pn, mc] = all_nspikes[pn, idx0:idx1].sum() / n_cells_per_mc

        return d_out


    def get_units_above_thresh(self, activity, thresh=.02):
        """
        Returns an array with the units above a certain threshold in activity.
        activity 
        """

        n_patterns = activity.shape[0]
        n_units = activity.shape[1]
#        activity_thresh = thresh * np.mean(activity)
#        activity_thresh = 1.
        activity_thresh = thresh * np.max(activity)
        print 'activity_thresh:', activity_thresh
        n_active_units_per_pattern = np.zeros(n_patterns)
        avg_activity_per_pattern = np.zeros((n_patterns, 2))
        for pn in xrange(n_patterns):
            idx = (activity[pn, :] > activity_thresh).nonzero()[0]
            n_active_units_per_pattern[pn] = idx.size
            avg_activity_per_pattern[pn, 0] = activity[pn ,idx].mean()
            avg_activity_per_pattern[pn, 1] = activity[pn ,idx].std()

            print 'n above thresh in pn %d: \t%d' % (pn, idx.size)
            print 'Average activity in pn %d: %.2f +- %.2f' % (pn, avg_activity_per_pattern[pn, 0], avg_activity_per_pattern[pn, 1])

        print 'Average firing rate for all patterns: %.2f +- %.2f' % (avg_activity_per_pattern[:, 0].mean() / self.params['t_sim'] * 1000., avg_activity_per_pattern[:, 1].mean() / self.params['t_sim'] * 1000.)
        n_active_patterns = np.zeros(n_units)
        for unit in xrange(n_units):
            idx = (activity[:, unit] > activity_thresh).nonzero()[0]
#            print 'unit %d is %d times above thresh' % (unit, idx.size), idx
            n_active_patterns[unit] = idx.size

        n_active_at_least_once = (n_active_patterns > 0).nonzero()[0].size
        n_active_more_than_once = (n_active_patterns > 1).nonzero()[0].size
        n_active_more_than_twice = (n_active_patterns > 2).nonzero()[0].size
        n_active_more_than_threetimes = (n_active_patterns > 3).nonzero()[0].size

        print 'Mean active units per pattern: %.2f +- %.2f ~ (%.2f +- %.2f percent)' % (n_active_units_per_pattern.mean(), n_active_units_per_pattern.std(), \
                n_active_units_per_pattern.mean() / float(n_units) * 100, n_active_units_per_pattern.std() / float(n_units) * 100.)
        print 'Mean number of active patterns per unit: %.2f +- %.2f ~ (%.2f +- %.2f percent)' % (n_active_patterns.mean(), n_active_patterns.std(), \
                n_active_patterns.mean() / float(n_patterns) * 100., n_active_patterns.std() / float(n_patterns) * 100.)
        print 'Number of units active at least once : %d ~ %.2f percent' % (n_active_at_least_once, n_active_at_least_once/ float(n_units) * 100.)
        print 'Number of units more active than once: %d ~ %.2f percent' % (n_active_more_than_once, n_active_more_than_once / float(n_units) * 100.)
        print 'Number of units more active than twice: %d ~ %.2f percent' % (n_active_more_than_twice, n_active_more_than_twice / float(n_units) * 100.)
        print 'Number of units more active than threetimes: %d ~ %.2f percent' % (n_active_more_than_threetimes, n_active_more_than_threetimes / float(n_units) * 100.)
        print 'Threshold is %.2e, absolute value: %.2f spikes' % (thresh, activity_thresh)






if __name__ == '__main__':

    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        param_tool = simulation_parameters.parameter_storage()

    params = param_tool.params
    print 'debug n_cells', params['n_cells']

    try:
        cell_type = sys.argv[2]
    except:
        print 'Missing cell_type argument'
        print 'Usage: python plot_activity_as_colormap.py FOLDER_NAME CELL_TYPE [PN_MAX]'
        exit(1)

    try:
        pn_max = int(sys.argv[3])
    except:
        print 'Plotting all patterns'
        pn_max = params['n_patterns']

#    pn_max = 4
    pylab.rcParams.update(plot_params)
    P = Plotter(params)

    if (cell_type == 'pyr' or cell_type == 'rsnp'):
        ok = raw_input('Cluster nspikes by minicolumn?\n')
        if ok.capitalize() == 'Y' or ok == '':
            P.plot(pn_max, cell_type, mc_clustering=True)
        else:
            P.plot(pn_max, cell_type)
    else:
        P.plot(pn_max, cell_type)
