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
from ExponentialFiltering import filter_spike_train

class Plotter(object):
    def __init__(self, params):

        self.params = params
        self.Merger = MergeSpikefiles.MergeSpikefiles(params)
        self.binwidth = 1 # [ms]
        self.n_timebins = self.params['t_sim'] / self.binwidth
        self.tau = 5 # [ms] time constant for exponential filtering of spike trains

    def plot(self, pattern_nr, cell_type, mc_clustering=False):
        self.cell_type = cell_type
        self.pattern_nr = pattern_nr

        n_cells = self.params['n_%s' % cell_type]
        gid_offset = self.params['%s_offset' % cell_type]
        d = np.zeros((n_cells, self.n_timebins))

        print "Celltype: %s pattern_nr: %d" % (cell_type, pattern_nr)
        nspikes_fn = self.params['%s_spikes_merged_fn_base' % cell_type] + '%d.dat' % (pattern_nr)

        spiketimes_fn = self.params['%s_spiketimes_merged_fn_base' % (cell_type)] + str(pattern_nr) + '.dat'

        if (not os.path.exists(spiketimes_fn)) or (not os.path.exists(nspikes_fn)):
            print 'DEBUG Merging ...\n\n'
            self.Merger.merge_nspike_files(self.params['%s_spike_fn_base' % (cell_type)], self.params['%s_spikes_merged_fn_base' % (cell_type)], pattern_nr)
            self.Merger.merge_spiketimes_files(self.params['%s_spiketimes_fn_base' % (cell_type)], self.params['%s_spiketimes_merged_fn_base' % (cell_type)], pattern_nr)
        print "Loading data ", nspikes_fn
        nspikes = np.loadtxt(nspikes_fn)
        print "Loading data ", spiketimes_fn
        spiketimes = np.loadtxt(spiketimes_fn)

        d_change = np.zeros((n_cells, self.n_timebins))
        for cell in xrange(n_cells):
            gid = cell + gid_offset
            idx = (nspikes[:, 0] == gid).nonzero()[0]
            if nspikes[idx, 1] > 1:
                print 'Cell %d spikes %d times' % (gid, nspikes[idx, 1])
                spike_idx = (spiketimes[:, 1] == gid).nonzero()[0]
                spikes = spiketimes[spike_idx, 0]
                spikes.sort()
                t_vec, filtered_trace = filter_spike_train(spikes, dt=self.binwidth, tau=self.tau, t_max=self.params['t_sim'])
                d[cell, :] = filtered_trace
#                d_change[cell, :] = np.abs(d[cell, 1:] - d[cell, :-1])
#                if d[cell, :].mean() != 0:
#                    d_change[cell, :] -= d.mean()
#                    d_change[cell, :] /= d.mean()


        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Activity for %s cells during pattern %d' % (cell_type.capitalize(), pattern_nr))
        ax.set_xlabel('Time [ms]')
        ax.set_ylabel('%s cell' % cell_type.capitalize())
        print "plotting ...."
        cax = ax.pcolormesh(d)#, cmap='binary')
        ax.set_xlim((0, d.shape[1]))
        ax.set_ylim((0, d.shape[0]))
        cbar = pylab.colorbar(cax)
        clabel = 'Exponentially filtered activity'
        cbar.set_label(clabel)
        fig_fn = self.params['figure_folder'] + '/' + 'filtered_%s_activity_pn%d.png' % (cell_type, pattern_nr)
        print 'Saving to:', fig_fn
        pylab.savefig(fig_fn)

#        fig = pylab.figure()
#        ax = fig.add_subplot(111)
#        ax.set_title('Change in activity for %s cells during pattern %d' % (cell_type.capitalize(), pattern_nr))
#        ax.set_xlabel('Time [ms]')
#        ax.set_ylabel('%s cell' % cell_type.capitalize())
#        print "plotting ...."
#        cax = ax.pcolormesh(d_change)#, cmap='binary')
#        ax.set_xlim((0, d_change.shape[1]))
#        ax.set_ylim((0, d_change.shape[0]))
#        cbar = pylab.colorbar(cax)
#        clabel = 'Change in exponentially filtered activity'
#        cbar.set_label(clabel)
#        fig_fn = self.params['figure_folder'] + '/' + 'filtered_%s_activity_change_pn%d.png' % (cell_type, pattern_nr)
#        print 'Saving to:', fig_fn
#        pylab.savefig(fig_fn)

        pylab.show()



def get_arguments():
    info_txt = 'Usage: python plot_activity_as_colormap.py FOLDER_NAME CELL_TYPE [PATTERN_NR]'
    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        print info_txt
        print 'Taking the parameters currently in simulation_parameters.py'
        param_tool = simulation_parameters.parameter_storage()

    params = param_tool.params
    print 'debug n_cells', params['n_cells']

    try:
        cell_type = sys.argv[2]
    except:
        print 'Missing cell_type argument'
        print info_txt
        exit(1)

    try:
        pn = int(sys.argv[3])
    except:
        print info_txt
        print 'Plotting pattern 0'
        pn = 0

    return params, cell_type, pn

if __name__ == '__main__':
    params, cell_type, pn = get_arguments()
    P = Plotter(params)
    P.plot(pn, cell_type)
