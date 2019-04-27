import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import matplotlib
#matplotlib.use('')
#matplotlib.use('Qt4Agg')
import simulation_parameters
import MergeSpikefiles
import SetOfCurvesPlotter
import numpy as np
import pylab
from FigureCreator import plot_params
plot_params['xtick.labelsize'] = 22
plot_params['ytick.labelsize'] = 22

plot_params['axes.labelsize'] = 32
plot_params['axes.titlesize'] = 32

def wta(d):
    """
    d : 1-D numpy.array e.g. of activity
    returns the array after the winner-take-all operation
    """
    winner = np.argmax(d)
    d_new = np.zeros(d.size)
    d_new[winner] = 1.0
    return d_new


class Plotter(object):
    def __init__(self, params, pn_max):

        self.params = params
        self.pn_max = pn_max
        self.Merger = MergeSpikefiles.MergeSpikefiles(params)
        self.cell_type = 'readout'


    def plot_activity_vs_morphstage(self, normed=True):
        """
        plots the activity of the two readout cells that should be active dring the respective pattern
        if none of the two is active the activity during that pattern is summed to a third cell indicating 'wrong' readout
        """
        if normed:
            activity = self.activity_norm
        else:
            activity = self.activity
        x_axis = self.params['rivalry_morph_stages']# = [.3, .4, .5, .6, .7]
        n_x = len(x_axis) + 2
        x_axis += [.0 , 1.]
        x_axis.sort()
        n_patterns = self.params['n_patterns_test_rivalry']# = 10     # number of different odor patterns to be checked with different concentrations
        y_data = np.zeros((n_patterns, n_x, 3))
        avg_curve = np.zeros((n_x, 2, 3))
        n_readout = activity.shape[1]
        print 'y_data.shape', y_data.shape
        print 'activity.shape', activity.shape

        pattern_idx = 0
        print 'n_patterns', n_patterns
        for pn in xrange(n_patterns):
            readout_0 = pn 
            readout_1 = (pn + 1) % n_readout
            for i_ in xrange(1, n_x - 1):
                y_data[pn, i_, 0] = activity[pattern_idx, readout_0]
                y_data[pn, i_, 1] = activity[pattern_idx, readout_1]
                y_data[pn, i_, 2] = (activity[pattern_idx, :].sum() - activity[pattern_idx, readout_0] - activity[pattern_idx, readout_1]) / (n_readout - 2)
                pattern_idx += 1

        # blue
        y_data[:, 0, 0] = 0.
        y_data[:, -1, 0] = 1.

        # red
        y_data[:, 0, 1] = 1.
        y_data[:, -1, 1] = .0

        y_data[:, 0, 2] = .0
        y_data[:, -1, 2] = .0

        pylab.rcParams.update(plot_params)
        fig = pylab.figure()
        ax = fig.add_subplot(111)
#        for pn in xrange(n_patterns):
#            ax.plot(x_axis, y_data[pn, :, 0], '-', color='b')
#            ax.plot(x_axis, y_data[pn, :, 1], '-', color='r')
#            ax.plot(x_axis, y_data[pn, :, 2], '-', color='k')

#        fig2 = pylab.figure()
#        ax2 = fig2.add_subplot(111)

        subset = []
        for pn in subset:
            ax.plot(x_axis, y_data[pn, :, 0], '-', color='b')
            ax.plot(x_axis, y_data[pn, :, 1], '-', color='r')
            ax.plot(x_axis, y_data[pn, :, 2], '-', color='k')


        for i_ in xrange(1, n_x - 1):
            avg_curve[i_, 0, 0] = y_data[:, i_, 0].mean()
            avg_curve[i_, 1, 0] = y_data[:, i_, 0].std() / np.sqrt(n_patterns)

            avg_curve[i_, 0, 1] = y_data[:, i_, 1].mean()
            avg_curve[i_, 1, 1] = y_data[:, i_, 1].std() / np.sqrt(n_patterns)

            avg_curve[i_, 0, 2] = y_data[:, i_, 2].mean()
            avg_curve[i_, 1, 2] = y_data[:, i_, 2].std() / np.sqrt(n_patterns)

        # BLUE 
        # left blue
        avg_curve[0, 0, 0] = 0.
        # right blue
        avg_curve[-1, 0, 0] = 59.1 
        avg_curve[-1, 1, 0] = 22.74 / np.sqrt(50)

        # RED
        # left red
        avg_curve[0, 0, 1] = 59.1 
        avg_curve[0, 1, 1] = 22.74 / np.sqrt(50)
        # right blue
        avg_curve[-1, 0, 1] = 0.

        # most left point 
#        avg_curve[0, 0, 0] = .0
#        avg_curve[-1, 1, 0] = .0

        ax.errorbar(x_axis, avg_curve[:, 0, 0], yerr=avg_curve[:, 1, 0], lw=5, color='b')
        ax.errorbar(x_axis, avg_curve[:, 0, 1], yerr=avg_curve[:, 1, 1], lw=5, color='r')
        ax.errorbar(x_axis, avg_curve[:, 0, 2], yerr=avg_curve[:, 1, 2], lw=5, color='k')

        ax.set_xlabel('Fraction of pattern B')
        ax.set_ylabel('Mean number of output spikes')

        output_fn = self.params['figure_folder'] + '/rivalry_morphing.pdf'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=300)

        output_fn = self.params['figure_folder'] + '/rivalry_morphing.png'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=300)



    def load_spikes(self):
        gid_offset = self.params['%s_offset' % self.cell_type]

        # NSPIKES
        n_patterns = self.pn_max
        list_of_winners = -1 * np.ones(n_patterns) # list of most active readout cells
        if self.params['concentration_invariance']:
            n_readout = self.params['n_patterns_test_conc_inv']
        else:
            n_readout = self.params['n_readout']
            n_patterns = self.pn_max

        print 'N_readout', n_readout
        self.activity = np.zeros((n_patterns, n_readout))
        self.activity_norm = np.zeros((n_patterns, n_readout))
        self.activity_wta = np.zeros((n_patterns, n_readout))
        missing_activity = []
        for pn in xrange(n_patterns):
            fn = self.params['%s_spikes_merged_fn_base' % (self.cell_type)] + str(pn) + '.dat'
            if (os.path.exists(fn) == False):
                print "Merging data for file: ", fn
                self.Merger.merge_nspike_files(self.params['%s_spike_fn_base' % (self.cell_type)], self.params['%s_spikes_merged_fn_base' % (self.cell_type)], pn)
            try: 
                print "loading data ", fn
                d = np.loadtxt(fn)
                if (d.ndim == 1):
                    n_cells = 1
                    index = (d[0] - gid_offset)
                    self.activity[pn, index] += d[1]

                else:
                    n_cells = d[:, 0].size
                    for i in xrange(n_cells):
                        index = (d[i, 0] - gid_offset)
                        self.activity[pn, index] = (d[i, 1] / self.params['t_sim']) * 1000.
                most_active_neuron = self.activity[pn, :].argmax() 
                list_of_winners[pn] = most_active_neuron
#                print 'most active neuron in pattern %d is %d' % (pn, most_active_neuron)

            except:
                print "Missing activity for pn %d from file %s" % (pn, fn)
                missing_activity.append(pn)

        print "Missing activity of the following %d patterns:" % len(missing_activity), missing_activity
        for pn in xrange(n_patterns):
            if (self.activity[pn, :].sum() > 0):
                self.activity_norm[pn, :] = self.activity[pn, :] / (self.activity[pn, :].sum())
            self.activity_wta[pn, :] = wta(self.activity[pn, :])

        output_fn = self.params['readout_activity_data']
        np.savetxt(output_fn, self.activity)
        output_fn = self.params['readout_activity_data_normalized']
        np.savetxt(output_fn, self.activity_norm)
        output_fn = self.params['readout_activity_data_wta']
        np.savetxt(output_fn, self.activity_wta)


    def plot_colormaps(self):

        title = "Average firing rate readout cells whole run"
        self.activity_rev = np.zeros(self.activity.shape)
        n_row = self.activity.shape[0]
        for row in xrange(self.activity.shape[0]):
            self.activity_rev[row, :] = self.activity[n_row - row - 1, :]


        pylab.rcParams.update(plot_params)

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(self.activity)
        ax.set_ylim((0, self.activity[:,0].size))
#        ax.set_xlim((0, 10))
        ax.set_xlim((0, self.activity[0, :].size))
        #ax.set_title(title)
        ax.set_xlabel("Readout cell")
        ax.set_ylabel("Pattern")
        cb = pylab.colorbar(cax)
        cb.ax.set_ylabel('Spike rate [Hz]')
        output_fn = self.params['%s_activity_cmap' % self.cell_type]
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=(300))

        title = "Normalized readout activity "
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(self.activity_norm)
        ax.set_ylim((0, self.activity_norm[:,0].size))
#        ax.set_xlim((0, 10))
        ax.set_xlim((0, self.activity_norm[0, :].size))
        ax.set_title(title)
        ax.set_ylabel('Pattern number')
        ax.set_xlabel('Readout cell')
        pylab.colorbar(cax)
        output_fn = self.params['%s_activity_cmap' % self.cell_type].rsplit('.png')[0] + '_normalized.png'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=300)

        title = "Winner-take-all readout activity"
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(self.activity_wta, cmap='binary')
        ax.set_ylim((0, self.activity_wta[:,0].size))
#        ax.set_xlim((0, 10))
        ax.set_xlim((0, self.activity_wta[0, :].size))
        ax.set_title(title)
        ax.set_ylabel('Pattern number')
        ax.set_xlabel('Readout cell')
        pylab.colorbar(cax)
        output_fn = self.params['%s_activity_cmap' % self.cell_type].rsplit('.png')[0] + '_wta.png'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=300)

if __name__ == '__main__':

#    folder = sys.argv[1]
    folder = 'Cluster_PatternRivalryMorphingPLWithSniff_wRsnpPyr3.0e-03_ORnoise0.00_nGlom40_nHC12_nMC30_vqOvrlp4_np350_FullSystem'
    params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
    param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    params = param_tool.params

    try:
        pn_max = int(sys.argv[2])
    except:
        print 'Plotting all patterns'
        pn_max = params['n_patterns']


    P = Plotter(params, pn_max)
    P.load_spikes()
#    P.plot_colormaps()
    P.plot_activity_vs_morphstage(normed=False)

    pylab.show()

