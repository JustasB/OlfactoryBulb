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
    def __init__(self, params, pn_max, abstract_readout):

        self.params = params
        self.Merger = MergeSpikefiles.MergeSpikefiles(params)

        cell_type = 'readout'
        gid_offset = params['%s_offset' % cell_type]


        # NSPIKES
        n_patterns = pn_max
        list_of_winners = -1 * np.ones(n_patterns) # list of most active readout cells
        if params['concentration_invariance']:
            n_readout = params['n_patterns_test_conc_inv']
        else:
            n_readout = params['n_readout']
            n_patterns = pn_max
        print 'N_readout', n_readout
        activity = np.zeros((n_readout, n_patterns))
        activity_norm = np.zeros((n_readout, n_patterns))
        activity_wta = np.zeros((n_readout, n_patterns))
        mean_activity = np.zeros(n_patterns)
        missing_activity = []
        for pn in xrange(n_patterns):
        #    print "Celltype: %s pattern_nr: %d" % (cell_type, pn)
            fn = params['%s_spikes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
            if (os.path.exists(fn) == False):
                print "Merging data for file: ", fn
                self.Merger.merge_nspike_files(params['%s_spike_fn_base' % (cell_type)], params['%s_spikes_merged_fn_base' % (cell_type)], pn)
            try: 
                print "loading data ", fn
                d = np.loadtxt(fn)
                if (d.ndim == 1):
                    n_cells = 1
                    index = (d[0] - gid_offset)
                    activity[index, pn] += d[1]

                else:
                    n_cells = d[:, 0].size
                    for i in xrange(n_cells):
                        index = (d[i, 0] - gid_offset)
                        activity[index, pn] = (d[i, 1] / params['t_sim']) * 1000.
                most_active_neuron = activity[:, pn].argmax() 
                list_of_winners[pn] = most_active_neuron
                mean_activity[pn] = activity[most_active_neuron, pn]
#                print 'most active neuron in pattern %d is %d' % (pn, most_active_neuron)

            except:
                print "Missing activity for pn %d from file %s" % (pn, fn)
                missing_activity.append(pn)

        print 'Mean activity of winner neurons:', mean_activity.mean(), mean_activity.std()
        print "Missing activity of the following %d patterns:" % len(missing_activity), missing_activity
        cnt_correct = 0
        incorrect = []
        incorrect_without_silent = []
        for pn in xrange(n_patterns):
            if (list_of_winners[pn] == np.argmax(abstract_readout[pn, :])):
                cnt_correct += 1
            else:
                incorrect.append(pn)
                if (pn not in missing_activity):
                    incorrect_without_silent.append(pn)
            if (activity[:, pn].sum() > 0):
                activity_norm[:, pn] = activity[:, pn] / (activity[:, pn].sum())
            activity_wta[:, pn] = wta(activity[:, pn])

        np.savetxt(params['incorrect_patterns_without_silent'], np.array(incorrect_without_silent), fmt='%d')
        np.savetxt(params['silent_patterns'], np.array(missing_activity), fmt='%d')
        #exit(1)

        print 'Missing patterns', missing_activity
        print "Number correctly recognized patterns: %d / %d = %.2f percent (%.2f percent when silent are not counted)" % (cnt_correct, n_patterns, (cnt_correct / float(n_patterns)) * 100., (cnt_correct / float(n_patterns - len(missing_activity))) * 100.)
        print "Incorrectly recognized patterns: ", len(incorrect), incorrect
        print "Incorrectly recognized patterns without silent: ", len(incorrect_without_silent), incorrect_without_silent

        output_fn = params['readout_activity_data']
        np.savetxt(output_fn, activity)
        output_fn = params['readout_activity_data_normalized']
        np.savetxt(output_fn, activity_norm)
        output_fn = params['readout_activity_data_wta']
        np.savetxt(output_fn, activity_wta)

        #print "List of winner readouts:", list_of_winners
        #print "List of winner readouts:", np.unique(list_of_winners)
        #winner_readouts = np.unique(list_of_winners)
        #print "Number of different winner readouts:", np.unique(list_of_winners).size
        #winners_corrected = list(np.unique(list_of_winners))
        #for i in xrange(len(missing_activity)):
        #    if (missing_activity[i] in winner_readouts):
        #        winners_corrected.remove(missing_activity[i])
        #print "Corrected list (without silent patterns):", winners_corrected
        #print "Num correct", len(winners_corrected)


        title = "Average firing rate readout cells whole run"
        activity = activity.transpose()
        activity_rev = np.zeros(activity.shape)
        n_row = activity.shape[0]
        for row in xrange(activity.shape[0]):
            activity_rev[row, :] = activity[n_row - row - 1, :]


        pylab.rcParams.update(plot_params)

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(activity)
        ax.set_ylim((0, activity[:,0].size))
#        ax.set_xlim((0, 10))
        ax.set_xlim((0, activity[0, :].size))
        #ax.set_title(title)
        ax.set_xlabel("Readout cell")
        ax.set_ylabel("Pattern")
        cb = pylab.colorbar(cax)
        cb.ax.set_ylabel('Spike rate [Hz]')
        output_fn = params['%s_activity_cmap' % cell_type]
        output_fn = 'Cluster_ORnoise0.0_nGlom40_nHC12_nMC30_vqOvrlp4_np50_fullSystem_ConcInv/Figures/readout_activity.pdf'
#        output_fn = 'Cluster_ConcInvTraining_nGlom40_nHC12_nMC30_vqOvrlp0_np50_OcOnly/Figures/readout_activity.pdf'
#        output_fn = 'Cluster_OcOcLearning_nGlom40_nHC12_nMC30_vqOvrlp4_np50_OcOnly/Figures/readout_activity.pdf'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=(300))

        title = "Normalized readout activity "
        activity_norm = activity_norm.transpose()
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(activity_norm)
        ax.set_ylim((0, activity_norm[:,0].size))
#        ax.set_xlim((0, 10))
        ax.set_xlim((0, activity_norm[0, :].size))
        ax.set_title(title)
        ax.set_ylabel('Pattern number')
        ax.set_xlabel('Readout cell')
        pylab.colorbar(cax)
        output_fn = params['%s_activity_cmap' % cell_type].rsplit('.png')[0] + '_normalized.pdf'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=300)

        title = "Winner-take-all readout activity"
        activity_wta = activity_wta.transpose()
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(activity_wta, cmap='binary')
        ax.set_ylim((0, activity_wta[:,0].size))
#        ax.set_xlim((0, 10))
        ax.set_xlim((0, activity_wta[0, :].size))
        ax.set_title(title)
        ax.set_ylabel('Pattern number')
        ax.set_xlabel('Readout cell')
        pylab.colorbar(cax)
        output_fn = params['%s_activity_cmap' % cell_type].rsplit('.png')[0] + '_wta.pdf'
        print "saving ....", output_fn
        pylab.savefig(output_fn, dpi=300)
        pylab.show()


    def plot_spiketimes(self):
        pass
        # SPIKETIMES
        #list_of_winners = np.ones(n_patterns) # list of most active readout cells
        #list_of_winners *= -1
        #t1 = 100
        #t2 = 200
        #activity_interval = np.zeros((params['n_readout'], n_patterns))
        #for pn in xrange(n_patterns):
        #    fn = params['%s_spiketimes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
        #    if (os.path.exists(fn) == False):
        #        Merger.merge_spiketimes_files(params['%s_spike_fn_base' % (cell_type)], params['%s_spikes_merged_fn_base' % (cell_type)], pn)
        #    print "loading data ", fn
        #    try:
        #        d = np.loadtxt(fn)
        #        n_spikes = d[:, 0].size
        #        for i in xrange(n_spikes):
        #            t = d[i, 0]
        #            if ((t > t1) and (t < t2)):
        #                index = d[i, 1] - gid_offset
        #                activity_interval[index, pn] += 1
        #        most_active_neuron = activity_interval[:, pn].argmax()
        #        list_of_winners[pn] = most_active_neuron
        #    except:
        #        pass
        #winner_readouts = np.unique(list_of_winners)
        #print "Within interval: List of winner readouts:", np.unique(list_of_winners)
        #print "Within interval: Number of different winner readouts:", np.unique(list_of_winners).size
        #activity_interval /= (t2 - t1) / 1000.
        #title = "Average number of spikes in a minicolumn within interval %d - %d" % (t1, t2)
        #activity_interval = activity_interval.transpose()
        #n_row = activity_interval.shape[0]
        #fig = pylab.figure()
        #ax = fig.add_subplot(111)
        #print "plotting ...."
        #cax = ax.pcolor(activity_interval)
        #ax.set_ylim((0, activity_interval[:,0].size))
        #ax.set_xlim((0, activity_interval[0, :].size))
        #ax.set_title(title)
        #pylab.colorbar(cax)
        #output_fn = params['%s_activity_interval_cmap' % cell_type]
        #print "saving ....", output_fn
        #pylab.savefig(output_fn)
        #print "saving ....", params['readout_activity_interval_data']
        #np.savetxt(params['readout_activity_interval_data'], activity_interval)
        #self.params['readout_activity_cmap'] =  '%s/readout_activity.png' % ( self.params['figure_folder'])
        #self.params['readout_activity_interval_cmap'] =  '%s/readout_activity_interval.png' % ( self.params['figure_folder'])
        #self.params['readout_activity_data'] =  '%s/readout_activity.dat' % ( self.params['figure_folder'])
        #self.params['readout_activity_interval_data'] =  '%s/readout_activity_interval.dat' % ( self.params['figure_folder'])
        #figure_params = {
        #    'figure.subplot.bottom': 0.10,
        #    'figure.subplot.hspace': 0.9,
        #    'figure.subplot.left': 0.125,
        #    'figure.subplot.right': 0.90,
        #    'figure.subplot.top': 0.90,
        #    'figure.subplot.wspace': 0.50,
        #    'legend.fontsize' : 6
        #    }
        #pylab.rcParams.update(figure_params)



if __name__ == '__main__':

    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        param_tool = simulation_parameters.parameter_storage()

    params = param_tool.params


    try:
        pn_max = int(sys.argv[2])
    except:
        print 'Plotting all patterns'
        pn_max = params['n_patterns']

#    params['concentration_invariance'] = 0
    print 'concentration_invariance:', params['concentration_invariance']
    try:
        fn = sys.argv[2]
        abstract_readout = np.loadtxt(fn)
    except:
        print 'Plotting all patterns'
        readout_activation = np.eye(pn_max)
        if params['concentration_invariance']:
            readout_activation = np.zeros((params['n_readout'], params['n_patterns']))
            row = 0 
            for readout_idx in xrange(params['n_patterns_test_conc_inv']):
                for conc_ in xrange(params['n_conc_check']):
                    readout_activation[row, readout_idx] = 1.
                    row += 1

#    pn_max = 4
    P = Plotter(params, pn_max, readout_activation)


