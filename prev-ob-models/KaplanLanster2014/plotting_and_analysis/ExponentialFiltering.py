import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import simulation_parameters
import numpy as np
import pylab
import MergeSpikefiles
import utils


def filter_spike_train(spikes, dt=1., tau=30., t_max=None):
    """
    spikes: list or array of spikes
    """
    if t_max == None:
        t_max = spikes[-1] + tau
    t_vec = np.arange(0, t_max, dt)
    y = np.zeros(t_vec.size)

    spike_idx = []
    for spike in spikes:
        spike_idx.append((t_vec < spike).nonzero()[0][-1])

    for i_, spike in enumerate(spikes):
        y[spike_idx[i_]:] += np.exp(-(t_vec[spike_idx[i_]:] - spike) / tau)

    return t_vec, y


class ExponentialFilter(object):

    def __init__(self, params, cell_type='pyr'):
        self.params = params

        self.cell_type = cell_type
        self.Merger = MergeSpikefiles.MergeSpikefiles(self.params)



    def plot_filtered_spikes(self, pn, pyr_gid):

        fn = self.params['%s_spiketimes_merged_fn_base' % self.cell_type] + '%d.dat' % pn
        if not os.path.exists(fn):
            self.merge_spike_files(pn, self.cell_type)

        all_spikes = np.loadtxt(fn) 
        # all_spikes column 0 -- spike time
        # all_spikes column 1 -- GID of cell
        if all_spikes.size == 0:
            print 'No spikes found in %s' % fn
            return
        cell_spike_idx = (all_spikes[:, 1] == pyr_gid).nonzero()[0]
        cell_spikes = all_spikes[cell_spike_idx, 0]
        t_vec, y0 = filter_spike_train(cell_spikes)

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.plot(t_vec, y0, label='Filtered spikes', ls='-', c='g', lw=3)
        return fig




    def plot_voltage(self, pn, pyr_gid, fig=None):

        fn = self.params['%s_volt_fn_base' % self.cell_type] + '%d_%d.v' % (pn, pyr_gid)
        if not os.path.exists(fn):
            print 'Voltage file not found:', fn
            return
        
        print 'Plotting', fn
        d = np.loadtxt(fn)
        dt_rec = d.size / params['t_sim']
        t = np.arange(0, params['t_sim'] + 1. / dt_rec, 1. / dt_rec)

        if fig == None:
            fig = pylab.figure()
            ax = fig.add_subplot(111)

        axes = fig.get_axes()
        ax2 = axes[0].twinx()
        ax2.plot(t, d)
        ax2.set_ylabel('Voltage [mV]')
        pylab.legend()


#        p



    def merge_spike_files(self, pn, cell_type):
        print 'Merging nspike file for %s pattern %d' % (cell_type, pn)
        self.Merger.merge_nspike_files(params['%s_spike_fn_base' % cell_type], params['%s_spikes_merged_fn_base' % cell_type], pn)
        print 'Merging spiketimes file for %s pattern %d' % (cell_type, pn)
        self.Merger.merge_spiketimes_files(params['%s_spiketimes_fn_base' % cell_type], params['%s_spiketimes_merged_fn_base' % cell_type], pn)



if __name__ == '__main__':
    info_txt = \
    """
    Usage:
        python ExponentialFiltering.py [FOLDER] [CELLTYPE] 
        or
        python ExponentialFiltering.py [FOLDER] [CELLTYPE] [PATTERN_NUMBER]
    """

    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params

    try:
        cell_type = sys.argv[2]
        assert (cell_type in params['cell_types']), 'Wrong arguments\n%s' % info_txt
    except:
        cell_type = 'pyr'

    print 'cell_type:', cell_type
    try:
        pn = int(sys.argv[3])
    except:
        print 'Using the default pattern number 0'
        pn = 0

    F = ExponentialFilter(params, cell_type=cell_type)
    spiking_cells, nspikes_sorted = utils.sort_gids_by_nspikes(params['%s_spikes_merged_fn_base' % cell_type] + '%d.dat' % pn, gid_offset=params['%s_offset' % cell_type])
    filter_fn = params['%s_volt_fn_base' % cell_type].rsplit('/')[-1] + '%d_' % (pn)
    recorded_gids = set(utils.get_recorded_gids(params['volt_folder'], filter_fn))
    candidate_gids = list(recorded_gids.intersection(spiking_cells))

    print 'spiking cells', len(candidate_gids), candidate_gids
    nspikes = np.loadtxt(params['%s_spikes_merged_fn_base' % cell_type] + '%d.dat' % pn)
    nspikes_candidates = np.zeros(len(candidate_gids))
    for i_, gid in enumerate(candidate_gids):
        idx = (nspikes[:, 0] == gid).nonzero()[0]
        nspikes_candidates[i_] = nspikes[idx, 1]
    idx_sorted = nspikes_candidates.argsort() 
    print nspikes_candidates[idx_sorted]
    print np.array(candidate_gids)[idx_sorted]
    gid = candidate_gids[0]
#    pyr_gid = 80992 #np.params['


    fig = F.plot_filtered_spikes(pn, gid)
    F.plot_voltage(pn, gid, fig=fig)

    pylab.show()
