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


def plot_raster(params, fn, ax, pn, title='', color='k', alpha=1.):
    print 'Loading Spikes from:', params['%s_spikes_merged_fn_base' % cell_type]
    if (os.path.exists(fn) == False):
        Merger = MergeSpikefiles.MergeSpikefiles(params)
        Merger.merge_spiketimes_files(params['%s_spiketimes_fn_base' % (cell_type)], params['%s_spiketimes_merged_fn_base' % (cell_type)], pn)

    print 'Loading ', fn
    data = np.loadtxt(fn)
    assert (data.size > 0), 'ERROR file %s has 0 size\nIf there was a problem when merging them, delete the empty one and rerun' % (fn)

    ax.plot(data[:,0], data[:,1], 'o', markersize=5, markeredgewidth=.0, color=color, alpha=alpha)
    ax.set_xlim((0, params['t_sim']))
    ax.set_title(title)
    ax.set_xlabel('Time [ms]')
#    ax.set_ylabel('Cell GID')

    ylabels = ax.get_yticklabels()
    yticks = ax.get_yticks()
    new_ylabels = []
    for i_, y in enumerate(yticks[0:]):
#    for i_, y in enumerate(yticks[1:]):
        new_ylabels.append('%d' % (y - params['%s_offset' % cell_type]))
        
    ax.set_ylim((-1 + params['%s_offset' % cell_type], params['n_%s' % cell_type] + params['%s_offset' % cell_type] + 1))
    if len(new_ylabels) > 0:
        ax.set_yticklabels(new_ylabels)

    xlabels = ax.get_xticklabels()
    xticks = ax.get_xticks()
    new_xlabels = ['']
    for i_, x in enumerate(xticks[1:-1]):
#    for i_, x in enumerate(xticks[1:]):
        new_xlabels.append('%d' % x)
    new_xlabels.append('')
    ax.set_xticklabels(new_xlabels)


def get_sniff_amplitude(x, tstart, tstop, T, t_shift, amp):
    f_x = 0
    if (x > tstart) and (x < tstop):
        f_x = (amp * (np.sin(x / (T) - t_shift))**2)
    return f_x


def plot_sniff_input(params, ax):
    if params['with_sniffing_input']:
        tstop = params['t_stop'] = 1200 # [ms]
        tstart = params['t_start'] = 200 # [ms]
        T = params['sniff_period'] = 80. # [ms]
        t_shift = params['t_shift_sniff'] = 40. # [ms]
    else:
        print 'This was run without sniffing input\nReturn None'
        return None

    times = np.arange(0, params['t_sim'], 5)

    ylim = ax.get_ylim()
    alpha_max = .2
    c = 'b'
    for t in times:
        f_x = get_sniff_amplitude(t, tstart, tstop, T, t_shift, 1.0)
#        print 'f_x', f_x
        ax.plot((t, t), (ylim[0], ylim[1]), lw=4, ls='-', c=c, alpha=f_x * alpha_max)




if __name__ == '__main__':

    info_txt = \
    """
    Usage:
        python plot_pattern_completion_rivalry.py [PATTERN_NUMBER]
    """
#        python plot_pattern_completion_rivalry.py [TRAINING_FOLDER] [TEST_FOLDER] [PATTERN_NUMBER_MIN] [PATTERN_NUMBER_MAX]
    assert (len(sys.argv) > 1), 'ERROR: pattern number not given\n' + info_txt
    pn_max = int(sys.argv[1])
    
    training_folder = 'Cluster_OcOcLearning_nGlom40_nHC12_nMC30_vqOvrlp4_np50_OcOnly/'
    plot_folder = 'Cluster_PatternCompletionTestPostLearningWithSniff_fOR0.50_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/'
    params_fn = os.path.abspath(plot_folder) + '/Parameters/simulation_parameters.json'
    param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    params = param_tool.params
    training_params_fn = os.path.abspath(training_folder) + '/Parameters/simulation_parameters.json'
    training_param_tool = simulation_parameters.parameter_storage(params_fn=training_params_fn)
    training_params = training_param_tool.params
    cell_type = 'readout'
#    cell_type = 'pyr'
#    cell_type = 'mit'

    for pn in xrange(pn_max):
        training_fn = training_params['%s_spiketimes_merged_fn_base' % cell_type] + str(pn) + '.dat'
        test_fn = params['%s_spiketimes_merged_fn_base' % cell_type] + str(pn) + '.dat'

        plot_params['figure.subplot.left'] = .11
        plot_params['figure.subplot.top'] = .92
        plot_params['figure.subplot.right'] = .98
        plot_params['xtick.labelsize'] = 24
        plot_params['ytick.labelsize'] = 24

        plot_params['axes.labelsize'] = 32
        plot_params['axes.titlesize'] = 32
        pylab.rcParams.update(plot_params)

        fig = pylab.figure()
        ax = fig.add_subplot(111)


        color_0 = '#A6A6A6'
        color_1 = 'b'
#        title = 'Pattern completion test pattern %d' % (pn)
#        title = 'MT spikes'
        title = '%s spikes ' % (cell_type.capitalize())
        plot_raster(training_params, training_fn, ax, pn, title=title, color=color_0, alpha=0.9)
        plot_raster(params, test_fn, ax, pn, title=title, color=color_1, alpha=1.)
#        plot_sniff_input(params, ax)

        output_fn = params['figure_folder'] + '/' + 'competion_raster_%s_%d.png' % (cell_type, pn)
        print 'Saving figure to', output_fn
        pylab.savefig(output_fn, dpi=(300))

    pylab.show()

