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
from mpl_toolkits.axes_grid1 import ImageGrid

def plot_raster_for_celltype(params, ax, cell_type, pn, title='', show=True, color='k', alpha=1.):
    print 'Loading Spikes from:', params['%s_spikes_merged_fn_base' % cell_type]

    fn = params['%s_spiketimes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
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
    ax.set_ylabel('Cell GID')

    ax.set_ylim((0 + params['%s_offset' % cell_type], params['n_%s' % cell_type] + params['%s_offset' % cell_type]))
    ylabels = ax.get_yticklabels()
    yticks = ax.get_yticks()
    new_ylabels = []
    for i_, y in enumerate(yticks[1:]):
        new_ylabels.append('%d' % (y - 72000))

#    ax.set_xticklabels

#    output_fn = params['figure_folder'] + '/' + 'rasterplot_%s_%d.png' % (cell_type, pn)
#    print 'Saving figure to', output_fn
#    pylab.savefig(output_fn, dpi=(400))

        
#    if len(new_ylabels) > 0:
#        ax.set_yticklabels(new_ylabels)
#    if cell_type == 'pyr':
#        for mc in xrange(params['n_mc'] * params['n_hc']):
#            idx0 = params['pyr_offset'] + mc * params['n_pyr_per_mc']
#            idx1 = params['pyr_offset'] + (mc + 1) * params['n_pyr_per_mc']
#            ax.plot((0, params['t_sim']), (idx0, idx0), ls='-', c='k', alpha=.2)
#            ax.plot((0, params['t_sim']), (idx1, idx1), ls='-', c='k', alpha=.2)


def plot_raster_grid(params, pn_min, pn_max, cell_type):

    fig = pylab.figure(figsize=(16, 8))
    n_fig_x = pn_max - pn_min
    grid = ImageGrid(fig, 111, nrows_ncols=(1, n_fig_x), axes_pad=0.1)#, aspect=(1, 2))

    w, h = 1. / (n_fig_x), 1.
    for i_, pn in enumerate(xrange(pn_min, pn_max)):
        plot_raster_for_celltype(params, grid[pn], cell_type, pn)
        print 'get axes', grid[pn].set_aspect(.1, 'box-forced')
        pos = [float(i_) / n_fig_x, .1, 1, h]
        grid[pn].set_position(pos)
        print 'pos', grid[pn].properties()['position']




if __name__ == '__main__':
    info_txt = \
    """
    Usage:
        python plot_response_curve.py [FOLDER] [CELLTYPE] [PATTERN_NUMBER_MAX]
        or
        python plot_response_curve.py [FOLDER] [CELLTYPE] [PATTERN_NUMBER_MIN] [PATTERN_NUMBER_MAX]
    """
    assert (len(sys.argv) > 3), 'ERROR: folder and cell_type not given\n' + info_txt
    folder = sys.argv[1]
    cell_type = sys.argv[2]
    pn_0 = int(sys.argv[3])
    try:
        pn_max = int(sys.argv[4])
    except:
        pn_max = pn_0 + 2
    print 'plotting from pn %d - %d' % (pn_0, pn_max)

    params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
    param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    params = param_tool.params


    plot_params['figure.subplot.left'] = .15
    plot_params['figure.subplot.top'] = .85
    pylab.rcParams.update(plot_params)

#    plot_raster_grid(params, pn_0, pn_max, cell_type)

#    figsize = (12, 9)
#    fig = pylab.figure(figsize=figsize)

    fig = pylab.figure()
    ax = fig.add_subplot(111)

#    pylab.subplots_adjust(left=.2)


#    colorlist = ["#0000FF", "#006600", "#FF0000", "#00FFFF", "#CC00FF", "#FFFF00", "#000000", \
#            "#00FF00", "#663300", "#FF3399", "#66CCCC", "#FFCC99", "#FFFFCC"]
    colorlist = ['k', '#0EC5FF', 'r', '#0AC800', 'c', 'm']
    for i_, pn in enumerate(range(pn_0, pn_max)):
        print 'Plotting raster for %s pattern %d ' % (cell_type, pn)
        c = cm.jet(i_/ float(pn_max - pn_0 - 1), 1)
        c = cm.brg(.5 * i_/ float(pn_max - pn_0), 1)
        c = list(c)
        alpha = (1. - .3 * i_ / float(pn_max - pn_0))
        c = colorlist[i_]
        alpha = 1.
        print 'alpha', alpha
#        plot_raster_for_celltype(params, ax, cell_type, pn, title='Pyramidal cell spikes for patterns %d - %d' % (pn_0, pn_max-1), \
#        plot_raster_for_celltype(params, ax, cell_type, pn, title='Pattern activities for \nincreasing concentration', \
        plot_raster_for_celltype(params, ax, cell_type, pn, title='%s spikes patterns %d - %d' % (cell_type.capitalize(), pn_0, pn_max), \
                show=False, color=c, alpha=alpha)

    pylab.show()
