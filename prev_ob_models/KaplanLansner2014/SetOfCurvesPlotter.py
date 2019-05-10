import numpy as np
import pylab
from FigureCreator import plot_params
import os
import sys


class SetOfCurvesPlotter(object):

    def __init__(self, params):
        self.params = params


    def plot_set_of_curves(self, pn=0, output_fn=None, cell_type='orn', ax=None):
        """
        Keyword arguments (optional):
        pn -- (int) to identify the parameter and spike files to load
        output_fn -- figure file name
        cell_type -- 'orn' or 'mit'
        fig = a pylab.subplot instance, if the class should paint in a certain subplot
        """

        self.cell_type = cell_type
        self.pattern_nr = pn
        if cell_type == 'orn':
            n_per_group = self.params['rel_orn_mit']
            shared_param_idx = 3
            x_idx = 2 # 1 = oor, 2 = c/Kd
        else:
            n_per_group = 1
            shared_param_idx = 4
            x_idx = 3

        self.gid_idx = 0 # index in parameter file which contains the cell gids

        x_data, y_data = self.get_xy_data(shared_param_idx, x_idx)
        y_avg, y_std = self.get_average_group_curve(x_data, y_data, n_per_group)

        colorlist = ["#0000FF", "#006600", "#FF0000", "#00FFFF", "#CC00FF", "#FFFF00", "#000000", \
                "#00FF00", "#663300", "#FF3399", "#66CCCC", "#FFCC99", "#FFFFCC"]
        
        if ax == None:
            pylab.rcParams.update(plot_params)
            fig = pylab.figure()
            ax = fig.add_subplot(111)
        color_group_idx = 0
        n_cells = len(x_data)
        lw = 1
        plots, labels = [], []
        for i_ in xrange(n_cells):
            if cell_type == 'orn':
                color_group_idx = i_ / n_per_group
                color = colorlist[color_group_idx % len(colorlist)]
                p, = ax.plot(x_data[i_], y_data[i_], c=color, lw=lw)
                plots.append(p)
                labels.append('%d' % color_group_idx)
            else:
                ax.plot(x_data[i_], y_data[i_], label='%d' % i_, lw=lw)
#            ax.set_xscale('log')

        for group in xrange(len(y_avg)):
            color = colorlist[group % len(colorlist)]
            ax.errorbar(x_data[0], y_avg[group], yerr=y_std[group], c=color, lw=4)

        ax.set_title('%s response curves' % cell_type.upper())

        ax.set_xlabel('Concentration [a.u.]')
        ax.set_ylabel('Output rate [Hz]')
#        ax.legend(plots, labels, {'legend.loc' : 'upper left', 'legend.fontsize':8})
        if cell_type == 'mit':
            ax.legend(plots, labels, loc='upper left')

        ax.set_ylim((0, ax.get_ylim()[1]))
        ax.set_xlim((1e-2, 1.))

        self.ax = ax
        if output_fn != None:
            pylab.savefig(output_fn)
        else:
            pylab.show()
        return x_data, y_data


    def get_average_group_curve(self, x_data, y_data, n_per_group):
        """
        A group is a set of ORNs projecting to the same MIT cell.
        Keyword arguments:
        x_data, y_data  --  the response curves of single ORNs in lists, obtained from get_xy_data
        n_per_group -- (int) should be params['rel_orn_mit']
        """

        n_total = len(y_data)
        n_groups = n_total / n_per_group
        x_dim = x_data[0].size
        y_avg = np.zeros((x_dim, n_groups))
        y_std = np.zeros((x_dim, n_groups))

        for group_idx in xrange(n_groups):
            i0 = group_idx * n_per_group
            i1 = (group_idx + 1) * n_per_group
            for x_idx in xrange(x_dim):
                y_avg[x_idx, group_idx] = y_data[i0:i1, x_idx].mean()
                y_std[x_idx, group_idx] = y_data[i0:i1, x_idx].std()

        y_std = y_std.transpose()
        y_avg = y_avg.transpose()
        print 'y avg std', len(y_std)
        return y_avg, y_std


    def annotate(self, info_txt):

        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        
        x = xlim[0] + (xlim[1] - xlim[0]) * .0005
        y = ylim[0] + (ylim[1] - ylim[0]) * .8
        self.ax.annotate(info_txt, (x, y))

    def get_xy_data(self, shared_param_idx, x_idx):


        print 'DEBUG', self.pattern_nr, type(self.pattern_nr), type(str(self.cell_type))
        param_fn = self.params['%s_params_fn_base' % (self.cell_type)] + '%d.dat' % self.pattern_nr
        print 'Loading cell parameters from:', param_fn
        import os
        print 'debug CHECK', os.path.exists(param_fn)
        cell_params = np.loadtxt(param_fn, skiprows=1)

        nspike_fn = self.params['%s_spikes_merged_fn_base' % self.cell_type] + '%d.dat' % self.pattern_nr
        print 'Loading nspike data from:', nspike_fn
        self.nspikes = np.loadtxt(nspike_fn)


        # get GIDs that have a certain parameter in common
        group_params = np.unique(cell_params[:, shared_param_idx])
        n_groups = group_params.size
        # the data to be plotted
        x_data = [[] for i in xrange(n_groups)]
        y_data = [[] for i in xrange(n_groups)]

        for i_, val in enumerate(group_params):
            group_idx = np.array(cell_params[:, shared_param_idx] == val).nonzero()[0]
            x_values_for_group = cell_params[group_idx, x_idx]
            gids = cell_params[group_idx, self.gid_idx]
            for j_, gid in enumerate(gids):
                f_out = self.get_nspikes(gid) / (self.params['t_sim'] / 1000.)
                y_data[i_].append(f_out)
                x_data[i_].append(x_values_for_group[j_])

        x_data = np.array(x_data)
        y_data = np.array(y_data)
        return x_data, y_data
        


    def get_nspikes(self, gid):
        """
        Assuming that GIDs are stored in column 0 in self.nspikes the number of spikes fired by cell gid is returned.
        """
        try:
            idx = self.nspikes[:, 0] == gid
            return self.nspikes[idx, 1][0]
        except:
            return 0

        


if __name__ == '__main__':

    if len(sys.argv) > 1:
        param_fn = sys.argv[1]
        if os.path.isdir(param_fn):
            param_fn += '/Parameters/simulation_parameters.json'
        import json
        f = file(param_fn, 'r')
        print 'Loading parameters from', param_fn
        params = json.load(f)

    else:
        import simulation_parameters
        param_tool = simulation_parameters.parameter_storage()
        params = param_tool.params


    sim_cnt = 0
    import MergeSpikefiles
    Merger = MergeSpikefiles.MergeSpikefiles(params)
    Merger.merge_ob_spiketimes_file(pattern=sim_cnt)
    Merger.merge_ob_nspike_files(pattern=sim_cnt)

    SOCP = SetOfCurvesPlotter(params)
    output_fn = params['figure_folder'] + '/ob_response_curve_%d.png' % sim_cnt
    SOCP.plot_set_of_curves(pn=sim_cnt, output_fn=output_fn, cell_type='mit')
    print 'Opening with ristretto: %s' % (output_fn)
    os.system('ristretto %s' % output_fn)
