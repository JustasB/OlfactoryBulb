"""
This script plots the MT cells in the MI-space (mds_ob_oc_output_fn) generated by MDSVQ.py
and colors all mitral cells projecting to one target Hypercolumn.
The mitral cell -> Hypercolumn mapping is stored in vq_ob_oc_output_fn

Can only be called after running CreateObOcConnections.py
"""

import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import simulation_parameters
import numpy as np
import pylab
import matplotlib
from mpl_toolkits.mplot3d import Axes3D


class Plotter(object):

    def __init__(self, params):
        self.params = params


    def plot_mt_cells_to_hc(self, hc=0):

        mt_pos = np.loadtxt(self.params['mds_ob_oc_output_fn'])
        prj_matrix = np.loadtxt(self.params['vq_ob_oc_output_fn'])
        fig = pylab.figure()
        ax = Axes3D(fig)

        mt_projecting_to_hc = prj_matrix[:, hc].nonzero()[0]
        print 'n mt_projecting_to_hc', mt_projecting_to_hc.size

        # map centroid numbers to different colors

        norm = matplotlib.mpl.colors.Normalize(vmin=0, vmax=1)
        m = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.RdBu)
        m.set_array(np.arange(0, 1, 0.01))
        colors = m.to_rgba(prj_matrix[:, hc])
        cax = ax.scatter(mt_pos[:,0], mt_pos[:,1], mt_pos[:,2], c=colors, marker='o', linewidth='5', edgecolor=colors)

        output_fig = self.params['figure_folder'] + '/prj_mt_hc%d_ovrlp%d.png' % (hc, self.params['vq_ob_oc_overlap'])
        print 'Saving to:', output_fig
        pylab.savefig(output_fig)
        pylab.show()



if __name__ == '__main__':
    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        param_tool = simulation_parameters.parameter_storage()

    params = param_tool.params
    P = Plotter(params)
    P.plot_mt_cells_to_hc(hc=2)

