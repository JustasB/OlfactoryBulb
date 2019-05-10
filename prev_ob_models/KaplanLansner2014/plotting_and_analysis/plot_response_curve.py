import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import simulation_parameters
import MergeSpikefiles
import SetOfCurvesPlotter

info_txt = \
"""
Usage:
    python plot_response_curve.py [FOLDER] [CELLTYPE]
"""
assert (len(sys.argv) > 2), 'ERROR: folder and cell_type not given\n' + info_txt
folder = sys.argv[1]
cell_type = sys.argv[2]

params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
params = param_tool.params

pn = 0
Merger = MergeSpikefiles.MergeSpikefiles(params)
if cell_type == 'mit':
    Merger.merge_ob_spiketimes_file(pattern=pn)
    Merger.merge_ob_nspike_files(pattern=pn)
elif cell_type == 'orn':
    Merger.merge_epth_spiketimes_file(pattern=pn)
    Merger.merge_epth_nspike_files(pattern=pn)
else:
    print 'ERROR: Invalid cell type given %s\n' % cell_type
    print 'Use lower case, e.g. orn or mit for cell_type'


sim_cnt = 0
SOCP = SetOfCurvesPlotter.SetOfCurvesPlotter(params)
if cell_type == 'mit':
    output_fn = params['figure_folder'] + '/ob_response_curve_%d.png' % sim_cnt
else:
    output_fn = params['figure_folder'] + '/epth_response_curve_%d.png' % sim_cnt
SOCP.plot_set_of_curves(pn=sim_cnt, output_fn=output_fn, cell_type=cell_type)
#import pylab
#pylab.show()
print 'Opening with ristretto: %s' % (output_fn)
os.system('ristretto %s' % output_fn)

# ------- Merge spike files
