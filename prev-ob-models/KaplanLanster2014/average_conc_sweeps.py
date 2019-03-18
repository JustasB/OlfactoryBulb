import os
import sys
import numpy as np
#import matplotlib
#matplotlib.use("WXAgg") 
import pylab
from scipy import interpolate
import MergeSpikefiles 


if (len(sys.argv) < 2):
    info_txt = 'Please give folder names to average over after the script name:'
    info_txt =+ '\nUsage:\n\tpython average_conc_sweeps.py [FOLDER_1] [FOLDER_2] [...]'
    print info_txt
    exit(1)
else:
    from_folders = sys.argv[1:]



if len(sys.argv) > 1:
    param_fn = os.path.abspath(from_folders[0]) + '/Parameters/simulation_parameters.json'
    import json
    f = file(param_fn, 'r')
    print 'Loading parameters from', param_fn
    params = json.load(f)

else:
    import simulation_parameters
    param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params


Merger = MergeSpikefiles.MergeSpikefiles(params)
cell_type = 'mit'
n_cells = params['n_%s' % cell_type]
offset = params['%s_offset' % cell_type]

#print "debug,", params['mit_spike_fn_base'].rsplit('/')
nspikes_fn_base = params['mit_spike_fn_base'].rsplit('/')[1] + '/' + params['mit_spike_fn_base'].rsplit('/')[2]
nspikes_merged_fn_base = params['mit_spikes_merged_fn_base'].rsplit('/')[1] + '/' + params['mit_spikes_merged_fn_base'].rsplit('/')[2]
print "Taking the following merged MT spike files :"
pn = 0
file_list = []
n_data = []
for folder in from_folders:
    fn = os.path.abspath(folder) + '/' + nspikes_merged_fn_base + str(pn) + '.dat'
    if not (os.path.exists(fn)):
        print "Merging to", fn
        Merger.merge_nspike_files(folder + '/' + nspikes_fn_base, os.path.abspath(folder) + '/' + nspikes_merged_fn_base, pn)

    data_new = np.zeros((n_cells, 2))
    data = np.loadtxt(fn)
    for i in xrange(data.shape[0]):
        gid = data[i, 0]
        index = gid - offset
        data_new[index, 1] = data[i, 1]

    if (folder[-1] == '/'):
        fn_new = folder + nspikes_merged_fn_base + 'full'
    else:
        fn_new = folder + '/' + nspikes_merged_fn_base + 'full'
    print 'Saving data to:', fn_new
    np.savetxt(fn_new, data_new)
    file_list.append(fn_new)
    # check if all files contain same amount of data
    d = np.loadtxt(fn)
    n_data.append(d.size)

#assert (np.unique(n_data) != 1), "Data files have different sizes! Not able to average, filelist %s" % file_list

n_data = n_data[0]
n_conc = params['n_mit_y'] # the concentration axis
n_x = params['n_mit_x']
n_files = len(from_folders)
all_data = np.zeros((n_x, n_conc, n_files))
sum_data = np.zeros((n_x, n_conc))

print "File_list:"
for file_id in xrange(len(file_list)):
    fn = file_list[file_id]
    print fn

y_scale = params['t_sim'] / 1000.
# (params["t_stop"]-params["t_start"]) / 1000.

# collect all_data 
for file_id in xrange(len(file_list)):
    # assume the data is sorted according to cell gids
    # gid array -> gor / curve_id, | conc
    fn = file_list[file_id]
    d = np.loadtxt(fn)
    assert ((n_x * n_conc) == d[:,0].size), "Data file %s has wrong number of data points, check if network_parameters contains right parameters!" % (fn)
    for conc in xrange(n_conc):
        sum_data[:, conc] += d[(conc * params['n_mit_x']):((conc + 1) * params['n_mit_x']), 1]
        all_data[:, conc, file_id] = d[(conc * params['n_mit_x']):((conc + 1) * params['n_mit_x']), 1]
#        print "debug all_data[:, %d, %d] :" % (conc, file_id), all_data[:, conc, file_id]
#        print "debug all_data[%d, %d,:] :" % (curve, conc), all_data[curve, conc, :]


# average over all files
average_curve = np.zeros((n_x, n_conc))
stds = np.zeros((n_x, n_conc)) # the standard deviations
y_max = 0.
for curve in xrange(params['n_mit_x']):
#    average_curve[curve, :] = sum_data[curve,:] / float(n_files)
    for conc in xrange(n_conc):
        stds[curve, conc] = all_data[curve, conc, :].std() / np.sqrt(n_files)
#        print "debug all_data[%d, %d,:] :" % (curve, conc), all_data[curve, conc, :]
#        print "debug std[%d, %d] = %f" % (curve, conc, stds[curve, conc])
        average_curve[curve, conc] = all_data[curve, conc, :].mean()
        average_curve[curve, conc] /= y_scale
    y_max = max(y_max, max(average_curve[curve, :]))
    

# load MT parameters to get the actual concentration of the stimulus
mit_params = np.loadtxt(params['mit_params_fn_base'] + str(0) + '.dat', skiprows=1)
x_axis = np.unique(mit_params[:, 3]) # concentration axis

# Plotting
def get_figsize(fig_width_pt):
    inches_per_pt = 1.0/72.0                # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

from FigureCreator import plot_params
pylab.rcParams.update(plot_params)
colorlist = ["#0000FF", "#006600", "#FF0000", "#00FFFF", "#CC00FF", "#FFCC00", "#000000", "#00FF00", "#663300", "#FF3399", "#66CCCC", "#FFCC99", "#666666"]

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')
#ax_2 = ax.twinx()

for curve in xrange(params['n_mit_x']):
#    ax.plot(x_axis, average_curve[curve,:], 
    ax.errorbar(x_axis, average_curve[curve,:], yerr=stds[curve,:], color=colorlist[curve%len(colorlist)], lw=3)

#output_fn = "ob_conc_sweep_average_cray_%d.dat" % params['date']
output_fn = params['other_folder'] + '/averaged_ob_conc_sweep.dat'
print "Saving data to:", output_fn
output_data = "# row 0: x_axis, row 1: mean curve 1, row 2: std, row 3: mean_curve2, row 4: std, ..."
output_data += "# n_x = %d\t n_curves = %d\n" % (n_conc, params['n_mit_x'])

for x in x_axis:
    output_data += "%.6e\t" % x
output_data += "\n"

#for x in xrange(average_curve[0,:].size):
#    output_data += "%.6e\t" % (average_curve[:, x].mean())
#output_data += "\n"

# print the average curves and the std of that curve in the next line
for curve in xrange(params['n_mit_x']):
    for x in xrange(n_conc):
        output_data += "%.6e\t" % (average_curve[curve, x])
    output_data += "\n"
    for x in xrange(n_conc):
        output_data += "%.6e\t" % (stds[curve, x])
    output_data += "\n"
output_data += "\n"

output_data_file = file(output_fn, 'w')
output_data_file.write(output_data)
output_data_file.close()

# spline curves
num_curves = params['n_mit_x']
color_cnt = 0
x_start = 2 # start fitting splines from this x-value upwards
x_stop = params['n_or']

ax.set_ylim((0, ax.get_ylim()[1]))
#ax.set_ylim((0, 35))

#ax_2.set_ylim((0, y_max / ((params["t_stop"]-params["t_start"]) / 1000.)))
#ax_2.set_ylim((0, y_max / (params["t_sim"] / 1000.)))
ax.set_title('Mitral cell response curves')
ax.set_xlim((x_axis.min(), x_axis.max()))
ax.set_xlabel("Concentration [a.u.]")
ax.set_ylabel("Output rate [Hz]")
#ax.set_xlabel("x")
#ax.set_ylabel("y1")
#ax_2.set_ylabel("y2")
#F = pylab.gcf()

#F.set_size_inches(get_figsize(800))
#output_fn = "ob_set_of_curves_%dngor.svg" % (params['n_gor'])
output_fn = params['figure_folder'] + '/' + 'averaged_ob_conc_sweep.png'
print 'Saving figure to:', output_fn
pylab.savefig(output_fn, dpi=(200))
output_fn = params['figure_folder'] + '/' + 'averaged_conc_sweep.pdf'
print 'Saving figure to:', output_fn
pylab.savefig(output_fn, dpi=(200))
pylab.show()


