"""
This script plots the result from cluster_odorant_space.py, 
the so-called F-Test for placing the ORs in the odorant space.
http://en.wikipedia.org/wiki/F-test

The F-value plotted here is the between-group-variability divided by 
the within-group-variability, where a group is the odorants assigned to the same cluster center i.e. an OR.


"""
import pylab
import numpy as np
import sys

if (len(sys.argv) < 2):
    fn = raw_input("Please enter data file to be plotted\n")
else:
    fn = sys.argv[1]

# --------------------------------------------------------------------------
import math
def get_figsize(fig_width_pt):
    inches_per_pt = 1.0/72.0                # Convert pt to inch
    golden_mean = (math.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

params2 = {'backend': 'eps',
          'axes.labelsize': 12,
          'text.fontsize': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'legend.pad': 0.2,     # empty space around the legend box
          'legend.fontsize': 12,
           'lines.markersize': 0,
          'font.size': 12,
          'path.simplify': False,
          'figure.figsize': get_figsize(800)}

def set_figsize(fig_width_pt):
    pylab.rcParams['figure.figsize'] = get_figsize(fig_width_pt)

pylab.rcParams.update(params2)

fig = pylab.figure()
ax = fig.add_subplot(111)
# --------------------------------------------------------------------------
data = pylab.loadtxt(fn)

variance = data[:, 1].mean()

#bgv = between group variability or variance
bgv = data[:, 3]
F = data[:, 5]



# average the data for each OR
col = 5
x_axis = np.unique(data[:, 0])
y_axis = np.zeros(x_axis.size)
y_std = np.zeros(x_axis.size)
for i_, OR in enumerate(x_axis):
    idx = (data[:, 0] == OR).nonzero()
#    data[:, col] == OR
#    print 'debug', OR, idx, data[:, col]
    y_axis[i_] = data[idx, col].mean()
    y_std[i_] = data[idx, col].std()

#ax.plot(x_axis, y_axis)
ax.errorbar(x_axis, y_axis, yerr=y_std)
ax.set_xlabel('Number of ORs put as centroid clusters')
ax.set_ylabel('F = explained / unexplained variance')
#ax.scatter(data[:,0], F, color='k')

#pylab.savefig(fn.rsplit(".")[0] + ".eps", dpi=200)


pylab.show()
