
import pylab
import numpy as np
import sys
import os
import matplotlib.cm as cm

fn = sys.argv[1]

d = np.loadtxt(fn)

x_index = 0
color_index = 2
y_index = 3

fig = pylab.figure(figsize=(12,8))
ax = fig.add_subplot(111)

n_hcs = np.unique(d[:, x_index])
vq_overlaps = np.unique(d[:, color_index]).tolist()
print 'vq_overlaps', vq_overlaps
n_colors = float(len(vq_overlaps))

x_values = []
y_values = []
colors = []

legend_txt = []
plots = {}
for hc in n_hcs:
    idx = (d[:, x_index] == hc).nonzero()[0]
    y_vals = d[idx, y_index]
    for i_, val in enumerate(y_vals):
        y_val = d[idx[i_], y_index]
        z = d[idx[i_], color_index]
        cidx = vq_overlaps.index(z)
        color = cm.jet(cidx / n_colors, 1)
        x_values.append(hc)
        y_values.append(y_val)
        p, = ax.plot(hc, y_val, 'o', c=color)
        if not plots.has_key(cidx):
            plots[cidx] = p
            legend_txt.append('vq_overlap=%d' % z)
#        print p, 


ax.legend(plots.values(), legend_txt, loc='lower left')
ax.set_xlabel('N_HC')
ax.set_ylabel('correct classifications')
#pylab.legend()
pylab.show()
