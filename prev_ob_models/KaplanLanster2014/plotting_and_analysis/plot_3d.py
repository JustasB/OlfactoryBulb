import sys
import numpy as np
import pylab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

fn = sys.argv[1]
d = np.loadtxt(fn)

#x_axis_idx = 2
#y_axis_idx = 3
#x_label = '$\\tau_{z_i}$'
#y_label = '$v_{stim}$'
x_axis_idx = 0
y_axis_idx = 1
#x_label = '$dx$'
#y_label = '$dv$'
x_label = 'x'
y_label = 'y'

z_axis_idx = 2# w_avg at the end

x_data = d[:, x_axis_idx]
y_data = d[:, y_axis_idx]
z_data = d[:, z_axis_idx]

#if z_axis_idx == 4:
#    z_label = '$w_{max}$'
if z_axis_idx == 10:
    z_label = '$w_{avg, end}$'
else:
    z_label = 'z-axis'

fig = pylab.figure()
ax = Axes3D(fig)

colored = False
#colored = True

if (colored):
    color_code_axis = z_axis_idx
    colorbar_label = '$w_{avg, end}$'
    code = d[:, color_code_axis]
    min_4d = np.min(code)
    max_4d = np.max(code)
    range_4d = float(max_4d - min_4d)
    norm = matplotlib.mpl.colors.Normalize(vmin=min_4d, vmax=max_4d)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.jet)
    m.set_array(np.arange(min_4d, max_4d, 0.01))
    colors = m.to_rgba(code)
    cax = ax.scatter(x_data, y_data, z_data, c=colors, marker='o', linewidth='5', edgecolor=colors)

else:
    cax = ax.scatter(x_data, y_data, z_data, marker='o')

ax.set_title(fn)
ax.set_xlabel(x_label, fontsize=24)
ax.set_ylabel(y_label, fontsize=24)
ax.set_zlabel(z_label, fontsize=24)

if colored: 
    cb = fig.colorbar(m, ax=ax, shrink=0.8)
#    cb.set_label('$t_{max}$ [ms]', fontsize=24)
    cb.set_label(colorbar_label, fontsize=24)

pylab.show()
