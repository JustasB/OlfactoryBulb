
import pylab
import numpy as np
import sys

if (len(sys.argv) < 2):
    fn = raw_input("Please enter data file to be plotted\n")
else:
    fn = sys.argv[1]

rcp = {'figure.subplot.wspace': 0.45, 
    'figure.subplot.hspace':.25, 
    'figure.subplot.top':.85, 
    'figure.subplot.left':.16, 
    'figure.subplot.rigt':.96, 
     'axes.labelsize' : 16,
    'label.fontsize': 14,
    'xtick.labelsize' : 12, 
    'ytick.labelsize' : 12, 
    'axes.titlesize'  : 16,
    'legend.fontsize': 9
        }
pylab.rcParams.update(rcp)

fig = pylab.figure()
fig.suptitle('Parameters of Gaussian fits to OR-distance distributions', fontsize=20)

ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)
data = pylab.loadtxt(fn)



ax1.scatter(data[:, 0], data[:, 1])
ax2.scatter(data[:, 0], data[:, 2])
ax3.scatter(data[:, 0], data[:, 3])
ax4.scatter(data[:, 0], data[:, 4])
ax5.scatter(data[:, 0], data[:, 5])
ax6.scatter(data[:, 0], data[:, 6])
ax7.scatter(data[:, 0], data[:, 7])
ax8.scatter(data[:, 0], data[:, 8])
ax9.scatter(data[:, 0], data[:, 9])


ax1.set_ylabel('$w_1$')
ax4.set_ylabel('$w_2$')
ax7.set_ylabel('$w_3$')
ax2.set_ylabel('$\mu_1$')
ax5.set_ylabel('$\mu_2$')
ax8.set_ylabel('$\mu_3$')
ax3.set_ylabel('$\sigma_1$')
ax6.set_ylabel('$\sigma_2$')
ax9.set_ylabel('$\sigma_3$')

#    fit_params_txt = '#n_OR\tw0\tmu1\tsigma1\tw1\tmu1\tsigma1\tw2\tmu2\tsigma2\n'

for ax in fig.axes:
    ax.set_xlim((0, 70))

ax7.set_xlabel('n_OR')
ax8.set_xlabel('n_OR')
ax9.set_xlabel('n_OR')

ax1.set_title('Weights')
ax2.set_title('Mean values')
ax3.set_title('Widths')


pylab.show()


