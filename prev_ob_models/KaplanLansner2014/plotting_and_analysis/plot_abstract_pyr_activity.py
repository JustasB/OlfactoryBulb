import pylab
import numpy as np
import sys
import os
from FigureCreator import plot_params


def apply_hc_wta(d, n_hc, n_mc):

    n_row = d.shape[0]
    d_out = np.zeros(d.shape)
    for i_ in xrange(n_row):
        for hc in xrange(n_hc):
            mc1 = hc * n_mc
            mc2 = (hc + 1) * n_mc
#            d_out[i_, mc1:mc2] = np.exp(d[i_, mc1:mc2] - d[i_, mc1:mc2].max())
            d_out[i_, mc1:mc2] = 0
            d_out[i_, np.argmax(d[i_, mc1:mc2]) + mc1] = 1.
#            summed_activity = d[i_, mc1:mc2].sum()
#            if (summed_activity > 1): # normalize
#            print 'debug hc: %d \tactivity: %.2e' % (hc, summed_activity)
#            d_out[i_, mc1:mc2] = np.exp(d[i_, mc1:mc2]) / summed_activity
#                d_out[i_, mc1:mc2] /= summed_activity
    return d_out

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print("Please enter data folder to be analysed after the script name")
        print("e.g.\npython analyse_data_folder.py data_today/ ")
        exit(1)
    else:
        fn = sys.argv[1]


    path = os.path.abspath(fn)

    print "loading data ...."
    data = np.loadtxt(path)

    if data.ndim == 1:
        data_new = np.zeros((2, data.shape[0]))
        data_new[0, :] = data
        data_new[1, :] = data
        data = data_new.copy()
        # if it's only one line:
        #d = np.loadtxt(path)
        #data = np.zeros((1, d.size))
        #for i in xrange(d.size):
        #    data[0, i] = d[i]


    pn_max = 10
    data = data[:pn_max, :]

    pylab.rcParams.update(plot_params)
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    print "plotting ...."
    cax = ax.pcolormesh(data)#, edgecolor='k', linewidths='1')
    pylab.ylim(0, data.shape[0])
    pylab.xlim(0, data.shape[1])
    cbar = pylab.colorbar(cax)
    clabel = 'Exponentially filtered activity'
    clabel = 'Abstract model activity'
    cbar.set_label(clabel)
#    title = 'Abstract model activity'
#    ax.set_title(title)
    ax.set_ylabel('Pattern number')
    ax.set_xlabel('Minicolumn')

    try:
        output_fn = sys.argv[2]
        print 'Saving to:', output_fn
        pylab.savefig(output_fn, dpi=300)
    except:
        pass

#    data_wta = apply_hc_wta(data, 12, 30)
#    pylab.rcParams.update(plot_params)
#    fig = pylab.figure()
#    ax = fig.add_subplot(111)
#    print "plotting ...."
#    cax = ax.pcolormesh(data_wta)#, edgecolor='k', linewidths='1')
#    pylab.ylim(0, data_wta.shape[0])
#    pylab.xlim(0, data_wta.shape[1])
#    cbar = pylab.colorbar(cax)
#    clabel = 'Exponentially filtered activity'
#    cbar.set_label(clabel)
#    title = 'Abstract model activity softmax'
#    ax.set_title(title)
#    ax.set_ylabel('Pattern number')
#    ax.set_xlabel('Minicolumn')

    pylab.show()
