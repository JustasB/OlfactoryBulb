import pylab
import numpy as np
import sys
import os
from FigureCreator import plot_params

def get_mean_max_min_median(d):
    print 'In row   mean    max     min     mean of smallest 10 percent   median'
    for i in xrange(d[:, 0].size):
        idx_sorted = d[i, :].argsort()
        lowest = d[i, idx_sorted[:10]]
        print i, d[i, :].mean(), d[i, :].max(), d[i, :].min(), lowest.mean(), np.median(d[i, :])
        

def get_units_above_thresh(d, thresh):

    if d.ndim == 1:
        above_thresh = (d > thresh).nonzero()[0]
        n_above_thresh = above_thresh.size
        n_zero = (d == 0).nonzero()[0].size
        if n_zero ==  d.size:
            n_all_zero = 1
    else:
        n_above_thresh = np.zeros(d[:, 0].size)
        n_all_zero = np.zeros(d[:, 0].size)
        for i in xrange(d[:, 0].size):
            above_thresh = (d[i, :] > thresh).nonzero()[0]
            n_above_thresh[i] = above_thresh.size
            n_zero = (d[i, :] == 0).nonzero()[0].size
            if n_zero ==  d[i, :].size:
                n_all_zero[i] = 1

#        print 'In row %d: %d units have value higher than %.2e' % (i, n_above_thresh[i], thresh), above_thresh
        print '%d rows have all elements filled with zero' % (n_all_zero[i].sum())
        print 'On average, per row %.2f +- %.2f units have higher values than %.2e, i.e. %.1f +- %.1f percent' % (n_above_thresh.mean(), n_above_thresh.std(), thresh, \
                n_above_thresh.mean() / d[0, :].size * 100., n_above_thresh.std() / d[0, :].size * 100.)
        print 'In total %d units have higher values than %.2e' % (n_above_thresh.sum(), thresh)

def plot_hist(d, n_bins=20, title=''):

    count, bins = np.histogram(d, bins=n_bins)
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.bar(bins[:-1], count, width=bins[1]-bins[0])
    ax.set_title(title)
#    pass


def plot_sorted_values_as_hist(d, title=''):

    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(d.size), d, width=1)
    ax.set_title(title)


def plot_wta(d):
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    print "plotting ...."
    data = np.zeros(d.shape)
    n_correct = 0
    for row in xrange(d.shape[0]):
        winner = d[row, :].argmax()
        data[row, winner] = 1
        if row == winner:
            n_correct += 1

    print 'n_correct: ', n_correct
    cax = ax.pcolormesh(data, cmap='binary')
    pylab.colorbar(cax)


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print("Please enter data folder to be analysed after the script name")
        print("e.g.\npython analyse_data_folder.py data_today/ ")
        exit(1)
    else:
        fn = sys.argv[1]


    path = os.path.abspath(fn)

    print "loading data ...."

    # if it's only one line:
    #d = np.loadtxt(path)
    #data = np.zeros((1, d.size))
    #for i in xrange(d.size):
    #    data[0, i] = d[i]

    try:
        data = np.loadtxt(path, delimiter=",")#.transpose()
    except:
        data = np.loadtxt(path)

    print 'debug data shape', data.shape
    if data.ndim == 1:
        data_new = np.zeros((2, data.shape[0]))
        data_new[0, :] = data
        data_new[1, :] = data
        data = data_new.copy()
    print 'debug data shape', data.shape

    #get_mean_max_min_median(data)
    thresh = 0.001
    get_units_above_thresh(data, thresh)
    thresh = 0.1
    get_units_above_thresh(data, thresh)
#    thresh = -1.5
#    get_units_above_thresh(data, thresh)

#    plot_wta(data)
#    plot_hist(data, n_bins=40)
    print 'Sum of all elements:', data.sum()
#    for OR in xrange(0, 5):
#        d = data[:, OR]
#        idx_sorted = np.argsort(d)
#        plot_sorted_values_as_hist(d[idx_sorted], title='OR %d' % OR)

    try: # if you want to take the log of the data to plot
        if sys.argv[2] == 'log':
            n_row = data[:, 0].size
            n_col = data[0, :].size
            log_data = np.zeros((n_row, n_col))
            for i in xrange(n_row):
                for j in xrange(n_col):
                    if data[i, j] > 0:
                        log_data[i, j] = np.log(data[i, j])
            fig = pylab.figure()
            ax = fig.add_subplot(111)
            ax.set_title('log data')
            cax = ax.pcolormesh(log_data)#, edgecolor='k', linewidths='1')
            pylab.ylim(0, data.shape[0])
            pylab.xlim(0, data.shape[1])
            pylab.colorbar(cax)
    except:
        pass

    #data_rev = np.zeros(data.shape)
    #n_row = data[:, 0].size - 1
    #for row in xrange(data[:, 0].size):
    #    data_rev[n_row - row, :] = data[row, :]


    pylab.rcParams.update(plot_params)
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    print "plotting ...."
    #cax = ax.imshow(data[:,:12])
    #cax = ax.pcolor(data, edgecolor='k', linewidths='1')

    cax = ax.pcolormesh(data)#, edgecolor='k', linewidths='1')
    #cax = ax.pcolor(data, cmap='binary')
    #cax = ax.pcolor(data, cmap='RdBu')

    #cax = ax.pcolor(log_data)#, edgecolor='k', linewidths='1')


    pylab.ylim(0, data.shape[0])
    pylab.xlim(0, data.shape[1])
    pylab.colorbar(cax)

    try: 
        plot_fn = sys.argv[2]
        print "saving ....", plot_fn
        pylab.savefig(plot_fn)
        title = plot_fn.rsplit('/')[-1]
    except:
        title = fn.rsplit('/')[-1]

    ax.set_title(title)
    pylab.show()
