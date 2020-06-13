import numpy
import pylab
import sys

fn = sys.argv[1]

title = fn

n_bins = 100
#d = numpy.loadtxt(fn, skiprows=1)
d = numpy.loadtxt(fn)


d_thresh = 1e-1
d_new = []

d = d.flatten()
#log_d = numpy.zeros(d.size)
#log_d = []
#for i in xrange(d.size):
#    if d[i] > 0:
#        log_d.append(numpy.log(d[i]))
#    if abs(d[i]) > d_thresh:
#        d_new.append(d[i])
#log_d = numpy.array(log_d)

# d_new
#d_new = numpy.array(d_new)
#fig = pylab.figure()
#x_max = d_new.max()
#bin_width = (d_new.max() - d_new.min()) / n_bins
#data, bins = numpy.histogram(d_new, n_bins)
#ax = fig.add_subplot(111)
#ax.bar(bins[:-1], data, width=bin_width)
#pylab.title("Bin width = %f" % (bin_width))
#pylab.ylabel("Count")
#pylab.title("w_thresh = %.1e" % d_thresh)


# log
#fig = pylab.figure()
#x_max = log_d.max()
#bin_width = (log_d.max() - log_d.min()) / n_bins
#data, bins = numpy.histogram(log_d, n_bins)
#ax = fig.add_subplot(111)
#ax.bar(bins[:-1], data, width=bin_width)
#pylab.title("Bin width = %f" % (bin_width))
#pylab.ylabel("Count")
#pylab.title("log(weights) (if w > 0)")

#pylab.show()


# d normal
fig = pylab.figure()
x_max = d.max()
data, bins = numpy.histogram(d, n_bins)
bin_width = bins[1] - bins[0]
ax = fig.add_subplot(111)
ax.bar(bins[:-1], data, width=bin_width)
pylab.title(fn)
pylab.ylabel("Count")
pylab.show()

