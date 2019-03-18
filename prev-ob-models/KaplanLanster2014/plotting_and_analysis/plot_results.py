import pylab
import numpy
import sys

if (len(sys.argv) < 2):
    fn = raw_input("Please enter data file to be plotted\n")
else:
    fn = sys.argv[1]

data = np.loadtxt(fn)

# if the first line contains crap use skiprows=1
#data = np.loadtxt(fn, skiprows=1)

fig = pylab.figure()
ax = fig.add_subplot(111) 

# if you want to use multiple figures in one, use
#ax1 = fig.add_subplot(211) 
#ax2 = fig.add_subplot(212) 
# and 

if (data.ndim == 1):
    x_axis = numpy.arange(data.size)
    ax.plot(x_axis, data)
else:

     
#    ax.errorbar(data[:,0], data[:,1], yerr=data[:, 2])
#    print 'mean y-value:', data[:, 1].mean()
    ax.plot(data[:, 0], data[:, 1], ls='-', lw=3, c='b')
#    ax.scatter(data[:,0], data[:,2])
#    ax.plot(data[:,3], data[:,6])

# saving:
# fig.savefig('output_figure.png') 

# otherwise nothing is shown
pylab.show()
