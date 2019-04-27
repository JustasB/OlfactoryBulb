import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import pylab
import numpy as np 
import sys
import os

import simulation_parameters # defines simulation parameters
param_tool = simulation_parameters.parameter_storage()
params = param_tool.params

def load_time_axis():

    # loading the time axis doesn't work for some unknown reason
    # --> re-create the time axis
    t = np.arange(0, params['t_sim'], params['time_step_rec'])
    return t



num_files = len(sys.argv) - 1

if (num_files == 0):
    fn = raw_input("Please enter data file to be plotted\n")

#print sys.argv[1:]

fig = pylab.figure()
for fn in sys.argv[1:]:
    data = np.loadtxt(fn, skiprows=1)

    print 'Plotting', fn
    if (data.ndim == 1):
        try:
            x_axis = load_time_axis()
            go = x_axis.shape == data.shape
            assert (go == True)
        except:
            x_axis = np.arange(data.size)
        pylab.plot(x_axis, data, label=fn, lw=2)
#        pylab.plot(x_axis, data)
    else:
        pylab.plot(data[:,0], data[:,1], lw=2, label=fn)#, yerr=data[:,3])
        #pylab.title("%s for n_rows=1000, averaged over all processes" % (title))
        pylab.xticks(size='large')
        pylab.yticks(size='large')
        #pylab.savefig("%s_nrows1000.eps" % (title))

if len(sys.argv[1:]) < 10:
    pylab.legend()
pylab.show()
