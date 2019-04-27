import sys, pickle

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *
from simset_odor import *

from moose_utils import *
from neuro_utils import * # has the dual exponential functions
from data_utils import *

from pylab import * # part of matplotlib that depends on numpy but not scipy
import scipy.optimize
import scipy.interpolate
import scipy.special # special function erf()
import scipy.stats
## for scipy.signal.deconvolve but it doesn't work
## if any of the respiration pulse coeff's are zero - see below.
#import scipy.signal

def compute_odor_params(delay_side, risetime_side, duration_side):
    """ the args are the number of SDs away, you want the corresponding param to be
    first compute delay, risetime and duration. risetime and delay are correlated.
    then get the double sigmoid params for the generated delay, risetime and duration."""
    ######### generate delay, risetime and duration
    x1 = delay_side
    x2 = risetime_side
    delay = x1*delay_sd + delay_mean
    #### delay and rise_time are correlated normal variables
    risetime = ( delay_risetime_correlation*x1 + \
        sqrt(1-delay_risetime_correlation**2)*x2 ) * risetime_sd + risetime_mean
    duration = duration_side*duration_sd + duration_mean
    
    return delay, risetime, duration

def allcombo_stimuli():

    ## fig with subpanels for each ORN combo
    fig = figure(facecolor='w')
    delaxes() # delete the main axes to accomodate subpanels

    row_num = 0
    rows = 2*2*2
    cols = rows
    for exc_delay in [-2,2]:
        for exc_rise in [-2,2]:
            for exc_duration in [-2,2]:
                for exc_amplitude in [1]:
                    col_num = 0
                    for inh_delay in [-2,2]:
                        for inh_rise in [-2,2]:
                            for inh_duration in [-2,2]:
                                for inh_amplitude in [1]:
                                    exc_delay, exc_risetime, exc_duration = \
                                        compute_odor_params(exc_delay, exc_rise, exc_duration)
                                    inh_delay, inh_risetime, inh_duration = \
                                        compute_odor_params(inh_delay, inh_rise, inh_duration)
                                    ## WARNING: need to add air response too
                                    ## WARNING the decay slope is taken ad hoc twice the risetime!
                                    exc_response = [(0,0),\
                                        (exc_delay,0),\
                                        (exc_delay+exc_risetime,exc_amplitude*FIRINGMEANA),\
                                        (exc_delay+exc_duration,exc_amplitude*FIRINGMEANA),\
                                        (exc_delay+exc_duration+2*exc_risetime,0),\
                                        (RESPIRATION,0)]
                                    inh_response = [(0,0),\
                                        (inh_delay,0),\
                                        (inh_delay+inh_risetime,-inh_amplitude*FIRINGMEANA),\
                                        (inh_delay+inh_duration,-inh_amplitude*FIRINGMEANA),\
                                        (inh_delay+inh_duration+2*inh_risetime,0),\
                                        (RESPIRATION,0)]
                                    ax = fig.add_subplot(rows,cols,row_num*cols+col_num+1)
                                    x,y = zip(*exc_response) # unpack and zip!
                                    ax.plot(x,y)
                                    x,y = zip(*inh_response) # unpack and zip!
                                    ax.plot(x,y)
                                    ax.set_xlim(0,RESPIRATION)
                                    col_num += 1
                    row_num += 1
                    print "Finished row",row_num

if __name__ == "__main__":
    allcombo_stimuli()
    show()
