import sys, pickle

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *
from simset_odor import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

from data_utils import * # has axes_labels()

## USAGE: python2.6 plot_firerates_odors.py <firefiles_name>
## EXAMPLE: python2.6 plot_firefiles.py ../firefiles/firefiles760.0/firetimes_rndpulse_glom_0_pulse_5_avgnum0.txt

#### We have dual exponential functions as impulse response / kernels
#### for respiration, odor A and odor B (each is different for every glomerulus)
#### These kernels are convolved with the respiration pulse
#### to generate the firing rate which is fed to a Poissonian generator.

RUNTIME = PULSE_RUNTIME + SETTLETIME

## time points for the firing rate which is read from a pickled file
firingtsteps = arange(0,RUNTIME+1e-10,FIRINGFILLDT)# include the last RUNTIME point also.
numt = len(firingtsteps)
extratime = arange(0,3*RESPIRATION+1e-10,FIRINGFILLDT)
randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)

if __name__ == "__main__":
    fname = sys.argv[1]
    f = open(fname,'r')
    lines = f.read() ## reads everything
    f.close()
    strtimes = lines.split()
    firetimes = [ float(strtime) for strtime in strtimes ]
    firetimes.sort()
    
    bindt = 25e-3
    if 'pulse' in fname:
        bintimes = arange(0.0,PULSE_RUNTIME,bindt)
        binfrate = array(plotBins(firetimes, int(round(PULSE_RUNTIME/bindt)),\
            PULSE_RUNTIME+SETTLETIME, SETTLETIME)) / NUM_ORN_FILES_PER_GLOM
    else:
        bintimes = arange(0.0,ODORRUNTIME-SETTLETIME,bindt)
        binfrate = array(plotBins(firetimes, int(round((ODORRUNTIME-SETTLETIME)/bindt)),\
            ODORRUNTIME, SETTLETIME)) / NUM_ORN_FILES_PER_GLOM

    fig = figure(facecolor='w')
    ax = fig.add_subplot(111)
    title(fname, fontsize=24)
    frateperglomList = []
    plot(bintimes[:len(binfrate)], binfrate, '-r', linewidth=2, marker=',')
    axes_labels(ax,'time (s)','ORN firing rate (Hz)',fontsize=24)
    
    show()
