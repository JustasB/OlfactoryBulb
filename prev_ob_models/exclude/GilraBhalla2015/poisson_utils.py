import string
from pylab import *
#from scipy.interpolate import interp1d # needs liblapack, not available on cluster nodes
# Have just extracted the Poisson generating routines from neuroutils.py
# so that I don't need to import scipy and thus it will run on all gj nodes.

def poissonTrain(runtime, firingrate, refractory):
    lastfiredtime = 0 # s
    ornstimvector = []
    if firingrate > 0:
        while lastfiredtime < runtime:
            # numpy.random.exponential
            isi = exponential(1/firingrate) # inter-spike interval 'isi' of a Poisson spike generator is exponentially distributed r*exp(-r*isi)
            if isi < refractory:
                continue ###### Do not fire
            else:
                #### FIRED! append time of firing
                lastfiredtime += isi
                if lastfiredtime < runtime:
                    ornstimvector.append(lastfiredtime)
    return ornstimvector

def poissonTrainVaryingRate(runtime, maxfiringrate, refractory, tvec, ratevec):
    """
    This uses the spike thinning algorithm outlined in Dayan and Abbott 2001.
    firing rate should never exceed maxfiringrate. I typically give twice the estimated max to be safe.
    tvec and ratevec are lists t and rate respectively,
    from which I obtain a rate(t) function by linear interpolation (numpy).
    """
    lastfiredtime = 0 # s
    ornstimvector = []
    # Do not want to increase dependency on scipy.interpolate
    # as not installed on all nodes on cluster -- use numpy
    #rate = interp1d(tvec, ratevec, kind='linear')
    while lastfiredtime < runtime:
        # numpy.random.exponential
        # inter-spike interval 'isi' of a homogeneous (const rate 'r') Poisson spike generator
        # is exponentially distributed r*exp(-r*isi)
        isi = exponential(1/maxfiringrate)
        if isi < refractory:
            continue ###### Do not fire
        else:
            lastfiredtime += isi
            if lastfiredtime > runtime:
                break
            # spike thinning by rate(t)/maxfiringrate
            # numpy.random.uniform
            if uniform(0,1) < interp(lastfiredtime,tvec,ratevec)/maxfiringrate:
                #### FIRED! append time of firing
                ornstimvector.append(lastfiredtime)
    return ornstimvector
