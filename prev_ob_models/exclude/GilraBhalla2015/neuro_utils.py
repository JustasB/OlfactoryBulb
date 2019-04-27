import string
from pylab import *
#from scipy.interpolate import interp1d # needs liblapack, not available on cluster nodes
from scipy.optimize import fsolve, fmin_tnc
import scipy.stats

#### This returns a difference of AlphaGaussians -- hybrid between alpha-fn and Gaussian.
#### Multiply Gaussian by (t/tp) so that it starts from zero.
def subtractedAlphaGaussians(t, peakp, peakm, tp, tm, sigmap, sigmam):
    return peakp * (t/tp)**(tp/2.0/sigmap) * np.exp(-(t-tp)**2/(2*sigmap**2)) \
        - peakm * (t/tm)**(tp/2.0/sigmam) * np.exp(-(t-tm)**2/(2*sigmam**2))
    #return peakp # just a const firing rate as odor response

def subtractedGaussians(t, peakp, peakm, tp, tm, sigmap, sigmam):
    return peakp * np.exp(-(t-tp)**2/(2*sigmap**2)) \
        - peakm * np.exp(-(t-tm)**2/(2*sigmam**2))
    #return peakp # just a const firing rate as odor response

def trunc_stdnormal():
    randomnum = standard_normal()
    if randomnum<0.0:
        randomnum = 0.0
    return randomnum

def gamma_shape_scale(mean,sd):
    theta = sd**2/mean
    return mean/theta, theta

############################### All about dual exponential functions ################################

def dualexpfunction(t, delay, tau1, tau2):
    if t<delay:
        return 0.0
    else:
        if tau1!=tau2:
            return ( math.exp(-(t-delay)/tau1)-math.exp(-(t-delay)/tau2) ) / (tau1-tau2)
        else:
            return ((t-delay)/tau1)*math.exp(1-(t-delay)/tau1)

def dualexp_norm(tau1, tau2):
    """
    dualexp_norm() returns the factor that must multiply the dualexpfunction() above to ensure that it has a peak of 1.0 .
    The peak of the dualexpfunction() is the multiplicative inverse of this factor returned by norm_calc().
    """
    if tau1==tau2: return 1.0
    else: return (tau2-tau1)/( (tau2/tau1)**(tau1/(tau1-tau2)) - (tau2/tau1)**(tau2/(tau1-tau2)) )

def dualexp_tpeak(tau1,tau2):
    if tau1==tau2: return tau1 # same as alpha function
    else: return tau1*tau2*log(tau1/tau2) / (tau1-tau2)

def invert_2exp(tau1,tau2,matchfraction,risephase=True):
    """ return the time t at which dualexpfunction equals matchfraction*peak.
    It returns the initial t during rising phase if risephase is True, else later t during falling phase.
    It is upto you to provide matchfraction that lies between 0.0 and 1.0 ."""
    peak = 1.0/dualexp_norm(tau1,tau2)
    tp = dualexp_tpeak(tau1,tau2)
    if risephase: init_t = matchfraction*tp
    else: init_t = tp/matchfraction
    matchvalue = matchfraction*peak
    return fsolve(lambda t: dualexpfunction(t,0.0,tau1,tau2)-matchvalue, init_t)[0]

def get_2exp_from_tpeak_duration(tpeak,duration):
    """ Doesn't work, fmin_tnc doesn't converge. Need to do bounded fmin, since tau1>0 and tau2>0.
    Supposed to return tau1 and tau2 for a 2exp, having tpeak and duration (FWHM). """
    eqns_func = lambda (tau1,tau2): \
        (invert_2exp(tau1,tau2,0.5,False) - invert_2exp(tau1,tau2,0.5,True) - duration)**2 + \
        (dualexp_tpeak(tau1,tau2) - tpeak)**2
    return fmin_tnc(eqns_func, (tpeak,duration), bounds=((0,+inf),(0,+inf)), approx_grad=True)

### Test DoG, DoAG, 2exp
#tpeak,duration = 400e-3,400e-3
#tau1,tau2 = get_2exp_from_tpeak_duration(tpeak,duration)[0]
#print 'tpeak=',dualexp_tpeak(tau1,tau2)
#plot(subtractedGaussians(arange(0,1.0,0.001),1.0,0.0,tpeak,10,duration/2.0,200e-3),'r-')
#plot(subtractedAlphaGaussians(arange(0,1.0,0.001),1.0,0.0,tpeak,10,duration/2.0,200e-3),'b-')
### 2exp doesn't work.
#plot(array([dualexpfunction(t,0,tau1,tau2) for t in arange(0,1.0,0.001)])*dualexp_norm(tau1,tau2),'g-')
#show()

#################################### All about double sigmoid functions #####################################

def dblsigmoid(t, center1, spread1, center2, spread2):
    return 1.0/(1+exp(-(t-center1)/spread1)) - 1.0/(1+exp(-(t-center2)/spread2))

def dblsigmoid_tpeak(center1, spread1, center2, spread2):
    derivative_dblsigmoid = ( lambda t: -exp(-(t-center1)/spread1)/spread1 / (1+exp(-(t-center1)/spread1))**2 + \
        exp(-(t-center2)/spread2)/spread2 / (1+exp(-(t-center2)/spread2))**2 )
    ## fsolve returns ndarray even for one value, return only the value
    return fsolve(derivative_dblsigmoid, (spread1+spread2)/2.0)[0]

def invert_dblsigmoid(center1, spread1, center2, spread2, matchfraction,risephase=True):
    """ return the time t at which dblsigmoid equals matchfraction*peak.
    It returns the initial t during rising phase if risephase is True, else later t during falling phase.
    It is upto you to provide matchfraction that lies between 0.0 and 1.0 ."""
    tp = dblsigmoid_tpeak(center1, spread1, center2, spread2)
    peak = dblsigmoid(tp, center1, spread1, center2, spread2)
    if risephase: init_t = center1
    else: init_t = center2
    matchvalue = matchfraction*peak
    ## fsolve returns ndarray even for one value, return only the value
    return fsolve(lambda t: dblsigmoid(t, center1, spread1, center2, spread2)-matchvalue, init_t)[0]

################################################################################################################

