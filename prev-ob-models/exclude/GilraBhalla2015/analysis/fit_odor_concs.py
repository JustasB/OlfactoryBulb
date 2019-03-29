# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO BE A CLONE OF MUKUND'S AND ADIL'S MATLAB ONE
## USAGE: python2.6 fit_odor_concs.py ../results/odor_morphs/2011-01-13_odormorph_SINGLES_JOINTS_PGS.pickle
## This is only for concentration series, only one of odorAon or odorBon should be True

from scipy import optimize
from scipy.special import * # has error function erf() and inverse erfinv()
from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from simset_odor import * # has NUMBINS
from networkConstants import * # has central_glom
from sim_utils import * # has rebin() to alter binsize

## use error function(x) for x>=0 (zero for x<0),
## OR use sigmoid(x) (non-zero for -ve x)
USE_ERF = False#True

errcut = 1e-20
iterationnum = 1

NUMBINS = 17 # I override the NUMBINS in simset_odor above, and I rebin() below
bin_width_time = 2.0*RESPIRATION/NUMBINS # smooths / overlapping bins. non-overlapping would be RESPIRATION/NUMBINS

## index of the odor that is on, in inputList
if odorAon:
    odoridx = 0
    sortreverse = False
else:
    odoridx = 1
    sortreverse = True

NUMMIX = len(inputList)
NUMWTS = NUMMIX-3 # remove the two pure odors and one pure air weights
NUMPARAMS = 2*NUMBINS+1*NUMWTS+1 # one pure + one air response, weights of A or B for conc series, max of output sigmoid
firstrun = True

# numbers of mitral to be fitted.
fitted_mitral_list = [2*central_glom+0, 2*central_glom+1]

log81 = math.log(81)

def constrain0to1(x):
    return math.exp(x)/(1+math.exp(x))

# define sigmoid which runs from (-0.5,0.1) to (+0.5,0.9)    
# Ideally the fitted sigmoid should be shifted by 0.5 i.e.
# exp((x-0.5)*log81)/(1+exp((x-0.5)*log81))
# This will overlap most of the linear part.
# But for fitting it doesn't matter, the fit routine will shift the parameters as required.
# But while plotting the internal response parameters, shift by 0.5 and plot -- see below
def outputsigmoid(x):
    if USE_ERF:
        if x<0: return 0
        else: return erf(x)
    else:
        try:
            return exp(x*log81)/(1+exp(x*log81)) # use numpy's exp
        except OverflowError as overflowerr:
            print overflowerr
            print x
            return 1.0
    
def inversesigmoid(x):
    if USE_ERF:
        if x<0: return x
        else: return erfinv(x)
    else:
        ## just to set initial values, value doesn't matter too much when x tends to 0
        if x>1e-200: return math.log(x/(1-x))
        else: return -5e2

def chisqfunc(params, ydata, errdata):
    R = params[0:NUMBINS]
    Rair = params[NUMBINS:2*NUMBINS]
    #### for the weights also, we use exactly what is done by Mukund and Adil in matlab:
    #### constrain weights to be between 0 and 1
    #### sort the weights to ensure monotonicity
    inputs = [ constrain0to1(x) for x in params[2*NUMBINS:(2*NUMBINS+NUMWTS)] ]
    ## important to put these in else along with sort(),
    ## weights saturate at 0.9 or so rather than at 1.0
    inputs.extend([0.0,1.0]) # for pure odors
    inputs.sort(reverse=sortreverse) # in place sort
    inputs.append(0.0) # for air - keep this after sort!
    #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
    sigmoidmax = ydata.max() + exp(params[2*NUMBINS+NUMWTS])

    global iterationnum
    chisqarray = [0.0]
    for i,input in enumerate(inputList):
        C = inputs[i]
        for bin in range(NUMBINS):
            Rmix = sigmoidmax*outputsigmoid( Rair[bin] + C*R[bin] )
            chisqarray.append( (ydata[i][bin] - Rmix)/errdata[i][bin] ) # divide by error to do chi-square fit
            
    chisqarray = array(chisqarray) / (ydata.size - NUMPARAMS) # normalized chi-sq to the number of dof
    if iterationnum%100==0: print 'iteration number =',iterationnum
    iterationnum += 1
    return chisqarray
    
def fit_morphs(filename, fitted_mitral):
    f = open(filename,'r')
    #### mitral_responses_list[avgnum][odornum][mitralnum][binnum]
    #### mitral_responses_binned_list[avgnum][odornum][mitralnum][binnum]
    mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
    f.close()

    ###################### Input conditioning
    mitral_responses_binned_list = \
        rebin(mitral_responses_list, NUMBINS, bin_width_time)
    #### very important to convert to numpy array,
    #### else where() below returns empty list.
    mitral_responses_binned_list = array(mitral_responses_binned_list)
    mitral_responses_mean = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    # since I fit the mean response, I must use standard error/deviation of the _mean_
    # = standard deviation of a repeat / sqrt(num of repeats).
    NUMAVGs = len(mitral_responses_binned_list)
    mitral_responses_se= mitral_responses_std/sqrt(NUMAVGs)
    # take the odor responses of the mitral to be fitted
    firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]
    firingbinserrList = mitral_responses_se[:,fitted_mitral]
    # put in a minimum error, else divide by zero problems.
    errmin = firingbinserrList[where(firingbinserrList>errcut)].min() # find the minimum error > errcut
    # numpy where(), replace by errmin, all those elements in firingbinsList which are less than errmin 
    firingbinserrList = where(firingbinserrList>errcut, firingbinserrList, errmin)

    ########################## Initial values for the parameters
    if firstrun:
        params0 = []
        spikesmax = firingbinsmeanList.max()
        if odorAon: R = firingbinsmeanList[-2] # odor A is last but one
        else: R = firingbinsmeanList[0] # odor B is first
        Rair = firingbinsmeanList[-1] # air response is last
        # The initial parameters are for odor 
        # the small value 0.001 should be put, else divide by zero errors in chi-sq!
        # extend(): Don't add the list as an element but add the elements of the list
        params0.extend([ ( inversesigmoid(0.998*R[i]/spikesmax+0.001) - \
            inversesigmoid(0.998*Rair[i]/spikesmax+0.001) )/log81 for i in range(NUMBINS) ])
        params0.extend([ inversesigmoid(0.998*Rair[i]/spikesmax+0.001)/log81 for i in range(NUMBINS) ]) # air is last
        # initial params for the air vector
        params0.extend([0.0]*NUMWTS) # weights of mixtures
        # argument for the exp in sigmoidmax as per Mukund and Adil.
        # -1 gives match for generated data, -3 went into local minimum.
        params0.append(-1)
        ##### pure odor concentrations are not parameters.
        ##### They are set to (CA=1,CB=0) and (CA=0,CB=1) and act as normalization.
        # take only the mixture values, not the start and end-1 points which are pure odors,
        # nor end point which is pure air
        for i,input in enumerate(inputList[1:-2]):
            # to constrain weights between 0 and 1, sigmoid is used,
            # so use inversesigmoid to set initial value
            input = input[odoridx]
            params0[2*NUMBINS+i] = inversesigmoid(input)
    else:
        f = open(sys.argv[1]+'_params','r')
        params0 = pickle.load(f)
        f.close()

    ###################################### Fitting
    #params = array(params0) ## uncomment this line and comment next fitting line if you just want to plot parameters.
    # args is a tuple! if only one element write (elem, )
    params = optimize.leastsq( chisqfunc, params0,
        args=(firingbinsmeanList, firingbinserrList), full_output=1, maxfev=50000 )
    print params
    params = params[0] # leastsq returns a whole tuple of stuff - errmsg etc.

    paramsfile = open(sys.argv[1]+'_params','w')
    pickle.dump(params, paramsfile)
    paramsfile.close()

    ############################## Calculate fitted responses and return them
    # Calculate sum of squares of the chisqarray
    chisqarraysq = [i**2 for i in chisqfunc(params, firingbinsmeanList, firingbinserrList)]
    chisq = reduce(lambda x, y: x+y, chisqarraysq)

    #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
    sigmoidmax = firingbinsmeanList.max() + math.exp(params[2*NUMBINS+NUMWTS])
    
    #### for the weights also, we use exactly what is done by Mukund and Adil in matlab:
    #### constrain weights to be between 0 and 1
    #### sort the weights to ensure monotonicity
    inputs = [ constrain0to1(x) for x in params[2*NUMBINS:(2*NUMBINS+NUMWTS)] ]
    inputs.extend([0.0,1.0])
    inputs.sort(reverse=sortreverse) # in place sort

    fitted_responses = []
    for inpnum,input in enumerate(inputList[:-1]):
        fitted_responses.append(\
            [ sigmoidmax*outputsigmoid( \
            inputs[inpnum]*params[i] + params[NUMBINS+i]\
            ) for i in range(NUMBINS) ] )
    fitted_responses.append( [ sigmoidmax*outputsigmoid( params[NUMBINS+i] )\
            for i in range(NUMBINS) ] )

    return (params,chisq,inputs,fitted_responses,firingbinsmeanList,firingbinserrList)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        print "Specify data file containing pickled list."
        sys.exit(1)

    for fitted_mitral in fitted_mitral_list:
        params,chisq,inputs,fitted_responses,firingbinsmeanList,firingbinserrList\
            = fit_morphs(filename, fitted_mitral)
        print "Mit",fitted_mitral,"chisq =",chisq

        ################# Plot simulated responses
        fig = figure(facecolor='w') # 'none' is transparent
        ax = fig.add_subplot(111)
        pure_i = 0
        brightness = 0.1
        for i,(inputA,inputB) in enumerate(inputList):
            ## The inputA/inputB acts to morph odor response from red/green to black (air)
            ## if not a pure odor/air, bring down its brightness.
            if odorAon: pureslist = [len(inputList)-2,len(inputList)-1]
            else: pureslist = [0,len(inputList)-1]
            if i not in pureslist:
                ax.errorbar(x=range(NUMBINS),y=firingbinsmeanList[i],yerr=firingbinserrList[i],\
                    color=( (1-brightness)+brightness*inputA, (1-brightness)+brightness*inputB, (1-brightness) ),\
                    marker='+',linestyle='solid')
            else:
                ax.errorbar(x=range(NUMBINS),y=firingbinsmeanList[i],yerr=firingbinserrList[i],\
                    color=(inputA,inputB,0),marker='+',linestyle='solid',\
                    label=['odor A','odor B','air'][pure_i])
                pure_i += 1

        ##################### Plot fitted responses.
        if odorAon:
            # RA + Rair
            ax.plot(fitted_responses[-2],\
                color=(1,0,1),marker='o',linestyle='dashed', linewidth=2.0, label='fit A')
        else:
            # RB + Rair
            ax.plot(fitted_responses[0],\
                color=(0,1,1),marker='o',linestyle='dashed', linewidth=2.0, label='fit B')
        # Rair
        ax.plot(fitted_responses[-1],\
            color=(0.5,0.5,0.5),marker='o',linestyle='dashed', linewidth=2.0, label='fit air')
        title('Mitral %d responses & linear-sigmoid fit'%fitted_mitral,fontsize=24 )
        axes_labels(ax,'respiratory phase bin','firing rate (Hz)',adjustpos=True)
        #ylim(ymin=-6, ymax=4)
        legend()

        ################### Linearity of weights plot
        print 'weights =',inputs

        fig2 = figure(facecolor='w') # none is transparent
        ax2 = fig2.add_subplot(111)
        maxerror = sqrt(sum(array([0.8,0.6,0.4,0.2])**2)/4.0) # max rms error
        actualweights = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        if odorAon:
            ax2.plot(actualweights,inputs,color=(1,0,0),marker='+',linestyle='solid',label='weight odor A')
            ax2.plot(actualweights,arange(0.0,1.01,0.2),color=(1,0,1),marker='+',linestyle='dashed',label='linear A')
            ## normalized score = 1 - norm-ed rms error
            score = 1 - sqrt( sum( (inputs[1:-1]-arange(0.2,0.81,0.2))**2 )/4.0 )/maxerror
        else:
            ax2.plot(actualweights,inputs,color=(0,1,0),marker='+',linestyle='solid',label='weight odor B')
            ax2.plot(actualweights,arange(1.0,-0.01,-0.2),color=(0,1,1),marker='+',linestyle='dashed',label='linear B')
            ## normalized score = 1 - norm-ed rms error
            score = 1 - sqrt( sum( (inputs[1:-1]-arange(0.8,0.19,-0.2))**2 )/4.0 )/maxerror
        ax2.plot(actualweights,arange(0.0,1.01,0.2),color=(1,0,1),marker='+',linestyle='dashed',label='linear')
        #title( 'chisquare normalized = '+str(chisq) )
        title( 'Linearity mitral %d: \nscore=%.2f'%(fitted_mitral,score), fontsize=24 )
        axes_labels(ax2,'weight','fitted weight',adjustpos=True)
        legend(loc='center right')
    
    show()
