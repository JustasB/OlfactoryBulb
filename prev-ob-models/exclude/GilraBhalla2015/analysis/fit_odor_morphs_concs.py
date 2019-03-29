# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO BE A CLONE OF MUKUND'S AND ADIL'S MATLAB ONE
## USAGE: python2.6 fit_odor_morphs.py <morphresultsfile.pickle> <concAresultsfile.pickle> <concBresultsfile.pickle>

from scipy import optimize
from scipy.special import * # has error function erf() and inverse erfinv()
from pylab import *
import pickle
import sys
import math

## use error function(x) for x>=0 (zero for x<0),
## OR use sigmoid(x) (non-zero for -ve x)
USE_ERF = False#True

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, inputList, GLOMS_ODOR, GLOMS_NIL
## we override inputList since we need three - for morph, concA, and concB
morphInputList = [ (0.0,1.0), (0.2,0.8), (0.4,0.6), (0.6,0.4), (0.8,0.2), (1.0,0.0), (0.0,0.0) ]
odorAInputList = [ (0.0,0.0), (0.2,0.0), (0.4,0.0), (0.6,0.0), (0.8,0.0), (1.0,0.0), (0.0,0.0) ]
odorBInputList = [ (0.0,1.0), (0.0,0.8), (0.0,0.6), (0.0,0.4), (0.0,0.2), (0.0,0.0), (0.0,0.0) ]

from simset_odor import * # has NUMBINS
from networkConstants import * # has central_glom
from sim_utils import * # has rebin() to alter binsize

errcut = 1e-20
iterationnum = 1

NUMBINS = 17 # I override the NUMBINS in simset_odor above, and I rebin() below
## overlapping bins - smooths results.
## for non-overlapping set to RESPIRATION/NUMBINS
bin_width_time = 2.0*RESPIRATION/NUMBINS

NUMMIX = len(inputList)
NUMWTS = NUMMIX-3 # remove the two pure odors and one pure air weights
## two pure + one air response, weights of A and B for mixtures, max of output sigmoid
NUMPARAMS = 3*NUMBINS+2*NUMWTS+1
firstrun = True

### numbers of mitrals to be fitted.
fitted_mitral_list = [2*central_glom+0, 2*central_glom+1]

## Fit type: 'lin' : linear or 'arb' : monotonic arbitrary
## if arbitrary fit_type, weights are also free params,
## if linear fit_type, weights are not free params.
fit_type = 'arb'

log81 = math.log(81)
 
def constrain0to1(x):
    try:
        return exp(x)/(1+exp(x)) # use numpy's exp
    except OverflowError as overflowerr:
        print overflowerr
        print x
        return 1.0

# define sigmoid which runs from (-0.5,0.1) to (+0.5,0.9)    
# Ideally the fitted sigmoid should be shifted by 0.5 i.e.
# exp((x-0.5)*log81)/(1+exp((x-0.5)*log81))
# This will overlap most of the linear part.
# But for fitting it doesn't matter,
# the fit routine will shift the parameters as required.
# But while plotting the internal response parameters,
# shift by 0.5 and plot -- see below
##### The slope of the sigmoid should be a parameter?
##### No, the individual activation functions R_odor(t) etc. will adjust.
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

def chisqfunc(params, morphdata, concAdata, concBdata):
    """
    Each of morphdata, concAdata and concBdata = (firingmeandata,firingerrdata)
    """
    RA = params[0:NUMBINS]
    RB = params[NUMBINS:2*NUMBINS]
    Rair = params[2*NUMBINS:3*NUMBINS]
    if fit_type == 'arb':
        #### for the weights also, we use exactly
        #### what is done by Mukund and Adil in matlab:
        #### constrain weights to be between 0 and 1
        #### sort the weights to ensure monotonicity
        inputsA = [ constrain0to1(x) for x in params[3*NUMBINS:(3*NUMBINS+NUMWTS)] ]
        ## important to put these in else along with sort(),
        ## weights saturate at 0.9 or so rather than at 1.0
        inputsA.extend([0.0,1.0]) # for pure odors
        inputsA.sort() # in place sort
        inputsA.append(0.0) # for air - keep this after sort!
        inputsB = [ constrain0to1(x) for x in params[(3*NUMBINS+NUMWTS):(3*NUMBINS+2*NUMWTS)] ]
        ## important to put these in else along with sort(),
        ## weights saturate at 0.9 or so rather than at 1.0
        inputsB.extend([0.0,1.0]) # for pure odors
        inputsB.sort(reverse=True) # weights of odor B need to be used in reverse
        inputsB.append(0.0) # for air - keep this after sort!
        #### Mukund and Adil constrained sigmoidmax > max firing rate response (note exp(x)>0.)
        sigmoidmax = array((morphdata[0],concAdata[0],concBdata[0])).max() + exp(params[3*NUMBINS+2*NUMWTS])
    else:
        ## *-operator unpacks the list which become args of zip()
        ## zip collects the i-th elements together of all the args.
        inputsA,inputsB = zip(*inputList) # keep the last (0,0) air input
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = array((morphdata[0],concAdata[0],concBdata[0])).max() + exp(params[3*NUMBINS])
        

    global iterationnum
    if iterationnum%100==0: print 'iteration number =',iterationnum
    #if iterationnum%100==0: print inputsA, inputsB
    chisqarray = [0.0]
    for i in range(len((morphInputList))):
        CA = inputsA[i]
        CB = inputsB[i]
        for bin in range(NUMBINS):
            Rmix = sigmoidmax*outputsigmoid( Rair[bin] + CA*RA[bin] + CB*RB[bin] )
            ## divide by error to do chi-square fit
            chisqarray.append( (morphdata[0][i][bin] - Rmix)/morphdata[1][i][bin] )

            Rmix = sigmoidmax*outputsigmoid( Rair[bin] + CA*RA[bin] )
            ## divide by error to do chi-square fit
            chisqarray.append( (concAdata[0][i][bin] - Rmix)/concAdata[1][i][bin] )

            Rmix = sigmoidmax*outputsigmoid( Rair[bin] + CB*RB[bin] )
            ## divide by error to do chi-square fit
            chisqarray.append( (concBdata[0][i][bin] - Rmix)/concBdata[1][i][bin] )
            
    ## normalized chi-sq to the number of dof
    chisqarray = array(chisqarray) / (3*morphdata[0].size - NUMPARAMS)
    iterationnum += 1
    return chisqarray

def read_morphfile(filename,fitted_mitral):
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
    ## since I fit the mean response,
    ## I must use standard error/deviation of the _mean_
    ## = standard deviation of a repeat / sqrt(num of repeats).
    NUMAVGs = len(mitral_responses_binned_list)
    mitral_responses_se= mitral_responses_std/sqrt(NUMAVGs)
    ## take the odor responses of the mitral to be fitted
    firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]
    firingbinserrList = mitral_responses_se[:,fitted_mitral]
    #### put in a minimum error, else divide by zero problems.
    ## find the minimum error > errcut
    errmin = firingbinserrList[where(firingbinserrList>errcut)].min()
    ## numpy where(), replace by errmin,
    ## all those elements in firingbinsList which are less than errmin 
    firingbinserrList = where(firingbinserrList>errcut, firingbinserrList, errmin)
    
    return (firingbinsmeanList,firingbinserrList)

def fit_morphs(filenamemorph, filenameconcA, filenameconcB, fitted_mitral):
    morphdata = read_morphfile(filenamemorph,fitted_mitral)
    concAdata = read_morphfile(filenameconcA,fitted_mitral)
    concBdata = read_morphfile(filenameconcB,fitted_mitral)
    ## compute initial params & return fits using only the morphdata
    ## however fitting in the chisq function uses morph and conc series data
    (firingbinsmeanList,firingbinserrList) = morphdata

    ########################## Initial values for the parameters
    if fit_type=='lin':
        params_filename = filenamemorph+'_morphconc_paramslin'+str(fitted_mitral)
    else:
        params_filename = filenamemorph+'_morphconc_params'+str(fitted_mitral)
    
    if firstrun:
        params0 = []
        spikesmax = firingbinsmeanList.max()
        RA = firingbinsmeanList[-2] # odor A is last but one
        RB = firingbinsmeanList[0] # odor B is first
        Rair = firingbinsmeanList[-1] # air response is last
        # The initial parameters are for odor A followed by odor B
        # the small value 0.001 should be put, else divide by zero errors in chi-sq!
        # extend(): Don't add the list as an element but add the elements of the list
        params0.extend([ ( inversesigmoid(0.998*RA[i]/spikesmax+0.001) - \
            inversesigmoid(0.998*Rair[i]/spikesmax+0.001) )/log81 for i in range(NUMBINS) ])
        params0.extend([ ( inversesigmoid(0.998*RB[i]/spikesmax+0.001) - \
            inversesigmoid(0.998*Rair[i]/spikesmax+0.001) )/log81 for i in range(NUMBINS) ])
        # initial params for the air vector # air is last
        params0.extend([ inversesigmoid(0.998*Rair[i]/spikesmax+0.001)/log81 for i in range(NUMBINS) ])
        if fit_type == 'arb':
            params0.extend([0.0]*2*NUMWTS) # weights of mixtures
        # argument for the exp in sigmoidmax as per Mukund and Adil.
        # -1 gives match for generated data, -3 went into local minimum.
        params0.append(-1)
        ##### pure odor concentrations are not parameters.
        ##### They are set to (CA=1,CB=0) and (CA=0,CB=1) and act as normalization.
        ## if arbitrary fit_type, weights are also free params,
        ## if linear fit_type, weights are not free params.
        if fit_type == 'arb':
            ## take only the mixture values, not the start and end-1 points which are pure odors,
            ## nor end point which is pure air
            for i,(inputA,inputB) in enumerate(inputList[1:-2]):
                # to constrain weights between 0 and 1, sigmoid is used,
                # so use inversesigmoid to set initial value for weights
                params0[3*NUMBINS+i] = inversesigmoid(inputA)
                params0[3*NUMBINS+NUMWTS+i] = inversesigmoid(inputB)
    else:
        f = open(params_filename,'r')
        params0,chisq = pickle.load(f)
        f.close()

    ###################################### Fitting
    #params = array(params0)
    ## uncomment above line and comment next fitting line if you just want to plot parameters.
    ## args is a tuple! if only one element write (elem, )
    params = optimize.leastsq( chisqfunc, params0,
        args=(morphdata, concAdata, concBdata), full_output=1 )
    print params
    params = params[0] # leastsq returns a whole tuple of stuff - errmsg etc.

    ## Calculate sum of squares of the chisqarray
    chisqarraysq = [i**2 for i in chisqfunc(params, morphdata, concAdata, concBdata)]
    chisq = reduce(lambda x, y: x+y, chisqarraysq)

    paramsfile = open(params_filename,'w')
    pickle.dump((params,chisq), paramsfile)
    paramsfile.close()

    ############################## Calculate fitted responses and return them
    
    if fit_type == 'arb':
        #### for the weights also, we use exactly what is done by Mukund and Adil in matlab:
        #### constrain weights to be between 0 and 1
        #### sort the weights to ensure monotonicity
        inputsA = [ constrain0to1(x) for x in params[3*NUMBINS:(3*NUMBINS+NUMWTS)] ]
        inputsA.extend([0.0,1.0])
        inputsA.sort() # in place sort
        inputsB = [ constrain0to1(x) for x in params[(3*NUMBINS+NUMWTS):(3*NUMBINS+2*NUMWTS)] ]
        inputsB.extend([0.0,1.0])
        inputsB.sort(reverse=True) # weights of odor B need to be used in reverse
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = array((morphdata[0],concAdata[0],concBdata[0])).max()\
            + math.exp(params[3*NUMBINS+2*NUMWTS])
    else:
        ## *-operator unpacks the list which become args of zip()
        ## zip collects the i-th elements together of all the args.
        inputsA,inputsB = zip(*(inputList[:-1])) # leave out the last (0,0) air input
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = array((morphdata[0],concAdata[0],concBdata[0])).max()\
            + math.exp(params[3*NUMBINS])

    fitted_responses = []
    for inpnum,(inputA,inputB) in enumerate(morphInputList[:-1]):
        fitted_responses.append(\
            [ sigmoidmax*outputsigmoid( \
            inputsA[inpnum]*params[i] + inputsB[inpnum]*params[NUMBINS+i] + params[2*NUMBINS+i]\
            ) for i in range(NUMBINS) ] )
    fitted_responses.append( [ sigmoidmax*outputsigmoid( params[2*NUMBINS+i] )\
        for i in range(NUMBINS) ] )

    return (params,chisq,inputsA,inputsB,fitted_responses,firingbinsmeanList,firingbinserrList)

if __name__ == "__main__":
    if len(sys.argv) > 3:
        filenamemorph = sys.argv[1]
        filenameconcA = sys.argv[2]
        filenameconcB = sys.argv[3]
    else:
        print "Specify data files containing pickled lists of morphs & conc series."
        sys.exit(1)

    for fitted_mitral in fitted_mitral_list:
        params,chisq,inputsA,inputsB,fitted_responses,firingbinsmeanList,firingbinserrList\
            = fit_morphs(filenamemorph, filenameconcA, filenameconcB, fitted_mitral)
        print "Mit",fitted_mitral,"chisq =",chisq

        ################# Plot simulated responses
        fig = figure(facecolor='w') # 'none' is transparent
        ax = fig.add_subplot(111)
        pure_i = 0
        brightness = 0.1
        for i,(inputA,inputB) in enumerate(morphInputList):
            ## The inputA acts to morph odor response from red to green color
            ## air response in black
            ## if not a pure odor/air, bring down its brightness.
            if i not in [0,len(inputList)-2,len(inputList)-1]:
                ax.errorbar(x=range(NUMBINS),y=firingbinsmeanList[i],yerr=firingbinserrList[i],\
                    color=( (1-brightness)+brightness*inputA, (1-brightness)+brightness*inputB, (1-brightness) ),\
                    marker='+',linestyle='solid', linewidth=2)
            else:
                ax.errorbar(x=range(NUMBINS),y=firingbinsmeanList[i],yerr=firingbinserrList[i],\
                    color=(inputA,inputB,0),marker='+',linestyle='solid',linewidth=2,\
                    label=['odor A','odor B','air'][pure_i])
                pure_i += 1

        ##################### Plot fitted responses.
        # RA + Rair
        ax.plot(fitted_responses[-2],\
            color=(1,0,1),marker='o',linestyle='dashed', linewidth=2.0, label='fit A')
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
        
        ################### Plot the activation (argument of sigmoid/erf)
        fig2 = figure(facecolor='w') # none is transparent
        ax2 = fig2.add_subplot(111)
        RA = params[0:NUMBINS]
        RB = params[NUMBINS:2*NUMBINS]
        Rair = params[2*NUMBINS:3*NUMBINS]
        ## Adil did not shift 'sigmoid argument' i.e. activation by 0.5,
        ## the shift is needed to get sigmoid(0)=0.1,
        ## doesn't matter for fitting though, Rodor(t), etc. accomodate this. 
        if not USE_ERF:
            RA = RA + 0.5
            RB = RB + 0.5
            Rair = Rair + 0.5
        plot(RA,color=(1,0,0),marker='+',linestyle='solid',linewidth=2,label='activation_A')
        plot(RB,color=(0,1,0),marker='+',linestyle='solid',linewidth=2,label='activation_B')
        plot(Rair,color=(0,0,0),marker='+',linestyle='solid',linewidth=2,label='activation_air')
        title( 'Activation of mitral %d '%(fitted_mitral,), fontsize=24 )
        axes_labels(ax2,'respiratory phase bin','firing rate (Hz)',adjustpos=True)

        ################### Linearity of weights plot
        print 'weightsA =',inputsA
        print 'weightsB =',inputsB

        fig2 = figure(facecolor='w') # none is transparent
        ax2 = fig2.add_subplot(111)
        actualweights = [ wts[0] for wts in morphInputList[:-1]]
        ax2.plot(actualweights,inputsA,color=(1,0,0),marker='+',\
            linestyle='solid',linewidth=2,label='weight odorA')
        ax2.plot(actualweights,inputsB,color=(0,1,0),marker='+',\
            linestyle='solid',linewidth=2,label='weight odorB')
        ax2.plot(actualweights,arange(0.0,1.01,0.2),color=(1,0,1),\
            marker='+',linestyle='dashed',linewidth=2,label='linear A')
        ax2.plot(actualweights,arange(1.0,-0.01,-0.2),color=(0,1,1),\
            marker='+',linestyle='dashed',linewidth=2,label='linear B')
        #title( 'chisquare normalized = '+str(chisq) )
        maxerror = sqrt(sum(array([0.8,0.6,0.4,0.2])**2)/4.0) # max rms error
        ## normalized score = 1 - norm-ed rms error
        scoreA = 1 - sqrt( sum( (inputsA[1:-1]-arange(0.2,0.81,0.2))**2 )/4.0 )/maxerror
        scoreB = 1 - sqrt( sum( (inputsB[1:-1]-arange(0.8,0.19,-0.2))**2 )/4.0 )/maxerror
        title( 'Linearity mitral %d: \nscoreA=%.2f, scoreB=%.2f'\
            %(fitted_mitral,scoreA,scoreB), fontsize=24 )
        axes_labels(ax2,'weight','fitted weight',adjustpos=True)
        legend(loc='center right')
    
    show()
