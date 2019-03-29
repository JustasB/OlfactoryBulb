# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO BE A CLONE OF MUKUND'S AND ADIL'S MATLAB ONE
## USAGE: python2.6 fit_odor_morphs.py ../results/odor_morphs/2011-01-13_odormorph_SINGLES_JOINTS_PGS.pickle [CHISQ_HIST] [SAVEFIG]

from scipy import optimize
from scipy.special import * # has error function erf() and inverse erfinv()
from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, inputList and pulseList, GLOMS_ODOR, GLOMS_NIL
from networkConstants import * # has central_glom
from sim_utils import * # has rebin() to alter binsize
from analysis_utils import * # has read_morphfile() and NUM_REBINS, etc.

## use error function(x) for x>=0 (zero for x<0),
## OR use sigmoid(x) (non-zero for -ve x)
USE_ERF = False#True

iterationnum = 1

## I don't use the NUMBINS in simset_odor.py, rather I rebin() with below NUM_REBINS
## Adil used 17 bins for a 1s rat respiration cycle.
## I'm using 9 bins to get the same binwidth, else there are oscillations ~ 35 Hz gamma?
NUM_REBINS = 9#17

NUMMIX = len(inputList)
## remove the two pure odors and one pure air weights
NUMWTS = NUMMIX-3
## two pure + one air response, weights of A and B for mixtures, max of output sigmoid
NUMPARAMS = 3*NUM_REBINS+2*NUMWTS+1
firstrun = False#True

### numbers of mitral to be fitted.
fitted_mitral_list = [2*central_glom+0, 2*central_glom+1]

## Fit type: 'lin' : linear or 'arb' : monotonic arbitrary
## if arbitrary fit_type, weights are also free params,
## if linear fit_type, weights are not free params.
## This param is passed to fit_morphs()
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

def chisqfunc(params, ydata, errdata, fit_type):
    RA = params[0:NUM_REBINS]
    RB = params[NUM_REBINS:2*NUM_REBINS]
    Rair = params[2*NUM_REBINS:3*NUM_REBINS]
    if fit_type == 'arb':
        #### for the weights also, we use exactly what is done by Mukund and Adil in matlab:
        #### constrain weights to be between 0 and 1
        #### sort the weights to ensure monotonicity
        inputsA = [ constrain0to1(x) for x in params[3*NUM_REBINS:(3*NUM_REBINS+NUMWTS)] ]
        ## important to put these in else along with sort(),
        ## weights saturate at 0.9 or so rather than at 1.0
        inputsA.extend([0.0,1.0]) # for pure odors
        inputsA.sort() # in place sort
        inputsA.append(0.0) # for air - keep this after sort!
        inputsB = [ constrain0to1(x) for x in params[(3*NUM_REBINS+NUMWTS):(3*NUM_REBINS+2*NUMWTS)] ]
        ## important to put these in else along with sort(),
        ## weights saturate at 0.9 or so rather than at 1.0
        inputsB.extend([0.0,1.0]) # for pure odors
        inputsB.sort(reverse=True) # weights of odor B need to be used in reverse
        inputsB.append(0.0) # for air - keep this after sort!
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = ydata.max() + exp(params[3*NUM_REBINS+2*NUMWTS])
    else:
        ## *-operator unpacks the list which become args of zip()
        ## zip collects the i-th elements together of all the args.
        inputsA,inputsB = zip(*inputList) # keep the last (0,0) air input
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = ydata.max() + exp(params[3*NUM_REBINS])
        

    global iterationnum
    if iterationnum%1000==0: print 'iteration number =',iterationnum
    #if iterationnum%100==0: print inputsA, inputsB
    chisqarray = [0.0]
    for i,(inputA,inputB) in enumerate(inputList):
        CA = inputsA[i]
        CB = inputsB[i]
        for bin in range(NUM_REBINS):
            Rmix = sigmoidmax*outputsigmoid( Rair[bin] + CA*RA[bin] + CB*RB[bin] )
            chisqarray.append( (ydata[i][bin] - Rmix)/errdata[i][bin] ) # divide by error to do chi-square fit
            
    ## not yet squared, so normalized 'chi' to sqrt of number of dof
    chisqarray = array(chisqarray) / sqrt(ydata.size-NUMPARAMS)
    iterationnum += 1
    return chisqarray

def fit_morphs(filename, fitted_mitral, fit_type='arb', refit=True):
    ## The model predicts the individual response not the mean.
    ## Hence below fitting uses standard deviation, not standard error of the mean.
    numavgs,firingbinsmeanList,firingbinserrList = read_morphfile(filename,fitted_mitral,NUM_REBINS)

    ########################## Initial values for the parameters
    if fit_type=='lin':
        params_filename = filename+'_paramslin'+str(fitted_mitral)
    else:
        params_filename = filename+'_params'+str(fitted_mitral)
    
    if firstrun or refit:
        params0 = []
        spikesmax = firingbinsmeanList.max()
        RA = firingbinsmeanList[-2] # odor A is last but one
        RB = firingbinsmeanList[0] # odor B is first
        Rair = firingbinsmeanList[-1] # air response is last
        # The initial parameters are for odor A followed by odor B
        # the small value 0.001 should be put, else divide by zero errors in chi-sq!
        # extend(): Don't add the list as an element but add the elements of the list
        params0.extend([ ( inversesigmoid(0.998*RA[i]/spikesmax+0.001) - \
            inversesigmoid(0.998*Rair[i]/spikesmax+0.001) )/log81 for i in range(NUM_REBINS) ])
        params0.extend([ ( inversesigmoid(0.998*RB[i]/spikesmax+0.001) - \
            inversesigmoid(0.998*Rair[i]/spikesmax+0.001) )/log81 for i in range(NUM_REBINS) ])
        # initial params for the air vector # air is last
        params0.extend([ inversesigmoid(0.998*Rair[i]/spikesmax+0.001)/log81 for i in range(NUM_REBINS) ])
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
                params0[3*NUM_REBINS+i] = inversesigmoid(inputA)
                params0[3*NUM_REBINS+NUMWTS+i] = inversesigmoid(inputB)
    else:
        f = open(params_filename,'r')
        params0,chisq = pickle.load(f)
        f.close()

    ###################################### Fitting
    if not refit:
        params = array(params0) ## only use params, do not fit again
    else:
        ## args is a tuple! if only one element write (elem, )
        params = optimize.leastsq( chisqfunc, params0,
            args=(firingbinsmeanList, firingbinserrList, fit_type), full_output=1, maxfev=50000)
        params = params[0] # leastsq returns a whole tuple of stuff - errmsg etc.

    ## Calculate sum of squares of the chisqarray
    chisqarraysq = [i**2 for i in chisqfunc(params, firingbinsmeanList, firingbinserrList, fit_type)]
    chisq = reduce(lambda x, y: x+y, chisqarraysq)

    if refit:
        paramsfile = open(params_filename,'w')
        pickle.dump((params,chisq), paramsfile)
        paramsfile.close()

    ############################## Calculate fitted responses and return them
    
    if fit_type == 'arb':
        #### for the weights also, we use exactly what is done by Mukund and Adil in matlab:
        #### constrain weights to be between 0 and 1
        #### sort the weights to ensure monotonicity
        inputsA = [ constrain0to1(x) for x in params[3*NUM_REBINS:(3*NUM_REBINS+NUMWTS)] ]
        inputsA.extend([0.0,1.0])
        inputsA.sort() # in place sort
        inputsB = [ constrain0to1(x) for x in params[(3*NUM_REBINS+NUMWTS):(3*NUM_REBINS+2*NUMWTS)] ]
        inputsB.extend([0.0,1.0])
        inputsB.sort(reverse=True) # weights of odor B need to be used in reverse
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = firingbinsmeanList.max() + math.exp(params[3*NUM_REBINS+2*NUMWTS])
    else:
        ## *-operator unpacks the list which become args of zip()
        ## zip collects the i-th elements together of all the args.
        inputsA,inputsB = zip(*(inputList[:-1])) # leave out the last (0,0) air input
        #### Mukund and Adil constrained sigmoidmax > ydatamax (note exp(x)>0.)
        sigmoidmax = firingbinsmeanList.max() + math.exp(params[3*NUM_REBINS])

    fitted_responses = []
    for inpnum,(inputA,inputB) in enumerate(inputList[:-1]):
        fitted_responses.append(\
            [ sigmoidmax*outputsigmoid( \
            inputsA[inpnum]*params[i] + inputsB[inpnum]*params[NUM_REBINS+i] + params[2*NUM_REBINS+i]\
            ) for i in range(NUM_REBINS) ] )
    fitted_responses.append( [ sigmoidmax*outputsigmoid( params[2*NUM_REBINS+i] )\
            for i in range(NUM_REBINS) ] )

    return (params,chisq,inputsA,inputsB,fitted_responses,numavgs,firingbinsmeanList,firingbinserrList)

def plot_example_onemit(ax1,ax2,fitted_mitral,mit_fit_params):
    bindt = RESPIRATION/float(NUM_REBINS)
    respiration2time = arange(RESPIRATION,2*RESPIRATION,bindt) + bindt/2.0

    params,chisq,inputsA,inputsB,fitted_responses,numavgs,firingbinsmeanList,firingbinserrList =\
        mit_fit_params
    print "Mit",fitted_mitral,"normalized chisq =",chisq
    brightness = 0.2
    num_morphs = len(inputList)-1
    for i,(inputA,inputB) in enumerate(inputList):
        ## The inputA acts to morph odor response from red to blue color
        ## air response in black
        ## if not a pure odor/air, bring down its brightness.
        if i==0: color,alpha = 'b',1.0
        elif i==num_morphs-1: color,alpha = 'r',1.0
        elif i==num_morphs: color,alpha = 'k',1.0
        else: color,alpha = (i/float(num_morphs),0,1.0-i/float(num_morphs)),brightness
        if i in [0,num_morphs-1,num_morphs]:
            simresponse = firingbinsmeanList[i]
            ## For the plots, show std error of the mean
            simerr = firingbinserrList[i]/sqrt(numavgs)
            ax1.fill_between(respiration2time,simresponse+simerr,simresponse-simerr,
                color=color,alpha=alpha*0.4,linewidth=0)
            ax1.plot(respiration2time,simresponse,\
                color=color,alpha=alpha,marker='.',markersize=marker_size,\
                linewidth=linewidth,clip_on=False)

    ##################### Plot fitted responses.
    ## RA + Rair
    line, = ax1.plot(respiration2time,fitted_responses[-2],\
        color='m',marker='+',markersize=marker_size,\
        linestyle='dashed', linewidth=linewidth, label='fit A',clip_on=False)
    line.set_dashes((3,1))
    ## RB + Rair
    line, = ax1.plot(respiration2time,fitted_responses[0],\
        color='c',marker='+',markersize=marker_size,\
        linestyle='dashed', linewidth=linewidth, label='fit B',clip_on=False)
    line.set_dashes((3,1))
    ## Rair
    line, = ax1.plot(respiration2time,fitted_responses[-1],\
        color=(0.5,0.5,0.5),marker='+',markersize=marker_size,\
        linestyle='dashed', linewidth=linewidth, label='fit air',clip_on=False)
    line.set_dashes((3,1))
    #title('Mitral %d responses & linear-sigmoid fit'%fitted_mitral,fontsize=24 )
    #axes_labels(ax,'respiratory phase bin','firing rate (Hz)',adjustpos=True)
    #ylim(ymin=-6, ymax=4)
    #legend()

    ################### Linearity of weights plot
    print 'weightsA =',inputsA
    print 'weightsB =',inputsB

    actualweights = [ wts[0] for wts in inputList[:-1]]
    ax2.plot(actualweights,arange(0.0,1.01,0.2),color='r',\
        marker='.',markersize=marker_size,clip_on=False,\
        linestyle='solid',linewidth=linewidth,label='linear A')
    ax2.plot(actualweights,arange(1.0,-0.01,-0.2),color='b',\
        marker='.',markersize=marker_size,clip_on=False,\
        linestyle='solid',linewidth=linewidth,label='linear B')
    line, = ax2.plot(actualweights,inputsA,color='m',clip_on=False,\
        marker='+',linestyle='dashed',markersize=marker_size,\
        linewidth=linewidth,label='weight odorA')
    line.set_dashes((3,1))
    line, = ax2.plot(actualweights,inputsB,color='c',clip_on=False,\
        marker='+',linestyle='dashed',markersize=marker_size,\
        linewidth=linewidth,label='weight odorB')
    line.set_dashes((3,1))
    #title( 'chisquare normalized = '+str(chisq) )
    maxerror = sqrt(sum(array([0.8,0.6,0.4,0.2])**2)/4.0) # max rms error
    ## normalized score = 1 - norm-ed rms error
    scoreA = 1 - sqrt( sum( (inputsA[1:-1]-arange(0.2,0.81,0.2))**2 )/4.0 )/maxerror
    scoreB = 1 - sqrt( sum( (inputsB[1:-1]-arange(0.8,0.19,-0.2))**2 )/4.0 )/maxerror
    #title( 'Linearity mitral %d: \nscoreA=%.2f, scoreB=%.2f'%(fitted_mitral,scoreA,scoreB), fontsize=24 )
    #axes_labels(ax2,'weight','fitted weight',adjustpos=True)
    #legend(loc='center right')

    ## beautify plots
    for ax in [ax1,ax2]:
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=False,y0min=True,xticksposn='bottom',yticksposn='left')

def plot_example_chisq():
    fig = figure(figsize=(columnwidth,columnwidth/2.0),dpi=300,facecolor='w') # 'none' is transparent
    if 'CHISQ_HIST' in sys.argv:
        axgrid = (2,3)
        ax5 = plt.subplot2grid(axgrid,(0,2),frameon=False)
        ax6 = plt.subplot2grid(axgrid,(1,2),frameon=False)
        import average_odor_morphs as chisq_hist
        ## chi-sq histograms for both non-lin and lin weights
        chisq_hist.plot_chisq_hist_paperfigure(ax5,ax6,'../results/odor_morphs'+dirextn)
    else: axgrid = (2,2)
    for fitted_mitral in fitted_mitral_list:
        mit_fit_params = fit_morphs(filename, fitted_mitral, fit_type=fit_type)

        ################# Plot simulated responses
        ax1 = plt.subplot2grid(axgrid,(fitted_mitral,0),frameon=False)
        #text(0.1,1.0,['A','C'][fitted_mitral],fontsize=label_fontsize,transform=ax1.transAxes)
        ax2 = plt.subplot2grid(axgrid,(fitted_mitral,1),frameon=False)
        #text(0.1,1.0,['B','D'][fitted_mitral],fontsize=label_fontsize,transform=ax2.transAxes)
        
        plot_example_onemit(ax1,ax2,fitted_mitral,mit_fit_params)

    fig.tight_layout()
    fig_clip_off(fig)
    subplots_adjust(left=0.1,top=0.92,bottom=0.15,wspace=0.5)
    fig.text(0.015,0.7,'firing rate (Hz)',fontsize=label_fontsize, rotation='vertical', transform=fig.transFigure)
    fig.text(0.15,0.025,'time (s)',fontsize=label_fontsize, transform=fig.transFigure)
    fig.text(0.35,0.7,'fitted weight',fontsize=label_fontsize, rotation='vertical', transform=fig.transFigure)
    fig.text(0.43,0.025,'ORN weight',fontsize=label_fontsize, transform=fig.transFigure)
    if 'SAVEFIG' in sys.argv:
        fig.savefig('../figures/sim_morphs.svg', bbox_inches='tight',dpi=fig.dpi)
        fig.savefig('../figures/sim_morphs.png', bbox_inches='tight',dpi=fig.dpi)

def plot_example(refit=True):
    fig = figure(figsize=(columnwidth*2./3.,columnwidth/2.0),dpi=300,facecolor='w') # 'none' is transparent
    axgrid = (2,2)
    for fitted_mitral in fitted_mitral_list:
        mit_fit_params = fit_morphs(filename, fitted_mitral, fit_type=fit_type, refit=refit)

        ################# Plot simulated responses
        ax1 = plt.subplot2grid(axgrid,(fitted_mitral,0),frameon=False)
        #text(0.1,1.0,['A','C'][fitted_mitral],fontsize=label_fontsize,transform=ax1.transAxes)
        ax2 = plt.subplot2grid(axgrid,(fitted_mitral,1),frameon=False)
        #text(0.1,1.0,['B','D'][fitted_mitral],fontsize=label_fontsize,transform=ax2.transAxes)
        
        plot_example_onemit(ax1,ax2,fitted_mitral,mit_fit_params)
        axes_labels(ax1,['','time (s)'][fitted_mitral],\
            ['firing rate (Hz)',''][fitted_mitral],xpad=0,ypad=0)
        axes_labels(ax2,['','ORN weight'][fitted_mitral],\
            ['fitted weight',''][fitted_mitral],xpad=0,ypad=0)
        ax1.yaxis.set_label_coords(-0.25,-0.3)
        ax2.yaxis.set_label_coords(-0.16,-0.3)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.2,wspace=0.4)
    fig_clip_off(fig)
    if 'SAVEFIG' in sys.argv:
        fig.savefig('../figures/sim_morphs.svg', bbox_inches='tight',dpi=fig.dpi)
        fig.savefig('../figures/sim_morphs.png', bbox_inches='tight',dpi=fig.dpi)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        post_pulses = filename.split('odor_morphs')[1]
        dirextn = post_pulses.split('/')[0]
        print 'directory extension is',dirextn
    else:
        print "Specify data file containing pickled list."
        sys.exit(1)

    ### old paper figure that shows example and chisq distribs
    #plot_example_chisq()
    ### OBSOLETE: PAPER FIGURE supplementary fig 1 that shows only example
    ### use average_odor_morphs.py to get the R2N distribution plot
    ##plot_example()
    ## PAPER FIGURE supplementary Fig 4
    ## use average_odor_morph.py to get the sqrt(R2N) distribution plot
    ## and it also calls plot_example_onemit() to plot one mitral example
    
    ## Plot the fitting result of the commandline data file, do not refit
    plot_example(refit=False)

    show()
