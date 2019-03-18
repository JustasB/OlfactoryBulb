# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO ROUGHLY FOLLOW PRIYANKA'S ANALYSIS
########## This variant does not use an air kernel, rather a constant air rate.
## USAGE1: python2.6 fit_odor_pulses.py <../results/odor_pulses/pulseresult.pickle> 191.0
## USAGE2: python2.6 fit_odor_pulses.py <../results/odor_pulses/pulseresult.pickle> 191.0 <../results/odor_morphs/morphresult.pickle>

########## OBSOLETE ######### OBSOLETE ########## OBSOLETE ######### OBSOLETE ######## OBSOLETE ########
####### Use fit_odor_pulses_constair_separateodors.py instead. Priyanka fits each odor separately for each mitral.

from scipy import optimize
from scipy import special
import scipy.interpolate
from scipy import signal
import scipy.io
from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, PULSE_RUNTIME
from networkConstants import * # has central_glom
from data_utils import * # has axes_labels()
from analysis_utils import * # has read_morphfile() and various constants
import generate_firerates_odors as gen_frate_odors # has respiration_function()

iterationnum = 1
sigmoid_op = False
SMOOTH_KERNELS = False

## Have overridden the kernel_time=1.5s (for generating ORN frates) in stimuliConstants.py,
## as Priyanka uses 2.0s kernels.
kernel_time = 2.0

numpulses = RANDOM_PULSE_NUMS*2
## rebin the spikes as below, irrespective of previous binning
pulserebindt = fitting_dt # 50ms as Priyanka uses
pulserebins = int(PULSE_RUNTIME/pulserebindt)
bin_width_time = 2*pulserebindt ## overlapping bins with this bin width
pulsetlist = arange(pulserebindt/2.0,PULSE_RUNTIME,pulserebindt)
## very IMPORTANT to have the same dt for pulse and kernel
## when doing discrete convolution!
kerneldt = pulserebindt # same as fitting_dt above
kernel_size = int(kernel_time/kerneldt)
kerneltlist = arange(0.0,kernel_time,kerneldt)[:kernel_size]
#NUMWTS = NUMMIX-3 # remove the two pure odors and one pure air weights
#NUMPARAMS = 3*NUMBINS+2*NUMWTS+1 # two pure + one air response,
# weights of A and B for mixtures, max of output sigmoid
firstrun = True

## If you want to plot in style of Priyanka
fill_below = True

## number of the mitral to be fitted.
fitted_mitrals = [2*central_glom,2*central_glom+1]

global iternum
iternum = 0

resp_pulses = gen_frate_odors.respirationPulses
st_idx = gen_frate_odors.startidx
numt = gen_frate_odors.numt

def responseconvolve(pulses, kernel):
    ## don't use mode='same', as it takes only the valid part of the convolution.
    ## always multiply convolution_result by dt
    ## to leave response invariant on decimating stim and kernel simultaneously.
    responsecalc = convolve( pulses, kernel, mode='full' ) * pulserebindt
    ## bin0 of the simulated response is at time pulserebindt
    ## so the fitted response should also be shifted by a pulserebindt
    return array( responsecalc[1:pulserebins+1] )
    
def sigmoid_fn(valarray, max, half_point):
    ## ensure that valarray is a numpy array
    ## else where() returns empty tuple, and below exp() statement won't work.
    valarray = array(valarray)
    if sigmoid_op:
        ## exp() - logistic function - not used
        #return max/(1.0+exp(-valarray+half_point)) # numpy's exp works for arrays.
        ## thresholded erf() - used
        ## erf(0.5)~=0.5; set half_point as center of sigmoid.
        valarray[where(valarray<0)[0]] = 0
        posindices = where(valarray<0)[0]
        valarray[posindices] = max * special.erf( 0.5*valarray[posindices]/half_point )
        return valarray
    else:
        ## where value in valarray < 0, set it to zero.
        valarray[where(valarray<0)[0]]=0
        return valarray

#def inv_sigmoid_fn(valarray, max, half_point):
#    ## ensure that valarray is a numpy array
#    ## else where() returns empty tuple, and below exp() statement won't work.
#    valarray = array(valarray)
#    if sigmoid_op:
#        ## numpy has overloaded operators, also numpy ln is overloaded for arrays.
#        return half_point - log(max/valarray-1)
#    else:
#        return valarray

def SGfilter2(kernelA,kernelB):
    ## savitsky-golay filter maintains high frequency signal, unlike moving average filter
    ## window_size must be odd, savitsky_golay copied from scipy cookbook online
    window_size = 11 # same as Priyanka
    kernelA = savitzky_golay(kernelA, window_size, order=4)
    kernelB = savitzky_golay(kernelB, window_size, order=4)
    return kernelA, kernelB

def chisqfunc(params, firingbinsmeanList, firingbinserrList, bgnd, \
        pulserebinList, ret_array=True, ORNkernels=()):
    kernelA = params[0:kernel_size]
    kernelB = params[kernel_size:2*kernel_size]
    if sigmoid_op:
        sigmoid_max = params[-2]
        sigmoid_threshold = params[-1]
    else:
        sigmoid_max = None
        sigmoid_threshold = None
    if SMOOTH_KERNELS:
        kernelA, kernelB = SGfilter2(kernelA, kernelB)
    ## These below should be numpy arrays,
    ## else + acts like append in python (not element-wise addition)!!!
    responsecalcA1 = sigmoid_fn(
        responseconvolve(pulserebinList[2],kernelA) + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcB1 = sigmoid_fn(
        responseconvolve(pulserebinList[3],kernelB) + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcA2 = sigmoid_fn(
        responseconvolve(pulserebinList[4],kernelA) + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcB2 = sigmoid_fn(
        responseconvolve(pulserebinList[5],kernelB) + bgnd,\
        sigmoid_max, sigmoid_threshold)
    ### don't use A+B to fit, just plot the result on A+B.
    #responsecalcAB = sigmoid_fn(
    #    responseconvolve(pulserebinList[6],kernelA) + responseconvolve(pulserebinList[7],kernelB) + \
    #    + bgnd, sigmoid_max, sigmoid_threshold)
    ## fitted response is also shifted by a pulserebindt, to compare index to index
    ## only compare after odor onset at PULSE_START
    starti = int(PULSE_START/pulserebindt)
    chisqarray = [ (val-firingbinsmeanList[1][i])/firingbinserrList[1][i]\
        for i,val in enumerate(responsecalcA1[starti:],start=starti) ]
    chisqarray.extend( [ (val-firingbinsmeanList[2][i])/firingbinserrList[2][i]\
        for i,val in enumerate(responsecalcB1[starti:],start=starti) ] )
    chisqarray.extend( [ (val-firingbinsmeanList[3][i])/firingbinserrList[3][i]\
        for i,val in enumerate(responsecalcA2[starti:],start=starti) ] )
    chisqarray.extend( [ (val-firingbinsmeanList[4][i])/firingbinserrList[4][i]\
        for i,val in enumerate(responsecalcB2[starti:],start=starti) ] )
    #chisqarray.extend( [ (val-firingbinsmeanList[5][i])/firingbinserrList[5][i]\
    #    for i,val in enumerate(responsecalcAB[starti:],start=starti) ] )
        
    global iternum
    iternum += 1
    ## normalize chi-sq to the number of dof
    num_dof = float(firingbinsmeanList[1][starti:].size*(numpulses-2) - params.size)
    if ret_array:
        if iternum % 10000 == 0:
            chisqarraysq = [i**2 for i in chisqarray]
            chisq = reduce(lambda x, y: x+y, chisqarraysq) / num_dof
            print chisq
        ## len(chisqarray) must be >= len(params):
        ## since leastsq() assumes it is solving M eqns in N unknowns where M>=N.
        return array(chisqarray) / sqrt(num_dof) # this must be squared and added for chi-sq
    else:
        chisq = reduce(lambda x, y: x+y**2, chisqarray) / num_dof
        if iternum % 10000 == 0:
            print chisq
        return chisq

def fit_pulses(filename,stimseed,match_morph,NOSHOW):
    ######### load in the mitral pulse responses
    f = open(filename,'r')
    #### mitral_responses_list[avgnum][pulsenum][mitralnum][binnum]
    #### mitral_responses_binned_list[avgnum][pulsenum][mitralnum][binnum]
    mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
    f.close()
    ########## load in the stimuli
    seednum = stimseed
    f = open('../generators/firerates/firerates_2sgm_'+str(seednum)+'.pickle','r')
    frateOdorList,fratePulseList,randomPulseList,\
        randomPulseStepsList,randomResponseList,kernels\
        = pickle.load(f)
    #del frateOdorList, fratePulseList
    pulseList = array(randomPulseList)
    pulseStepsList = array(randomPulseStepsList)
    ## randomResponseList[glomnum][frate,...]
    ORNfrateList = array(randomResponseList)[central_glom]
    f.close()

    ##-------------------------- rebin the responses and pulses ------------------------------
    ### for overlapping pulses (uncomment in analysis_utils.py)
    #mitral_responses_binned_list = \
    #    rebin_pulses(mitral_responses_list, int(PULSE_RUNTIME/pulserebindt), PULSE_RUNTIME, 0.0, bin_width_time)
    mitral_responses_binned_list = \
        rebin_pulses(mitral_responses_list, int(PULSE_RUNTIME/pulserebindt), PULSE_RUNTIME, 0.0)
    ## don't just decimate, rather take average of each bin
    decimatefactor = int(round(pulserebindt/FIRINGFILLDT)) # ensure that these are multiples
    ## binning the pulses: average of every decimatefactor number of values
    pulserebinList = array( [ array([ pulseList[pulsenum][i:i+decimatefactor].mean() \
            for i in range(0,len(pulseList[pulsenum]),decimatefactor) ]) \
        for pulsenum in range(len(pulseList)) ] )
    ## decimating the pulses, i.e. taking one value after every decimatefactor
    pulserebinList = array( [ pulseList[pulsenum][::decimatefactor] \
        for pulsenum in range(len(pulseList)) ] )

    ## very important to convert to numpy array(),
    ## else where() below does not find any element satisfying condition.
    mitral_responses_binned_list = array(mitral_responses_binned_list)
    numtrials = len(mitral_responses_binned_list)
    mitral_responses_mean = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    ## since I fit the mean response, I must use standard error/deviation of the _mean_
    ## = standard deviation of a repeat / sqrt(num of repeats).
    ## NO! The model is for an individual response, not for the mean response!
    #mitral_responses_se = mitral_responses_std/sqrt(numtrials)

    fittedkernels = []
    chisqs = []
    global sigmoid_op
    ## full fitting data for both mitrals
    mdict_full = []
    fitted_responses_full = []
    if not NOSHOW:
        fig = figure(facecolor='w') # 'none' is transparent
        ax_kernelsA = plt.subplot2grid((6,3),(0,2),rowspan=3)
        ax_kernelsB = plt.subplot2grid((6,3),(3,2),rowspan=3)
    for fitted_mitral in fitted_mitrals:
        ## take the odor responses of the mitral to be fitted
        firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]
        firingbinserrList = mitral_responses_std[:,fitted_mitral]
        ## slightly different from Priyanka -- not taking full air buffer, allowing to settle
        bgnd = mean( [response[int((SETTLETIME+PULSE_AIR_BUFFER/2.0)/pulserebindt):\
                int((SETTLETIME+PULSE_AIR_BUFFER)/pulserebindt)] \
            for response in firingbinsmeanList] )
        print "bgnd =",bgnd
        ### force the errors to be 1.2*sqrt(mean) which is what Adil observes.
        ### The model seems to have constant noise independent of the mean firing rate.
        ### standard error of the mean, hence /sqrt(numtrials)
        #firingbinserrList = 1.2*firingbinsmeanList/sqrt(numtrials)

        errcut = 5e-1
        #errcut = 1.0/pulserebindt/sqrt(numtrials) # errcut = 3.33Hz!
        #errcut = 10.0 #Hz
        ## put in a minimum error, else divide by zero problems.
        ## find the minimum error >= errcut
        largeerrors = firingbinserrList[where(firingbinserrList>errcut)]
        if largeerrors!=[]: errmin = largeerrors.min()
        else: errmin = errcut
        print firingbinserrList
        print 'least error = ',errmin
        errmin = errcut
        ## numpy where(), replace by errmin,
        ## all those elements in firingbinsList which are less than errmin 
        firingbinserrList = where(firingbinserrList>errcut, firingbinserrList, errmin)

        #---------------------------- fitting ------------------------------------------------

        decimatefactor = int(kerneldt/FIRINGFILLDT)
        ## kernels are only kernel_size long.
        ## numpy append() is like python's extend() - makes flat not nested lists
        ## I take only kernelsize length after decimating the ORN kernel 
        ## by taking only every int(kerneldt/FIRINGFILLDT) element
        ## since convolve_result is always multiplied by dt,
        ## the response is invariant of decimating the stimulus and the kernel simultaneously.
        ## NOTE: We are using ORN kernel as mitral kernel, typical mitral
        ## firing rates are 5 times larger. So put in factor of 5!!!
        ## If ORN kernel is too long, take only kernel_size, if too short, pad with zeroes.
        MitralORNfiringratio = 5.0
        linORNkernelA = array(kernels[1][::decimatefactor][:kernel_size])*MitralORNfiringratio
        if len(linORNkernelA)<kernel_size:
            linORNkernelA = append(linORNkernelA,zeros(kernel_size-len(linORNkernelA)))
        linORNkernelB = array(kernels[2][::decimatefactor][:kernel_size])*MitralORNfiringratio
        if len(linORNkernelB)<kernel_size:
            linORNkernelB = append(linORNkernelB,zeros(kernel_size-len(linORNkernelB)))
        ## Initial values for the parameters
        if firstrun:
            ORNkernels = (linORNkernelA, linORNkernelB)
            params = linORNkernelA # ORN kernelA is initial kernelA
            params = append(params,linORNkernelB) # ORN kernelB is initial kernelB
            if sigmoid_op:
                ## I don't have a steepness param,\
                ## else there will be redundancy in kernel scaling vs steepness.
                params = append(params, [120.0,60.0]) # 120Hz max and 60Hz threshold for sigmoid
        else:
            f = open(sys.argv[1]+'_params'+str(fitted_mitral),'r')
            sigmoid_op,chisq,params,discard_firingmeans,discard_firingsd,discard_fitresponses = pickle.load(f)
            f.close()
            ORNkernels = ()

        ## export the data for use in matlab in Priyanka's fitting function
        mdict={'firingbinsmeanList':firingbinsmeanList,
            'firingbinserrList':firingbinserrList, 'bgnd':bgnd, 'FIRINGFILLDT':FIRINGFILLDT,
            'pulseList':pulserebinList,'pulseStepsList':pulseStepsList,
            'start_kernels':ORNkernels,'pulserebindt':pulserebindt,'pulsetlist':pulsetlist}
        scipy.io.savemat('mit'+str(fitted_mitral)+'.mat', mdict=mdict)
        mdict_full.append(mdict)

        ## comment next fitting line if you just want to plot parameters.
        ## args is a tuple! if only one element write (elem, )

        #########
        ### fmin_powell is faster than fmin (simplex) but slower than leastsq
        ### However it does not enter into the very bad local minima of leastsq
        ### WARNING! But often it is worse at fitting than leastsq()
        #iternum = 0
        #params = optimize.fmin_powell( chisqfunc, params, \
        #    args=(firingbinsmeanList, firingbinserrList, pulserebinList, False),\
        #    full_output=1, maxfun=50000, maxiter=50000, ftol = 1 )
        #chisq = params[1]
        #print "The normalized chisq for mitral",fitted_mitral,"after fmin_powell =",chisq
        #params = params[0]

        #########
        ## leastsq() uses a modified version of Levenberg-Marquardt algorithm
        ## it optimizes M equations in N unknowns, where M>=N.
        ## Hence chisqfunc must return M numbers in an array
        ## which must be >= the number of params N in params0, 
        ## else: 'TypeError: Improper input parameters.'
        iternum = 0
        params = optimize.leastsq( chisqfunc, params,\
            args=(firingbinsmeanList, firingbinserrList, bgnd, pulserebinList, True, ORNkernels),\
            full_output=1, maxfev=100000, ftol=1e-100, xtol=1e-100 )
        print params[3] # print the status message
        params = params[0] # take only fitted params; leastsq returns a tuple with errmsg, etc.
        ### Calculate sum of squares of the chisqarray
        chisqarraysq = [i**2 for i in chisqfunc(params,\
            firingbinsmeanList, firingbinserrList, bgnd, pulserebinList, True, ORNkernels)]
        #chisqnumpy = array(chisqarraysq)
        #pos = where(chisqnumpy>100)[0]
        #lenf = len(firingbinsmeanList[0])
        #print "The [(pulsenum,index),...] where chisq>100 are",[ (i/lenf,i%lenf) for i in pos]
        chisq = reduce(lambda x, y: x+y, chisqarraysq)
        print "The normalized chisq for mitral",fitted_mitral,"after leastsq =",chisq
        
        ##########
        ### minimize is a wrapper for all the various fitting algorithms!
        #iternum = 0
        #fit_method = "Nelder-Mead"#"CG"#"Powell"#"Nelder-Mead"
        #params,info = optimize.minimize( chisqfunc, params,\
        #    args=(firingbinsmeanList, firingbinserrList, pulserebinList, False, ORNkernels),\
        #    method = fit_method, full_output=True,\
        #    options={'maxfev':100000, 'ftol':1e-100, 'xtol':1e-100, 'maxfun':500000, 'maxiter':500000} )
        #print info['message']
        #chisq = info['fun']
        #print "The normalized chisq for mitral",fitted_mitral,"after",fit_method,"=",chisq

        kernelA = params[0:kernel_size]
        kernelB = params[kernel_size:2*kernel_size]
        ## irrespective of smooth kernels flag, filter the output kernels for display
        kernelA, kernelB = SGfilter2(kernelA, kernelB)

        if sigmoid_op:
            sigmoid_max = params[-2]
            sigmoid_threshold = params[-1]
        else:
            sigmoid_max = None
            sigmoid_threshold = None

        #chisq = chisqfunc(params,firingbinsmeanList, firingbinserrList, pulserebinList)

        fittedkernels.append( (kernelA,kernelB) )
        chisqs.append(chisq)

        fitted_responses = []
        for pulsenum in range(1,6):
            if pulsenum in [1,3]:
                response = sigmoid_fn(responseconvolve(pulserebinList[pulsenum+1],kernelA) + bgnd,\
                    sigmoid_max, sigmoid_threshold)
            elif pulsenum in [2,4]:
                response = sigmoid_fn(responseconvolve(pulserebinList[pulsenum+1],kernelB) + bgnd,\
                    sigmoid_max, sigmoid_threshold)
            else:
                response = responseconvolve(pulserebinList[pulsenum+1],kernelA)
                response += responseconvolve(pulserebinList[pulsenum+2],kernelB)
                response += bgnd
                response = sigmoid_fn(response, sigmoid_max, sigmoid_threshold)
            fitted_responses.append(response)
        fitted_responses_full.append(fitted_responses)

        paramsfile = open(filename+'_params'+str(fitted_mitral),'w')
        pickle.dump((sigmoid_op,chisq,params,firingbinsmeanList,firingbinserrList,fitted_responses), paramsfile)
        paramsfile.close()


        ##---------------------------- kernel via the m-sequence method -----------------------------------------
        ## Doesn't work! One reason is that it doesn't take care of thresholding of neural response!

        ## only compare for PULSE_DURATION from PULSE_START
        starti = int(PULSE_START/pulserebindt)
        endi = int((PULSE_START+PULSE_DURATION)/pulserebindt)
        stimulus = pulserebinList[2][starti:endi]
        response = firingbinsmeanList[1][starti:endi] # response to odor A first pulse
        kernelA_mseq0 = circular_correlate(stimulus,response) # from data_utils.py
        stimulus = pulserebinList[4][starti:endi]
        response = firingbinsmeanList[3][starti:endi] # response to odor A first pulse
        kernelA_mseq1 = circular_correlate(stimulus,response) # from data_utils.py
        kernelA_mseq = (kernelA_mseq0+kernelA_mseq1)/2.0
        #figure()
        #pulsepoints = arange(0,PULSE_DURATION,pulserebindt)
        #plot(pulsepoints,kernelA_mseq)
        #plot(pulsepoints,stimulus)
        #plot(pulsepoints,response)


        #---------------------------- plot kernels and pulse responses -----------------------------------------

        if not NOSHOW:
            ## fig defined before the mitnum loop
            ############################## plot the kernels
            text(0.6,0.3,'arb units', fontsize=label_fontsize, rotation='vertical')

            ax_kernelsA.plot(kerneltlist, kernelA, color=(1-fitted_mitral,fitted_mitral,0), marker=',',\
                linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
            #ax.plot(kerneltlist, linORNkernelA, color=(1,0,1), marker=',',\
            #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
            ax_kernelsB.plot(kerneltlist, kernelB, color=(fitted_mitral,1-fitted_mitral,0), marker=',',\
                linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
            #ax.plot(kerneltlist, linORNkernelB, color=(0,1,1), marker=',',\
            #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
            ax_kernelsA.set_yticks([0])
            ax_kernelsB.set_yticklabels(['0'])
            ax_kernelsA.set_xlim(0,kernel_time)
            ax_kernelsB.set_xlim(0,kernel_time)
            biglegend(ax=ax_kernelsA)
            biglegend(ax=ax_kernelsB)
            ax_kernelsB.set_xticks([0,2.0])
            ax_kernelsB.set_xticklabels(['0','2.0'])
            axes_off(ax_kernelsA,x=True,y=False)
            axes_labels(ax_kernelsA,'','',adjustpos=False)
            axes_labels(ax_kernelsB,'time (s)','',adjustpos=False)

            ############################### plot the responses and the fits
            text(-0.125,0.25,'firing rate (Hz)', fontsize=label_fontsize, rotation='vertical')
            for pulseiter, pulsenum in enumerate([1,2,5]):#range(1,numpulses): # similar to Priyanka, only 3 plots
                sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
                ## similar to Priyanka: only 3 plots, skip 2 pulse responses
                ax = plt.subplot2grid((6,3),(2*pulseiter,fitted_mitral),rowspan=2)
                ################# random pulses and ORN responses
                xpulse = int(max(firingbinsmeanList[pulsenum]+firingbinserrList[pulsenum]/sqrt(9))+1)
                randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
                if pulsenum in [1,2,3,4]: # odor A or B
                    if fill_below:
                        if pulsenum in [1,3]: col = (1,0,0)
                        if pulsenum in [2,4]: col = (0,1,0)
                        fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+1]+
                            pulseList[0]+ORN_BGND_MEAN,linewidth=0,color=col,alpha=0.4)
                    else:
                        plot(randompulsetime,2*pulseList[pulsenum+1]+pulseList[0]+ORN_BGND_MEAN,
                            linewidth=linewidth,label='Air+Odor waveform')
                        plot(randompulsetime,ORNfrateList[pulsenum+1]+ORNfrateList[0]+ORN_BGND_MEAN,
                            linewidth=linewidth,label='Receptors response')
                elif pulsenum in [5]: # odor A & odor B
                    if fill_below:
                        fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+1],\
                            color=(1,0,0),linewidth=0,alpha=0.4)
                        fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+2],\
                            color=(0,1,0),linewidth=0,alpha=0.4)
                    else:
                        plot(randompulsetime,pulseList[0]+ORN_BGND_MEAN,\
                            linewidth=linewidth, label='Air waveform')
                        plot(randompulsetime,2*pulseList[pulsenum+1],\
                            linewidth=linewidth, label='Odor A waveform')
                        plot(randompulsetime,2*pulseList[pulsenum+2],\
                            linewidth=linewidth, label='Odor B waveform')
                        plot(randompulsetime,ORNfrateList[pulsenum+2]+ORNfrateList[pulsenum+1]\
                            +ORNfrateList[0]+ORN_BGND_MEAN, linewidth=linewidth, label='Receptors response')
                ################### Plot the simulated responses
                ## smooth the simulated response
                ## windowsize=5 and SD=0.65 are defaults from matlab's smoothts() for gaussian smoothing
                Gwindow = signal.gaussian(5,0.65)
                ## help from http://www.scipy.org/Cookbook/SignalSmooth
                simresponse = convolve(Gwindow/Gwindow.sum(),firingbinsmeanList[pulsenum],mode='same')
                if fill_below:
                    ## numpy array, hence adds element by element
                    fill_between(pulsetlist,
                        simresponse+firingbinserrList[pulsenum]/sqrt(9),
                        simresponse-firingbinserrList[pulsenum]/sqrt(9),
                        color=(0.7,0.7,0.7))
                    plot(pulsetlist,simresponse,linewidth=linewidth,color=(0.3,0.3,0.3))
                else:
                    errorbar(pulsetlist,y=simresponse,\
                        yerr=firingbinserrList[pulsenum]/sqrt(9),\
                        color=(pulsenum/float(numpulses),1-pulsenum/float(numpulses),sister_ratio),\
                        marker='+',linestyle='solid',linewidth=linewidth,label='Mitral response')
                ################## Plot the fitted responses
                starti = int(PULSE_START/pulserebindt)
                if pulsenum in [1,3]: titlestr = 'Odor A random pulse'
                elif pulsenum in [2,4]: titlestr = 'Odor B random pulse'
                else: titlestr = 'Odor A and odor B random pulse'
                response = fitted_responses[pulsenum-1]
                if fill_below:
                    plot(pulsetlist[starti:],response[starti:],linestyle='solid',
                        linewidth=linewidth,color=(0,0,0),
                        label=['air','A','B','A','B','A&B'][pulsenum])
                    lgd = legend()
                    for k in lgd.get_lines():
                        k.set_linewidth(0)
                    lgd.draw_frame(False)
                    ltext  = lgd.get_texts()
                    for l in ltext:
                        l.set_fontsize(label_fontsize)
                    #fig.setp(ltext, fontsize=20)
                else:
                    plot(pulsetlist[starti:],response[starti:],
                        marker='o',linestyle='dashed',
                        linewidth=linewidth, label='Linear fit')
                    #title('Odor morphs w/ smooth fit',fontsize=24 )
                    #title( 'pulsenum = '+str(pulsenum)+', chisquare normalized = '+str(chisq) )
                    title(titlestr, fontsize=label_fontsize)
                    lgd = legend()
                    ltext  = lgd.get_texts()
                    for l in ltext:
                        l.set_fontsize(label_fontsize)
                    ax.set_ylim(0,xpulse)
                ax.set_yticks([xpulse])
                ax.set_yticklabels([str(xpulse)])
                for label in ax.get_yticklabels():
                    label.set_fontsize(label_fontsize)
                ax.set_xlim(0,9)
                ax.set_xticks([])
            ax.set_xticks([0,9])
            ax.set_xticklabels(['0','9'])    
            axes_labels(ax,'time (s)','',adjustpos=False)

    if not NOSHOW:
        show()
    
    return fittedkernels,chisqs,fitted_responses_full,mdict_full

if __name__ == "__main__":
    if len(sys.argv) > 2:
        ######### should we match mitral morph responses using kernel
        if '--match_morph' in sys.argv:
            match_morph = True
            mmfilenameindex = sys.argv.index('--match_morph')+1
            morph_filename = sys.argv[mmfilenameindex]
        else: match_morph = False
        if 'NOSHOW' in sys.argv: NOSHOW = True
        else: NOSHOW = False
        fit_pulses(sys.argv[1],sys.argv[2],match_morph,NOSHOW)
    else:
        print "At least specify data file containing pickled mitral responses, and ORN frate seed."
        sys.exit(1)
