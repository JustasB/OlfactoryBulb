# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO ROUGHLY FOLLOW PRIYANKA'S ANALYSIS
#### OBSOLETE -- Priyanka does not fit air kernel! Use fit_odor_pulses_constair.py instead
#### This program also has no rebinning with 50ms bin and no padding of init kernel length with zeros.
## USAGE1: python2.6 fit_odor_pulses.py <../results/odor_pulses/pulseresult.pickle> 191.0
## USAGE2: python2.6 fit_odor_pulses.py <../results/odor_pulses/pulseresult.pickle> 191.0 <../results/odor_morphs/morphresult.pickle>

########## OBSOLETE ######### OBSOLETE ########## OBSOLETE ######### OBSOLETE ######## OBSOLETE ########
####### Use fit_odor_pulses_constair.py instead. Priyanka does not fit air kernel, nor A+B as in _constair.

from scipy import optimize
from scipy import special
import scipy.interpolate
from pylab import *
import pickle
import sys
import math

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, PULSE_RUNTIME, pulsebins, pulsebindt
from networkConstants import * # has central_glom
from data_utils import * # has axes_labels()
from analysis_utils import * # has read_morphfile() and various constants
import generate_firerates_odors as gen_frate_odors # has respiration_function()

iterationnum = 1
sigmoid = False
FitScaledORNkernels = False#True
SMOOTH_KERNELS = False

numpulses = RANDOM_PULSE_NUMS*2
pulsetlist = arange(pulsebindt/2.0,PULSE_RUNTIME,pulsebindt)
## very IMPORTANT to have the same dt for pulse and kernel
## when doing discrete convolution!
kerneldt = pulsebindt
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

## for matching morphs
match_morph = False
## for the morphs, I don't use the NUMBINS in simset_odor.py,
## rather I rebin() with an overlap
NUM_REBINS = 20
resp_pulses = gen_frate_odors.respirationPulses
st_idx = gen_frate_odors.startidx
numt = gen_frate_odors.numt

def responseconvolve(pulses, kernel):
    ## pulsebindt is typically unchanged, pulse_time can change
    ## changes number of pulsebins in stimuliConstants.py
    decimatefactor = int(round(pulsebindt/FIRINGFILLDT)) # ensure that these are multiples
    pulses = pulses[::decimatefactor]
    ## don't use mode='same', as it takes only the valid part of the convolution.
    ## always multiply convolution_result by dt
    ## to leave response invariant on decimating stim and kernel simultaneously.
    responsecalc = convolve( pulses, kernel, mode='full' ) * pulsebindt
    ## bin0 of the simulated response is at time pulsebindt
    ## so the fitted response should also be shifted by a pulsebindt
    return array( responsecalc[1:pulsebins+1] )
    
def sigmoid_fn(valarray, max, half_point):
    ## ensure that valarray is a numpy array
    ## else where() returns empty tuple, and below exp() statement won't work.
    valarray = array(valarray)
    if sigmoid:
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
#    if sigmoid:
#        ## numpy has overloaded operators, also numpy ln is overloaded for arrays.
#        return half_point - log(max/valarray-1)
#    else:
#        return valarray

def SGfilter3(kernelR,kernelA,kernelB):
    ## savitsky-golay filter maintains high frequency signal, unlike moving average filter
    ## window_size must be odd, savitsky_golay copied from scipy cookbook online
    if int(kernel_size/2)%2==0: window_size = int(kernel_size/2)+1
    else: window_size = int(kernel_size/2)
    kernelR = savitzky_golay(kernelR, window_size, order=5)
    kernelA = savitzky_golay(kernelA, window_size, order=5)
    kernelB = savitzky_golay(kernelB, window_size, order=5)
    return kernelR, kernelA, kernelB

def chisqfunc(params, firingbinsmeanList, firingbinserrList,\
        pulseList, ret_array=True, ORNkernels=()):
    if not FitScaledORNkernels:
        kernelR = params[0:kernel_size]
        kernelA = params[kernel_size:2*kernel_size]
        kernelB = params[2*kernel_size:3*kernel_size]
        bgnd = params[3*kernel_size:3*kernel_size+1]
    else:
        scaleR = params[0]
        if scaleR<0: scaleR = 0.0
        scaleA = params[1]
        if scaleA<0: scaleA = 0.0
        scaleB = params[2]
        if scaleB<0: scaleB = 0.0
        kernelR = scaleR*ORNkernels[0]
        kernelA = scaleA*ORNkernels[1]
        kernelB = scaleB*ORNkernels[2]
        bgnd = params[3]
    if sigmoid:
        sigmoid_max = params[-2]
        sigmoid_threshold = params[-1]
    else:
        sigmoid_max = None
        sigmoid_threshold = None
    if SMOOTH_KERNELS:
        kernelR, kernelA, kernelB = SGfilter3(kernelR, kernelA, kernelB)
    ## These below should be numpy arrays,
    ## else + acts like append in python (not element-wise addition)!!!
    pulseair = responseconvolve(pulseList[0],kernelR)
    responsecalcair = sigmoid_fn(
        responseconvolve(pulseList[1],kernelR) + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcA1 = sigmoid_fn(
        responseconvolve(pulseList[2],kernelA) + pulseair + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcB1 = sigmoid_fn(
        responseconvolve(pulseList[3],kernelB) + pulseair + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcA2 = sigmoid_fn(
        responseconvolve(pulseList[4],kernelA) + pulseair + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcB2 = sigmoid_fn(
        responseconvolve(pulseList[5],kernelB) + pulseair + bgnd,\
        sigmoid_max, sigmoid_threshold)
    responsecalcAB = sigmoid_fn(
        responseconvolve(pulseList[6],kernelA) + responseconvolve(pulseList[7],kernelB) + \
        + pulseair + bgnd, sigmoid_max, sigmoid_threshold)
    ## fitted response is also shifted by a pulsebindt, to compare index to index
    chisqarray = [ (val-firingbinsmeanList[0][i])/firingbinserrList[0][i]\
        for i,val in enumerate(responsecalcair) ]
    chisqarray.extend( [ (val-firingbinsmeanList[1][i])/firingbinserrList[1][i]\
        for i,val in enumerate(responsecalcA1) ] )
    chisqarray.extend( [ (val-firingbinsmeanList[2][i])/firingbinserrList[2][i]\
        for i,val in enumerate(responsecalcB1) ] )
    chisqarray.extend( [ (val-firingbinsmeanList[3][i])/firingbinserrList[3][i]\
        for i,val in enumerate(responsecalcA2) ] )
    chisqarray.extend( [ (val-firingbinsmeanList[4][i])/firingbinserrList[4][i]\
        for i,val in enumerate(responsecalcB2) ] )
    chisqarray.extend( [ (val-firingbinsmeanList[5][i])/firingbinserrList[5][i]\
        for i,val in enumerate(responsecalcAB) ] )
        
    global iternum
    iternum += 1
    ## normalize chi-sq to the number of dof
    num_dof = float(firingbinsmeanList.size - params.size)
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

def match_morphs_from_kernels(kernels,bgnd):
    kernelR = kernels[0]
    kernelA = kernels[1]
    kernelB = kernels[2]
    ## resp_pulses is peak-normalized to 1; the on-off random pulses also go to 1.
    ## Thus, no factor in response for scaling of stimulus
    ## However the pulse, response, kernels were using dt=kerneldt=pulsebindt
    ## But the kernel is a binned kernel, now interpolate to have dt = FIRINGFILLDT
    ## FIRINGFILLDT is the dt for respirationPulses and odor ORN frates
    interpol_type = 'cubic' # can have 'cubic', 'linear', etc.
    kernelR_fn = scipy.interpolate.interp1d(kerneltlist, kernelR, kind=interpol_type)
    kernelA_fn = scipy.interpolate.interp1d(kerneltlist, kernelA, kind=interpol_type)
    kernelB_fn = scipy.interpolate.interp1d(kerneltlist, kernelB, kind=interpol_type)
    ## make a finer dt for kernel = FIRINGFILLDT,
    ## need to subtract kerneldt from endpoint to keep within interpolation function range
    kerneltlist_finer = arange(0.0,kernel_time-kerneldt,FIRINGFILLDT)
    kernelR_full = [kernelR_fn(t) for t in kerneltlist_finer]
    kernelA_full = [kernelA_fn(t) for t in kerneltlist_finer]
    kernelB_full = [kernelB_fn(t) for t in kerneltlist_finer]
    ## always multiply convolution_result by dt to leave response invariant
    ## on decimating/un-decimating stim and kernel simultaneously.
    respirationPulses = gen_frate_odors.respirationPulses
    responseR = convolve(respirationPulses,kernelR_full,mode='full') * FIRINGFILLDT
    responseR = responseR[st_idx:st_idx+numt]
    #responseR = responseR[:numt]
    responseA = convolve(respirationPulses,kernelA_full,mode='full') * FIRINGFILLDT
    responseA = responseA[st_idx:st_idx+numt]
    #responseA = responseA[:numt]
    responseB = convolve(respirationPulses,kernelB_full,mode='full') * FIRINGFILLDT
    responseB = responseB[st_idx:st_idx+numt]
    #responseB = responseB[:numt]

    matched_responses = []
    for inpnum,(inputA,inputB) in enumerate(inputList[:-1]):
        ## ideally I should use the fitted input weights for this morph
        matched_response = inputA*responseA + inputB*responseB + responseR + bgnd
        ## clip values going negative
        matched_response[ where(matched_response<0)[0] ] = 0
        matched_responses.append( matched_response )
    ## clip values going negative
    responseR[ where(responseR<0)[0] ] = 0
    matched_responses.append( responseR +bgnd )
    return matched_responses

if __name__ == "__main__":
    if len(sys.argv) > 2:
        ######### load in the mitral pulse responses
        f = open(sys.argv[1],'r')
        #### mitral_responses_list[avgnum][pulsenum][mitralnum][binnum]
        #### mitral_responses_binned_list[avgnum][pulsenum][mitralnum][binnum]
        mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
        f.close()
        ########## load in the stimuli
        seednum = sys.argv[2]
        f = open('../generators/firerates/firerates_2sgm_'+seednum+'.pickle','r')
        frateOdorList,fratePulseList,randomPulseList,\
            randomPulseStepsList,randomResponseList,kernels\
            = pickle.load(f)
        #del frateOdorList, fratePulseList
        pulseList = array(randomPulseList)
        pulseStepsList = array(randomPulseStepsList)
        ## randomResponseList[glomnum][frate,...]
        ORNfrateList = array(randomResponseList)[central_glom]
        f.close()
        ######### should we match mitral morph responses using kernel
        if len(sys.argv)>3:
            match_morph = True
            morph_filename = sys.argv[3]
    else:
        print "Specify data file containing pickled mitral responses, and ORN frate seed."
        sys.exit(1)
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
    
    for fitted_mitral in fitted_mitrals:
        ## take the odor responses of the mitral to be fitted
        firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]
        firingbinserrList = mitral_responses_std[:,fitted_mitral]

        ### force the errors to be 1.2*sqrt(mean) which is what Adil observes.
        ### The model seems to have constant noise independent of the mean firing rate.
        ### standard error of the mean, hence /sqrt(numtrials)
        #firingbinserrList = 1.2*firingbinsmeanList/sqrt(numtrials)

        errcut = 1e-20
        #errcut = 1.0/pulsebindt/sqrt(numtrials) # errcut = 3.33Hz!
        #errcut = 10.0 #Hz
        ## put in a minimum error, else divide by zero problems.
        ## find the minimum error >= errcut
        largeerrors = firingbinserrList[where(firingbinserrList>errcut)]
        if largeerrors!=[]: errmin = largeerrors.min()
        else: errmin = errcut
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
        MitralORNfiringratio = 5.0
        linORNkernelR = array(kernels[0][::decimatefactor][:kernel_size])*MitralORNfiringratio
        linORNkernelA = array(kernels[1][::decimatefactor][:kernel_size])*MitralORNfiringratio
        linORNkernelB = array(kernels[2][::decimatefactor][:kernel_size])*MitralORNfiringratio
        ## Initial values for the parameters
        if firstrun:
            ORNkernels = (linORNkernelR, linORNkernelA, linORNkernelB)
            if not FitScaledORNkernels:
                params = linORNkernelR # ORN kernelR is initial kernelR
                params = append(params,linORNkernelA) # ORN kernelA is initial kernelA
                params = append(params,linORNkernelB) # ORN kernelB is initial kernelB
                ## constant background mitral firing w/o air and odor due to ORN bgnd.
                params = append(params,ORN_BGND_MEAN*MitralORNfiringratio)
            else:
                ## initial scaling of the ORN kernels to fit mitral kernels
                params = array([1.0,1.0,1.0])
                params = append(params,0.0) # mitral background firing rate
            #figure()
            #plot(linORNkernelR,color=(0,0,0))
            #plot(linORNkernelA,color=(1,0,0))
            #plot(linORNkernelB,color=(0,1,0))
            #figure()
            #plot(responseconvolve(pulseList[5],linORNkernelB))
            ##show()
            ##sys.exit()
            if sigmoid:
                ## I don't have a steepness param,\
                ## else there will be redundancy in kernel scaling vs steepness.
                params = append(params, [120.0,60.0]) # 120Hz max and 60Hz threshold for sigmoid
        else:
            f = open(sys.argv[1]+'_params'+str(fitted_mitral),'r')
            sigmoid,chisq,params = pickle.load(f)
            f.close()
            ORNkernels = ()

        ## comment next fitting line if you just want to plot parameters.
        ## args is a tuple! if only one element write (elem, )

        #########
        ### fmin_powell is faster than fmin (simplex) but slower than leastsq
        ### However it does not enter into the very bad local minima of leastsq
        ### WARNING! But often it is worse at fitting than leastsq()
        #iternum = 0
        #params = optimize.fmin_powell( chisqfunc, params, \
        #    args=(firingbinsmeanList, firingbinserrList, pulseList, False),\
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
            args=(firingbinsmeanList, firingbinserrList, pulseList, True, ORNkernels),\
            full_output=1, maxfev=100000, ftol=1e-100, xtol=1e-100 )
        print params[3] # print the status message
        params = params[0] # take only fitted params; leastsq returns a tuple with errmsg, etc.
        ### Calculate sum of squares of the chisqarray
        chisqarraysq = [i**2 for i in chisqfunc(params,\
            firingbinsmeanList, firingbinserrList, pulseList, True, ORNkernels)]
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
        #    args=(firingbinsmeanList, firingbinserrList, pulseList, False, ORNkernels),\
        #    method = fit_method, full_output=True,\
        #    options={'maxfev':100000, 'ftol':1e-100, 'xtol':1e-100, 'maxfun':500000, 'maxiter':500000} )
        #print info['message']
        #chisq = info['fun']
        #print "The normalized chisq for mitral",fitted_mitral,"after",fit_method,"=",chisq

        if not FitScaledORNkernels:
            kernelR = params[0:kernel_size]
            kernelA = params[kernel_size:2*kernel_size]
            kernelB = params[2*kernel_size:3*kernel_size]
            bgnd = params[3*kernel_size:3*kernel_size+1]
        else:
            kernelR = ORNkernels[0]*params[0]
            kernelA = ORNkernels[1]*params[1]
            kernelB = ORNkernels[2]*params[2]
            print params[0],params[1],params[2]
            bgnd = params[3]
        if SMOOTH_KERNELS:
            kernelR, kernelA, kernelB = SGfilter3(kernelR, kernelA, kernelB)

        if sigmoid:
            sigmoid_max = params[-2]
            sigmoid_threshold = params[-1]
        else:
            sigmoid_max = None
            sigmoid_threshold = None
        pulseair = responseconvolve(pulseList[0],kernelR)
        
        #chisq = chisqfunc(params,firingbinsmeanList, firingbinserrList, pulseList)

        paramsfile = open(sys.argv[1]+'_params'+str(fitted_mitral),'w')
        pickle.dump((sigmoid,chisq,params), paramsfile)
        paramsfile.close()


#---------------------------- matching morphs using kernels -----------------------------------


        if match_morph:
            fitkernels = (kernelR,kernelA,kernelB)
            matched_responses = match_morphs_from_kernels(fitkernels,bgnd)

#---------------------------- plot matched responses -------------------------------------------

            mfiringbinsmeanList,mfiringbinserrList = \
                read_morphfile(morph_filename,fitted_mitral,NUM_REBINS)
            bindt = RESPIRATION/float(NUM_REBINS)
            responsetlist = arange( SETTLETIME+RESPIRATION+bindt, \
                SETTLETIME+2*RESPIRATION+bindt/2.0, bindt)
            figure()
            firingtsteps = gen_frate_odors.firingtsteps
            plot(firingtsteps, matched_responses[0], color=(0,1,0))
            plot(firingtsteps, matched_responses[-2], color=(1,0,0))
            plot(firingtsteps, matched_responses[-1], color=(0,0,0))
            plot(responsetlist, mfiringbinsmeanList[0], color=(0,1,0.5))
            plot(responsetlist, mfiringbinsmeanList[-2], color=(1,0,0.5))
            plot(responsetlist, mfiringbinsmeanList[-1], color=(0,0,0.5))
            #plot(firingtsteps, frateOdorList[central_glom][0], color=(0,1,0.5))
            #plot(firingtsteps, frateOdorList[central_glom][-2], color=(1,0,0.5))
            #plot(firingtsteps, frateOdorList[central_glom][-1], color=(0,0,0.5))
            #plot(firingtsteps, gen_frate_odors.respirationPulses[st_idx:st_idx+numt], color=(1,1,0.5))

#---------------------------- plot kernels and pulse responses -----------------------------------------


        ############################## plot the kernels
        fig = figure(facecolor='w') # 'none' is transparent
        ## A super axes to set a common ylabel
        bigAxes = axes(frameon=False) # hide frame
        xticks([]) # don't want to see any ticks on this axis
        yticks([])
        text(-0.125,0.3,'arb units', fontsize=24, rotation='vertical')

        #ax = fig.add_subplot(111)
        ax = fig.add_subplot(3,1,1)
        ax.set_yticks([0,75])
        ax.set_yticklabels(['0','75'])
        for label in ax.get_yticklabels():
            label.set_fontsize(24)
        axes_off(ax,x=True,y=False)
        ax.plot(kerneltlist, kernelA, color=(1,0,0), marker=',',\
            linestyle='solid',linewidth=2,label='OdorA kernel')
        ax.plot(kerneltlist, linORNkernelA, color=(1,0,1), marker=',',\
            linestyle='solid',linewidth=2,label='5 * ORN kernel')
        #ax.set_ylim(-10,25)
        #ax.set_xlim(0,2)
        biglegend()

        ax = fig.add_subplot(3,1,2)
        ax.set_yticks([-50,0,50])
        ax.set_yticklabels(['-50','0','50'])    
        for label in ax.get_yticklabels():
            label.set_fontsize(24)
        axes_off(ax,x=True,y=False)
        ax.plot(kerneltlist, kernelB, color=(0,1,0), marker=',',\
            linestyle='solid',linewidth=2,label='OdorB kernel')
        ax.plot(kerneltlist, linORNkernelB, color=(0,1,1), marker=',',\
            linestyle='solid',linewidth=2,label='5 * ORN kernel')
        #ax.set_xlim(0,2)
        #ax.set_ylim(-20,20)
        biglegend()

        ax = fig.add_subplot(3,1,3)
        ax.set_yticks([0,50])
        ax.set_yticklabels(['0','50'])    
        for label in ax.get_yticklabels():
            label.set_fontsize(24)
        ax.set_xticks([0,2])
        ax.set_xticklabels(['0','2'])    
        ax.plot(kerneltlist, kernelR, color=(0,0,0), marker=',',\
            linestyle='solid',linewidth=2,label='Air kernel')
        ax.plot(kerneltlist, linORNkernelR, color=(0,0,1), marker=',',\
            linestyle='solid',linewidth=2,label='5 * ORN kernel')
        axes_labels(ax,'time (s)','',adjustpos=False)
        #ax.set_xlim(0,2)
        #ax.set_ylim(-10,20)
        biglegend()
        fig.suptitle('Mitral '+str(fitted_mitral)+\
            ': odor and air kernels',fontsize=24)

        ############################### plot the responses and the fits
        fig = figure(facecolor='w')
        ## A super axes to set a common ylabel
        bigAxes = axes(frameon=False) # hide frame
        xticks([]) # don't want to see any ticks on this axis
        yticks([])
        #ylabel('firing rate (Hz)', fontsize=24)
        text(-0.125,0.25,'firing rate (Hz)', fontsize=28, rotation='vertical')
        for pulsenum in range(numpulses):
            sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            ax = fig.add_subplot(6,1,pulsenum+1)
            ################# random pulses and ORN responses
            xpulse = 50
            randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
            if pulsenum in [0]: # air response does not have air pedestal
                if fill_below:
                    fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+1]+
                        ORN_BGND_MEAN,linewidth=0,color=(0,0.9,0.9),alpha=0.4)
                else:
                    plot(randompulsetime,2*pulseList[pulsenum+1]+ORN_BGND_MEAN,\
                        linewidth=2, label='Air waveform')
                    plot(randompulsetime,ORNfrateList[pulsenum+1]+ORN_BGND_MEAN,\
                        linewidth=2, label='Receptors response')
            elif pulsenum in [1,2,3,4]: # odor A or B
                if fill_below:
                    if pulsenum in [1,3]: col = (1,0,0)
                    if pulsenum in [2,4]: col = (0,1,0)
                    fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+1]+
                        pulseList[0]+ORN_BGND_MEAN,linewidth=0,color=col,alpha=0.4)
                else:
                    plot(randompulsetime,2*pulseList[pulsenum+1]+pulseList[0]+ORN_BGND_MEAN,
                        linewidth=2,label='Air+Odor waveform')
                    plot(randompulsetime,ORNfrateList[pulsenum+1]+ORNfrateList[0]+ORN_BGND_MEAN,
                        linewidth=2,label='Receptors response')
            elif pulsenum in [5]: # odor A & odor B
                if fill_below:
                    fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+1],\
                        color=(1,0,0),linewidth=0,alpha=0.4)
                    fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+2],\
                        color=(0,1,0),linewidth=0,alpha=0.4)
                else:
                    plot(randompulsetime,pulseList[0]+ORN_BGND_MEAN,\
                        linewidth=2, label='Air waveform')
                    plot(randompulsetime,2*pulseList[pulsenum+1],\
                        linewidth=2, label='Odor A waveform')
                    plot(randompulsetime,2*pulseList[pulsenum+2],\
                        linewidth=2, label='Odor B waveform')
                    plot(randompulsetime,ORNfrateList[pulsenum+2]+ORNfrateList[pulsenum+1]\
                        +ORNfrateList[0]+ORN_BGND_MEAN, linewidth=2, label='Receptors response')
            ################### Plot the simulated responses
            if fill_below:
                ## numpy array, hence adds element by element
                fill_between(pulsetlist,
                    firingbinsmeanList[pulsenum]+firingbinserrList[pulsenum],
                    firingbinsmeanList[pulsenum]-firingbinserrList[pulsenum],
                    color=(0.7,0.7,0.7))
                plot(pulsetlist,firingbinsmeanList[pulsenum],linewidth=2,color=(0.3,0.3,0.3))
            else:
                errorbar(pulsetlist,y=firingbinsmeanList[pulsenum],\
                    yerr=firingbinserrList[pulsenum],\
                    color=(pulsenum/float(numpulses),1-pulsenum/float(numpulses),sister_ratio),\
                    marker='+',linestyle='solid',linewidth=2,label='Mitral response')
            ################## Plot the fitted responses
            if pulsenum == 0:
                response = sigmoid_fn(responseconvolve(pulseList[pulsenum+1],kernelR) + bgnd,\
                    sigmoid_max, sigmoid_threshold)
                titlestr = 'Air random pulse'
            elif pulsenum in [1,3]:
                response = sigmoid_fn(responseconvolve(pulseList[pulsenum+1],kernelA) + pulseair + bgnd,\
                    sigmoid_max, sigmoid_threshold)
                titlestr = 'Odor A random pulse'
            elif pulsenum in [2,4]:
                response = sigmoid_fn(responseconvolve(pulseList[pulsenum+1],kernelB) + pulseair + bgnd,\
                    sigmoid_max, sigmoid_threshold)
                titlestr = 'Odor B random pulse'
            else:
                response = responseconvolve(pulseList[pulsenum+1],kernelA) + pulseair
                response += responseconvolve(pulseList[pulsenum+2],kernelB)
                response += bgnd
                response = sigmoid_fn(response, sigmoid_max, sigmoid_threshold)
                titlestr = 'Odor A and odor B random pulse'
            if fill_below:
                plot(pulsetlist,response,linestyle='solid',linewidth=2.0,color=(0,0,0),
                    label=['air','A','B','A','B','A&B'][pulsenum])
                lgd = legend()
                for k in lgd.get_lines():
                    k.set_linewidth(0)
                lgd.draw_frame(False)
                ltext  = lgd.get_texts()
                for l in ltext:
                    l.set_fontsize(20)
                #fig.setp(ltext, fontsize=20)
            else:
                plot(pulsetlist,response,marker='o',linestyle='dashed',
                    linewidth=2.0, label='Linear fit')
                #title('Odor morphs w/ smooth fit',fontsize=24 )
                #title( 'pulsenum = '+str(pulsenum)+', chisquare normalized = '+str(chisq) )
                title(titlestr, fontsize=24)
                lgd = legend()
                ltext  = lgd.get_texts()
                for l in ltext:
                    l.set_fontsize(8)
                ax.set_ylim(0,xpulse)
            ax.set_yticks([xpulse])
            ax.set_yticklabels([str(xpulse)])
            for label in ax.get_yticklabels():
                label.set_fontsize(24)
            ax.set_xlim(0,9)
            ax.set_xticks([])
        ax.set_xticks([0,9])
        ax.set_xticklabels(['0','9'])    
        axes_labels(ax,'time (s)','',adjustpos=False)

    show()
