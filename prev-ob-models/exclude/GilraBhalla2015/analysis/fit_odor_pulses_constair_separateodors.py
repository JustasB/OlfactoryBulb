# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO ROUGHLY FOLLOW PRIYANKA'S ANALYSIS
########## This variant does not use an air kernel, rather a constant air rate.
## USAGE1: python2.6 fit_odor_pulses_constair_separateodors.py <../results/odor_pulses/pulseresult.pickle> <stimseed>
## USAGE2: python2.6 fit_odor_pulses_constair_separateodors.py <../results/odor_pulses/pulseresult.pickle> <stimseed> [NOSHOW] [SAVEFIG]
## USAGE3: python2.6 fit_odor_pulses_constair_separateodors.py <../results/odor_pulses/pulseresult.pickle> <stimseed> match_resp
## USAGE3: python2.6 fit_odor_pulses_constair_separateodors.py <../results/odor_pulses/pulseresult.pickle> <stimseed> mit[0|1] INDIVIDUAL [match_resp]
## USAGE4: python2.6 fit_odor_pulses_constair_separateodors.py <../results/odor_pulses/pulseresult.pickle> <stimseed> TEST
## example: python2.6 fit_odor_pulses_constair_separateodors.py ../results/odor_pulses/odorpulses_netseed770.0_stimseed770.0_SINGLES_JOINTS_PGS_LAT_NOVARMIT_numgloms3_directed0.01.pickle 770.0 SAVEFIG match_resp
## example: python2.6 fit_odor_pulses_constair_separateodors.py ../results/odor_pulses/odorpulses_netseed864.0_stimseed864.0_SINGLES_JOINTS_PGS_NOLAT_NOVARMIT_numgloms3_directed0.01.pickle 864.0 TEST

from scipy import optimize
from scipy import special
import scipy.interpolate
from scipy import signal
import scipy.io
from pylab import *
import pickle
import sys
import math
import copy as cp

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has SETTLETIME, PULSE_RUNTIME
from networkConstants import * # has central_glom
from data_utils import * # has axes_labels()
from analysis_utils import * # has load_morph...(), predict_respresp() and set_min_err()

iterationnum = 1
## Since each odor is fitted for each mitral separately,
## there cannot be a sigmoid function, as the maxfrate and half-point need to be common.
#sigmoid_op = False
SMOOTH_KERNELS = False

## Have overridden the kernel_time=1.5s (for generating ORN frates) in stimuliConstants.py,
## as Priyanka uses 2.0s kernels.
kernel_time = 2.0

numpulses = RANDOM_PULSE_NUMS*2
## rebin the spikes as below, irrespective of previous binning
pulserebindt = 50e-3#fitting_dt # 50ms as Priyanka uses
pulserebins = int(PULSE_RUNTIME/pulserebindt)
#bin_width_time = 2*pulserebindt ## overlapping bins with this bin width
pulsetlist = arange(pulserebindt/2.0,PULSE_RUNTIME,pulserebindt)

## convolutiondt is to bin the flow rate / pressure sensor waveform for pulses:
## namely pulseList for pulse response fits 
## above pulserebindt/convolutiondt should be an integer for rebinning in responseconvolve()
## For respiration convolution, I always use 10ms directly in predict_respresp()
convolutiondt = 50e-3 # small dt to bin the narrow inspiration and pulseList dynamics

## Can have different kerneldt than pulserebindt leading to smoother kernel
## yet taking into account finer features of pulse and response.
## This is not same as smoothing kernel inside the chi-sq calc function.
## In the latter, number of params remain the same: fitting function will 
## adjust the raw values to get jagged kernel even after smoothing.
## cf. pulserebindt for responses above and 
## cf. fitting_dt for generated ORN kernels (in stimuliConstants.py)
kerneldt = 50e-3 # 50ms as Priyanka uses 
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

## how many respiration pulses to fit, starting from the last.
## 1 to NUM_RESPS (in stimuliConstants.py)
numrespsfit = 1
## respdt is to bin respiratory responses
respdt = 50e-3 # different from pulserebindt and kerneldt above

## whether to use chisq or leastsq for fitting
usechisq = False

def responseconvolve(pulses, kernel, \
        convdt=convolutiondt, bindt=pulserebindt, numbins=pulserebins):
    ## Assume pulses are sent in binned at convdt
    ## kernel is sent in at kerneldt
    ## response has to be return-ed at bindt resolution.
    ## bindt/convdt MUST be an integer.
    ## convolution is performed at convdt,
    ## which is less than or equal to bindt and kerneldt.
    ## if kerneldt and convdt are not the same,
    ## interpolate the kernel to convdt.
    ## very IMPORTANT to have the same dt for pulse and kernel
    ## when doing discrete convolution!
    ## Convolution is multiplied by dt to leave result invariant of dt.
    if convdt==kerneldt:
        kernel_interpol = kernel
    else:
        kernel_interpol_fn = scipy.interpolate.interp1d(kerneltlist,kernel)
        kernel_interpol = [ kernel_interpol_fn(t) \
            for t in arange(0.0,kernel_time-kerneldt+1e-3,convdt) ]
    ## don't use mode='same', as it takes only the valid part of the convolution.
    ## always multiply convolution_result by dt
    ## to leave response invariant on decimating stim and kernel simultaneously.
    responsecalc = convolve( pulses, kernel_interpol, mode='full' ) * convdt
    #responsecalc = myconvolve( pulses, kernel_interpol ) * convdt
    decimate_bins = int(bindt/convdt)
    responsecalc = responsecalc[::decimate_bins]
    ## bin0 of the simulated response is at time bindt
    ## so the fitted response should also be shifted by a bindt
    ## NO! if you run in TEST mode, you need to start from 0, else kernel is shifted.
    return array( responsecalc[0:numbins] )
    
def rectify(valarray):
    ## ensure that valarray is a numpy array
    ## else where() returns empty tuple, and below exp() statement won't work.
    valarray = array(valarray)
    ## where value in valarray < 0, set it to zero.
    valarray[where(valarray<0)[0]]=0
    return valarray

def predict_respresp(kernelA,kernelB,numresps=numrespsfit):
    respconvdt = 10e-3
    ## respirationtime starts in value from RESPIRATION-SETTLETIME.
    ## As it is sliced from extratime, so remove that offset.
    ## ORN stimulus is generated by first pre-chopping than convolving:
    ## See my comments in odor_response() in generate_firerates_odors.py
    respirationtime = gen_frate_odors.extratime[gen_frate_odors.startidx:-1] \
                    - (RESPIRATION-SETTLETIME)
    respPulses = gen_frate_odors.respirationPulses[gen_frate_odors.startidx:-1]
    startrespbin = int( (SETTLETIME+(NUM_RESPS-numresps)*RESPIRATION) / respdt )
    ## decimate respPulses to respconvdt
    bins = int(respconvdt/FIRINGFILLDT)
    respPulses = respPulses[::bins]
    ## decimate respirationtime to respdt
    bins = int(respdt/FIRINGFILLDT)
    respirationtime = respirationtime[::bins]
    numbins = len(respirationtime)
    morph_responseA = responseconvolve(respPulses,kernelA,respconvdt,respdt,numbins)
    morph_responseB = responseconvolve(respPulses,kernelB,respconvdt,respdt,numbins)

    ## report midpoint of bin as the time point for the bin value
    ## returned first list starts in value from
    ## SETTLETIME + (NUM_RESPS-numresps)*RESPIRATION + respdt/2.0
    return respirationtime[startrespbin:]+respdt/2.0, \
        morph_responseA[startrespbin:],morph_responseB[startrespbin:]

def SGfilter2(kernel):
    ## savitsky-golay filter maintains high frequency signal, unlike moving average filter
    ## window_size must be odd, savitsky_golay copied from scipy cookbook online
    window_size = 11 # same as Priyanka
    kernel = savitzky_golay(kernel, window_size, order=4) # same as Priyanka
    return kernel

def chisqfunc(params, firingbinsmeanList, firingbinserrList, bgnd, \
        pulserebinList, ret_array=True):
    kernel = params[0:kernel_size]
    ## smoothing kernels here does not help, as number of params remain same
    ## and fitting routine can always modify raw values to nullify effect of smoothing.
    if SMOOTH_KERNELS:
        kernel = SGfilter2(kernel)
    starti = int(PULSE_START/pulserebindt)
    chisqarray = []    
    for pulsenum,pulse in enumerate(pulserebinList):
        ## These below should be numpy arrays,
        ## else + acts like append in python (not element-wise addition)!!!
        responsecalc = \
            rectify(responseconvolve(pulse,kernel) + bgnd)
        ## Priyanka minimizes least sq unlike Adil who minimizes chi-sq
        ## With chi-sq, there is 'what min err to choose' issue.
        ## If usechisq, then do least squares, else chi-sq
        if not usechisq:
            chisqarray.extend( [ (val-firingbinsmeanList[pulsenum][i])\
                for i,val in enumerate(responsecalc[starti:],start=starti) ] )
        else:
            chisqarray.extend( [ (val-firingbinsmeanList[pulsenum][i])/firingbinserrList[pulsenum][i]\
                for i,val in enumerate(responsecalc[starti:],start=starti) ] )
        
    global iternum
    iternum += 1
    ## normalize chi-sq to the number of dof
    num_dof = float(firingbinsmeanList[0][starti:].size*len(pulserebinList) - params.size)
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

def test_fitpulses(filename,stimseed):
    ########## load in the stimuli
    pulseList,pulseStepsList,ORNfrateList,ORNrespList,genORNkernels \
        = read_pulsestimuli(stimseed)
    ## decimate pulseList by taking every convolutiondt/FIRINGFILLDT bins
    ## and decimate ORNkernels by taking every kerneldt/FIRINGFILLDT=50/1=50 bins
    pulserebinList,linORNkernelA,linORNkernelB = \
            decimate(pulseList,convolutiondt,genORNkernels,kerneldt,kernel_size)
    ## ORN kernelA/B is initial kernelA/B
    param_kernels = array([linORNkernelA,linORNkernelB])/5.0
    ORNkernels = cp.deepcopy(param_kernels) # deepcopy else python just copies the reference

    #figure()
    #plot(gen_frate_odors.extratime,gen_frate_odors.respirationPulses)
    #figure()
    #plot(gen_frate_odors.respirationPulses[::int(convolutiondt/FIRINGFILLDT)])

    #figure()
    #for respnum,respresp in enumerate(ORNrespList[:-1]):
        #plot(respresp-ORNrespList[-1],label=str(respnum))
    #legend()
    
    fig = figure()
    frate_air = ORNfrateList[0] # rectangle pulse
    ## odor A has pulsenums 2 and 4, odor B has 3 and 5
    for odornum,mitodor_pulsenums in enumerate(array([[2,4],[3,5]])):
        frates = []
        for pulsenum in mitodor_pulsenums:
            frate = ORNfrateList[pulsenum] + frate_air + ORN_BGND_MEAN
            decimatefactor = int(pulserebindt/FIRINGFILLDT)
            ## decimate the frate, but from middle of the bins
            frate = array(frate[decimatefactor/2::decimatefactor])
            frates.append(frate)

        bgnd = mean( [response[int((SETTLETIME+PULSE_AIR_BUFFER*3.0/4.0)/pulserebindt):\
                int((SETTLETIME+PULSE_AIR_BUFFER)/pulserebindt)] \
            for response in frates] )
        print "bgnd =",bgnd

        iternum = 0
        params = optimize.leastsq( chisqfunc, param_kernels[odornum],\
            args=(frates, None, bgnd, pulserebinList[mitodor_pulsenums], True),\
            full_output=1, maxfev=100000, ftol=1e-100, xtol=1e-100 )
        print params[3] # print the status message
        param_kernels[odornum] = params[0] # take only fitted params; leastsq returns a tuple with errmsg, etc.
        ### Calculate sum of squares of the chisqarray
        chisqarraysq = [i**2 for i in chisqfunc(param_kernels[odornum],\
            frates, None, bgnd, pulserebinList[mitodor_pulsenums], True)]
        #chisqnumpy = array(chisqarraysq)
        #pos = where(chisqnumpy>100)[0]
        #lenf = len(firingbinsmeanList[0])
        #print "The [(pulsenum,index),...] where chisq>100 are",[ (i/lenf,i%lenf) for i in pos]
        chisq = reduce(lambda x, y: x+y, chisqarraysq)
        if usechisq: chisqstr = 'chisq'
        else: chisqstr = 'leastsq' 
        print "The normalized",chisqstr,"for odornum",odornum,"=",chisq
        
        ## pulse response
        kernel = param_kernels[odornum]
        response_unrect = responseconvolve( \
            pulserebinList[mitodor_pulsenums[0]],kernel) + bgnd
        response = rectify(response_unrect)
        orig_response = rectify(responseconvolve( \
            pulserebinList[mitodor_pulsenums[0]],ORNkernels[odornum]) + bgnd)
            
        ## respiratory response
        respirationtime,morph_responseA,morph_responseB = \
                predict_respresp(param_kernels[0],param_kernels[1])
        morph_responses = (morph_responseA,morph_responseB)

        fig.add_subplot(2,3,odornum*3+1)
        plot(frates[0],'r')
        plot(response_unrect,'b')
        plot(orig_response,'g')
        fig.add_subplot(2,3,odornum*3+2)
        plot(ORNkernels[odornum],'r')
        plot(kernel,'b')
        fig.add_subplot(2,3,odornum*3+3)
        ## decimate respPulses to respdt
        bins = int(respdt/FIRINGFILLDT)
        ## subtract the resp air responses to get just the odor response
        respresp = ORNrespList[(5,0)[odornum]][::bins] - ORNrespList[6][::bins]
        respresp = respresp[int((SETTLETIME+(NUM_RESPS-numrespsfit)*RESPIRATION)/respdt):\
                int((SETTLETIME+NUM_RESPS*RESPIRATION)/respdt)]
        plot(respirationtime,respresp,'r')
        plot(respirationtime,morph_responses[odornum],'b')
    
    show()

class fit_plot_pulses():
    def __init__(self,filename,stimseed,match_resp,NOSHOW,SAVEFIG=False):
        self.filename = filename
        self.stimseed = stimseed
        self.match_resp = match_resp
        self.NOSHOW = NOSHOW
        self.SAVEFIG = SAVEFIG

    def fit_pulses(self,filename,stimseed,match_resp,NOSHOW,SAVEFIG=False,dirextn='',dirextn4stim=True):
        ######### load in the mitral pulse responses
        #### mitral_responses_list[avgnum][pulsenum][mitralnum][binnum]
        #### mitral_responses_binned_list[avgnum][pulsenum][mitralnum][binnum]
        mitral_responses_list, mitral_responses_binned_list = read_pulsefile(filename)
        ########## load in the stimuli
        if dirextn4stim: stimdir = '../generators/firerates'+dirextn
        else: stimdir = '../generators/firerates'
        pulseList,pulseStepsList,ORNfrateList,_,genORNkernels \
            = read_pulsestimuli(stimseed,stimdir)
        self.pulseList = pulseList

        ##-------------------------- rebin the responses and pulses ------------------------------
        ## rebin sim responses to pulserebindt=50ms, then take mean
        numtrials,mitral_responses_mean,mitral_responses_std = \
                rebin_mean(mitral_responses_list,pulserebindt)
        ## for pulseList, decimate by taking every convolutiondt/FIRINGFILLDT bins
        ## for ORNkernels, decimate by taking every kerneldt/FIRINGFILLDT=50/1=50 bins
        pulserebinList,linORNkernelA,linORNkernelB = \
                decimate(pulseList,convolutiondt,genORNkernels,kerneldt,kernel_size)

        ## Need to boot up the plotting, before the loop over the two fitted_mitral-s (central sisters)
        if not NOSHOW:
            fig = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
            ax_kernelsA = plt.subplot2grid((12,3),(0,2),rowspan=3,frameon=False)
            text(0.1,1.0,'G', fontsize=label_fontsize, transform = ax_kernelsA.transAxes)
            ax_kernelsB = plt.subplot2grid((12,3),(3,2),rowspan=3,frameon=False)
            text(0.1,1.0,'H', fontsize=label_fontsize, transform = ax_kernelsB.transAxes)
            if match_resp:
                ax_respA = plt.subplot2grid((12,3),(6,2),rowspan=3,frameon=False)
                text(0.1,1.0,'I', fontsize=label_fontsize, transform = ax_respA.transAxes)
                ax_respB = plt.subplot2grid((12,3),(9,2),rowspan=3,frameon=False)
                text(0.1,1.0,'J', fontsize=label_fontsize, transform = ax_respB.transAxes)

        ## full fitting data for both mitrals
        mdict_full = []
        fittedkernels = []
        chisqs = []
        fitted_responses_full = []
        fitted_responses_unrect_full = []
        for fitted_mitral in fitted_mitrals:
            ## take the odor responses of the mitral to be fitted
            firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]

            ## The model predicts the individual response not the mean.
            ## Hence below fitting uses standard deviation, not standard error of the mean.
            firingbinserrList = mitral_responses_std[:,fitted_mitral]
            ### My model seems to have constant noise independent of the mean firing rate.
            ### TEST: force the SDs to be 1.2*sqrt(mean) which is what Adil observes:
            #firingbinserrList = 1.2*firingbinsmeanList

            ## If using chisq for fitting, set a minimum error (for zero-valued frate bins with zero error)
            if usechisq:
                ## each mitral has its own min err, so setting individually for each fitted_mitral
                #errmin = 1.0/pulserebindt/sqrt(numtrials) # errmin = 3.33Hz!
                #errmin = 10.0 #Hz
                #firingbinserrList = set_min_err(firingbinserrList,errmin)
                firingbinserrList = set_min_err(firingbinserrList) # default errmin = min non-zero SD

            self.numtrials = numtrials

            ## slightly different from Priyanka -- not taking full air buffer, allowing to settle
            bgnd = mean( [response[int((SETTLETIME+PULSE_AIR_BUFFER*7.0/8.0)/pulserebindt):\
                    int((SETTLETIME+PULSE_AIR_BUFFER)/pulserebindt)] \
                for response in firingbinsmeanList[1:]] ) ## air response pulsenum=0 doesn't have air pedestal
            print "bgnd =",bgnd

            for pulsenum,response in enumerate(firingbinsmeanList[1:]):
                starti = int(PULSE_START/pulserebindt)
                endi = len(response)
                avgfrate = sum(response[starti:])/float(endi-starti)
                print "Average firing rate for odor",['A','B','A','B','A+B'][pulsenum],\
                    'for mitnum',fitted_mitral,'is',avgfrate

            #---------------------------- fitting ------------------------------------------------

            ## Initial values for the parameters
            if firstrun:
                ## ORN kernelA/B is initial kernelA/B
                param_kernels = [linORNkernelA,linORNkernelB]
                ### can also use all zeros as init kernels, works for reasonably positive responses
                #param_kernels = [[0.0]*len(linORNkernelA),[0.0]*len(linORNkernelB)]
                ORNkernels = (linORNkernelA, linORNkernelB)
            else:
                f = open(sys.argv[1]+'_params'+str(fitted_mitral),'r')
                chisqs,param_kernels,discard_firingmeans,\
                        discard_firingsd,discard_fitresponses = pickle.load(f)
                f.close()
                ORNkernels = ()

            ## export the data for use in matlab in Priyanka's fitting function
            mdict={'firingbinsmeanList':firingbinsmeanList,
                'firingbinserrList':firingbinserrList, 'bgnd':bgnd, 'FIRINGFILLDT':FIRINGFILLDT,
                'pulseList':pulserebinList,'pulseStepsList':pulseStepsList,
                'start_kernels':ORNkernels,'pulserebindt':pulserebindt,'pulsetlist':pulsetlist,
                'start_i':int(PULSE_START/pulserebindt)}
            scipy.io.savemat('mit'+str(fitted_mitral)+'.mat', mdict=mdict)
            mdict_full.append(mdict)

            chisqs_mit = []
            ## fit each odor A and B separately, but fitting the corresponding 2 random pulses for each
            for odornum,mitodor_pulsenums in enumerate(array([[1,3],[2,4]])):
                ## comment fitting line below if you just want to plot parameters.
                ## args is a tuple! if only one element write (elem, )

                #########
                ### fmin_powell is faster than fmin (simplex) but slower than leastsq
                ### However it does not enter into the very bad local minima of leastsq
                ### WARNING! But often it is worse at fitting than leastsq()
                #iternum = 0
                #params = optimize.fmin_powell( chisqfunc, param_kernels[odornum], \
                #    args=(firingbinsmeanList[mitodor_pulsenums], firingbinserrList[mitodor_pulsenums],\
                #    bgnd, pulserebinList[mitodor_pulsenums+1], False),\
                #    full_output=1, maxfun=50000, maxiter=50000, ftol = 1 )
                #chisq = params[1]
                #print "The normalized chisq for mitral",fitted_mitral,"after fmin_powell =",chisq
                #param_kernels[odornum] = params[0]

                #########
                ## leastsq() uses a modified version of Levenberg-Marquardt algorithm
                ## it optimizes M equations in N unknowns, where M>=N.
                ## Hence chisqfunc must return M numbers in an array
                ## which must be >= the number of params N in params0, 
                ## else: 'TypeError: Improper input parameters.'
                iternum = 0
                params = optimize.leastsq( chisqfunc, param_kernels[odornum],\
                    args=(firingbinsmeanList[mitodor_pulsenums], firingbinserrList[mitodor_pulsenums],\
                    bgnd, pulserebinList[mitodor_pulsenums+1], True),\
                    full_output=1, maxfev=100000, ftol=1e-20, xtol=1e-100 )
                print params[3] # print the status message
                param_kernels[odornum] = params[0] # take only fitted params; leastsq returns a tuple with errmsg, etc.
                ### Calculate sum of squares of the chisqarray
                chisqarraysq = [i**2 for i in chisqfunc(param_kernels[odornum],\
                    firingbinsmeanList[mitodor_pulsenums], firingbinserrList[mitodor_pulsenums],\
                    bgnd, pulserebinList[mitodor_pulsenums+1], True)]
                #chisqnumpy = array(chisqarraysq)
                #pos = where(chisqnumpy>100)[0]
                #lenf = len(firingbinsmeanList[0])
                #print "The [(pulsenum,index),...] where chisq>100 are",[ (i/lenf,i%lenf) for i in pos]
                chisq = reduce(lambda x, y: x+y, chisqarraysq)
                if usechisq: chisqstr = 'chisq'
                else: chisqstr = 'leastsq' 
                print "The normalized",chisqstr,"for mitral",fitted_mitral,"for odornum",odornum,"=",chisq
                
                ##########
                ### minimize is a wrapper for all the various fitting algorithms!
                ### Unfortunately, it's only in scipy 0.11 version, not on Ubuntu 12.04 (has 0.9).
                #iternum = 0
                #fit_method = "Powell"#"CG"#"Powell"#"Nelder-Mead"
                #params = optimize.minimize( chisqfunc, param_kernels[odornum],\
                #    args=(firingbinsmeanList[mitodor_pulsenums], firingbinserrList[mitodor_pulsenums],\
                #    bgnd, pulserebinList[mitodor_pulsenums+1], False),\
                #    method = fit_method,\
                #    options={'maxfev':100000, 'ftol':1e-100, 'xtol':1e-100, 'maxfun':500000, 'maxiter':500000} )
                #print params.message
                #param_kernels[odornum] = params.x
                #chisq = params.fun # to check how to obtain chisq
                #print "The normalized chisq for mitral",fitted_mitral,"for odornum",odornum,"after leastsq =",chisq
                
                chisqs_mit.append(chisq)

            kernelA,kernelB = param_kernels

            fittedkernels.append( (kernelA,kernelB) )
            chisqs.append(chisqs_mit)

            ## construct rectified and unrectified fitted pulse responses from above fitted kernels
            fitted_responses = []
            fitted_responses_unrect = []
            for pulsenum in range(1,6):
                if pulsenum in [1,3]:
                    response_unrect = responseconvolve(pulserebinList[pulsenum+1],kernelA) + bgnd
                    response = rectify(response_unrect)
                elif pulsenum in [2,4]:
                    response_unrect = responseconvolve(pulserebinList[pulsenum+1],kernelB) + bgnd
                    response = rectify(response_unrect)
                else:
                    response_unrect = responseconvolve(pulserebinList[pulsenum+1],kernelA)
                    response_unrect += responseconvolve(pulserebinList[pulsenum+2],kernelB)
                    response_unrect += bgnd
                    response = rectify(response_unrect)
                fitted_responses_unrect.append(response_unrect)
                fitted_responses.append(response)
            fitted_responses_unrect_full.append(fitted_responses_unrect)
            fitted_responses_full.append(fitted_responses)

            paramsfile = open(filename+'_params'+str(fitted_mitral),'w')
            pickle.dump( \
                (chisqs_mit,(kernelA,kernelB),bgnd,firingbinsmeanList,firingbinserrList,fitted_responses), \
                paramsfile )
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

            ## ----------------------- smooth kernels only for display of prediction & fit --------------------------

            ## irrespective of smooth kernels flag, filter the output kernels for display
            ## IMPORTANT: SMOOTHED KERNEL IS USED FOR BELOW FOR PREDICTION & FIT
            ## ONLY FOR DISPLAY PURPOSES!!!!
            ## In saving _params file above and in fit_odor_pulses_.....py, I use the original kernels.
            smooth_kernelA = SGfilter2(append([0],kernelA))
            smooth_kernelB = SGfilter2(append([0],kernelB))

            ##---------------------------- predict respiring odor response for plotting only ------------------------
            
            if match_resp and not NOSHOW:
                ## Important to first convolve, then pre-truncate the flow rate trace
                ## Else the first resp cycle does not match.
                respirationtime,morph_responseA,morph_responseB = \
                    predict_respresp(kernelA,kernelB)
                ## plot simulated resp responses
                ## Load the mitral responses also
                morph_numavgs,morph_mitral_responses_avg,morph_mitral_responses_std = \
                    load_morph_responses_from_pulsefilename(filename,fitted_mitral,respdt,numresps=numrespsfit)
                morph_color = ['r','b'][fitted_mitral]
                ## Need to subtract the resp air response, not the const 'bgnd'
                airresponse = morph_mitral_responses_avg[6]
                simresponse = morph_mitral_responses_avg[5] - airresponse
                simerr = morph_mitral_responses_std[5]/sqrt(morph_numavgs)
                ax_respA.fill_between(respirationtime, simresponse+simerr, simresponse-simerr,
                    color=morph_color,alpha=0.4,linewidth=0)
                ax_respA.plot(respirationtime,simresponse,color='k',linewidth=linewidth)
                simresponse = morph_mitral_responses_avg[0] - airresponse
                simerr = morph_mitral_responses_std[0]/sqrt(morph_numavgs)
                ax_respB.fill_between(respirationtime, simresponse+simerr, simresponse-simerr,
                    color=morph_color,alpha=0.4,linewidth=0)
                ax_respB.plot(respirationtime,simresponse,color='k',linewidth=linewidth)
                ## plot predicted resp responses
                morph_color = ['r','b'][fitted_mitral]
                ax_respA.plot(respirationtime,morph_responseA,color=morph_color,linewidth=linewidth)
                ax_respB.plot(respirationtime,morph_responseB,color=morph_color,linewidth=linewidth)
                
            ##---------------------------- plot kernels and pulse responses -----------------------------------------

            if not NOSHOW:
                ## fig defined before the mitnum loop
                ############################## plot the kernels

                ## kernel actually starts from kerneldt/2.0, but add a time point a bin earlier,
                ## as while smoothing Savistsky-Golay above, I set kernel[-kerneldt/2.0]=0.0
                ## Ideally I should have set kernel[0]=0, but SG filter cannot handle non-uniform time-series.
                plot_kerneltlist = append([-kerneldt/2.0],kerneltlist+kerneldt/2.0)
                ax_kernelsA.plot(plot_kerneltlist, smooth_kernelA, color=['r','b'][fitted_mitral], marker=',',\
                    linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
                #ax.plot(kerneltlist, linORNkernelA, color=(1,0,1), marker=',',\
                #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
                ax_kernelsB.plot(plot_kerneltlist, smooth_kernelB, color=['r','b'][fitted_mitral], marker=',',\
                    linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
                #ax.plot(kerneltlist, linORNkernelB, color=(0,1,1), marker=',',\
                #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')

                ############################### plot the responses and the fits
                for pulseiter, pulsenum in enumerate([1,2,5]):#range(1,numpulses): # similar to Priyanka, only 3 plots
                    sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
                    ## similar to Priyanka: only 3 plots, skip 2 pulse responses
                    ax = plt.subplot2grid((12,3),(4*pulseiter,fitted_mitral),rowspan=4,frameon=False)
                    panel_label = ['A','B','C','D','E','F'][fitted_mitral*3+pulseiter]
                    ## ax.transAxes ensures relative to axes size, rather than to data units.
                    text(0.1,1.0,panel_label, fontsize=label_fontsize, transform = ax.transAxes)
                    ################# random pulses and ORN responses
                    xpulse = int(max(firingbinsmeanList[pulsenum]+firingbinserrList[pulsenum]/sqrt(numtrials))+1)
                    randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
                    if pulsenum in [1,2,3,4]: # odor A or B
                        if fill_below:
                            if pulsenum in [1,3]: col = 'm'
                            if pulsenum in [2,4]: col = 'c'
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
                                color='m',linewidth=0,alpha=0.4)
                            fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+2],\
                                color='c',linewidth=0,alpha=0.4)
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
                            simresponse+firingbinserrList[pulsenum]/sqrt(numtrials),
                            simresponse-firingbinserrList[pulsenum]/sqrt(numtrials),
                            color=(0.7,0.7,0.7))
                        plot(pulsetlist,simresponse,linewidth=linewidth,color=(0.3,0.3,0.3))
                    else:
                        errorbar(pulsetlist,y=simresponse,\
                            yerr=firingbinserrList[pulsenum]/sqrt(numtrials),\
                            color=(pulsenum/float(numpulses),1-pulsenum/float(numpulses),sister_ratio),\
                            marker='+',linestyle='solid',linewidth=linewidth,label='Mitral response')
                    ################## Plot the fitted responses
                    starti = int(PULSE_START/pulserebindt)
                    if pulsenum in [1,3]: titlestr = 'Odor A random pulse'
                    elif pulsenum in [2,4]: titlestr = 'Odor B random pulse'
                    else: titlestr = 'Odor A and odor B random pulse'
                    response_unrect = fitted_responses_unrect[pulsenum-1]
                    if fill_below:
                        plot(pulsetlist[starti:],response_unrect[starti:],
                            linestyle='solid',linewidth=linewidth,color=['r','b'][fitted_mitral],
                            label=['air','A','B','A','B','A&B'][pulsenum])
                        #lgd = legend()
                        #for k in lgd.get_lines():
                        #    k.set_linewidth(0)
                        #lgd.draw_frame(False)
                        #ltext  = lgd.get_texts()
                        #for l in ltext:
                        #    l.set_fontsize(label_fontsize)
                        ##fig.setp(ltext, fontsize=20)
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
                    ax.set_yticks([0,xpulse])
                    ax.set_yticklabels(['0',str(xpulse)])
                    ax.set_xlim(0,9)
                    ax.set_xticks([])
                    axes_labels(ax,'','',adjustpos=False) # sets default label_fontsize
                    ax.get_xaxis().set_ticks_position('none')
                    ax.get_yaxis().set_ticks_position('left')
                    xmin, xmax = ax.get_xaxis().get_view_interval()
                    ymin, ymax = ax.get_yaxis().get_view_interval()
                    ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
                    if fitted_mitral==0 and pulseiter==1:
                        axes_labels(ax,'','firing rate (Hz)',adjustpos=False)
                ax.get_xaxis().set_ticks_position('bottom')
                ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=axes_linewidth))
                ax.set_xticks([0,9])
                ax.set_xticklabels(['0','9'])    
                axes_labels(ax,'time (s)','',adjustpos=False)

        if not NOSHOW:
            ## Below text is wrt to last axes drawn. Tried setting it to bigAxes, but that doesn't work.
            ## transform = ax.transAxes sets the test position as axes units and not data units.
            text(1.04,3.5,'arb units', fontsize=label_fontsize, rotation='vertical', transform = ax.transAxes)
            #text(-0.25,1.0,'firing rate (Hz)', fontsize=label_fontsize, rotation='vertical', transform = ax.transAxes)

            ## beautify panels for kernels
            ax_kernelsA.set_yticks([0])
            ax_kernelsB.set_yticks([0])
            ax_kernelsA.set_xlim(0,kernel_time)
            ax_kernelsB.set_xlim(0,kernel_time)
            #biglegend(ax=ax_kernelsA)
            #biglegend(ax=ax_kernelsB)
            ax_kernelsB.set_xticks([0,kernel_time])
            #ax_kernelsB.set_xticklabels(['0','2'])
            axes_off(ax_kernelsA,x=True,y=False)
            ## turn on/off the side axes (after turning off the frame above):
            ## http://www.shocksolution.com/2011/08/removing-an-axis-or-both-axes-from-a-matplotlib-plot/
            ax_kernelsA.get_xaxis().set_ticks_position('none')
            ax_kernelsA.get_yaxis().set_ticks_position('left')
            ax_kernelsB.get_xaxis().set_ticks_position('bottom')
            ax_kernelsB.get_yaxis().set_ticks_position('left')
            ## need to add the axes lines that I want, after deleting full frame.
            xmin, xmax = ax_kernelsA.get_xaxis().get_view_interval()
            ymin, ymax = ax_kernelsA.get_yaxis().get_view_interval()
            ax_kernelsA.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
            xmin, xmax = ax_kernelsB.get_xaxis().get_view_interval()
            ymin, ymax = ax_kernelsB.get_yaxis().get_view_interval()
            ax_kernelsB.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
            ax_kernelsB.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=axes_linewidth))
            ## axes_labels() also sets label sizes to default label_fontsize
            axes_labels(ax_kernelsA,'','',adjustpos=False)
            axes_labels(ax_kernelsB,'','',adjustpos=False)
            
            if match_resp:
                ## beautify plot respA
                #ax_respA.set_xlim(SETTLETIME+(NUM_RESPS-numrespsfit)*RESPIRATION, \
                #        SETTLETIME+NUM_RESPS*RESPIRATION)
                ax_respA.get_xaxis().set_ticks_position('none')
                ax_respA.get_yaxis().set_ticks_position('left')
                ## need to add the axes lines that I want, after deleting full frame.
                xmin, xmax = ax_respA.get_xaxis().get_view_interval()
                ymin, ymax = ax_respA.get_yaxis().get_view_interval()
                ax_respA.set_xticks([])
                ymin,ymax = int(floor(ymin)),int(ceil(ymax))
                ax_respA.set_yticks([0,ymax])
                ax_respA.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
                axes_labels(ax_respA,'','') # sets default fontsize too
                ## beautify plot respB
                #ax_respB.set_xlim(SETTLETIME+(NUM_RESPS-numrespsfit)*RESPIRATION, \
                #        SETTLETIME+NUM_RESPS*RESPIRATION)
                ax_respB.get_xaxis().set_ticks_position('bottom')
                ax_respB.get_yaxis().set_ticks_position('left')
                ## need to add the axes lines that I want, after deleting full frame.
                xmin, xmax = ax_respB.get_xaxis().get_view_interval()
                ymin, ymax = ax_respB.get_yaxis().get_view_interval()
                ymin,ymax = int(floor(ymin)),int(ceil(ymax))
                ax_respB.set_xticks([xmin,xmax])
                ax_respB.set_yticks([0,ymax])
                ax_respB.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
                ax_respB.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=axes_linewidth))
                axes_labels(ax_respB,'time(s)','') # sets default fontsize too

            fig.tight_layout()
            subplots_adjust(top=0.94)
            if SAVEFIG:
                fig.savefig('../figures/linearity_example_'+str(stimseed)+'.svg',dpi=fig.dpi)
                fig.savefig('../figures/linearity_example_'+str(stimseed)+'.png',dpi=fig.dpi)

        self.fittedkernels = fittedkernels
        self.chisqs = chisqs
        self.fitted_responses_full = fitted_responses_full
        self.fitted_responses_unrect_full = fitted_responses_unrect_full
        self.mdict_full = mdict_full
        return fittedkernels,chisqs,fitted_responses_full,mdict_full

    ## OBSOLETE PAPER FIGURE 4: linearity example, plot for one mitnum (fitted_mitral)...
    def plot_pulses_mitnum_papersubfigure(self,fitted_mitral):
        ## load the simulated and fitted variables
        match_resp = self.match_resp
        (kernelA,kernelB) = self.fittedkernels[fitted_mitral]
        fitted_responses_unrect = self.fitted_responses_unrect_full[fitted_mitral]
        ## mdict has below variables, insert them into the local namespace
        ##mdict={'firingbinsmeanList':firingbinsmeanList,
        ##    'firingbinserrList':firingbinserrList, 'bgnd':bgnd, 'FIRINGFILLDT':FIRINGFILLDT,
        ##    'pulseList':pulserebinList,'pulseStepsList':pulseStepsList,
        ##    'start_kernels':ORNkernels,'pulserebindt':pulserebindt,'pulsetlist':pulsetlist,
        ##    'start_i':int(PULSE_START/pulserebindt)}
        ## mdict_full has mdict-s for mit0 and mit1.
        mdict = self.mdict_full[fitted_mitral]
        firingbinsmeanList = mdict['firingbinsmeanList']
        firingbinserrList = mdict['firingbinserrList']
        pulseStepsList = mdict['pulseStepsList']
        pulseList = self.pulseList # mdict['pulseList'] is actually pulserebinList
        ## careful, updating locals() doesn't insert variables into the namespace! see:
        ## http://stackoverflow.com/a/8028772
        #locals().update(self.mdict_full[fitted_mitral])
        numtrials = self.numtrials
        
        fig = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
        ax_kernelsA = plt.subplot2grid((3,6),(0,0),rowspan=1,colspan=2,frameon=False)
        #text(0.1,1.0,'f', fontweight='bold', fontsize=label_fontsize, transform = ax_kernelsA.transAxes)
        ax_kernelsB = plt.subplot2grid((3,6),(1,0),rowspan=1,colspan=2,frameon=False)
        #text(0.1,1.0,'g', fontweight='bold', fontsize=label_fontsize, transform = ax_kernelsB.transAxes)
        if match_resp:
            ax_resp = plt.subplot2grid((3,6),(2,0),rowspan=1,colspan=2,frameon=False)
            #text(0.1,1.0,'k', fontweight='bold', fontsize=label_fontsize, transform = ax_resp.transAxes)

        ## ----------------------- smooth kernels only for display of prediction & fit --------------------------

        ## irrespective of smooth kernels flag, filter the output kernels for display
        ## IMPORTANT: SMOOTHED KERNEL IS USED FOR BELOW FOR PREDICTION & FIT
        ## ONLY FOR DISPLAY PURPOSES!!!!
        ## In saving _params file above and in fit_odor_pulses_.....py, I use the original kernels.
        smooth_kernelA = SGfilter2(append([0],kernelA))
        smooth_kernelB = SGfilter2(append([0],kernelB))
        
        if match_resp:
            ## Important to first convolve, then pre-truncate the flow rate trace
            ## Else the first resp cycle does not match.
            respirationtime,morph_responseA,morph_responseB = \
                predict_respresp(kernelA,kernelB)
            ## plot simulated resp responses
            ## Load the mitral responses also
            morph_numavgs,morph_mitral_responses_avg,morph_mitral_responses_std = \
                load_morph_responses_from_pulsefilename(\
                self.filename,fitted_mitral,respdt,numresps=numrespsfit)
            ## Need to subtract the resp air response, not the const 'bgnd'
            airresponse = morph_mitral_responses_avg[6]
            simresponseA = morph_mitral_responses_avg[5] - airresponse
            simerrA = morph_mitral_responses_std[5]/sqrt(morph_numavgs)
            ax_resp.fill_between(respirationtime, simresponseA+simerrA, simresponseA-simerrA,
                color='r',alpha=0.4,linewidth=0)
            ax_resp.plot(respirationtime,simresponseA,color='m',linewidth=linewidth)
            simresponseB = morph_mitral_responses_avg[0] - airresponse
            simerrB = morph_mitral_responses_std[0]/sqrt(morph_numavgs)
            ax_resp.fill_between(respirationtime, simresponseB+simerrB, simresponseB-simerrB,
                color='b',alpha=0.4,linewidth=0)
            ax_resp.plot(respirationtime,simresponseB,color='c',linewidth=linewidth)
            ## plot predicted resp responses
            ax_resp.plot(respirationtime,morph_responseA,color='r',linewidth=linewidth)
            ax_resp.plot(respirationtime,morph_responseB,color='b',linewidth=linewidth)
            add_scalebar(ax_resp,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizex=0.1,labelx='0.1 s',bbox_to_anchor=[0.8,-0.1],bbox_transform=ax_resp.transAxes)
            add_scalebar(ax_resp,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizey=15,labely='15 Hz',bbox_to_anchor=[0.7,0.85],bbox_transform=ax_resp.transAxes)

            ## separate figure for resp plots of fitted mitral: odor A and odor B
            fig_resp = figure(figsize=(columnwidth/3.,linfig_height/3.),dpi=300,facecolor='w') # 'none' is transparent
            ax_resp2 = fig_resp.add_subplot(111)
            ax_resp2.fill_between(respirationtime, simresponseA+simerrA, simresponseA-simerrA,
                color='r',alpha=0.4,linewidth=0)
            ax_resp2.plot(respirationtime,simresponseA,color='m',linewidth=linewidth)
            ax_resp2.fill_between(respirationtime, simresponseB+simerrB, simresponseB-simerrB,
                color='b',alpha=0.4,linewidth=0)
            ax_resp2.plot(respirationtime,simresponseB,color='c',linewidth=linewidth)
            ## plot predicted resp responses
            ax_resp2.plot(respirationtime,morph_responseA,color='r',linewidth=linewidth)
            ax_resp2.plot(respirationtime,morph_responseB,color='b',linewidth=linewidth)
            add_scalebar(ax_resp2,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizex=0.1,labelx='0.1 s',\
                bbox_to_anchor=[0.8,-0.1],bbox_transform=ax_resp.transAxes)
            xmin,xmax,ymin,ymax = \
                beautify_plot(ax_resp2,x0min=False,y0min=False,\
                xticksposn='none',yticksposn='left',drawxaxis=False)
            ax_resp2.set_yticks([0,ymax])
            axes_labels(ax_resp2,"","Hz",ypad=-4)
            fig_resp.tight_layout()
            subplots_adjust(bottom=0.15)
            if self.SAVEFIG:
                fig_resp.savefig('../figures/resp_fit_'+str(self.stimseed)+\
                    '_mit'+str(fitted_mitral)+'.svg',dpi=fig.dpi)
                fig_resp.savefig('../figures/resp_fit_'+str(self.stimseed)+\
                    '_mit'+str(fitted_mitral)+'.png',dpi=fig.dpi)

        ##---------------------------- plot kernels and pulse responses -----------------------------------------

        ############################## plot the kernels

        ## kernel actually starts from kerneldt/2.0, but add a time point a bin earlier,
        ## as while smoothing Savistsky-Golay above, I set kernel[-kerneldt/2.0]=0.0
        ## Ideally I should have set kernel[0]=0, but SG filter cannot handle non-uniform time-series.
        plot_kerneltlist = append([-kerneldt/2.0],kerneltlist+kerneldt/2.0)
        ax_kernelsA.plot(plot_kerneltlist, smooth_kernelA, color='r', marker=',',\
            linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
        #ax.plot(kerneltlist, linORNkernelA, color=(1,0,1), marker=',',\
        #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
        ax_kernelsB.plot(plot_kerneltlist, smooth_kernelB, color='b', marker=',',\
            linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
        #ax.plot(kerneltlist, linORNkernelB, color=(0,1,1), marker=',',\
        #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
        kmin = min(min(smooth_kernelA),min(smooth_kernelB),0) # 0 should be shown, just in case min>0.
        kmax = max(max(smooth_kernelA),max(smooth_kernelB))
        kextra = (kmax-kmin)*0.02 # to avoid clipping in y
        ax_kernelsA.set_ylim(kmin-kextra,kmax+kextra)
        ax_kernelsB.set_ylim(kmin-kextra,kmax+kextra)

        ############################### plot the responses and the fits
        figure(fig.number) # set current figure, as subplot2grid() acts on current figure.
        for pulseiter, pulsenum in enumerate([1,2,5]):#range(1,numpulses): # similar to Priyanka, only 3 plots
            sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            ## similar to Priyanka: only 3 plots, skip 2 pulse responses
            ax = plt.subplot2grid((3,6),(pulseiter,2),rowspan=1,colspan=4,frameon=False)
            panel_label = ['h','i','j'][pulseiter]
            ## ax.transAxes ensures relative to axes size, rather than to data units.
            #text(0.05, 1.0, panel_label, fontweight='bold', fontsize=label_fontsize, transform = ax.transAxes)
            ################# random pulses and ORN responses
            xpulse = int(max(firingbinsmeanList[pulsenum]+firingbinserrList[pulsenum]/sqrt(numtrials))+1)
            randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
            if pulsenum in [1,2,3,4]: # odor A or B
                if fill_below:
                    if pulsenum in [1,3]: col = 'm'
                    if pulsenum in [2,4]: col = 'c'
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
                        color='m',linewidth=0,alpha=0.4)
                    fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+2],\
                        color='c',linewidth=0,alpha=0.4)
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
                    simresponse+firingbinserrList[pulsenum]/sqrt(numtrials),
                    simresponse-firingbinserrList[pulsenum]/sqrt(numtrials),
                    color=(0.7,0.7,0.7))
                plot(pulsetlist,simresponse,linewidth=linewidth,color=(0.3,0.3,0.3))
            else:
                errorbar(pulsetlist,y=simresponse,\
                    yerr=firingbinserrList[pulsenum]/sqrt(numtrials),\
                    color=(pulsenum/float(numpulses),1-pulsenum/float(numpulses),sister_ratio),\
                    marker='+',linestyle='solid',linewidth=linewidth,label='Mitral response')
            ################## Plot the fitted responses
            starti = int(PULSE_START/pulserebindt)
            if pulsenum in [1,3]: titlestr = 'Odor A random pulse'
            elif pulsenum in [2,4]: titlestr = 'Odor B random pulse'
            else: titlestr = 'Odor A and odor B random pulse'
            response_unrect = fitted_responses_unrect[pulsenum-1]
            if fill_below:
                plot(pulsetlist[starti:],response_unrect[starti:],
                    linestyle='solid',linewidth=linewidth,color=['r','b',(1,0,1)][pulseiter],
                    label=['air','A','B','A','B','A&B'][pulsenum])
                #lgd = legend()
                #for k in lgd.get_lines():
                #    k.set_linewidth(0)
                #lgd.draw_frame(False)
                #ltext  = lgd.get_texts()
                #for l in ltext:
                #    l.set_fontsize(label_fontsize)
                ##fig.setp(ltext, fontsize=20)
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
            ax.set_yticks([0,xpulse])
            ax.set_yticklabels(['0',str(xpulse)])
            ax.set_xlim(0,9)
            ax.set_xticks([])
            ax.get_xaxis().set_ticks_position('none')
            ax.get_yaxis().set_ticks_position('left')
            xmin, xmax = ax.get_xaxis().get_view_interval()
            ymin, ymax = ax.get_yaxis().get_view_interval()
            ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
            if pulseiter==1:
                axes_labels(ax,'','firing rate (Hz)',adjustpos=False)
                add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
                    sizex=1,labelx='1 s',bbox_to_anchor=[0.6,-0.4],bbox_transform=ax.transAxes)
            else:
                axes_labels(ax,'','',adjustpos=False) # sets default label_fontsize
        ax.get_xaxis().set_ticks_position('none')
        axes_labels(ax,'','',adjustpos=False)

        ## beautify panels for kernels
        ax_kernelsA.set_yticks([0])
        ax_kernelsB.set_yticks([0])
        ax_kernelsA.set_xlim(0,kernel_time)
        ax_kernelsB.set_xlim(0,kernel_time)
        ax_kernelsA.set_xticks([])
        ax_kernelsB.set_xticks([])
        axes_off(ax_kernelsA,x=True,y=False)
        ## turn on/off the side axes (after turning off the frame above):
        ## http://www.shocksolution.com/2011/08/removing-an-axis-or-both-axes-from-a-matplotlib-plot/
        ax_kernelsA.get_xaxis().set_ticks_position('none')
        ax_kernelsA.get_yaxis().set_ticks_position('left')
        ax_kernelsB.get_xaxis().set_ticks_position('none')
        ax_kernelsB.get_yaxis().set_ticks_position('left')
        add_scalebar(ax_kernelsA,matchx=False,matchy=False,hidex=True,hidey=False,\
            sizex=0.5,labelx='0.5 s',sizey=100,labely='arb',\
            bbox_to_anchor=[0.8,-0.6],bbox_transform=ax_kernelsA.transAxes)
        ## need to add the axes lines that I want, after deleting full frame.
        xmin, xmax = ax_kernelsA.get_xaxis().get_view_interval()
        ymin, ymax = ax_kernelsA.get_yaxis().get_view_interval()
        ax_kernelsA.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
        xmin, xmax = ax_kernelsB.get_xaxis().get_view_interval()
        ymin, ymax = ax_kernelsB.get_yaxis().get_view_interval()
        ax_kernelsB.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
        ## axes_labels() also sets label sizes to default label_fontsize
        axes_labels(ax_kernelsA,'','',adjustpos=False)
        axes_labels(ax_kernelsB,'','',adjustpos=False)
            
        if match_resp:
            ## beautify plot resp
            #ax_resp.set_xlim(SETTLETIME+(NUM_RESPS-numrespsfit)*RESPIRATION, \
            #        SETTLETIME+NUM_RESPS*RESPIRATION)
            ax_resp.get_xaxis().set_ticks_position('none')
            ax_resp.get_yaxis().set_ticks_position('left')
            ## need to add the axes lines that I want, after deleting full frame.
            xmin, xmax = ax_resp.get_xaxis().get_view_interval()
            ymin, ymax = ax_resp.get_yaxis().get_view_interval()
            ax_resp.set_xticks([])
            ymin,ymax = int(floor(ymin)),int(ceil(ymax))
            ax_resp.set_yticks([0,ymax])
            ax_resp.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
            axes_labels(ax_resp,'','') # sets default fontsize too

        fig.tight_layout()
        ## do not set clip off, since kernels begin from -25ms due to binning issues.
        #fig_clip_off(fig)
        subplots_adjust(top=0.94)
        if self.SAVEFIG:
            fig.savefig('../figures/linearity_example_'+str(self.stimseed)+\
                '_mit'+str(fitted_mitral)+'.svg',dpi=fig.dpi)
            fig.savefig('../figures/linearity_example_'+str(self.stimseed)+\
                '_mit'+str(fitted_mitral)+'.png',dpi=fig.dpi)

    ## PAPER FIGURE 5 & SUPPL FIG 1: linearity example, plot for one mitnum (fitted_mitral)...
    def plot_pulses_mitnum_papersubfigure_v2(self,fitted_mitral):
        ## load the simulated and fitted variables
        match_resp = self.match_resp
        (kernelA,kernelB) = self.fittedkernels[fitted_mitral]
        fitted_responses_unrect = self.fitted_responses_unrect_full[fitted_mitral]
        ## mdict has below variables, insert them into the local namespace
        ##mdict={'firingbinsmeanList':firingbinsmeanList,
        ##    'firingbinserrList':firingbinserrList, 'bgnd':bgnd, 'FIRINGFILLDT':FIRINGFILLDT,
        ##    'pulseList':pulserebinList,'pulseStepsList':pulseStepsList,
        ##    'start_kernels':ORNkernels,'pulserebindt':pulserebindt,'pulsetlist':pulsetlist,
        ##    'start_i':int(PULSE_START/pulserebindt)}
        ## mdict_full has mdict-s for mit0 and mit1.
        mdict = self.mdict_full[fitted_mitral]
        firingbinsmeanList = mdict['firingbinsmeanList']
        firingbinserrList = mdict['firingbinserrList']
        pulseStepsList = mdict['pulseStepsList']
        pulseList = self.pulseList # mdict['pulseList'] is actually pulserebinList
        ## careful, updating locals() doesn't insert variables into the namespace! see:
        ## http://stackoverflow.com/a/8028772
        #locals().update(self.mdict_full[fitted_mitral])
        numtrials = self.numtrials
        
        fig = figure(figsize=(columnwidth,linfig_height/2.),dpi=300,facecolor='w') # prediction figure
        ax_kernelsA = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1,frameon=False)
        ax_kernelsB = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1,frameon=False)
        fig1 = figure(figsize=(columnwidth,linfig_height/2.),dpi=300,facecolor='w') # fits figure
        ax_kernelsA1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1,frameon=False)
        ax_kernelsB1 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1,frameon=False)

        ## ----------------------- smooth kernels only for display of prediction & fit --------------------------

        ## irrespective of smooth kernels flag, filter the output kernels for display
        ## IMPORTANT: SMOOTHED KERNEL IS USED FOR BELOW FOR PREDICTION & FIT
        ## ONLY FOR DISPLAY PURPOSES!!!!
        ## In saving _params file above and in fit_odor_pulses_.....py, I use the original kernels.
        smooth_kernelA = SGfilter2(append([0],kernelA))
        smooth_kernelB = SGfilter2(append([0],kernelB))
        
        if match_resp:
            ## Important to first convolve, then pre-truncate the flow rate trace
            ## Else the first resp cycle does not match.
            respirationtime,morph_responseA,morph_responseB = \
                predict_respresp(kernelA,kernelB)
            ## plot simulated resp responses
            ## Load the mitral responses also
            morph_numavgs,morph_mitral_responses_avg,morph_mitral_responses_std = \
                load_morph_responses_from_pulsefilename(\
                self.filename,fitted_mitral,respdt,numresps=numrespsfit)
            ## Need to subtract the resp air response, not the const 'bgnd'
            airresponse = morph_mitral_responses_avg[6]
            simresponseA = morph_mitral_responses_avg[5] - airresponse
            simerrA = morph_mitral_responses_std[5]/sqrt(morph_numavgs)
            simresponseB = morph_mitral_responses_avg[0] - airresponse
            simerrB = morph_mitral_responses_std[0]/sqrt(morph_numavgs)

            ## separate figure for resp plots of fitted mitral: odor A and odor B
            fig_resp = figure(figsize=(columnwidth/3.,linfig_height/2.),dpi=300,facecolor='w') # 'none' is transparent
            ax_resp2 = fig_resp.add_subplot(2,1,1)
            ax_resp3 = fig_resp.add_subplot(2,1,2)
            ax_resp2.fill_between(respirationtime, simresponseA+simerrA, simresponseA-simerrA,
                color='r',alpha=0.4,linewidth=0)
            ax_resp2.plot(respirationtime,simresponseA,color='m',linewidth=linewidth)
            ax_resp3.fill_between(respirationtime, simresponseB+simerrB, simresponseB-simerrB,
                color='b',alpha=0.4,linewidth=0)
            ax_resp3.plot(respirationtime,simresponseB,color='c',linewidth=linewidth)
            ## plot predicted resp responses
            ax_resp2.plot(respirationtime,morph_responseA,color='r',linewidth=linewidth)
            ax_resp3.plot(respirationtime,morph_responseB,color='b',linewidth=linewidth)
            add_scalebar(ax_resp2,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizex=0.1,labelx='0.1 s',\
                bbox_to_anchor=[0.8,-0.25],bbox_transform=ax_resp2.transAxes)
            for i,ax_resp_iter in enumerate([ax_resp2,ax_resp3]):
                xmin,xmax,ymin,ymax = \
                    beautify_plot(ax_resp_iter,x0min=False,y0min=False,\
                    xticksposn='none',yticksposn='left',drawxaxis=False)
                if i==0: ax_resp_iter.set_yticks([ymin-5,0,ymax])
                else: ax_resp_iter.set_yticks([ymin,0,ymax])
            ax_resp3.set_xticks([])
            axes_labels(ax_resp2,"","Hz",ypad=-4)
            ax_resp2.yaxis.set_label_coords(-0.4,-0.2)
            fig_resp.tight_layout()
            fig_resp.subplots_adjust(left=0.35,hspace=0.3)
            if self.SAVEFIG:
                fig_resp.savefig('../figures/resp_fit_'+str(self.stimseed)+\
                    '_mit'+str(fitted_mitral)+'.svg',dpi=fig.dpi)
                fig_resp.savefig('../figures/resp_fit_'+str(self.stimseed)+\
                    '_mit'+str(fitted_mitral)+'.png',dpi=fig.dpi)

        ##---------------------------- plot kernels and pulse responses -----------------------------------------

        ############################## plot the kernels

        ## kernel actually starts from kerneldt/2.0, but add a time point a bin earlier,
        ## as while smoothing Savistsky-Golay above, I set kernel[-kerneldt/2.0]=0.0
        ## Ideally I should have set kernel[0]=0, but SG filter cannot handle non-uniform time-series.
        plot_kerneltlist = append([-kerneldt/2.0],kerneltlist+kerneldt/2.0)
        for ax_kernel in [ax_kernelsA,ax_kernelsA1]:
            ax_kernel.plot(plot_kerneltlist, smooth_kernelA, color='r', marker=',',\
                linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
            #ax.plot(kerneltlist, linORNkernelA, color=(1,0,1), marker=',',\
            #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
        for ax_kernel in [ax_kernelsB,ax_kernelsB1]:
            ax_kernel.plot(plot_kerneltlist, smooth_kernelB, color='b', marker=',',\
                linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))
            #ax.plot(kerneltlist, linORNkernelB, color=(0,1,1), marker=',',\
            #    linestyle='solid',linewidth=linewidth,label='5 * ORN kernel')
        kmin = min(min(smooth_kernelA),min(smooth_kernelB),0) # 0 should be shown, just in case min>0.
        kmax = max(max(smooth_kernelA),max(smooth_kernelB))
        kextra = (kmax-kmin)*0.02 # to avoid clipping in y
        for ax_kernel in [ax_kernelsA,ax_kernelsA1,ax_kernelsB,ax_kernelsB1]:
            ax_kernel.set_ylim(kmin-kextra,kmax+kextra)

        ############################### plot the responses and the fits
        ## similar to Priyanka: only 3 plots, skip 2 pulse responses
        for pulseiter, pulsenum in enumerate([1,2,5]):#range(1,numpulses): # similar to Priyanka, only 3 plots
            if pulsenum in [1,2]:
                figure(fig1.number) # set fitting figure, as subplot2grid() acts on current figure.
                ax = plt.subplot2grid((2,3),(pulseiter,1),rowspan=1,colspan=2,frameon=False)
            else:
                figure(fig.number) # set prediction figure, as subplot2grid() acts on current figure.
                ax = plt.subplot2grid((2,3),(0,1),rowspan=2,colspan=2,frameon=False)
            sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            panel_label = ['h','i','j'][pulseiter]
            ## ax.transAxes ensures relative to axes size, rather than to data units.
            #text(0.05, 1.0, panel_label, fontweight='bold', fontsize=label_fontsize, transform = ax.transAxes)
            ################# random pulses and ORN responses
            xpulse = int(round(max(firingbinsmeanList[pulsenum]+firingbinserrList[pulsenum]/sqrt(numtrials))+1))
            randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
            if pulsenum in [1,2,3,4]: # odor A or B
                if fill_below:
                    if pulsenum in [1,3]: col = 'm'
                    if pulsenum in [2,4]: col = 'c'
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
                        color='m',linewidth=0,alpha=0.4)
                    fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+2],\
                        color='c',linewidth=0,alpha=0.4)
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
                    simresponse+firingbinserrList[pulsenum]/sqrt(numtrials),
                    simresponse-firingbinserrList[pulsenum]/sqrt(numtrials),
                    color=(0.7,0.7,0.7))
                plot(pulsetlist,simresponse,linewidth=linewidth,color=(0.3,0.3,0.3))
            else:
                errorbar(pulsetlist,y=simresponse,\
                    yerr=firingbinserrList[pulsenum]/sqrt(numtrials),\
                    color=(pulsenum/float(numpulses),1-pulsenum/float(numpulses),sister_ratio),\
                    marker='+',linestyle='solid',linewidth=linewidth,label='Mitral response')
            ################## Plot the fitted responses
            starti = int(PULSE_START/pulserebindt)
            if pulsenum in [1,3]: titlestr = 'Odor A random pulse'
            elif pulsenum in [2,4]: titlestr = 'Odor B random pulse'
            else: titlestr = 'Odor A and odor B random pulse'
            response_unrect = fitted_responses_unrect[pulsenum-1]
            if fill_below:
                plot(pulsetlist[starti:],response_unrect[starti:],
                    linestyle='solid',linewidth=linewidth,color=['r','b',(1,0,1)][pulseiter],
                    label=['air','A','B','A','B','A&B'][pulsenum])
                #lgd = legend()
                #for k in lgd.get_lines():
                #    k.set_linewidth(0)
                #lgd.draw_frame(False)
                #ltext  = lgd.get_texts()
                #for l in ltext:
                #    l.set_fontsize(label_fontsize)
                ##fig.setp(ltext, fontsize=20)
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
                ## below line should be outside this else: block,
                ## but now I've used it in two figures and all alignment will get spoiled
                ## so leaving it as is!!
                ax.set_ylim(0,xpulse)
            ax.set_yticks([0,xpulse])
            ax.set_yticklabels(['0',str(xpulse)])
            ax.set_xlim(0,9)
            ax.set_xticks([])
            ax.get_xaxis().set_ticks_position('none')
            ax.get_yaxis().set_ticks_position('left')
            xmin, xmax = ax.get_xaxis().get_view_interval()
            ymin, ymax = ax.get_yaxis().get_view_interval()
            ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
            if pulseiter==1:
                axes_labels(ax,'','firing rate (Hz)',adjustpos=False)
                ax.yaxis.set_label_coords(-0.12,1.1)
                add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
                    sizex=1,labelx='1 s',bbox_to_anchor=[1.0,0.65],bbox_transform=ax.transAxes)
            elif pulseiter==2:
                axes_labels(ax,'','firing rate (Hz)',adjustpos=False,ypad=2)
                add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
                    sizex=1,labelx='1 s',bbox_to_anchor=[0.6,-0.1],bbox_transform=ax.transAxes)
            else:
                axes_labels(ax,'','',adjustpos=False) # sets default label_fontsize
        ax.get_xaxis().set_ticks_position('none')

        ## beautify panels for kernels
        for ax_kernel in [ax_kernelsA,ax_kernelsA1,ax_kernelsB,ax_kernelsB1]:
            ax_kernel.set_yticks([0])
            ax_kernel.set_xlim(0,kernel_time)
            ax_kernel.set_xticks([])
            axes_off(ax_kernel,x=True,y=False)
            ## turn on/off the side axes (after turning off the frame above):
            ## http://www.shocksolution.com/2011/08/removing-an-axis-or-both-axes-from-a-matplotlib-plot/
            ax_kernel.get_xaxis().set_ticks_position('none')
            ax_kernel.get_yaxis().set_ticks_position('left')
            xmin, xmax = ax_kernel.get_xaxis().get_view_interval()
            ymin, ymax = ax_kernel.get_yaxis().get_view_interval()
            ax_kernel.add_artist(Line2D((xmin, xmin), (ymin, ymax),\
                color='black', linewidth=axes_linewidth))
            ## need to add the axes lines that I want, after deleting full frame.
            ## axes_labels() also sets label sizes to default label_fontsize
            axes_labels(ax_kernel,'','',adjustpos=False)
        for ax_kernel in [ax_kernelsA,ax_kernelsA1]:
            add_scalebar(ax_kernel,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizex=0.5,labelx='0.5 s',sizey=100,labely='arb',\
                bbox_to_anchor=[1.0,-0.6],bbox_transform=ax_kernel.transAxes)

        fig.tight_layout()
        fig1.tight_layout()
        ## do not set clip off, since kernels begin from -25ms due to binning issues.
        #fig_clip_off(fig)
        fig.subplots_adjust(top=0.9,bottom=0.1) # prediction fig: Fig 5
        fig1.subplots_adjust(left=0.05,top=0.94,bottom=0.1,right=0.95,wspace=0.6,hspace=0.1) # fitting fig : Fig S1
        if self.SAVEFIG:
            fig.savefig('../figures/linearity_example_prediction_'+str(self.stimseed)+\
                '_mit'+str(fitted_mitral)+'.svg',dpi=fig.dpi)
            fig.savefig('../figures/linearity_example_prediction_'+str(self.stimseed)+\
                '_mit'+str(fitted_mitral)+'.png',dpi=fig.dpi)
            fig1.savefig('../figures/linearity_example_fits_'+str(self.stimseed)+\
                '_mit'+str(fitted_mitral)+'.svg',dpi=fig.dpi)
            fig1.savefig('../figures/linearity_example_fits_'+str(self.stimseed)+\
                '_mit'+str(fitted_mitral)+'.png',dpi=fig.dpi)

    ## plot for one odornum 5 / 0 for both sisters: sis mit 0 vs sis mit 1
    def plot_kernels_odornum_papersubfigure(self,odornum):
        ## load the simulated and fitted variables
        match_resp = self.match_resp
        fig = figure(figsize=(columnwidth/3.0,linfig_height/3.0),dpi=300,facecolor='w') # 'none' is transparent
        ax_kernels = fig.add_subplot(111)
        if match_resp:
            fig_resp = figure(figsize=(columnwidth/3.0,linfig_height/3.0),dpi=300,facecolor='w') # 'none' is transparent
            ax_resp = fig_resp.add_subplot(111)

        for fitted_mitral in [central_glom,central_glom+1]:
            (kernelA,kernelB) = self.fittedkernels[fitted_mitral]
        
            ## ----------------------- smooth kernels only for display of prediction & fit --------------------------

            ## irrespective of smooth kernels flag, filter the output kernels for display
            ## IMPORTANT: SMOOTHED KERNEL IS USED FOR BELOW FOR PREDICTION & FIT
            ## ONLY FOR DISPLAY PURPOSES!!!!
            ## In saving _params file above and in fit_odor_pulses_.....py, I use the original kernels.
            smooth_kernelA = SGfilter2(append([0],kernelA))
            smooth_kernelB = SGfilter2(append([0],kernelB))
            if odornum==5: smooth_kernel = smooth_kernelA
            else: smooth_kernel = smooth_kernelB
            
            ############################## plot the kernels

            ## kernel actually starts from kerneldt/2.0, but add a time point a bin earlier,
            ## as while smoothing Savistsky-Golay above, I set kernel[-kerneldt/2.0]=0.0
            ## Ideally I should have set kernel[0]=0, but SG filter cannot handle non-uniform time-series.
            plot_kerneltlist = append([-kerneldt/2.0],kerneltlist+kerneldt/2.0)
            ax_kernels.plot(plot_kerneltlist, smooth_kernel, color=['r','b'][fitted_mitral], \
                marker=['s','o'][fitted_mitral], ms=marker_size/2.0, \
                linestyle='solid',linewidth=linewidth,label='mit '+str(fitted_mitral))

            if match_resp:
                ## Important to first convolve, then pre-truncate the flow rate trace
                ## Else the first resp cycle does not match.
                respirationtime,morph_responseA,morph_responseB = \
                    predict_respresp(kernelA,kernelB)
                if odornum==5: morph_response = morph_responseA
                else: morph_response = morph_responseB
                ## plot simulated resp responses
                ## Load the mitral responses also
                morph_numavgs,morph_mitral_responses_avg,morph_mitral_responses_std = \
                    load_morph_responses_from_pulsefilename(\
                    self.filename,fitted_mitral,respdt,numresps=numrespsfit)
                ## Need to subtract the resp air response, not the const 'bgnd'
                airresponse = morph_mitral_responses_avg[6]
                simresponse = morph_mitral_responses_avg[odornum] - airresponse
                simerr = morph_mitral_responses_std[odornum]/sqrt(morph_numavgs)
                ax_resp.fill_between(respirationtime, simresponse+simerr, simresponse-simerr,
                    color=['r','b'][fitted_mitral],alpha=0.4,linewidth=0)
                ax_resp.plot(respirationtime,simresponse,color=['m','c'][fitted_mitral],\
                    marker=['s','o'][fitted_mitral], ms=marker_size/2.0, linewidth=linewidth)
                ## plot predicted resp responses
                ax_resp.plot(respirationtime,morph_response,color=['r','b'][fitted_mitral],\
                    marker=['+','x'][fitted_mitral], ms=marker_size, linewidth=linewidth)

        ## beautify panels for kernels
        beautify_plot(ax_kernels,x0min=True,y0min=False,xticksposn='none',\
                yticksposn='left',drawxaxis=True)
        ax_kernels.set_yticks([0])
        ax_kernels.set_xlim(0,kernel_time)
        ax_kernels.set_xticks([])
        add_scalebar(ax_kernels,matchx=False,matchy=False,hidex=True,hidey=False,\
            sizex=0.5,labelx='0.5 s',sizey=100,labely='arb',\
            bbox_to_anchor=[0.9,0.75],bbox_transform=ax_kernels.transAxes)
        ## axes_labels() also sets label sizes to default label_fontsize
        axes_labels(ax_kernels,'','',adjustpos=False)
        fig.tight_layout()
        ## do not set clip off, since kernels begin from -25ms due to binning issues.
        #fig_clip_off(fig)
        fig.subplots_adjust(top=0.7)
        if self.SAVEFIG:
            fig.savefig('../figures/decorr/decorr_kernels_'+\
                str(self.stimseed)+'.svg',dpi=fig.dpi)
            fig.savefig('../figures/decorr/decorr_kernels_'+\
                str(self.stimseed)+'.png',dpi=fig.dpi)
            
        if match_resp:
            ## beautify plot resp
            add_scalebar(ax_resp,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizex=0.1,labelx='0.1 s',sizey=15,labely='15Hz',\
                bbox_to_anchor=[0.9,0.75],bbox_transform=ax_resp.transAxes)
            xmin,xmax,ymin,ymax = beautify_plot(ax_resp,x0min=False,y0min=False,
                    xticksposn='none',yticksposn='left',drawxaxis=False)
            ymin,ymax = int(floor(ymin)),int(ceil(ymax))
            ax_resp.set_yticks([0,ymax])
            axes_labels(ax_resp,'','') # sets default fontsize too

            fig_resp.tight_layout()
            ## do not set clip off, since kernels begin from -25ms due to binning issues.
            #fig_clip_off(fig_resp)
            fig_resp.subplots_adjust(top=0.7)
            if self.SAVEFIG:
                fig_resp.savefig('../figures/decorr/decorr_resps_'+\
                    str(self.stimseed)+'.svg',dpi=fig.dpi)
                fig_resp.savefig('../figures/decorr/decorr_resps_'+\
                    str(self.stimseed)+'.png',dpi=fig.dpi)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        ######### should we match mitral morph responses using kernel
        if 'match_resp' in sys.argv:
            match_resp = True
            ## autocalculate the morphfilename
            ## by replacing 'morphs' with 'pulses' in pulses filename
            #mmfilenameindex = sys.argv.index('--match_resp')+1
            #morph_filename = sys.argv[mmfilenameindex]
        else: match_resp = False
        if 'NOSHOW' in sys.argv: NOSHOW = True
        else: NOSHOW = False
        if 'SAVEFIG' in sys.argv: SAVEFIG = True
        else: SAVEFIG = False
        if 'TEST' in sys.argv:
            test_fitpulses(sys.argv[1],sys.argv[2])
        else:
            filename = sys.argv[1]
            worker = fit_plot_pulses(filename,sys.argv[2],match_resp,NOSHOW,SAVEFIG)
            post_pulses = filename.split('odor_pulses')[1]
            dirextn = post_pulses.split('/')[0]
            print 'directory extension is',dirextn
            worker.fit_pulses(sys.argv[1],sys.argv[2],match_resp,NOSHOW,SAVEFIG,dirextn,dirextn4stim=False)
            if not NOSHOW:
                ## Below plots pulse fits and predictions for a given mitral
                worker.plot_pulses_mitnum_papersubfigure_v2(0) # set mitnum 0/1
                ## Below plots kernels and resp fits for a given odor for both mits
                worker.plot_kernels_odornum_papersubfigure(5) # set 5/0 for odor A/B
                show()
    else:
        print "At least specify data file containing pickled mitral responses, and ORN frate seed."
        sys.exit(1)
