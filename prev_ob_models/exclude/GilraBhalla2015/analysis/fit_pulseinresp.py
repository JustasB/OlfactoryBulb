# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO ROUGHLY FOLLOW PRIYANKA'S ANALYSIS
########## This variant does not use an air kernel, rather a constant air rate.
## USAGE1: python2.6 fit_pulseinresp.py <../results/odor_morphs/pulseinresp_result.pickle> <stimseed>

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

RUNTIME = (NUM_RESPS+1)*RESPIRATION
STARTTIME = RESPIRATION

numpulses = RANDOM_PULSE_NUMS*2
## rebin the spikes as below, irrespective of previous binning
pulserebindt = 50e-3#fitting_dt # 50ms as Priyanka uses
pulserebins = int(RUNTIME/pulserebindt)
#bin_width_time = 2*pulserebindt ## overlapping bins with this bin width
pulsetlist = arange(pulserebindt/2.0,RUNTIME,pulserebindt)

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
kerneldt = 50e-3 # 50 ms as Priyanka uses
kernel_time = 0.5 # s 
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

def SGfilter2(kernel):
    ## savitsky-golay filter maintains high frequency signal, unlike moving average filter
    ## window_size must be odd, savitsky_golay copied from scipy cookbook online
    window_size = 11 # same as Priyanka
    kernel = savitzky_golay(kernel, window_size, order=4) # same as Priyanka
    return kernel

def chisqfunc(params, firingbinsmean, firingbinserr, pulse, ret_array=True):
    kernel = params[0:kernel_size]
    ## smoothing kernels here does not help, as number of params remain same
    ## and fitting routine can always modify raw values to nullify effect of smoothing.
    if SMOOTH_KERNELS:
        kernel = SGfilter2(kernel)
    starti = int(STARTTIME/pulserebindt)
    chisqarray = []    
    ## These below should be numpy arrays,
    ## else + acts like append in python (not element-wise addition)!!!
    responsecalc = \
        rectify(responseconvolve(pulse,kernel))
    ## Priyanka minimizes least sq unlike Adil who minimizes chi-sq
    ## With chi-sq, there is 'what min err to choose' issue.
    ## If usechisq, then do least squares, else chi-sq
    if not usechisq:
        chisqarray.extend( [ (val-firingbinsmean[i])\
            for i,val in enumerate(responsecalc[starti:],start=starti) ] )
    else:
        chisqarray.extend( [ (val-firingbinsmean[i])/firingbinserr[i]\
            for i,val in enumerate(responsecalc[starti:],start=starti) ] )
        
    global iternum
    iternum += 1
    ## normalize chi-sq to the number of dof
    num_dof = float(firingbinsmean[starti:].size - params.size)
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

class fit_plot_pulseinresp():
    def __init__(self,filename,stimseed):
        self.filename = filename
        self.stimseed = stimseed

    def fit_pulses(self,dirextn=''):
        ######### load in the mitral pulse responses
        #### mitral_responses_list[avgnum][odornum][mitralnum][binnum]
        #### mitral_responses_binned_list[avgnum][odornum][mitralnum][binnum]
        mitral_responses_list, mitral_responses_binned_list = read_pulsefile(self.filename)
        ########## load in the stimuli
        pulseInRespList,genORNkernels \
            = read_pulseinresp_stimuli(stimseed,'../generators/firerates'+dirextn)
        self.pulseInRespList = pulseInRespList

        ##-------------------------- rebin the responses and pulses ------------------------------
        ## rebin sim responses to pulserebindt=50ms, then take mean
        numtrials,mitral_responses_mean,mitral_responses_std = \
                rebin_mean(mitral_responses_list,pulserebindt,RUNTIME)
        ## for pulseList, decimate by taking every convolutiondt/FIRINGFILLDT bins
        ## for ORNkernels, decimate by taking every kerneldt/FIRINGFILLDT=50/1=50 bins
        pulserebinList,linORNkernelA,linORNkernelB = \
                decimate([pulseInRespList[central_glom][0]],convolutiondt,genORNkernels,kerneldt,kernel_size)

        ## Need to boot up the plotting, before the loop over the two fitted_mitral-s (central sisters)
        if not NOSHOW:
            fig = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
            ax_kernelsA = plt.subplot2grid((12,3),(0,2),rowspan=3,frameon=False)
            text(0.1,1.0,'G', fontsize=label_fontsize, transform = ax_kernelsA.transAxes)
            ax_kernelsB = plt.subplot2grid((12,3),(3,2),rowspan=3,frameon=False)
            text(0.1,1.0,'H', fontsize=label_fontsize, transform = ax_kernelsB.transAxes)

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

            air_bgnd = firingbinsmeanList[0]

            #---------------------------- fitting ------------------------------------------------

            ## Initial values for the parameters
            if firstrun:
                ## ORN kernelA/B is initial kernelA/B
                param_kernels = [linORNkernelA,linORNkernelB]
                ### can also use all zeros as init kernels, works for reasonably positive responses
                #param_kernels = [[0.0]*len(linORNkernelA),[0.0]*len(linORNkernelB)]
                ORNkernels = (linORNkernelA, linORNkernelB)
            else:
                print 'saving not implemented'
                sys.exit()
                f = open(sys.argv[1]+'_params'+str(fitted_mitral),'r')
                chisqs,param_kernels,discard_firingmeans,\
                        discard_firingsd,discard_fitresponses = pickle.load(f)
                f.close()
                ORNkernels = ()

            chisqs_mit = []
            ## fit each odor A and B separately, but fitting the corresponding 2 random pulses for each
            for odornum in [0,1]:
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
                    args=(firingbinsmeanList[odornum+1]-air_bgnd, firingbinserrList[odornum+1],\
                    pulserebinList[0], True),\
                    full_output=1, maxfev=100000, ftol=1e-20, xtol=1e-100 )
                print params[3] # print the status message
                param_kernels[odornum] = params[0] # take only fitted params; leastsq returns a tuple with errmsg, etc.
                ### Calculate sum of squares of the chisqarray
                chisqarraysq = [i**2 for i in chisqfunc(param_kernels[odornum],\
                    firingbinsmeanList[odornum+1]-air_bgnd, firingbinserrList[odornum+1],\
                    pulserebinList[0], True)]
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
            response_unrect = responseconvolve(pulserebinList[0],kernelA) + air_bgnd
            response = rectify(response_unrect)
            fitted_responses_unrect.append(response_unrect)
            fitted_responses.append(response)
            response_unrect = responseconvolve(pulserebinList[0],kernelB) + air_bgnd
            response = rectify(response_unrect)
            fitted_responses_unrect.append(response_unrect)
            fitted_responses.append(response)

            ## ----------------------- smooth kernels only for display of prediction & fit --------------------------

            ## irrespective of smooth kernels flag, filter the output kernels for display
            ## IMPORTANT: SMOOTHED KERNEL IS USED FOR BELOW FOR PREDICTION & FIT
            ## ONLY FOR DISPLAY PURPOSES!!!!
            ## In saving _params file above and in fit_odor_pulses_.....py, I use the original kernels.
            smooth_kernelA = SGfilter2(append([0],kernelA))
            smooth_kernelB = SGfilter2(append([0],kernelB))
                
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
                for odoriter, odornum in enumerate([1,2]):#range(1,numpulses): # similar to Priyanka, only 3 plots
                    sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
                    ## similar to Priyanka: only 3 plots, skip 2 pulse responses
                    ax = plt.subplot2grid((12,3),(4*odoriter,fitted_mitral),rowspan=4,frameon=False)
                    ################# random pulses and ORN responses
                    pulseinresptime = arange(0,RUNTIME+1e-10,FIRINGFILLDT)
                    ################### Plot the simulated responses
                    ## smooth the simulated response
                    ## windowsize=5 and SD=0.65 are defaults from matlab's smoothts() for gaussian smoothing
                    Gwindow = signal.gaussian(5,0.65)
                    ## help from http://www.scipy.org/Cookbook/SignalSmooth
                    simresponse = convolve(Gwindow/Gwindow.sum(),firingbinsmeanList[odornum],mode='same')
                    ## numpy array, hence adds element by element
                    fill_between(pulsetlist,
                        simresponse+firingbinserrList[odornum]/sqrt(numtrials),
                        simresponse-firingbinserrList[odornum]/sqrt(numtrials),
                        color=(0.7,0.7,0.7))
                    plot(pulsetlist,simresponse,linewidth=linewidth,color=(0.3,0.3,0.3))
                    ################## Plot the fitted responses
                    starti = int(STARTTIME/pulserebindt)
                    response_unrect = fitted_responses_unrect[odornum-1]
                    plot(pulsetlist[starti:],response_unrect[starti:],
                        linestyle='solid',linewidth=linewidth,color=['r','b'][fitted_mitral])
                    axes_labels(ax,'','',adjustpos=False) # sets default label_fontsize
                    ax.get_xaxis().set_ticks_position('none')
                    ax.get_yaxis().set_ticks_position('left')
                    xmin, xmax = ax.get_xaxis().get_view_interval()
                    ymin, ymax = ax.get_yaxis().get_view_interval()
                    ax.add_artist(Line2D((xmin, xmin), (ymin, ymax), color='black', linewidth=axes_linewidth))
                ax.get_xaxis().set_ticks_position('bottom')
                ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=axes_linewidth))
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
            
            fig.tight_layout()
            subplots_adjust(top=0.94)

        self.fittedkernels = fittedkernels
        self.chisqs = chisqs
        self.fitted_responses_full = fitted_responses_full
        self.fitted_responses_unrect_full = fitted_responses_unrect_full
        self.mdict_full = mdict_full
        return fittedkernels,chisqs,fitted_responses_full,mdict_full

    ## plot for one mitnum (fitted_mitral)...
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
        text(0.05,1.0,'f', fontweight='bold', fontsize=label_fontsize, transform = ax_kernelsA.transAxes)
        ax_kernelsB = plt.subplot2grid((3,6),(1,0),rowspan=1,colspan=2,frameon=False)
        text(0.05,1.0,'g', fontweight='bold', fontsize=label_fontsize, transform = ax_kernelsB.transAxes)
        if match_resp:
            ax_resp = plt.subplot2grid((3,6),(2,0),rowspan=1,colspan=2,frameon=False)
            text(0.05,1.0,'k', fontweight='bold', fontsize=label_fontsize, transform = ax_resp.transAxes)

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
            simresponse = morph_mitral_responses_avg[5] - airresponse
            simerr = morph_mitral_responses_std[5]/sqrt(morph_numavgs)
            ax_resp.fill_between(respirationtime, simresponse+simerr, simresponse-simerr,
                color='r',alpha=0.4,linewidth=0)
            ax_resp.plot(respirationtime,simresponse,color='m',linewidth=linewidth)
            simresponse = morph_mitral_responses_avg[0] - airresponse
            simerr = morph_mitral_responses_std[0]/sqrt(morph_numavgs)
            ax_resp.fill_between(respirationtime, simresponse+simerr, simresponse-simerr,
                color='b',alpha=0.4,linewidth=0)
            ax_resp.plot(respirationtime,simresponse,color='c',linewidth=linewidth)
            ## plot predicted resp responses
            ax_resp.plot(respirationtime,morph_responseA,color='r',linewidth=linewidth)
            ax_resp.plot(respirationtime,morph_responseB,color='b',linewidth=linewidth)
            add_scalebar(ax_resp,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizex=0.1,labelx='0.1 s',bbox_to_anchor=[0.8,-0.1],bbox_transform=ax_resp.transAxes)
            add_scalebar(ax_resp,matchx=False,matchy=False,hidex=True,hidey=False,\
                sizey=15,labely='15 Hz',bbox_to_anchor=[0.7,0.85],bbox_transform=ax_resp.transAxes)

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
        for pulseiter, pulsenum in enumerate([1,2,5]):#range(1,numpulses): # similar to Priyanka, only 3 plots
            sister_ratio = 0#(mitnum%MIT_SISTERS)/float(MIT_SISTERS)
            ## similar to Priyanka: only 3 plots, skip 2 pulse responses
            ax = plt.subplot2grid((3,6),(pulseiter,2),rowspan=1,colspan=4,frameon=False)
            panel_label = ['h','i','j'][pulseiter]
            ## ax.transAxes ensures relative to axes size, rather than to data units.
            text(0.05, 1.0, panel_label, fontweight='bold', fontsize=label_fontsize, transform = ax.transAxes)
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
            sizex=0.5,labelx='0.5 s',bbox_to_anchor=[0.8,-0.6],bbox_transform=ax_kernelsA.transAxes)
        add_scalebar(ax_kernelsA,matchx=False,matchy=False,hidex=True,hidey=False,\
            sizey=100,labely='arb',bbox_to_anchor=[0.5,-0.4],bbox_transform=ax_kernelsA.transAxes)
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
    NOSHOW = False
    if len(sys.argv) > 2:
        filename = sys.argv[1]
        stimseed = sys.argv[2]
        worker = fit_plot_pulseinresp(filename,stimseed)
        post_pulses = filename.split('odor_morphs')[1]
        dirextn = post_pulses.split('/')[0]
        print 'directory extension is',dirextn
        worker.fit_pulses(dirextn)
        show()
    else:
        print "At least specify data file containing pickled mitral responses, and ORN frate seed."
        sys.exit(1)
