# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: python2.6 average_odor_pulses.py [SAVEFIG]

import os,sys,math,string
import os.path
import pickle
import subprocess
cwd = os.getcwd() # current working directory

sys.path.extend(["..","../networks","../generators","../simulations"])

from stimuliConstants import * # has RESPIRATION
from sim_utils import * # has rebin() and imports data_utils.py for axes_off()
from calc_corrs import * # has calc_corrs()
## BE CAREFUL: use _constair_separateodors below, others are obsolete.
from fit_odor_pulses_constair_separateodors import * # has fit_pulses() and kerneltlist
import average_odor_morphs as corr_utils # has get_filename() for morphs

from pylab import * # part of matplotlib that depends on numpy but not scipy
from scipy import stats
from scipy import signal # for gaussian smoothing

IN_VIVO = True
directed = True
FRAC_DIRECTED = 0.01 ## overrides the one in networkConstants.py (came in via sim_utils.py)
## Below two overide the variables that came in from stimuliConstantsMinimal.py via stimuliConstants.py
NONLINEAR_ORNS = False
NONLINEAR_TYPE = 'P' # P for primary glom non-linear, L for lateral gloms non-linear

#fullfilename = '../results/odor_morphs/morphs_random'
#if NONLINEAR_ORNS: fullfilename += 'NL'+NONLINEAR_TYPE
#fullfilename += '.pickle'
#fullfile = open(fullfilename,'r')
#morphs = pickle.load(fullfile)
#fullfile.close()

PLOTRESP_NUM = 1 # whether to plot 2 respiration cycles or 1
NUMBINS = 5
BIN_WIDTH_TIME = RESPIRATION/NUMBINS
bindt = RESPIRATION/float(NUMBINS)
## I take the last PLOTRESP_NUM of respiration cycles out of NUM_RESPS simulated
responsetlist = arange( SETTLETIME+(NUM_RESPS-PLOTRESP_NUM)*RESPIRATION+bindt/2.0, \
    SETTLETIME+NUM_RESPS*RESPIRATION, RESPIRATION/NUMBINS )

min_frate_cutoff = 0.25 # Hz

salient = False#True
if salient:
    stim_seeds = [-1,-19]#range(-1,-37,-1)#[-1,-2,-3,-4,-5,-6,-7,-8,-9]
    num_gloms_list = [3]
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ## in order,below options are:
    ## all cells; no lat; no joints, varyRMP; no PGs; no singles + no joints, only mitrals
    #inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)), \
    #    (2,(False,True,False,False,False)), (3,(False,False,False,False,True)), \
    #    (4,(False,False,False,True,False)), (5,(True,True,False,False,False)), \
    #    (6,(True,True,False,True,False))]
    inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)) ]
else:
    #stim_seeds = [157.0,160.0,190.0,191.0,212.0,441.0]
    #num_gloms_list = [5,2]
    stim_seeds = arange(750.0,800.0,1.0)#[157.0,160.0,190.0,191.0]
    num_gloms_list = [3]
    ## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
    ## in order,below options are: all cells; no lat; no joints, varyRMP; no PGs; only mitrals
    inh_options = [ (0,(False,False,False,False,False)), (1,(False,False,True,False,False)), \
        (2,(True,False,False,False,False)), (3,(True,True,False,False,True)), \
        (6,(False,False,False,True,False)), (8,(True,True,False,True,False))]
net_seeds = [200.0]

def get_filename(netseed,stimseed,inh,ngi,stimi,neti,\
        nl_orns=NONLINEAR_ORNS,nl_type=NONLINEAR_TYPE,\
        resultsdir='../results/odor_pulses'):
    ### read filename from the output file of automated run
    ## construct the filename
    if inh[0]: singles_str = '_NOSINGLES'
    else: singles_str = '_SINGLES'
    if inh[1]: joints_str = '_NOJOINTS'
    else: joints_str = '_JOINTS'
    if inh[3]: pgs_str = '_NOPGS'
    else: pgs_str = '_PGS'
    if inh[2]: lat_str = '_NOLAT'
    else: lat_str = '_LAT'
    if inh[4]: varmitstr = '_VARMIT'
    else: varmitstr = '_NOVARMIT'
    ## stable enough that time tags are not needed
    filename = resultsdir+'/odorpulses'+\
        '_netseed'+str(netseed)+'_stimseed'+str(stimseed)
    if nl_orns: filename += '_NL'+nl_type
    filename += singles_str+joints_str+pgs_str+lat_str+varmitstr+\
        '_numgloms'+str(num_gloms_list[ngi])
    if directed: filename += '_directed'+str(FRAC_DIRECTED)
    filename += '.pickle'
    return filename, (singles_str, joints_str, pgs_str, lat_str, varmitstr)

def plot_scaled_kernels():
    for neti,netseed in enumerate(net_seeds):
        numsubfigs_rows = len(stim_seeds)*len(num_gloms_list)
        numsubfigs_cols = len(inh_options)*2 # two mitral sisters
        fig = figure(facecolor='w')
        ax = fig.add_subplot(111)
        figtext(0.1,0.94,'OdorA & OdorB kernels x1 and x2 netseed='+str(netseed),fontsize=20)
        delaxes()

        if not salient: net_seeds = [stimseed]
        for stimi,stimseed in enumerate(stim_seeds):
            for ngi,num_gloms in enumerate(num_gloms_list):
                ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
                for plotyi,(inhi,inh) in enumerate(inh_options):

                    ## add_subplot(rows,cols,fignum), fignum starts from 1 not 0
                    ## fignum fills first row to end, then next row
                    ## one subfig for mit0, one for mit1
                    ax0 = fig.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                        (stimi+ngi)*numsubfigs_cols+plotyi*2+1) # *2 for two mit sisters
                    ax1 = fig.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                        (stimi+ngi)*numsubfigs_cols+plotyi*2+2) # *2 for two mit sisters

                    all_chisqs = []
                    for dualnum,dualstimseed in enumerate( (stimseed,stimseed-len(salient_seed_glom_mapA)) ):
                        filename, switch_strs \
                            = get_filename(netseed,dualstimseed,inh,ngi,stimi,neti)
                        switches_str = string.join(switch_strs,'')
                        print filename

                        ## read in the spiketimes/binned responses for this run
                        kernels,chisqs,fitted_responses_full,mdict_full = \
                            fit_pulses(filename,dualstimseed,False,NOSHOW=True,SAVEFIG=False)
                        all_chisqs.append(chisqs)

                        ax0.plot(kerneltlist,kernels[0][0],color=(1,0,dualnum),linewidth=2)
                        ax0.plot(kerneltlist,kernels[0][1],color=(0,1,dualnum),linewidth=2)
                        ax1.plot(kerneltlist,kernels[1][0],color=(1,0,dualnum),linewidth=2)
                        ax1.plot(kerneltlist,kernels[1][1],color=(0,1,dualnum),linewidth=2)
                        if dualnum==0: ## double the kernels to compare with x2 kernels
                            ax0.plot(kerneltlist,kernels[0][0]*2,color=(1,0,1),linewidth=2,linestyle=':')
                            ax0.plot(kerneltlist,kernels[0][1]*2,color=(0,1,1),linewidth=2,linestyle=':')
                            ax1.plot(kerneltlist,kernels[1][0]*2,color=(1,0,1),linewidth=2,linestyle=':')
                            ax1.plot(kerneltlist,kernels[1][1]*2,color=(0,1,1),linewidth=2,linestyle=':')

                    ax0.set_title(
                        'mit0_'+str(num_gloms)+'G_'+\
                        '%1.2f,%1.2f'%(all_chisqs[0][0],all_chisqs[1][0])+switches_str,fontsize=8)
                    ax1.set_title(
                        'mit1_'+str(num_gloms)+'G_'+\
                        '%1.2f,%1.2f'%(all_chisqs[0][1],all_chisqs[1][1])+switches_str,fontsize=8)

## put whichever netseed, stimseed, mit and odor you want plotted separately a la Priyanka
def plot_scaled_kernels_special():
    netseed = 200.0
    stimseed = -19
    ngi = 0
    mitnum = 1 # 0 or 1
    odornum = 0 # 0 for odorA, 1 for odorB
    # use 1 or 3 for odorA, 2 or 4 for odorB
    if odornum==0: pulsenum = 1
    else: pulsenum = 2

    fig = figure(facecolor='w')
    ## set main axes off
    bigAxes = axes(frameon=False) # hide frame
    xticks([]) # don't want to see any ticks on this axis
    yticks([])
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh = (False,False,False,False,False)
    ax = fig.add_subplot(1,3,1) # for kernels

    all_chisqs = []
    for dualnum,dualstimseed in enumerate( (stimseed,stimseed-len(salient_seed_glom_mapA)) ):
        filename, switch_strs \
            = get_filename(netseed,dualstimseed,inh,ngi,None,None,None)
        switches_str = string.join(switch_strs,'')
        print filename

        ## read in the spiketimes/binned responses for this run
        ##mdict={'firingbinsmeanList':firingbinsmeanList,
        ##    'firingbinserrList':firingbinserrList, 'bgnd':bgnd, 'FIRINGFILLDT':FIRINGFILLDT,
        ##    'pulseList':pulserebinList,'pulseStepsList':pulseStepsList,
        ##    'start_kernels':ORNkernels,'pulserebindt':pulserebindt,'pulsetlist':pulsetlist}
        ## _full means data for both mitrals
        kernels,chisqs,fitted_responses_full,mdict_full = \
            fit_pulses(filename,dualstimseed,False,NOSHOW=True)
        fitted_responses = fitted_responses_full[mitnum]
        mdict = mdict_full[mitnum]
        all_chisqs.append(chisqs)
        FIRINGFILLDT=mdict['FIRINGFILLDT']
        pulserebindt=mdict['pulserebindt']
        pulsetlist=mdict['pulsetlist']
        pulseStepsList = mdict['pulseStepsList']
        firingbinsmeanList = mdict['firingbinsmeanList']
        firingbinserrList = mdict['firingbinserrList']

        ## plot kernel
        if dualnum==0: scale = 1.0 
        else: scale = 0.5 ## halve the kernels to compare with 1x kernel
        ax.plot(kerneltlist,kernels[mitnum][odornum]*scale,\
            color=(1-dualnum,0,dualnum),linewidth=4)
        
        ax2 = fig.add_subplot(1,3,dualnum+2) # for pulse, response and fit

        ################# random pulses and ORN responses
        xpulse = 50
        randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)
        fill_between(randompulsetime,xpulse*pulseStepsList[pulsenum+1],
            linewidth=0,color=(0.9,0.7,0.7),alpha=1.0)#facecolor=(0.9,0.7,0.7),edgecolor=(0.9,0.7,0.7))
        ################### Plot the simulated responses
        ## smooth the simulated response
        ## windowsize=5 and SD=0.69 are defaults from matlab's smoothts() for gaussian smoothing
        Gwindow = signal.gaussian(5,0.69)
        ## help from http://www.scipy.org/Cookbook/SignalSmooth
        simresponse = convolve(Gwindow/Gwindow.sum(),firingbinsmeanList[pulsenum],mode='same')
        ## numpy array, hence adds element by element
        fill_between(pulsetlist,
            simresponse+firingbinserrList[pulsenum]/sqrt(9),
            simresponse-firingbinserrList[pulsenum]/sqrt(9),
            color=(0.6,0.9,0.6))
        plot(pulsetlist,simresponse,linewidth=2,color=(0.3,0.3,0.3))
        ################## Plot the fitted responses
        starti = int(PULSE_START/pulserebindt)
        response = fitted_responses[pulsenum-1]
        plot(pulsetlist[starti:],response[starti:],linestyle='solid',
            linewidth=2.0,color=(1-dualnum,0,dualnum))

        ax2.set_ylim(-xpulse/20,xpulse)
        ax2.set_yticks([0,xpulse])
        ax2.set_yticklabels(['0',str(xpulse)])
        ax2.set_xlim(0,8.5)
        ax2.set_xticks([0,8.5])
        ax2.set_xticklabels(['0','8.5'])
        axes_labels(ax2,'time (s)','',adjustpos=False,fontsize=34)
        ax2.set_title('firing rate (Hz)',fontsize=34)
        
    ax.set_yticks([0])
    ax.set_xlim(0,2.0)
    ax.set_xticks([0,2])
    ax.set_xticklabels(['0','2'])
    axes_labels(ax,'time(s)','amplitude (arb)',adjustpos=False,fontsize=34)

def plot_kernels(graph=True):
    """ Plots kernels for each result file if graph=True, else just fits each result file """
    for stimi,stimseed in enumerate(stim_seeds):
        if graph:
            numsubfigs_rows = len(net_seeds)*len(num_gloms_list)
            numsubfigs_cols = len(inh_options)*2 # two mitral sisters
            fig = figure(facecolor='w')
            ax = fig.add_subplot(111)
            figtext(0.1,0.94,'OdorA & OdorB kernels stimseed='+str(stimseed),fontsize=20)
            delaxes()

        if not salient: net_seeds = [stimseed]
        for neti,netseed in enumerate(net_seeds):
            for ngi,num_gloms in enumerate(num_gloms_list):
                ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
                for plotyi,(inhi,inh) in enumerate(inh_options):

                    if graph:
                        ## add_subplot(rows,cols,fignum), fignum starts from 1 not 0
                        ## fignum fills first row to end, then next row
                        ## one subfig for mit0, one for mit1
                        ax0 = fig.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                            (neti+ngi)*numsubfigs_cols+plotyi*2+1) # *2 for two mit sisters
                        ax1 = fig.add_subplot(numsubfigs_rows,numsubfigs_cols,\
                            (neti+ngi)*numsubfigs_cols+plotyi*2+2) # *2 for two mit sisters

                    filename, switch_strs \
                        = get_filename(netseed,stimseed,inh,ngi,stimi,neti)
                    switches_str = string.join(switch_strs,'')
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    print filename

                    ## fit the responses for this result file
                    kernels,chisqs,fitted_responses_full,mdict_full = \
                        fit_pulses(filename,stimseed,False,NOSHOW=True,SAVEFIG=False)

                    if graph:
                        ax0.plot(kerneltlist,kernels[0][0],color=(1,0,0),linewidth=2)
                        ax0.plot(kerneltlist,kernels[0][1],color=(0,1,0),linewidth=2)
                        ax1.plot(kerneltlist,kernels[1][0],color=(1,0,0),linewidth=2)
                        ax1.plot(kerneltlist,kernels[1][1],color=(0,1,0),linewidth=2)

                        ax0.set_title(
                            'mit0_'+str(num_gloms)+'G_'+\
                            '%1.2f,%1.2f'%chisqs[0]+switches_str,fontsize=8)
                        ax1.set_title(
                            'mit1_'+str(num_gloms)+'G_'+\
                            '%1.2f,%1.2f'%chisqs[1]+switches_str,fontsize=8)

def residual2noise(fitted_response,sim_mean,sim_std,starti):
    residual = mean( [(val-sim_mean[i])**2\
                for i,val in enumerate(fitted_response[starti:],start=starti)] )
    noise = mean( [val**2\
                for i,val in enumerate(sim_std[starti:],start=starti)] )
    signal_mean = mean(sim_mean[starti:])
    #signal_mean = sim_mean[starti] ### HACK presently - but should be this way
    signal = mean( [(val-signal_mean)**2\
                for i,val in enumerate(sim_mean[starti:],start=starti)] )
    signal2noise = sqrt(signal/noise)
    signal2residual = sqrt(signal/residual)
    return signal2noise,signal2residual

def residual2noise_resp(fitted_response,sim_odor_mean,sim_air_mean,sim_std,starti):
    sim_mean = sim_odor_mean - sim_air_mean
    ## uncomment below code to remove those points that have 0 odor frate
    #zeroval_indices = where(sim_odor_mean==0)[0]
    ### Where the val of the odor+air simulated response is zero,
    ### set bins of all arrays to zero, thus they don't contribute to S/N/R.
    #for idx in zeroval_indices:
    #    fitted_response[idx]=0
    #    sim_mean[idx]=0
    #    sim_std[idx]=0
    return residual2noise(fitted_response,sim_mean,sim_std,starti)

def resp_SRN(filename,fitted_mitral,kernelA,kernelB):
    ## goodness of fits: Signal-Residual-Noise for respiratory responses
    ## below functions are from fit_odor_pulses_constair_separateodors.py
    ## minimum error is set by load_morph...() to avoid divide by zero errors
    morph_numavgs,morph_mitral_responses_avg,morph_mitral_responses_std = \
        load_morph_responses_from_pulsefilename(filename,fitted_mitral,respdt,numresps=1)
    respirationtime,morph_responseA,morph_responseB = \
        predict_respresp(kernelA,kernelB,numresps=1)
    ## subtract air response from mitral odor response
    airmean = morph_mitral_responses_avg[6]
    sim_responseA = morph_mitral_responses_avg[5]
    sim_responseB = morph_mitral_responses_avg[0]
    ## normalize the predicted response to the simulated response
    #normA = abs(max(sim_responseA-airmean))/abs(max(morph_responseA))
    #normB = abs(max(sim_responseB-airmean))/abs(max(morph_responseB))
    #morph_responseA /= normA
    #morph_responseB /= normB
    signal2noiseA,signal2residualA = residual2noise_resp( \
        morph_responseA,sim_responseA,airmean,\
        morph_mitral_responses_std[5],0 ) # start from 0th bin
    signal2noiseB,signal2residualB = residual2noise_resp( \
        morph_responseB,sim_responseB,airmean,\
        morph_mitral_responses_std[0],0 ) # start from 0th bin
    return signal2noiseA,signal2residualA,signal2noiseB,signal2residualB,\
            morph_mitral_responses_avg,morph_mitral_responses_std,\
            morph_responseA,morph_responseB

def get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,inh,\
        nl_orns,nl_type,dirextn,resp_fits=True,dirextn4stim=True):
    chisqs = []
    signal2residuals = []
    signal2noises = []
    AandB_signal2residuals = []
    AandB_signal2noises = []
    resp_signal2residuals = []
    resp_signal2noises = []
    n_accept = 0
    n_total = 0
    for stimi,stimseed in enumerate(stim_seeds):
        if not salient: net_seeds = [stimseed]
        for neti,netseed in enumerate(net_seeds):
            for ngi,num_gloms in enumerate(num_gloms_list):

                filename, switch_strs \
                    = get_filename(netseed,stimseed,inh,ngi,stimi,neti,\
                        nl_orns,nl_type,resultsdir='../results/odor_pulses'+dirextn)
                switches_str = string.join(switch_strs,'')
                ## if the result file for these seeds & tweaks doesn't exist,
                ## then carry on to the next.
                print filename
                if not os.path.exists(filename): continue
                ## If the fitted params file does not exist, create it (them).
                if not os.path.exists(filename+'_params0'):
                    ## fit the responses for this result file
                    fitter = fit_plot_pulses(filename,stimseed,False,NOSHOW=True,SAVEFIG=False)
                    kernels,chisqs,fitted_responses_full,mdict_full = \
                        fitter.fit_pulses(filename,stimseed,False,NOSHOW=True,SAVEFIG=False,\
                            dirextn=dirextn,dirextn4stim=dirextn4stim)

                for fitted_mitral in [0,1]:
                    f = open(filename+'_params'+str(fitted_mitral),'r')
                    ## firingbinsmeanList,firingbinserrList [pulsenum][binnum] are for this mitral
                    ## fitted responses [pulsenum-1][binnum] are for this mitral (note no air pulses).
                    chisqs_mit,kernels,bgnd,firingbinsmeanList,firingbinserrList,fitted_responses \
                        = pickle.load(f)
                    f.close()
                    n_total += 1 # total number of mitral cells
                    kernelA,kernelB = kernels
                    
                    ## Priyanka's definitions, following Geffen et al, 2009
                    ## for the last A&B pulse train index 5
                    ## only compare after odor onset at PULSE_START
                    starti = int(PULSE_START/pulserebindt)
                    ## If mean firing rate is close to zero for any odour, the fits are terrible,
                    ## even if the chi-sq comes out low! So discard low firing cells to any odour A or B.
                    if mean(firingbinsmeanList[1][starti:])>min_frate_cutoff \
                            and mean(firingbinsmeanList[2][starti:])>min_frate_cutoff \
                            and mean(firingbinsmeanList[3][starti:])>min_frate_cutoff \
                            and mean(firingbinsmeanList[4][starti:])>min_frate_cutoff:
                        chisqs.append(chisqs_mit[0])
                        chisqs.append(chisqs_mit[1])
                        for pulsenum in [1,2,3,4]:
                            ## Returns with sqrt() taken for S2N and S2R.
                            signal2noise,signal2residual = residual2noise( \
                                fitted_responses[pulsenum-1],firingbinsmeanList[pulsenum],\
                                firingbinserrList[pulsenum],starti )
                            signal2noises.append(signal2noise)
                            signal2residuals.append(signal2residual)
                        signal2noise,signal2residual = residual2noise( \
                            fitted_responses[4],firingbinsmeanList[5],firingbinserrList[5],starti )
                        AandB_signal2noises.append(signal2noise)
                        AandB_signal2residuals.append(signal2residual)
                        n_accept += 1
                        print stimseed,chisqs_mit,'S2R =',signal2residual,'S2N =',signal2noise
                        ## Below are for AandB summation fits, just calculated above
                        if signal2residual<0.5:
                            print 'poor pulses A+B prediction',filename
                        elif signal2residual/signal2noise > 2.0:
                            print 'good pulses A+B prediction, ratio =',signal2residual/signal2noise,filename
                            
                        ## goodness of fits for respiratory responses
                        if resp_fits:
                            signal2noiseA,signal2residualA,signal2noiseB,signal2residualB,\
                                morph_mitral_responses_avg,morph_mitral_responses_std,\
                                morph_responseA,morph_responseB = \
                                    resp_SRN(filename,fitted_mitral,kernelA,kernelB)
                            resp_signal2noises.append(signal2noiseA)
                            resp_signal2residuals.append(signal2residualA)
                            resp_signal2noises.append(signal2noiseB)
                            resp_signal2residuals.append(signal2residualB)
                            #print 'resp ',signal2residualA,'/',signal2noiseA,\
                            #    signal2residualB,'/',signal2noiseB
                            if signal2residualA/signal2noiseA > 1.5:
                                print 'good resp odour A: sqrt(S2R/S2N) = ',signal2residualA,'/',signal2noiseA
                            elif signal2residualA/signal2noiseB > 1.5:
                                print 'good resp odour B: sqrt(S2R/S2N) = ',signal2residualB,'/',signal2noiseB
                            elif signal2residualA<0.5 or signal2residualB<0.5:
                                print 'bad resp odour A: sqrt(S2R/S2N) = ',signal2residualA,'/',signal2noiseA,\
                                    'resp odour B: sqrt(S2R/S2N) = ',signal2residualB,'/',signal2noiseB
                    else:
                        pass
                        #print 'low reponse, discarded',filename
                        #print mean(firingbinsmeanList[1][starti:])
                        #print mean(firingbinsmeanList[2][starti:])

    print "Number of mitral cells accepted =",n_accept,"out of",n_total

    return chisqs, array(signal2residuals), array(signal2noises), \
        array(AandB_signal2residuals), array(AandB_signal2noises), \
        array(resp_signal2residuals), array(resp_signal2noises)


def plot_all_stats():
    """ Plot chi-sq histogram, and residual vs noise for individual pulse fits,
    both odors pulse prediction, and respiratory prediction. """
    fig = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [ ('',0,(False,False,False,False,False),'lat inh'), \
        ('',1,(False,False,True,False,False),'no lat inh') ]
    for ploti,(dirextn,inhi,inh,inhstr) in enumerate(inh_options):

        chisqs, signal2residuals, signal2noises, \
            AandB_signal2residuals, AandB_signal2noises, \
            resp_signal2residuals, resp_signal2noises = \
                get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,\
                    inh,NONLINEAR_ORNS,NONLINEAR_TYPE,dirextn)

        ax1 = fig.add_subplot(2,4,4*ploti+1,frameon=False)
        ax2 = fig.add_subplot(2,4,4*ploti+2,frameon=False)
        ax3 = fig.add_subplot(2,4,4*ploti+3,frameon=False)
        ax4 = fig.add_subplot(2,4,4*ploti+4,frameon=False)
        signal2residuals = array(signal2residuals)
        #ax1.hist(chisqs,20,histtype='step',linewidth=linewidth,label=inhstr,color='k')
        ax2.scatter(signal2noises,signal2residuals,s=linewidth,\
            label=inhstr,color='k',marker='.',alpha=0.7)
        ax3.scatter(AandB_signal2noises,AandB_signal2residuals,s=linewidth,\
            label=inhstr,color='k',marker='.',alpha=0.7)
        ax4.scatter(resp_signal2noises,resp_signal2residuals,s=linewidth,\
            label=inhstr,color='k',marker='.',alpha=0.7)
        
        ## beautify plots
        for axnum,ax in enumerate([ax1,ax2,ax3,ax4]):
            panel_label = ['A','B','C','D','E','F','G','H'][ploti*4+axnum]
            ## ax.transAxes ensures relative to axes size, rather than to data units.
            text(0.15, 1.0, panel_label, fontsize=label_fontsize, transform=ax.transAxes)
            ax.get_yaxis().set_ticks_position('left')
            ax.get_xaxis().set_ticks_position('bottom')
            xmin, xmax = ax.get_xaxis().get_view_interval()
            ymin, ymax = ax.get_yaxis().get_view_interval()
            ax.set_xlim([0,xmax])
            ax.set_ylim([0,ymax])
            ax.set_xticks([0,xmax])
            ax.set_yticks([0,ymax])
            ax.add_artist(Line2D((0, 0), (0, ymax), color='black', linewidth=axes_linewidth))
            ax.add_artist(Line2D((0, xmax), (0, 0), color='black', linewidth=axes_linewidth))
            minax = min(xmax,ymax)
            if ax!=ax1:
                ax.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)
            ## axes_labels() sets sizes of tick labels too.
            if ploti==1:
                if ax==ax1:
                    axes_labels(ax,'chi-sq','',adjustpos=False)
                else:
                    axes_labels(ax,'','',adjustpos=False)
            else:
                axes_labels(ax,'','',adjustpos=False)
    fig.tight_layout()
    subplots_adjust(top=0.95,wspace=1.0)
    fig.text(0.01,0.8,'count',fontsize=label_fontsize, rotation='vertical', transform=fig.transFigure)
    fig.text(0.2,0.92,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
            rotation='vertical', transform=fig.transFigure)
    fig.text(0.5,0.02,'$\sqrt{signal/noise}$',fontsize=label_fontsize, transform=fig.transFigure)
    fig.savefig('../figures/linpulses_stats.svg', bbox_inches='tight',dpi=fig.dpi)
    fig.savefig('../figures/linpulses_stats.png', bbox_inches='tight',dpi=fig.dpi)

######## PAPER SUPPLEMENTARY FIGURE
def plot_1x2x_stats_papersupplfigure():
    """ Plot residual vs noise for individual pulse fits,
    both odors pulse prediction, and respiratory prediction. """
    
    fig = figure(figsize=(columnwidth*2/3.,linfig_height*2/3),dpi=300,facecolor='w')
    hist_bins = arange(0.0,2.01,0.05)

    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,3)
    ax3 = fig.add_subplot(2,2,2)
    ax4 = fig.add_subplot(2,2,4)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    ## inh_options = [(dirextn,inhi,inh,inhstr,(axfits,axpreds),resp_fits_on),...]
    inh_options = [ \
        ('',0,(False,False,False,False,False),'all',(ax1,ax2),False),
        ('_aug15_2x',1,(False,False,False,False,False),'all',(ax3,ax4),False) ]
    yR2Nmax = 0.
    for (dirextn,inhi,inh,inhstr,(axfits,axpreds),resp_fits_on) in inh_options:

        ## sqrt() already taken for S2Rs and S2Ns -- see above get_pulse_goodnessfits()
        chisqs, signal2residuals, signal2noises, \
            AandB_signal2residuals, AandB_signal2noises, \
            resp_signal2residuals, resp_signal2noises = \
                get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,inh,\
                False,'P',dirextn,resp_fits=resp_fits_on,dirextn4stim=False)

        ## alpha a in the marker's facecolor (r,g,b,a) doesn't work
        axfits.hist(clip(signal2noises/signal2residuals,0,2),\
            hist_bins,normed=True,edgecolor='b',facecolor='b')
        _,y1=axfits.get_ylim()
        axpreds.hist(clip(AandB_signal2noises/AandB_signal2residuals,0,2),\
            hist_bins,normed=True,edgecolor='b',facecolor='b')
        _,y2=axpreds.get_ylim()
        yR2Nmax = max((yR2Nmax,y1,y2))
        print "Mean of signal2noises =",mean(array(signal2noises)**2)
        print "Mean of signal2residuals =",mean(array(signal2residuals)**2)
        print "Mean of A+B signal2noises =",mean(array(AandB_signal2noises)**2)
        print "Mean of A+B signal2residuals =",mean(array(AandB_signal2residuals)**2)

    ## beautify plots
    for axnum,ax in enumerate([ax1,ax2,ax3,ax4]):
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
        ax.set_xlim(0,2)
        ax.set_xticks([0,1,2])
        ax.set_ylim(0,yR2Nmax)
        ax.set_yticks([0,yR2Nmax])
        if axnum in [0,2]: ax.set_xticklabels(['','',''])
        else: ax.set_xticklabels(['0','1','2'])
        ## axes_labels() sets sizes of tick labels too.
        if axnum==1:
            axes_labels(ax,"$\sqrt{residual/noise}$","")
            ax.xaxis.set_label_coords(1.2,-0.35)
        elif axnum==0:
            axes_labels(ax,"","prob. density")
            ax.yaxis.set_label_coords(-0.3,-0.3)
        else: axes_labels(ax,'','',adjustpos=False) # just to set fontsize
    fig.tight_layout()
    fig_clip_off(fig)
    #fig.text(-0.02,0.8,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
    #    rotation='vertical',transform=fig.transFigure)
    fig.subplots_adjust(top=0.95,left=0.17,bottom=0.23,wspace=0.35,hspace=0.4)
    if SAVEFIG:
        fig.savefig('../figures/linpulses_stats_1x2x.svg',dpi=fig.dpi)
        fig.savefig('../figures/linpulses_stats_1x2x.png',dpi=fig.dpi)

######## PAPER FIGURE 6
def plot_mitioNL_stats_paperfigure():
    """ Plot residual vs noise for individual pulse fits,
    both odors pulse prediction, and respiratory prediction. """
    
    fig = figure(figsize=(columnwidth*2,linfig_height*2/3),dpi=300,facecolor='w')
    hist_bins = arange(0.0,2.01,0.05)

    ax1 = fig.add_subplot(2,5,1)
    ax2 = fig.add_subplot(2,5,6)
    ax3 = fig.add_subplot(2,5,2)
    ax4 = fig.add_subplot(2,5,7)
    ax5 = fig.add_subplot(2,5,3)
    ax6 = fig.add_subplot(2,5,8)
    ax7 = fig.add_subplot(2,5,4)
    ax8 = fig.add_subplot(2,5,9)
    ax9 = fig.add_subplot(2,5,5)
    ax10 = fig.add_subplot(2,5,10)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    ## inh_options = [(dirextn,inhi,inh,inhstr,(axfits,axpreds),resp_fits_on),...]
    inh_options = [ \
        ('',0,(False,False,False,False,False),'all',(ax1,ax2),False),
        ('',0,(False,False,False,True,False),'all',(ax3,ax4),False),
        ('',0,(True,True,False,False,False),'all',(ax5,ax6),False),
        ('_aug17_2_PGmod',0,(False,False,False,False,False),'all',(ax7,ax8),False),
        ('_aug17_4_PGmod',1,(False,False,False,False,False),'all',(ax9,ax10),False) ]
    maxy1 = 0
    maxy2 = 0
    for (dirextn,inhi,inh,inhstr,(axfits,axpreds),resp_fits_on) in inh_options:

        ## sqrt() already taken for S2Rs and S2Ns -- see above get_pulse_goodnessfits()
        chisqs, signal2residuals, signal2noises, \
            AandB_signal2residuals, AandB_signal2noises, \
            resp_signal2residuals, resp_signal2noises = \
                get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,inh,\
                False,'P',dirextn,resp_fits=resp_fits_on,dirextn4stim=False)

        ## alpha a in the marker's facecolor (r,g,b,a) doesn't work
        axfits.hist(clip(signal2noises/signal2residuals,0,2),\
            hist_bins,normed=True,edgecolor='b',facecolor='b')
        axpreds.hist(clip(AandB_signal2noises/AandB_signal2residuals,0,2),\
            hist_bins,normed=True,edgecolor='b',facecolor='b')
        print "Mean of signal2noises =",mean(array(signal2noises)**2)
        print "Mean of signal2residuals =",mean(array(signal2residuals)**2)
        print "Mean of A+B signal2noises =",mean(array(AandB_signal2noises)**2)
        print "Mean of A+B signal2residuals =",mean(array(AandB_signal2residuals)**2)

    ## beautify plots
    for axnum,ax in enumerate([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10]):
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
        ## find the ymax across all plots in top row and in bottom row
        if axnum % 2 == 0: maxy1 = max(maxy1,ymax)
        else: maxy2 = max(maxy2,ymax)
        ax.set_xlim(0,2)
        ax.set_xticks([0,1,2])
        if axnum in [0,2,4,6,8]: ax.set_xticklabels(['','',''])
        else: ax.set_xticklabels(['0','1','2+'])
        ## axes_labels() sets sizes of tick labels too.
        if axnum==5:
            axes_labels(ax,"$\sqrt{residual/noise}$","",xpad=1)
            #ax.xaxis.set_label_coords(1.2,-0.35)
        elif axnum==0:
            axes_labels(ax,"","prob. density")
            ax.yaxis.set_label_coords(-0.25,-0.2)
        else: axes_labels(ax,'','',adjustpos=False) # just to set fontsize
    for ax in [ax1,ax3,ax5,ax7,ax9]:
        ax.set_ylim([0,maxy1])
        ax.set_yticks([0,maxy1])
    for ax in [ax2,ax4,ax6,ax8,ax10]:
        ax.set_ylim([0,maxy2])
        ax.set_yticks([0,maxy2])
    #fig.tight_layout()
    fig_clip_off(fig)
    #fig.text(-0.02,0.8,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
    #    rotation='vertical',transform=fig.transFigure)
    fig.subplots_adjust(top=0.95,left=0.07,right=0.95,bottom=0.23,wspace=0.25,hspace=0.4)
    if SAVEFIG:
        fig.savefig('../figures/linpulses_stats_mitioNL.svg',dpi=fig.dpi)
        fig.savefig('../figures/linpulses_stats_mitioNL.png',dpi=fig.dpi)

######## PAPER FIGURE 5
def plot_all_stats_paperfigure():
    """ Plot residual vs noise for individual pulse fits,
    both odors pulse prediction, and respiratory prediction. """
    
    fig = figure(figsize=(columnwidth,linfig_height*2/3),dpi=300,facecolor='w')
    hist_bins = arange(0.0,2.01,0.05)

    import scipy.io
    ## The first two columns of Score_grand and Score_twoodor are
    ## S2N and S2R respectively, without sqrt(). So be careful to take sqrt().
    fit_goodnesses = scipy.io.loadmat('priyanka_goodness_of_fits.mat',struct_as_record=True)
    Priyanka_fit_odor = fit_goodnesses['Score_grand']
    Priyanka_pred_binary = fit_goodnesses['Score_twoodor']
    R2N_odor = sqrt(Priyanka_fit_odor[:,0]/Priyanka_fit_odor[:,1])
    R2N_binary = sqrt(Priyanka_pred_binary[:,0]/Priyanka_pred_binary[:,1])

    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,4)
    ## histogram of residual/noise
    expcol = (0.55,0.55,0.22) ## olive colour from http://cloford.com/resources/colours/500col.htm
    ax1.hist(clip(R2N_odor,0,2),hist_bins,normed=True,edgecolor=expcol,facecolor=expcol)
    ax2.hist(clip(R2N_binary,0,2),hist_bins,normed=True,edgecolor=expcol,facecolor=expcol)

    ax3 = fig.add_subplot(2,3,2)
    ax4 = fig.add_subplot(2,3,5)
    ax5 = fig.add_subplot(2,3,3)
    ax6 = fig.add_subplot(2,3,6)
    #ax7 = fig.add_subplot(2,4,4)
    #ax8 = fig.add_subplot(2,4,8)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    ## inh_options = [(dirextn,inhi,inh,inhstr,(axfits,axpreds),resp_fits_on),...]
    inh_options = [ \
        ('',0,(False,False,False,False,False),'all',(ax3,ax4),True),
        ('_excinh_jun17',1,(False,False,False,False,False),'all',(ax5,ax6),False) ]
        #('',2,(False,False,True,False,False),'all',(ax7,ax8),False),
        #('_0.33x_may14',3,(False,False,True,False,False),'all',(ax7,ax8),False) ]
    for (dirextn,inhi,inh,inhstr,(axfits,axpreds),resp_fits_on) in inh_options:

        ## sqrt() already taken for S2Rs and S2Ns -- see above get_pulse_goodnessfits()
        chisqs, signal2residuals, signal2noises, \
            AandB_signal2residuals, AandB_signal2noises, \
            resp_signal2residuals, resp_signal2noises = \
                get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,inh,\
                False,'P',dirextn,resp_fits=resp_fits_on)

        ## alpha a in the marker's facecolor (r,g,b,a) doesn't work,
        ## only in edgecolor -- looks slightly weird
        if inhi==3: color = (1,0,0,0.5)
        else: color = 'b'
        axfits.hist(clip(signal2noises/signal2residuals,0,2),\
            hist_bins,normed=True,edgecolor=color,facecolor=color)
        axpreds.hist(clip(AandB_signal2noises/AandB_signal2residuals,0,2),\
            hist_bins,normed=True,edgecolor=color,facecolor=color)
        print "Mean of signal2noises =",mean(array(signal2noises)**2)
        print "Mean of signal2residuals =",mean(array(signal2residuals)**2)
        print "Mean of A+B signal2noises =",mean(array(AandB_signal2noises)**2)
        print "Mean of A+B signal2residuals =",mean(array(AandB_signal2residuals)**2)

        ## respiratory / freely breathing predictions
        if inhi==0:
            fig_resp = figure(figsize=(columnwidth/3.0,linfig_height/2.0),dpi=300,facecolor='w')
            ax = fig_resp.add_subplot(111)
            ax.hist(clip(resp_signal2noises/resp_signal2residuals,0,2),\
                hist_bins,normed=True,edgecolor='b',facecolor='b')
            xmin,xmax,ymin,ymax = \
                beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
            ## axes_labels() sets sizes of tick labels too.
            axes_labels(ax,'$\sqrt{residual/noise}$','prob. density',adjustpos=False,xpad=1,ypad=4)
            ## set from Priyanka's resp fits in Priyanka_expmnt_resp_fits(),
            ## we need same ymax to compare
            ymax = 3.5
            ax.set_xlim([0,2])
            ax.set_ylim([0,ymax])
            ax.set_xticks([0,1,2])
            ax.set_yticks([0,ymax])
            fig_clip_off(fig_resp)
            fig_resp.tight_layout()
            fig_resp.subplots_adjust(right=0.85)
            if SAVEFIG:
                fig_resp.savefig('../figures/linpulses_respstats.svg',dpi=fig_resp.dpi)
                fig_resp.savefig('../figures/linpulses_respstats.png',dpi=fig_resp.dpi)

    ## beautify plots
    maxy1 = 0.0
    maxy2 = 0.0
    for axnum,ax in enumerate([ax1,ax2,ax3,ax4,ax5,ax6]):
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
        if axnum % 2 == 0: maxy1 = max(maxy1,ymax) # get max of all plots
        else: maxy2 = max(maxy2,ymax) # get max of all plots
    for axnum,ax in enumerate([ax1,ax2,ax3,ax4,ax5,ax6]):
        if axnum % 2 == 0:
            ax.set_ylim([0,maxy1])
            ax.set_yticks([0,maxy1])
        else:
            ax.set_ylim([0,maxy2])
            ax.set_yticks([0,maxy2])
        ax.set_xlim(0,2)
        ax.set_xticks([0,1,2])
        if axnum in [0,2,4,6]: ax.set_xticklabels(['','',''])
        else: ax.set_xticklabels(['0','1','2'])
        ## axes_labels() sets sizes of tick labels too.
        if axnum==3:
            axes_labels(ax,"$\sqrt{residual/noise}$","",xpad=2)
            #ax.xaxis.set_label_coords(1.25,-0.35)
        elif axnum==0:
            axes_labels(ax,"","prob. density")
            ax.yaxis.set_label_coords(-0.2,-0.3)
        else: axes_labels(ax,'','',adjustpos=False) # just to set fontsize
    fig.tight_layout()
    fig_clip_off(fig)
    #fig.text(-0.02,0.8,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
    #    rotation='vertical',transform=fig.transFigure)
    fig.subplots_adjust(top=0.95,left=0.1,bottom=0.25,wspace=0.35,hspace=0.4)
    if SAVEFIG:
        fig.savefig('../figures/linpulses_stats.svg',dpi=fig.dpi)
        fig.savefig('../figures/linpulses_stats.png',dpi=fig.dpi)    

    return

    ################# Below this is the older method of plotting S2R vs S2N
    
    fig = figure(figsize=(columnwidth/2.0,linfig_height*2/3),dpi=300,facecolor='w')
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    (dirextn,inhi,inh,inhstr) = \
        ('',0,(False,False,False,False,False),'all')

    ## sqrt() already taken for S2Rs and S2Ns -- see above get_pulse_goodnessfits()
    chisqs, signal2residuals, signal2noises, \
        AandB_signal2residuals, AandB_signal2noises, \
        resp_signal2residuals, resp_signal2noises = \
            get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,inh,\
            False,'P',dirextn,resp_fits=False)

    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    signal2residuals = array(signal2residuals)
    ## alpha a in the marker's facecolor (r,g,b,a) doesn't work
    ax1.scatter(signal2noises,signal2residuals,s=marker_size,edgecolor='k',\
        label=inhstr,marker='o',alpha=0.7,facecolor=(0.7,0.7,0.7,0.3),linewidth=linewidth)
    ax2.scatter(AandB_signal2noises,AandB_signal2residuals,s=marker_size,edgecolor='k',\
        label=inhstr,marker='o',alpha=0.7,facecolor=(0.7,0.7,0.7,0.3),linewidth=linewidth)
    print "Mean of signal2noises =",mean(array(signal2noises)**2)
    print "Mean of signal2residuals =",mean(array(signal2residuals)**2)
    print "Mean of A+B signal2noises =",mean(array(AandB_signal2noises)**2)
    print "Mean of A+B signal2residuals =",mean(array(AandB_signal2residuals)**2)

    xmax_list = []
    ymax_list = []
    ## beautify plots
    for axnum,ax in enumerate([ax1,ax2]):
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
        xmax_list.append(xmax)
        ymax_list.append(ymax)
    xmax = max(xmax_list)
    ymax = max(ymax_list)
    for axnum,ax in enumerate([ax1,ax2]):
        ax.set_xlim([0,xmax])
        ax.set_ylim([0,ymax])
        if axnum==0: ax.set_yticks([0,ymax])
        else: ax.set_yticks([])
        ax.set_xticks([0,xmax])
        minax = min(xmax,ymax)
        ax.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)
        ## axes_labels() sets sizes of tick labels too.
        if axnum==0:
            axes_labels(ax,"$\sqrt{signal/noise}$","$\sqrt{signal/residual}$",ypad=-4)
            ax.xaxis.set_label_coords(1.25,-0.15)
        else: axes_labels(ax,'','',adjustpos=False) # just to set fontsize
    fig.tight_layout()
    fig_clip_off(fig)
    subplots_adjust(left=0.3)
    #fig.text(-0.02,0.8,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
    #    rotation='vertical',transform=fig.transFigure)
    fig.subplots_adjust(top=0.95,left=0.2,bottom=0.25,wspace=0.35,hspace=0.35)
    fig.savefig('../figures/linpulses_stats.svg',dpi=fig.dpi)
    fig.savefig('../figures/linpulses_stats.png',dpi=fig.dpi)
    
    return

    ## separate figure for respiratory fits -- older method of plotting S2R vs S2N
    fig = figure(figsize=(columnwidth/3.0,linfig_height/2.0),dpi=300,facecolor='w') # 'none' is transparent
    ax = fig.add_subplot(111)
    ax.scatter(resp_signal2noises,resp_signal2residuals,s=marker_size,edgecolor='k',\
    label=inhstr,marker='o',alpha=0.7,facecolor=(0.7,0.7,0.7,0.3),linewidth=linewidth)
    xmin,xmax,ymin,ymax = \
        beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
    ## axes_labels() sets sizes of tick labels too.
    axes_labels(ax,'$\sqrt{signal/noise}$','$\sqrt{signal/residual}$',adjustpos=False)
    ax.set_xlim([0,xmax])
    ax.set_ylim([0,ymax])
    ax.set_xticks([0,xmax])
    ax.set_yticks([0,ymax])
    minax = min(xmax,ymax)
    ax.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)
    fig.tight_layout()
    subplots_adjust(top=0.84)
    fig_clip_off(fig)
    fig.savefig('../figures/linpulses_respstats.svg',dpi=fig.dpi)
    fig.savefig('../figures/linpulses_respstats.png',dpi=fig.dpi)

def plot_lin_contribs_paperfigure():
    """ Plot residual vs noise for summation A+B predictions.
    Then plot residual vs noise of fits and predictions together
    for no lat, no singles+no PGs, nonlinear-primary.
    """
    fig1 = figure(figsize=(columnwidth,linfig_height*1.2),dpi=300,facecolor='w')
    ax1 = plt.subplot2grid((3,3),(0,1)) # full
    ax2 = plt.subplot2grid((3,3),(1,1)) # no lat
    ax3 = plt.subplot2grid((3,3),(1,2)) # no lat, 0.5x
    ax4 = plt.subplot2grid((3,3),(2,1)) # no self
    ax5 = plt.subplot2grid((3,3),(0,2)) # non-lin
    #ax6 = plt.subplot2grid((2,3),(1,2))
    ## inh = (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [
        ('',0,(False,False,False,False,False),'all',False,None,ax1), \
        ('',1,(False,False,True,False,False),'no lat-inh',False,None,ax2), \
        ('_mar12_0.5xcentralfrate',\
                2,(False,False,True,False,False),'no lat-inh',False,None,ax3), \
        ('',3,(True,False,False,True,False),'no self-inh',False,None,ax4), \
        ('',4,(False,False,False,False,False),'non-lin ORNs',True,'P',ax5) ]
    R2Ns_AB_all = [ [] for i in range(len(inh_options)) ]
    xmaxlist = []
    ymaxlist = []
    for ploti,(dirextn,inhi,inh,inhstr,nl_orns,nl_type,ax) in enumerate(inh_options):

        chisqs, signal2residuals, signal2noises, \
            AandB_signal2residuals, AandB_signal2noises, \
            resp_signal2residuals, resp_signal2noises = \
                get_pulse_goodnessfits(stim_seeds,net_seeds,num_gloms_list,inh,\
                    nl_orns,nl_type,dirextn,resp_fits=False)

        signal2residuals = array(signal2residuals)
        ax.scatter(AandB_signal2noises,AandB_signal2residuals,s=marker_size,edgecolor='k',\
            label=inhstr,marker='o',alpha=0.7,facecolor=(0.7,0.7,0.7,0.3),linewidth=linewidth)
        ## beautify plots
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
        xmaxlist.append(xmax)
        ymaxlist.append(ymax)

    xmax = max(xmaxlist)
    ymax = max(ymaxlist)
    axlist = [ax1,ax2,ax3,ax4,ax5]
    for axnum,ax in enumerate(axlist):
        ax.set_xlim([0,xmax])
        if axnum in [2,3]: ax.set_xticks([0,xmax])
        else: ax.set_xticks([])
        ax.set_ylim([0,ymax])
        if axnum in [0,1,3]: ax.set_yticks([0,ymax])
        else: ax.set_yticks([])
        minax = min(xmax,ymax)
        ax.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)
        #panel_label = ['a','b','c','d','e','f'][axnum]
        ### ax.transAxes ensures relative to axes size, rather than to data units.
        #ax.text(0.15, 1.0, panel_label, fontsize=label_fontsize, transform=ax.transAxes)

    fig1.tight_layout()
    fig_clip_off(fig1)
    fig1.subplots_adjust(top=0.95,left=0.13,bottom=0.17,wspace=0.4,hspace=0.4)
    fig1.text(0.31,0.65,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
            rotation='vertical', transform=fig1.transFigure)
    fig1.text(0.6,0.025,'$\sqrt{signal/noise}$',fontsize=label_fontsize,transform=fig1.transFigure)
    fig1.savefig('../figures/lin_contribs.svg',dpi=fig1.dpi)
    fig1.savefig('../figures/lin_contribs.png',dpi=fig1.dpi)

def plot_lin_contribs():
    """ Plot residual vs noise for individual pulse fits and summation A+B predictions.
    Then plot residual vs noise of fits and predictions together
    for no lat, no grans, no singles, no PGs, nonlinear-primary.
    Plot residual/noise of fits and predictions vs firing rate. """
    fig1 = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
    ax1 = plt.subplot2grid((2,3),(0,0),frameon=False)
    ax2 = plt.subplot2grid((2,3),(1,0),frameon=False)
    ax3 = plt.subplot2grid((2,3),(0,1),frameon=False)
    ax4 = plt.subplot2grid((2,3),(1,1),frameon=False)
    ax5 = plt.subplot2grid((2,3),(0,2),frameon=False)
    ax6 = plt.subplot2grid((2,3),(1,2),frameon=False)
    fig2 = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
    ax7 = plt.subplot2grid((1,4),(0,0),frameon=False)
    ax8 = plt.subplot2grid((1,4),(0,1),frameon=False)
    ax9 = plt.subplot2grid((1,4),(0,2),frameon=False)
    ax10 = plt.subplot2grid((1,4),(0,3),frameon=False)
    fig3 = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
    fig4 = figure(figsize=(columnwidth,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
    ax11 = plt.subplot2grid((1,1),(0,0),frameon=False)
    #plot_fits_vs_corr(ax1,ax2) # these 2 panels go into another figure
    ## inh = (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [
        (0,(False,False,False,False,False),'all',False,None,ax1), \
        (1,(False,False,True,False,False),'no lat',False,None,ax2), \
        (2,(True,True,False,False,False),'no granules',False,None,ax4), \
        (3,(True,False,False,False,False),'no singles',False,None,ax3), \
        (4,(False,False,False,True,False),'no PGs',False,None,ax5),
        (5,(False,False,False,False,False),'non-linear',True,'P',ax6) ]
    #R2Ns_AB_all = [[]]*len(inh_options) # doesn't work, creates copies of the same list!!
    R2Ns_AB_all = [ [] for i in range(len(inh_options)) ]
    R2Ns_resp_all = [ [] for i in range(len(inh_options)) ]
    fig3plotnum = 1
    for ploti,(inhi,inh,inhstr,nl_orns,nl_type,ax) in enumerate(inh_options):
        chisqs = []
        signal2residuals = []
        signal2noises = []
        AandB_signal2residuals = []
        AandB_signal2noises = []
        resp_signal2residuals = []
        resp_signal2noises = []
        fratemeans,R2Ns,fratemeans_AB,R2Ns_AB = [],[],[],[]
        n_accept = 0
        for stimi,stimseed in enumerate(stim_seeds):
            if not salient: net_seeds = [stimseed]
            for neti,netseed in enumerate(net_seeds):
                for ngi,num_gloms in enumerate([3]):

                    filename, switch_strs \
                        = get_filename(netseed,stimseed,inh,ngi,stimi,neti,nl_orns,nl_type)
                    switches_str = string.join(switch_strs,'')
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    #print filename
                    ## If the fitted params file does not exist, create it (them).
                    if not os.path.exists(filename+'_params0'):
                        print filename
                        ## fit the responses for this result file
                        kernels,chisqs,fitted_responses_full,mdict_full = \
                            fit_pulses(filename,stimseed,False,NOSHOW=True,SAVEFIG=False)

                    for fitted_mitral in [0,1]:
                        f = open(filename+'_params'+str(fitted_mitral),'r')
                        ## firingbinsmeanList,firingbinserrList [pulsenum][binnum] are for this mitral
                        ## fitted responses [pulsenum-1][binnum] are for this mitral (note no air pulses).
                        chisqs_mit,kernels,bgnd,firingbinsmeanList,firingbinserrList,fitted_responses \
                            = pickle.load(f)
                        f.close()
                        kernelA,kernelB = kernels
                        
                        ## Priyanka's definitions, following Geffen et al, 2009
                        ## for the last A&B pulse train index 5
                        ## only compare after odor onset at PULSE_START
                        starti = int(PULSE_START/pulserebindt)
                        ## If mean firing rate is close to zero for any odour, the fits are terrible,
                        ## even if the chi-sq comes out low! So discard low firing cells to any odour A or B.
                        discard = False
                        if mean(firingbinsmeanList[1][starti:])>1 and mean(firingbinsmeanList[2][starti:])>1 \
                                and mean(firingbinsmeanList[3][starti:])>1 and mean(firingbinsmeanList[4][starti:])>1:
                            #chisqs.append(chisqs_mit[0])
                            #chisqs.append(chisqs_mit[1])
                            R2N_A,R2N_B = 0.0,0.0
                            ## Save S2R and S2N initially, if not discarded in any pulse, keep.
                            S2R_temp = []
                            S2N_temp = []
                            for pulsenum in [1,2,3,4]:
                                signal2noise,signal2residual = residual2noise( \
                                    fitted_responses[pulsenum-1],firingbinsmeanList[pulsenum],\
                                    firingbinserrList[pulsenum],starti )
                                ## discard those with too high S2R
                                ## (usually zero frate but large air initially)
                                if signal2residual>5:
                                    print 'discarded S2R =',signal2residual,'for',filename
                                    discard = True
                                else:
                                    S2N_temp.append(signal2noise)
                                    S2R_temp.append(signal2residual)
                                    if pulsenum in [1,3]: R2N_A += signal2noise/signal2residual
                                    else:  R2N_B += signal2noise/signal2residual
                            if not discard:
                                n_accept += 1
                                ## A mitral is only accepted if both odor responses are reasonable
                                ## S2N and S2R are noted for all its pulse fits
                                signal2noises.extend(S2N_temp)
                                signal2residuals.extend(S2R_temp)
                                AandB_signal2noise,AandB_signal2residual = residual2noise( \
                                    fitted_responses[4],firingbinsmeanList[5],firingbinserrList[5],starti )
                                AandB_signal2noises.append(AandB_signal2noise)
                                AandB_signal2residuals.append(AandB_signal2residual)
                                R2N_AB = AandB_signal2noise/AandB_signal2residual
                                ## Below are for AandB summation fits, just calculated above
                                if AandB_signal2residual<0.1:
                                    print 'poor fit',filename
                                elif 1.0/R2N_AB > 2.5:
                                    print 'good fit ratio =',signal2residual/signal2noise,filename

                                ## Below are needed for plotting R2N vs mean frate
                                fratemeanA = (mean(firingbinsmeanList[1][starti:])+\
                                    mean(firingbinsmeanList[3][starti:]))/2.0
                                fratemeanB = (mean(firingbinsmeanList[2][starti:])+\
                                    mean(firingbinsmeanList[4][starti:]))/2.0
                                fratemeans.append(fratemeanA)
                                R2Ns.append(R2N_A/2.0)
                                fratemeans.append(fratemeanB)
                                R2Ns.append(R2N_B/2.0)
                                fratemeans_AB.append(mean(firingbinsmeanList[5][starti:]))
                                ## R2Ns_AB_all stores nan-s also (see below),
                                ## whereas R2Ns_AB does not (same size as fratemeans),
                                R2Ns_AB.append(R2N_AB)
                                R2Ns_AB_all[inhi].append(R2N_AB)

                                ## Below are for plotting goodness of fit vs linear contribs
                                if inhi!=5: ## non-linear NLP is not yet done for morphs
                                    signal2noiseA,signal2residualA,signal2noiseB,signal2residualB,\
                                        morph_mitral_responses_avg,morph_mitral_responses_std,\
                                        morph_responseA,morph_responseB =\
                                        resp_SRN(filename,fitted_mitral,kernelA,kernelB)
                                    resp_signal2noises.append(signal2noiseA)
                                    resp_signal2noises.append(signal2noiseB)
                                    resp_signal2residuals.append(signal2residualA)
                                    resp_signal2residuals.append(signal2residualB)
                                    R2N_resp = ( signal2noiseA/signal2residualA + \
                                            signal2noiseB/signal2residualB ) / 2.0
                                    R2Ns_resp_all[inhi].append( R2N_resp )

                                ## plot the good/bad predictions (conditions below) for the default network
                                if inhi==0 and (R2N_AB>1.0 or R2N_resp>1.0):
                                    ax_AB = fig3.add_subplot(10,10,fig3plotnum,frameon=False)
                                    ax_AB.plot(firingbinsmeanList[5],'-k',linewidth=0.25)
                                    ax_AB.plot(fitted_responses[4],'-r',linewidth=0.25)
                                    beautify_plot(ax_AB)
                                    ax_AB.set_xticklabels([])
                                    ax_AB.set_yticklabels([])
                                    #ax_AB.set_title('%1.2f'%(R2N_AB),fontsize=4)
                                    ## ax.transAxes ensures relative to axes size, rather than to data units.
                                    ax_AB.text(0.15, 1.0, '%1.2f'%(R2N_AB), fontsize=3, transform=ax_AB.transAxes)
                                    fig3plotnum += 1
                                    
                                    ax_AB = fig3.add_subplot(10,10,fig3plotnum,frameon=False)
                                    airmean = morph_mitral_responses_avg[6]
                                    ax_AB.plot(morph_mitral_responses_avg[5]-airmean,'-r',linewidth=0.25)
                                    ax_AB.plot(morph_responseA,'-m',linewidth=0.25)
                                    ax_AB.plot(morph_mitral_responses_avg[0]-airmean,'-b',linewidth=0.25)
                                    ax_AB.plot(morph_responseB,'-c',linewidth=0.25)
                                    beautify_plot(ax_AB)
                                    ax_AB.set_xticklabels([])
                                    ax_AB.set_yticklabels([])
                                    #ax_AB.set_title('%1.2f'%(R2N_resp),fontsize=4)
                                    ## ax.transAxes ensures relative to axes size, rather than to data units.
                                    ax_AB.text(0.15, 1.0, '%1.2f'%(R2N_resp), fontsize=3, transform=ax_AB.transAxes)
                                    fig3plotnum += 1
                                    
                                    print 'R2N_AB = %1.2f, R2N_respA = %1.2f,%1.2f, R2N_respB = %1.2f,%1.2f'\
                                        %(R2N_AB,signal2noiseA,signal2residualA,signal2noiseB,signal2residualB),\
                                        filename

                            else: # discarded as signal2residual is too high
                                R2Ns_AB_all[inhi].append(nan)
                                R2Ns_resp_all[inhi].append(nan)
                        else: # discarded as mean firing rate(s) are too low.
                            R2Ns_AB_all[inhi].append(nan)
                            R2Ns_resp_all[inhi].append(nan)
                            #print 'low reponse, discarded',filename
                            #print mean(firingbinsmeanList[1][starti:])
                            #print mean(firingbinsmeanList[2][starti:])

        signal2residuals = array(signal2residuals)
        ax.scatter(signal2noises,signal2residuals,\
            label=inhstr,color=(0.3,0.3,0.3),marker='o',s=marker_size,\
            facecolors='none',linewidth=linewidth)
        ax.scatter(AandB_signal2noises,AandB_signal2residuals,\
            label=inhstr,color='k',marker='x',s=marker_size,linewidth=linewidth)
        print "Number of mitral cells accepted =",n_accept
        
        if inhi==0:
            ## plot for R2N vs frate
            ax7.scatter(fratemeans,R2Ns,\
                label=inhstr,color=(0.3,0.3,0.3),marker='o',s=marker_size,\
                facecolors='none',linewidth=linewidth)
            ax7.scatter(fratemeans_AB,R2Ns_AB,\
                label=inhstr,color='k',marker='x',s=marker_size,linewidth=linewidth)
            beautify_plot(ax7)
            axes_labels(ax7,'firing rate (Hz)','$\sqrt{residual/noise}$',adjustpos=False)
            ## plot for R2N_AB vs R2N_resp
            ax8.scatter(R2Ns_AB_all[inhi],R2Ns_resp_all[inhi],\
                label=inhstr,color='k',marker='x',s=marker_size,linewidth=linewidth)
            xmin,xmax,ymin,ymax = beautify_plot(ax8)
            minax = min(xmax,ymax)
            ax8.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)
            axes_labels(ax8,'$\sqrt{R/N} for AB$','',adjustpos=False)
            
            ## plot S2R vs S2N for default network
            ax11.scatter(resp_signal2residuals,resp_signal2noises,\
                label=inhstr,color='k',marker='x',s=marker_size,\
                linewidth=linewidth)            
            axes_labels(ax11,'$\sqrt{signal/residual}','$\sqrt{signal/noise}$',adjustpos=False)
            xmin,xmax,ymin,ymax = beautify_plot(ax11)
            minax = min(xmax,ymax)
            ax11.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)

        ## beautify plots
        xmin,xmax,ymin,ymax = beautify_plot(ax)
        minax = min(xmax,ymax)
        ax.plot([0,minax],[0,minax],linestyle='dashed',color='k',alpha=0.5,linewidth=linewidth)

    ## plot for R2N (for AB and resp) vs lin tweaks
    for R2Ns_all,ax in ((R2Ns_AB_all,ax9),(R2Ns_resp_all,ax10)):
        maxR2N = max(R2Ns_all[0])
        for lincontrib_line in zip(*R2Ns_all[0:3]):
            color = (lincontrib_line[0]/maxR2N,0,1-lincontrib_line[0]/maxR2N)
            nanincolor = [ isnan(col) for col in color ]
            if True not in nanincolor:
                if lincontrib_line[0]>1.0: color='r'
                else: color='b'
                ax.plot(lincontrib_line,'-x',ms=marker_size,color=color,linewidth=linewidth)
        beautify_plot(ax)
        ax.set_xticks(range(3))
        ax.set_xticklabels(['all','no lat','no granules'])
        fig2.autofmt_xdate()
        axes_labels(ax,'','',adjustpos=False)

    axlists = [ [ax1,ax2,ax3,ax4,ax5,ax6] , [ax7,ax8,ax9,ax10] ]
    for axlist in axlists:
        for axnum,ax in enumerate(axlist):
            panel_label = ['A','B','C','D','E','F'][axnum]
            ## ax.transAxes ensures relative to axes size, rather than to data units.
            ax.text(0.15, 1.0, panel_label, fontsize=label_fontsize, transform = ax.transAxes)

    fig1.tight_layout()
    fig1.subplots_adjust(top=0.95,left=0.13,bottom=0.17,wspace=0.4,hspace=0.4)
    fig1.text(0.01,0.7,'$\sqrt{signal/residual}$',fontsize=label_fontsize,\
            rotation='vertical', transform=fig1.transFigure)
    fig1.text(0.45,0.025,'$\sqrt{signal/noise}$',fontsize=label_fontsize,transform=fig1.transFigure)
    fig1.savefig('../figures/lin_contribs.svg', bbox_inches='tight',dpi=fig1.dpi)
    fig1.savefig('../figures/lin_contribs.png', bbox_inches='tight',dpi=fig1.dpi)

    fig2.tight_layout()
    fig2.subplots_adjust(top=0.95,wspace=0.4)
    fig2.text(0.3,0.7,'$\sqrt{R/N} for resp$',fontsize=label_fontsize,\
            rotation='vertical', transform=fig2.transFigure)
    fig2.savefig('../figures/lin_mech.svg', bbox_inches='tight',dpi=fig2.dpi)
    fig2.savefig('../figures/lin_mech.png', bbox_inches='tight',dpi=fig2.dpi)

    fig3.tight_layout()
    
    fig4.tight_layout()

def plot_fits_vs_corr(ax1,ax2):
    """ Plot avg-chisq vs correlation for each sisters pair - odor combination.
    Plot residual2noise vs correlation for each sisters pair - odor pair combination. """
    ## ax1 and ax2 are now passed in for another figure
    #fig = figure(figsize=(columnwidth/3.0,linfig_height),dpi=300,facecolor='w') # 'none' is transparent
    #ax1 = fig.add_subplot(2,1,1,frameon=False)
    #ax2 = fig.add_subplot(2,1,2,frameon=False)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [ (0,(False,False,False,False,False),'lat inh') ]
    for ploti,(inhi,inh,inhstr) in enumerate(inh_options):
        corrs = []
        corrsA = []
        corrsB = []
        chisqsA = []
        chisqsB = []
        residuals2noises = []
        residuals2noisesAandB = []
        n_accept_chisq = 0
        n_accept_ratio = 0
        for stimi,stimseed in enumerate(stim_seeds):
            if not salient: net_seeds = [stimseed]
            for neti,netseed in enumerate(net_seeds):
                for ngi,num_gloms in enumerate(num_gloms_list):

                    filename, switch_strs \
                        = get_filename(netseed,stimseed,inh,ngi,stimi,neti)
                    switches_str = string.join(switch_strs,'')
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(filename): continue
                    #print filename

                    ## morphs results file to calculate correlation
                    corr_filename, corr_switch_strs \
                        = corr_utils.get_filename(netseed,stimseed,inh,num_gloms,\
                        None,None,None,directed,FRAC_DIRECTED)
                    ## if the result file for these seeds & tweaks doesn't exist,
                    ## then carry on to the next.
                    if not os.path.exists(corr_filename): continue

                    ## calc phase corr-s for printing on the title
                    (air_corr,odorA_corr,odorB_corr),\
                        (tcorrlist,airxcorrgram,odorAxcorrgram,odorBxcorrgram),\
                        Dfrates = \
                            calc_corrs(corr_filename, norm_str="overall", \
                            numbins=corr_utils.NUMBINS, bin_width_time=corr_utils.BIN_WIDTH_TIME,
                            printinfo=False)

                    ratio = 0
                    ratioAandB = 0
                    chisqA,chisqB = 0,0
                    for fitted_mitral in [0,1]:
                        f = open(filename+'_params'+str(fitted_mitral),'r')
                        ## firingbinsmeanList,firingbinserrList [pulsenum][binnum] are for this mitral
                        ## fitted responses [pulsenum-1][binnum] are for this mitral (note no air pulses).
                        chisqs_mit,kernels,bgnd,firingbinsmeanList,firingbinserrList,fitted_responses \
                            = pickle.load(f)
                        f.close()
                        
                        ## Priyanka's definitions, following Geffen et al, 2009
                        ## for the last A&B pulse train index 5
                        ## only compare after odor onset at PULSE_START
                        starti = int(PULSE_START/pulserebindt)
                        ## If mean firing rate is close to zero for any odour, the fits are terrible,
                        ## even if the chi-sq comes out low! So discard low firing cells to any odour A or B.
                        if mean(firingbinsmeanList[1][starti:])>1 and mean(firingbinsmeanList[2][starti:])>1 \
                                and mean(firingbinsmeanList[3][starti:])>1 and mean(firingbinsmeanList[4][starti:])>1:
                            for pulsenum in [1,2,3,4]:
                                signal2noise,signal2residual = residual2noise( \
                                    fitted_responses[pulsenum-1],firingbinsmeanList[pulsenum],\
                                    firingbinserrList[pulsenum],starti )
                                ratio += signal2noise/signal2residual
                            signal2noise,signal2residual = residual2noise( \
                                fitted_responses[4],firingbinsmeanList[5],firingbinserrList[5],starti )
                            ratioAandB += signal2noise/signal2residual
                        else:
                            ratio += nan
                            ratioAandB += nan
                        ## chisqs
                        if mean(firingbinsmeanList[1][starti:])>0.6: chisqA += chisqs_mit[0]
                        else: chisqA += nan
                        if mean(firingbinsmeanList[2][starti:])>0.6: chisqB += chisqs_mit[0]
                        else: chisqB += nan

                    if not isnan(ratio) and not isnan(odorA_corr) and not isnan(odorB_corr):
                        n_accept_ratio += 1
                        corrs.append(mean(odorA_corr,odorB_corr)) ## mean correlation
                        residuals2noises.append(ratio/8.0) ## 4 pulses per sister
                        residuals2noisesAandB.append(ratioAandB/2.0) ## 1 pulse per sister
                    if not isnan(chisqA) and not isnan(odorA_corr):
                        n_accept_chisq += 1
                        corrsA.append(odorA_corr)
                        chisqsA.append(chisqA/2.0) # chisq of odor A, mean over 2 sister mitrals
                    if not isnan(chisqB) and not isnan(odorB_corr):
                        n_accept_chisq += 1
                        corrsB.append(odorB_corr)
                        chisqsB.append(chisqB/2.0) # chisq of odor B, mean over 2 sister mitrals

        ## for plotting fit vs corr, odorA vs odorB need not be distinguished
        ## set clip_on=False, to allow points that are just on the edge to not get half-cut.
        ax1.scatter(corrsA,chisqsA,marker='.',s=marker_size,color='k',alpha=0.7,\
            linewidth=linewidth,clip_on=False)
        ax1.scatter(corrsB,chisqsB,marker='.',s=marker_size,color='k',alpha=0.7,\
            linewidth=linewidth,clip_on=False)
        print "Number of sister-pair--odor combos accepted for chisq =",n_accept_chisq
        ax2.scatter(corrs,residuals2noises,marker='.',s=marker_size,color='k',alpha=0.7,\
            linewidth=linewidth,clip_on=False)
        ax2.scatter(corrs,residuals2noisesAandB,marker='x',s=marker_size,color='k',alpha=0.7,\
            linewidth=linewidth,clip_on=False)
        print "Number of sister-pair--odor-pair combos accepted for residual2noise =",n_accept_ratio
    xmin, xmax = ax1.get_xaxis().get_view_interval()
    ymin, ymax = ax1.get_yaxis().get_view_interval()
    ax1.set_xlim(-1,1)
    ax1.set_xticks([-1,0,1])
    ax1.set_xticklabels(['','',''])
    ax1.set_ylim(0,ymax)
    ax1.set_yticks([0,ymax])
    ## turn on/off the side axes (after turning off the frame above):
    ## http://www.shocksolution.com/2011/08/removing-an-axis-or-both-axes-from-a-matplotlib-plot/
    ax1.get_xaxis().set_ticks_position('bottom')
    ax1.get_yaxis().set_ticks_position('left')
    ax1.add_artist(Line2D((-1, -1), (0, ymax), color='black', linewidth=linewidth))
    ax1.add_artist(Line2D((-1, 1), (0, 0), color='black', linewidth=linewidth))
    axes_labels(ax1,'','chi-sq',fontsize=label_fontsize)

    xmin, xmax = ax2.get_xaxis().get_view_interval()
    ymin, ymax = ax2.get_yaxis().get_view_interval()
    ax2.get_xaxis().set_ticks_position('bottom')
    ax2.get_yaxis().set_ticks_position('left')
    ax2.set_xlim(-1,1)
    ax2.set_xticks([-1,0,1])
    ax2.set_yticks([0,1,ymax])
    ax2.add_artist(Line2D((-1, -1), (0, ymax), color='black', linewidth=linewidth))
    ax2.add_artist(Line2D((-1, 1), (0, 0), color='black', linewidth=linewidth))
    axes_labels(ax2,'correlation','$\sqrt{residual/noise}$',fontsize=label_fontsize)
    ## now called as part of larger figure, hence not needed
    #fig.tight_layout()
    #fig.savefig('../figures/fits_vs_corr.svg', bbox_inches='tight',dpi=fig.dpi)
    #fig.savefig('../figures/fits_vs_corr.png', bbox_inches='tight',dpi=fig.dpi)

def plot_kernelscorr_vs_conc_papersupplfigure():
    #import mpl_toolkits.axisartist as AA
    fig_kc = figure(figsize=(columnwidth/3.0,linfig_height*2/3.),dpi=300,facecolor='w')
    ### This was to use new_fixed_axis(), but I couldn't change its tick properties.
    #ax1 = AA.Subplot(fig_kc,2,1,1)
    #fig_kc.add_subplot(ax1)
    #ax2 = AA.Subplot(fig_kc,2,1,2)
    #fig_kc.add_subplot(ax2)
    ax1 = fig_kc.add_subplot(2,1,1)
    ax2 = fig_kc.add_subplot(2,1,2)
    ## inh =  (no_singles,no_joints,no_lat,no_PGs,varyRMP)
    inh_options = [
        ('',0,(False,False,False,False,False),'1x conc',False,None), \
        ('_aug15_2x',1,(False,False,False,False,False),'2x_conc',False,None) ]
    kernel_corrs_1x2x = []
    kernel_corrs_AvsB = []
    for fitted_mitral in [0,1]:
        for stimi,stimseed in enumerate(stim_seeds):
            if not salient: net_seeds = [stimseed]
            for neti,netseed in enumerate(net_seeds):
                for ngi,num_gloms in enumerate(num_gloms_list):
                    ## inh = (no_singles,no_joints,no_lat,no_PGs,varyRMP)
                    kernelAs = []
                    kernelBs = []
                    for ploti,(dirextn,inhi,inh,inhstr,nl_orns,nl_type) in enumerate(inh_options):
                        filename, switch_strs \
                            = get_filename(netseed,stimseed,inh,ngi,stimi,neti,nl_orns,nl_type,\
                                    '../results/odor_pulses'+dirextn)
                        ## if the result file for these seeds & tweaks doesn't exist,
                        ## then carry on to the next.
                        if not os.path.exists(filename): continue
                        ## If the fitted params file does not exist, create it (them).
                        if not os.path.exists(filename+'_params0'):
                            ## fit the responses for this result file
                            fitter = fit_plot_pulses(filename,stimseed,False,NOSHOW=True,SAVEFIG=False)
                            kernels,chisqs,fitted_responses_full,mdict_full = \
                                fitter.fit_pulses(filename,stimseed,False,NOSHOW=True,SAVEFIG=False,dirextn=dirextn)
                        f = open(filename+'_params'+str(fitted_mitral),'r')

                        ## firingbinsmeanList,firingbinserrList [pulsenum][binnum] are for this mitral
                        ## fitted responses [pulsenum-1][binnum] are for this mitral (note no air pulses).
                        chisqs_mit,kernels,bgnd,firingbinsmeanList,firingbinserrList,fitted_responses \
                            = pickle.load(f)
                        f.close()
                        kernelAs.append(kernels[0])
                        kernelBs.append(kernels[1])
                        ## Control: corr between odor kernels A and B
                        kernel_corrs_AvsB.append(stats.pearsonr(kernelAs[-1],kernelBs[-1])[0])
                    if len(kernelAs) == 2 and len(kernelBs) == 2:
                        print "Comparing",filename
                        kernel_corrs_1x2x.append(stats.pearsonr(kernelAs[0],kernelAs[1])[0])
                        kernel_corrs_1x2x.append(stats.pearsonr(kernelBs[0],kernelBs[1])[0])

    ax1.hist(kernel_corrs_1x2x,20,(-1,1),normed=True,edgecolor='b',facecolor='b')
    _,y1 = ax1.get_ylim()
    ax2.hist(kernel_corrs_AvsB,20,(-1,1),normed=True,edgecolor='b',facecolor='b')
    _,y2 = ax2.get_ylim()
    ycorrmax = max(y1,y2)
    for ax in [ax1,ax2]:
        ##beautify_plot() doesn't work with AxisArtist
        xmin,xmax,ymin,ymax = \
            beautify_plot(ax,x0min=False,y0min=True,\
                drawyaxis=False,xticksposn='bottom',yticksposn='none')
        #ax.axis["left", "top", "right"].set_visible(False) ## for AxisArtist
        ax.set_xlim([-1,1])
        ax.set_ylim([0,ycorrmax])
        ax.set_xticks([-1,0,1])
        ax.set_yticks([])
        ### make an new axis along the second axis axis (y-axis) which passes
        ### throught x=0.
        #ax.axis["x=0"] = ax.new_floating_axis(nth_coord=1, value=0,
        #                 axis_direction="top")
        ### new_fixed_axis doesn't even use pixel coords for offset
        ### calculate pixel values for (0,0) and (-1,0) and get offset
        #offpixels = ax.transData.transform([(0,0)])[0] - \
        #            ax.transData.transform([(-1,0)])[0]
        ## CAUTION: Resizing window or modifying the figure
        ## will screw things up
        ## make new (right-side) yaxis, but wth some offset
        #ax.axis["left2"] = ax.new_fixed_axis(loc="left",
        #                    offset=(32.37,0))
        ######## Could not change ticklabelsizes for new_fixed_axis(), so draw lines!
        ax.add_artist(Line2D((0, 0), (0, ycorrmax),\
            color='black', linewidth=axes_linewidth)) # y-axis in center
        ax.add_artist(Line2D((-0.05, 0.05), (ycorrmax, ycorrmax),\
            color='black', linewidth=axes_linewidth)) # tick at the top
        ax.text(-0.2, ycorrmax-0.2, str(int(ycorrmax)), fontsize=label_fontsize, transform=ax.transData)

    ax1.set_xticklabels([])
    ## axes_labels() sets sizes of tick labels too.
    axes_labels(ax2,'kernels corr.','',adjustpos=False,fontsize=label_fontsize,xpad=2)
    #axes_labels(ax1,"","count")
    #ax1.yaxis.set_label_coords(-0.2,-0.4)    
    fig_clip_off(fig_kc)
    #fig_kc.tight_layout()
    fig_kc.subplots_adjust(top=0.95,left=0.1,bottom=0.23,wspace=0.25,hspace=0.4)
    if SAVEFIG:
        fig_kc.savefig('../figures/linpulses_kernelscorr.svg',dpi=fig_kc.dpi)
        fig_kc.savefig('../figures/linpulses_kernelscorr.png',dpi=fig_kc.dpi)

def Priyanka_expmnt_resp_fits():
    ## from Priyanka_resp_fits_xy_measurements_CALCULATED.ods measured using imagej
    Priyanka_sqrtR2Ns = [0.6434834828,0.8345060122,1.0775429321,0.9184755103,0.5617743816,\
        0.9951396655,0.9524784277,0.8871217591,1.0978010014,0.9840340849,0.9774409973,\
        0.8435189615,0.8145737706,0.7257902014,0.51483604,0.3591452397,0.4256076386,\
        0.4895118111,0.7157424823,0.7781404064,1.2538144314,1.6401905231,0.907039066]
    fig_resp = figure(figsize=(columnwidth/3.0,linfig_height/2.0),dpi=300,facecolor='w')
    ax = fig_resp.add_subplot(111)
    hist_bins = arange(0.0,2.01,0.05)
    ax.hist(clip(Priyanka_sqrtR2Ns,0,2),\
        hist_bins,normed=True,edgecolor='b',facecolor='b')
    xmin,xmax,ymin,ymax = \
        beautify_plot(ax,x0min=True,y0min=True,xticksposn='bottom',yticksposn='left')
    ## axes_labels() sets sizes of tick labels too.
    axes_labels(ax,'$\sqrt{residual/noise}$','prob. density',adjustpos=False,xpad=1,ypad=4)
    ax.set_xlim([0,2])
    ax.set_ylim([0,ymax])
    ax.set_xticks([0,1,2])
    ax.set_yticks([0,ymax])
    fig_clip_off(fig_resp)
    fig_resp.tight_layout()
    fig_resp.subplots_adjust(right=0.85)
    if SAVEFIG:
        fig_resp.savefig('../figures/Priyanka_respstats.svg',dpi=fig_resp.dpi)
        fig_resp.savefig('../figures/Priyanka_respstats.png',dpi=fig_resp.dpi)


if __name__ == "__main__":
    if 'SAVEFIG' in sys.argv: SAVEFIG = True
    else: SAVEFIG = False
    #plot_kernels(graph=False)
    ## below two functions fit the results file and
    ## create the _params files if they don't exist
    ## plot all statistics using below function:
    #plot_all_stats() # older with chi-sq, and nolat.

    ## PAPER FIGURE: plot the sqrt(S/R) vs sqrt(S/N) for fits, A+Bs, resp.: OLD
    ## PAPER FIGURE 5: Now this plots R2N histogram for fits and A+Bs.
    plot_all_stats_paperfigure()
    ## PAPER FIGURE 6: non-linear mitral i-o curve gives non-linear random pulse preds
    #plot_mitioNL_stats_paperfigure()

    ## plot goodness of fits with parts of the network removed / tweaked
    ## No longer a paper figure -- use average_odor_scaledpulses_noregress.py for Suppl. Fig.
    #plot_lin_contribs_paperfigure() # old paper figure -- now use average_odor_scaledpulses_noregress.py
    
    ## PAPER SUPPL FIGURE 2: 1x and 2x fits and preds of random pulses
    #plot_1x2x_stats_papersupplfigure()
    ## PAPER SUPPL FIGURE 2: Correlation of kernels 1x vs 2x and odor A vs odor B (control).
    #plot_kernelscorr_vs_conc_papersupplfigure()
    
    ## PAPER SUPPL FIGURE 3: resp preds of Priyanka using half-exhalation (the best ones)
    #Priyanka_expmnt_resp_fits()
    
    ## Priyanka style example plots:
    #plot_scaled_kernels()
    #plot_scaled_kernels_special()
    show()
