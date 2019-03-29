# -*- coding: utf-8 -*-

########## THIS FITTING PROGRAM IS MEANT TO FIT sinusoids to 'mitral responses to sinusoids'!
## USAGE: python2.6 fit_odor_morphs.py ../results/odor_morphs/2011-01-13_odormorph_SINGLES_JOINTS_PGS.pickle

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

iterationnum = 0
## amplitude[sinnum], phase[sinnum] and DC offset are the params
NUMPARAMS = 2*num_sins+1

## I don't use the NUMBINS in simset_odor.py, rather I rebin()
bindt = 5e-3 #s
NUM_REBINS = int(SIN_RUNTIME/bindt)

### numbers of mitral to be fitted.
fitted_mitral_list = [2*central_glom+0, 2*central_glom+1]

FIT_LAT_SINS = True

## those sims for which const rate at central glom and sinusoids at lateral glom
if FIT_LAT_SINS:
    filelist = [
        #(5.0,'../results/odor_sins/2012_02_14_15_36_sins_SINGLES_JOINTS_NOPGS_numgloms2.pickle') # 20 to 60 Hz
        #(5.0,'../results/odor_sins/2012_02_14_17_36_sins_SINGLES_JOINTS_NOPGS_numgloms2.pickle') # 1 to 15 Hz # 5 Hz central
        #(5.0,'../results/odor_sins/2012_02_14_19_55_sins_SINGLES_JOINTS_NOPGS_numgloms2.pickle') # 1 to 15 Hz # 9 Hz central
        #(5.0,'../results/odor_sins/2012_02_15_08_20_sins_SINGLES_JOINTS_NOPGS_numgloms2.pickle') # 1 to 15 Hz # 3 Hz central
        (3.0,'../results/odor_sins/2012_02_15_21_54_sins_SINGLES_JOINTS_NOPGS_numgloms2.pickle') # 1 to 15 Hz # 3 Hz central
    ]
## 
else:
    ## 30 trials, only mitrals
    #filelist = [
    #(1.0,'../results/odor_sins/2012_02_02_15_32_sins_NOSINGLES_NOJOINTS_NOPGS_NOLAT_numgloms2.pickle'),
    #(2.0,'../results/odor_sins/2012_02_02_17_10_sins_NOSINGLES_NOJOINTS_NOPGS_NOLAT_numgloms2.pickle'),
    #(3.0,'../results/odor_sins/2012_02_02_19_18_sins_NOSINGLES_NOJOINTS_NOPGS_NOLAT_numgloms2.pickle')
    #]

    ## 40 trials, higher frequencies, only mitrals
    filelist = [
    (1.0,'../results/odor_sins/2012_02_04_13_47_sins_NOSINGLES_NOJOINTS_NOPGS_NOLAT_numgloms2.pickle')
    ]
    ## 40 trials, higher frequencies, only mitrals, 10x longer time
    filelist = [
    (1.0,'../results/odor_sins/2012_02_04_17_51_sins_NOSINGLES_NOJOINTS_NOPGS_NOLAT_numgloms2.pickle')
    ]

    ## 30 trials, mitrals + spines + singles + PGs
    #filelist = [
    #(1.0,'../results/odor_sins/2012_02_02_16_17_sins_SINGLES_NOJOINTS_PGS_NOLAT_numgloms2.pickle'),
    #(2.0,'../results/odor_sins/2012_02_02_17_54_sins_SINGLES_NOJOINTS_PGS_NOLAT_numgloms2.pickle'),
    #(3.0,'../results/odor_sins/2012_02_02_20_08_sins_SINGLES_NOJOINTS_PGS_NOLAT_numgloms2.pickle')
    #]

    ## 40 trials, higher frequencies, mitrals + spines + singles + PGs:
    #filelist = [
    #(1.0,'../results/odor_sins/2012_02_04_14_00_sins_SINGLES_NOJOINTS_PGS_NOLAT_numgloms2.pickle')
    #]

def chisqfunc(params, mitnum, ydata, errdata):
    ampl = params[0:num_sins]
    phase = params[num_sins:2*num_sins]
    DC = params[-1]

    global iterationnum
    if iterationnum%100==0: print 'iteration number =',iterationnum
    chisqarray = [0.0]
    for sinnum,f in enumerate(sine_frequencies):
        ## Leave the first cycle of lowest frequency out for transient settling
        ## Take the first cycle after leaving above time out
        startcyclenum = 1
        startbin = int(startcyclenum/float(f)/bindt)
        ## ydata[sinnum][binnum], similar for errdata
        data = ydata[sinnum]
        error = errdata[sinnum]
        omegabindt = 2*pi*f*bindt
        for binnum in range(startbin,NUM_REBINS):
            ## ampl must be positive, sign appears via phase; phase modulo 2pi
            Rmodel = DC + abs(ampl[sinnum]) * sin( omegabindt*binnum + (phase[sinnum]%(2*pi)) )
            if Rmodel<0.0: Rmodel=0.0 # threshold if below zero
            ## divide by error to do chi-square fit
            chisqarray.append( (data[binnum] - Rmodel)/error[binnum] )
            
    ## not yet squared, so normalized 'chi' to sqrt of number of dof
    ## ydata[sinnum][binnum]
    chisqarray = array(chisqarray) / sqrt(ydata.size-NUMPARAMS)
    iterationnum += 1
    return chisqarray

def fit_sins(filename, fitted_mitral):
    f = open(filename,'r')
    mitral_responses_list = pickle.load(f)
    f.close()
    ## mitral_responses_list[avgnum][sinnum][mitnum][spikenum]

    mitral_responses_binned_list = \
        rebin_pulses(mitral_responses_list, NUM_REBINS, SIN_RUNTIME, 0.0)
    numavgs = len(mitral_responses_list)
    mitral_responses_mean = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    ## take only the responses of the mitral to be fitted
    firingbinsmeanList = mitral_responses_mean[:,fitted_mitral,:]
    firingbinserrList = mitral_responses_std[:,fitted_mitral,:]/sqrt(numavgs)
    ## amplitude of sine wave, phase shift and DC offset
    params0 = [0.0]*num_sins+[0.0]*num_sins+[0.0]
    
    ## put in a minimum error, else divide by zero problems, or NaN value params and fits
    ## find the minimum error >= errcut
    largeerrors = firingbinserrList[where(firingbinserrList>errcut)]
    if largeerrors is not (): errmin = largeerrors.min()
    else: errmin = errcut
    ## numpy where(), replace by errmin,
    ## all those elements in firingbinsList which are less than errmin 
    firingbinserrList = where(firingbinserrList>errcut, firingbinserrList, errmin)

    ###################################### Fitting
    params = optimize.leastsq( chisqfunc, params0,
        args=(fitted_mitral, firingbinsmeanList, firingbinserrList),
        full_output=1, maxfev=10000)
    print params[3]
    params = params[0] # leastsq returns a whole tuple of stuff - errmsg etc.
    print "ampl[sinnum]+phase[sinnum]+DC =",params

    ## Calculate sum of squares of the chisqarray
    chisqarraysq = [i**2 for i in 
        chisqfunc(params, fitted_mitral, firingbinsmeanList, firingbinserrList)]
    chisq = reduce(lambda x, y: x+y, chisqarraysq)

    ############################## Calculate fitted responses and return them
    
    DC_fit = params[-1]
    ampl_fit = abs(params[0:num_sins])
    phase_fit = params[num_sins:2*num_sins] % (2*pi)
    fitted_responses = [ [ \
                DC_fit + ampl_fit[sinnum] * sin( 2*pi*t*f + phase_fit[sinnum] ) \
            for t in arange(0.0, SIN_RUNTIME, bindt) ] \
        for sinnum,f in enumerate(sine_frequencies) ]

    return (params,chisq,fitted_responses,firingbinsmeanList,firingbinserrList)

if __name__ == "__main__":
    #if len(sys.argv) > 3:
        #filename = sys.argv[1]
        #ampl = float(sys.argv[2])
        #DC = float(sys.argv[3])
    #else:
        #print "Specify responses data filename, sine amplitude, DC."
        #sys.exit(1)

    for fitted_mitral in fitted_mitral_list:
        mainfig = figure(facecolor='w')
        mainax = mainfig.add_subplot(111)
        title('Mitral '+str(fitted_mitral)+' frequency response',fontsize=24)
        mainfig2 = figure(facecolor='w')
        mainax2 = mainfig2.add_subplot(111)
        title('Mitral '+str(fitted_mitral)+' phase response',fontsize=24)
        paramsList = []
        for ampl,filename in filelist:
            params,chisq,fitted_responses,firingbinsmeanList,firingbinserrList\
                = fit_sins(filename, fitted_mitral)
            print "Mit",fitted_mitral,"normalized chisq =",chisq
            paramsList.append((ampl,params))

            ################# Plot simulated and fitted responses
            if fitted_mitral != 0: continue
            for sinnum in range(num_sins):
                fig = figure(facecolor='w')
                ax = fig.add_subplot(3,1,2)
                sincolor = (sinnum+1) / float(num_sins)
                ## mean + error (lighter/whiter shade than mean below)
                ax.plot(range(NUM_REBINS),\
                    firingbinsmeanList[sinnum]+firingbinserrList[sinnum],\
                    color=(0,(1-sincolor)*0.25+0.75,sincolor*0.25+0.75),\
                    marker='+',linestyle='solid', linewidth=2)
                ## mean
                ax.plot(range(NUM_REBINS),firingbinsmeanList[sinnum],\
                    color=(0,1-sincolor,sincolor),\
                    marker='+',linestyle='solid', linewidth=2)
                ## fitted
                ax.plot(range(NUM_REBINS),fitted_responses[sinnum],\
                    color=(1,1-sincolor,sincolor),\
                    marker='x',linestyle='solid', linewidth=2)
                titlestr = 'Mitral %d response & sinusoid f=%f fit'\
                    %(fitted_mitral,sine_frequencies[sinnum])
                title(titlestr, fontsize=24)
                axes_labels(ax,'respiratory phase bin','firing rate (Hz)',adjustpos=True)
        
            ################# Plot frequency and phase responses
            mainax.plot(sine_frequencies,abs(params[0:num_sins])/float(ampl),label=str(ampl)+'Hz ORN')
            mainax2.plot(sine_frequencies,(params[0:num_sins]%(2*pi))/pi*180,label=str(ampl)+'Hz ORN')
        
        axes_labels(mainax,'input frequency (Hz)','stimulus normalized output',adjustpos=True)
        mainax.legend()
        axes_labels(mainax2,'input frequency (Hz)','output phase (degrees)',adjustpos=True)
        mainax2.legend()
    
    show()
