# -*- coding: utf-8 -*-

import math, os
import pickle

from scipy import optimize
from pylab import *

## USAGE: python2.6 calc_differential_kernels.py

sys.path.extend(["..","../networks","../generators","../simulations"])

from networkConstants import * # has central_glom
from data_utils import *
from stimuliConstants import *

## study the main mitral that mit2 of glom1 (sideglom) connects to
## typically this is 0 i.e. mit0 of glom0 (mainglom)
fitted_mitral = DIRECTED_CONNS[2]

paramsnum = 12
STAtime = 1.0 # s

if paramsnum == 1:
    ## filelist with 5ms NOISEDT
    NOISEDT = 5e-3 # s
    NUMWHITETRAINS = 6
    rate_seednum = 191.0
    filelist = [
    '2011_06_22_20_00_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',
    None,
    None,
    '2011_06_23_16_21_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle',
    '2011_06_23_15_17_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype3.pickle',
    '2011_06_23_14_20_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype4.pickle',
    '2011_06_23_13_27_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype5.pickle'
    ]
elif paramsnum == 2:
    ## filelist with 25ms NOISEDT
    NOISEDT = 25e-3 # s
    NUMWHITETRAINS = 6
    rate_seednum = 191.0
    filelist = [
    '2011_06_23_17_49_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',
    None,
    None,
    None,
    None,
    None,
    '2011_06_23_18_16_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype5.pickle'
    ]
elif paramsnum == 3:
    ## filelist with 5ms NOISEDT
    NOISEDT = 5e-3 # s
    NUMWHITETRAINS = 240
    rate_seednum = 191.0
    filelist = [
    '2011_06_24_12_02_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',
    '2011_06_24_02_19_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype0.pickle',
    None,
    None,
    None,
    None,
    '2011_06_23_20_52_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype5.pickle'
    ]
elif paramsnum == 4:
    ## filelist with 25ms NOISEDT
    NOISEDT = 25e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 191.0
    filelist = [
    None,
    None,
    None,
    [(191.0,'2011_06_24_15_36_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    (181.0,'2011_06_25_13_44_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    (171.0,'2011_06_27_20_15_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    (161.0,'2011_06_28_00_31_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    (151.0,'2011_06_28_12_13_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle')],
    [(191.0,'2011_06_24_17_22_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype3.pickle'),
    (181.0,'2011_06_25_16_02_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype3.pickle')],
    None,
    [(191.0,'2011_06_24_13_35_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype5.pickle'),
    (181.0,'2011_06_25_19_13_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype5.pickle')]
    ]
elif paramsnum == 5:
    ## filelist with 50ms NOISEDT
    NOISEDT = 50e-3 # s
    NUMWHITETRAINS = 250
    filelist = [
    '2011_06_24_20_20_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',
    None,
    None,
    '2011_06_24_18_50_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle',
    None,
    None,
    None
    ]
elif paramsnum == 6:
    ## with extra noise sd=meanA, not sqrt(meanA)
    ## you'll need to rename the rate file
    ## filelist with 25ms NOISEDT
    NOISEDT = 25e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 191.0
    filelist = [
    '2011_06_25_02_32_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',
    None,
    '2011_06_25_00_46_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype1.pickle',
    '2011_06_24_23_25_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle',
    None,
    None,
    None
    ]
elif paramsnum == 7:
    ## usual sqrt(meanA) noise
    ## filelist with 25ms NOISEDT
    NOISEDT = 25e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 171.0

    ## without any interneurons
    #'2011_06_29_13_27_whitenoise_NOSINGLES_NOJOINTS_NOPGS_numgloms2_simtype-1.pickle'
    ## without PGs
    #'2011_06_29_12_51_whitenoise_SINGLES_JOINTS_NOPGS_numgloms2_simtype-1.pickle'
    ## all present
    #'2011_06_29_11_37_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle'
    ## only joints, 0.025 directed granules, mit2 to mit0
    #'2011_06_29_18_04_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'
    ## only joints, 0.05 directed granules, mit2 and mit3 to mit0
    #'2011_06_29_19_46_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'
    ## all present for averaging:
    #[(151.0,'2011_06_28_16_30_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    #(161.0,'2011_06_28_18_43_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    #(171.0,'2011_06_28_20_20_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle')],

    filelist = [
    '2011_06_29_13_27_whitenoise_NOSINGLES_NOJOINTS_NOPGS_numgloms2_simtype-1.pickle',
    None,
    None,
    [(151.0,'2011_06_28_16_30_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    (161.0,'2011_06_28_18_43_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle'),
    (171.0,'2011_06_28_20_20_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle')],
    '2011_06_28_22_44_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype3.pickle',
    '2011_06_29_00_36_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype4.pickle',
    None
    ]
elif paramsnum == 8:
    ## with extra noise sd=meanA, not sqrt(meanA)
    ## filelist with 100ms NOISEDT
    ## only joints, 0.05 directed granules, mit2 and mit3 to mit0, dt = 100ms
    #[(171.0,'2011_06_30_00_38_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'),
    #(181.0,'2011_06_30_12_56_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'),
    #(191.0,'2011_06_30_14_24_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle')],
    ## without interneurons, 0.05 directed granules, mit2 and mit3 to mit0, dt = 100ms
    #'2011_06_30_10_22_whitenoise_NOSINGLES_NOJOINTS_NOPGS_numgloms2_simtype2.pickle'
    NOISEDT = 100e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 181.0
    filelist = [
    None,
    None,
    None,
    [(171.0,'2011_06_30_00_38_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'),
    (181.0,'2011_06_30_12_56_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'),
    (191.0,'2011_06_30_14_24_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle')],
    None,
    None,
    None
    ]
elif paramsnum == 9:
    ## with extra noise sd=meanA, not sqrt(meanA)
    ## filelist with 50ms NOISEDT
    NOISEDT = 100e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 171.0
    filelist = [
    None,
    None,
    None,
    #[(191.0,'2011_06_30_23_58_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle'),
    #(171.0,'2011_07_01_11_42_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle')],
    #'2011_07_01_14_27_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle',
    #'2011_07_01_17_17_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle',
    [(191.0,'2011_07_02_01_03_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype2.pickle')],
    None,
    ## 50ms
    #[(191.0,'2011_07_01_23_14_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype4.pickle')],
    [(191.0,'2011_07_02_00_49_whitenoise_NOSINGLES_JOINTS_NOPGS_numgloms2_simtype4.pickle')],
    None
    ]
elif paramsnum == 10:
    ## filelist with 5ms NOISEDT
    NOISEDT = 5e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 191.0
    noisemean = 5.0
    ## with extra noise sd=meanA, not sqrt(meanA)
    noisesd = sqrt(noisemean)
    
    ## Each element of filelist could be just a filenamestr,
    ## or [(rateseed,filenamestr),...] or
    ## [(rateseed,filenamestr,noisemean,noisesd),...]
    filelist = [
    #[(191.0,'2011_07_26_18_53_whitenoise_NOSINGLES_NOJOINTS_NOPGS_numgloms2_simtype-1.pickle',2.0,sqrt(2.0)),
    #(191.0,'2011_07_26_19_25_whitenoise_NOSINGLES_NOJOINTS_NOPGS_numgloms2_simtype-1.pickle',5.0,sqrt(5.0)),
    #(191.0,'2011_07_26_19_50_whitenoise_NOSINGLES_NOJOINTS_NOPGS_numgloms2_simtype-1.pickle',8.0,sqrt(8.0))],
    [(191.0,'2011_07_27_17_23_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',2.0,sqrt(2.0)),
    (191.0,'2011_07_27_18_41_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',5.0,sqrt(5.0)),
    (191.0,'2011_07_28_00_08_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',8.0,sqrt(8.0))],
    None,
    #[(191.0,'2011_07_27_12_13_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype0.pickle',2.0,sqrt(2.0)),
    #(191.0,'2011_07_27_10_50_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype0.pickle',5.0,sqrt(5.0)),
    #(191.0,'2011_07_27_09_28_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype0.pickle',8.0,sqrt(8.0))],
    None,
    #[(191.0,'2011_07_27_13_19_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype2.pickle',2.0,sqrt(2.0))],
    None,
    None,
    None,
    None
    ]
elif paramsnum == 11:
    ## filelist with 25ms NOISEDT
    NOISEDT = 25e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 191.0
    noisemean = 5.0
    ## with extra noise sd=meanA, not sqrt(meanA)
    noisesd = sqrt(noisemean)
    
    ## Each element of filelist could be just a filenamestr,
    ## or [(rateseed,filenamestr),...] or
    ## [(rateseed,filenamestr,noisemean,noisesd),...]
    filelist = [
    [(191.0,'2011_07_28_19_56_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',2.0,sqrt(2.0)),
    (191.0,'2011_07_28_16_58_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',5.0,sqrt(5.0)),
    (191.0,'2011_07_28_12_51_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',8.0,sqrt(8.0))],
    None,
    None,
    None,
    None,
    None,
    None
    ]
elif paramsnum == 12:
    ## filelist with 5ms NOISEDT
    NOISEDT = 5e-3 # s
    NUMWHITETRAINS = 250
    rate_seednum = 191.0
    noisemean = 5.0
    ## with extra noise sd=meanA, not sqrt(meanA)
    noisesd = sqrt(0.1*noisemean/NOISEDT)
    
    ## Each element of filelist could be just a filenamestr,
    ## or [(rateseed,filenamestr),...] or
    ## [(rateseed,filenamestr,noisemean,noisesd),...]
    filelist = [
    [(191.0,'2011_08_01_20_27_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',2.0,sqrt(0.1*2.0/NOISEDT)),
    (191.0,'2011_08_02_10_36_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',5.0,sqrt(0.1*5.0/NOISEDT)),
    (191.0,'/2011_08_02_19_50_whitenoise_SINGLES_JOINTS_PGS_numgloms2_simtype-1.pickle',8.0,sqrt(0.1*8.0/NOISEDT))],
    None,
    None,
    None,
    None,
    None,
    None
    ]

STAtimelist = arange(-STAtime+NOISEDT,0.0+1e-10,NOISEDT)
kerneltimelist = arange(0.0,STAtime,NOISEDT)
lenSTA = (STAtime/NOISEDT)

def get_avg_STA(responseList,frateList):
    STAtotal = zeros(lenSTA)
    totalspikes = 0
    numtrials = 0
    ## same mean and sd of noise, but different firing rate realization
    for trainnum,frate in enumerate(frateList):
        ## same white noise firing rate realization, but different ORN spike times
        for trialnum,trialresponse in enumerate(responseList[:,trainnum]):
            numspikes,STAsum = calc_STA(frate,trialresponse,NOISEDT,STAtime)
            STAtotal += STAsum
            totalspikes += numspikes
            numtrials += 1
    ## only those response spikes are counted which occur after STAtime
    responserate = totalspikes/float(numtrials)/(PULSE_RUNTIME-STAtime)
    if totalspikes == 0:
        print "No target spikes"
        STA = STAtotal # essentially zeros if no target spikes
    else:
        STA = STAtotal/totalspikes
    return STA, responserate

def load_files(seednum, noisedt, numtrains, fn, fileidx):
    fratefile = '../generators/firerates/firerates_whitenoise_seed'\
        +str(seednum)+'_dt'+str(noisedt)+'_trains'+str(numtrains)+'.pickle'
    f = open(fratefile,'r')
    frateResponseList = pickle.load(f)
    f.close()
    print "Loaded the stimulus", fratefile
    ## frateResponseList[glomnum][trainnum]

    filenamefull = '../results/odor_whitenoise/'+fn
    f = open(filenamefull,'r')
    responseList = pickle.load(f)
    f.close()
    print "Loaded responses file",filenamefull
    ## responseList[trainnum][trialnum][mitralnum][spikenum]
    
    ## since last dimension has variable length, I need dtype=object below
    ## Now responseList[trialnum][trainnum][mitralnum] is a numpy array() of lists
    ## index as responseList[trialnum,trainnum,mitralnum][spikenum]
    responseList = array(responseList,dtype=object)

    if fileidx == 0:
        ## exc STA - input is to mainglom, study effect on mainmit
        frateList = frateResponseList[central_glom]
    else:
        ## inh STA - input is to sideglom, study effect on mainmit
        frateList = frateResponseList[central_glom+1]

    return frateList, responseList

def chisqfunc(params, STAdata_rev, ratedata):
    chisq = []
    kernel = params[0:lenSTA]
    kinteg = sum(kernel)*NOISEDT
    for i,(ratemean,ratesd,responserate) in enumerate(ratedata):
        diffs = (kernel+kinteg*ratemean**2)*(ratesd**2)/responserate - STAdata_rev[i]
        chisq.extend( [diff**2 for diff in diffs] )
    numdatapoints = len(STAdata_rev)*len(STAdata_rev[0])
    return array(chisq)/float(numdatapoints-lenSTA) # normalize to number of dof

if __name__ == "__main__":
    kernelslist = []
    fig = figure(facecolor='w') # 'none' is transparent
    ax = fig.add_subplot(111)
    fig2 = figure(facecolor='w') # 'none' is transparent
    ax2 = fig2.add_subplot(111)
    
    numfiles = len(filelist)-1
    for fileidx,filename in enumerate(filelist):
        if filename is None: continue
        if type(filename).__name__ == 'list':
            ## compute average STA for different seeds, but same mean, sd of noise
            if len(filename[0])==2:
                STAtotal = zeros(lenSTA)
                responserate_total = 0.0
                for seednum,fn in filename:
                    frateList,responseList = \
                        load_files(seednum, NOISEDT, NUMWHITETRAINS, fn, fileidx)
                    STA,responserate = get_avg_STA(responseList[:,:,fitted_mitral],frateList)
                    STAtotal += STA
                    responserate_total += responserate
                numSTAs = float(len(filename))
                STA = STAtotal/numSTAs
                responserate = responserate_total/numSTAs
                STAdata_rev = [ STA[::-1] ] # reversed STA
                ratedata = [(noisemean,noisesd,responserate)]
            ## compute average STA for diff seeds, diff means, diff sd-s of noise
            else:
                STAlist = []
                STAdata_rev = []
                ratedata = []
                for seednum,fn,meanHz,sdHz in filename:
                    frateList,responseList = \
                        load_files(seednum, NOISEDT, NUMWHITETRAINS, fn, fileidx)
                    STA,responserate = get_avg_STA(responseList[:,:,fitted_mitral],frateList)
                    STAlist.append( STA )
                    STAdata_rev.append( STA[::-1] ) # reversed STA
                    ratedata.append((meanHz,sdHz,responserate))

        else:
            frateList,responseList = \
                load_files(rate_seednum, NOISEDT, NUMWHITETRAINS, filename, fileidx)
            STA,responserate = get_avg_STA(responseList[:,:,fitted_mitral],frateList)
            STAdata_rev = [ STA[::-1] ] # reversed STA
            ratedata = [(noisemean,noisesd,responserate)]

        #### avg kernel from the reversed, calibrated STA
        kernel = zeros(lenSTA)
        for i,(noisemean,noisesd,responserate) in enumerate(ratedata):
            print i,noisemean, noisesd, responserate
            STArev = array(STAdata_rev[i])
            STArev = STArev - STArev.mean()
            thiskernel = STArev*responserate/noisesd**2
            kernel += thiskernel
            ax2.plot(kerneltimelist,thiskernel,label=str(i))
        kernel /= float(len(ratedata))

        #### fit THE KERNEL
        ## but I'm not using the fitted kernel, as it doesn't remove the offset :(
        ## initial params are an estimate of the kernel
        params0 = kernel
        ## You need more data points that params to be fitted!!
        ## That means at least two STA for different stim firing rates should be there.
        #if len(ratedata)>1:
        #    ## args is a tuple! if only one element write (elem, )
        #    params = optimize.leastsq( chisqfunc, params0,
        #        args=(STAdata_rev, ratedata), full_output=1 )
        #    print params[3]
        #    params = params[0] # leastsq returns a whole tuple of stuff - errmsg etc.
        #    kernel = params[0:lenSTA]
        kernelslist.append((fileidx,kernel))

        if fileidx == 0:
            ## exc STA - input is to mainglom, study effect on mainmit
            color = (0,0,0)
            legend = 'exc'
        else:
            ## inh STA - input is to sideglom, study effect on mainmit
            color = ( fileidx/float(numfiles), 1-fileidx/float(numfiles), 0 )
            legend = 'inh @ '+str(varied_mainrate[fileidx-1])+'Hz'
        ax.plot(STAtimelist,STA,color=color,label=legend)
        ax2.plot(kerneltimelist,kernel,color=color,label=legend)
 
    biglegend("upper left",ax)
    biglegend("upper right",ax2)
    axes_labels(ax,'time (s)','mean rate of 400 ORNs (Hz)')
    axes_labels(ax2,'time (s)','mitral kernel')
    ax.set_title('exc and inh STAs', fontsize=24)
    ax2.set_title('exc and inh kernels', fontsize=24)
       
    filename = '../results/odor_whitenoise/kernels.pickle'
    kernelsfile = open(filename,'w')
    pickle.dump(kernelslist, kernelsfile)
    kernelsfile.close()
    print "wrote",filename
    
    show()
