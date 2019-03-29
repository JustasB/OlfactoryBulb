# -*- coding: utf-8 -*-
"""
This script allows frequency analysis and oscillation detection on LFP recordings, it is based on OpenElectrophy 0.2
which can be downloaded here:

http://neuralensemble.org/trac/OpenElectrophy
"""

import distutils.version
import numpy

if distutils.version.LooseVersion(numpy.__version__)<'1.4':
    print "Numpy version checked"
    from numpy import setmember1d as in1d
    
import sys
# Path to OpenElectrophy 0.2 must be provided here
sys.path=["/home/nfourcau/OpenElectrophy_svn/OpenElectrophy/trunk/"]+sys.path

try:
    import OpenElectrophy as oe
    test_oe= (distutils.version.LooseVersion(oe.__version__)<'0.3') and (distutils.version.LooseVersion(oe.__version__)>='0.2')
    assert test_oe, 'Bad version of OpenElectrophy {}'.format( OE2.__version__ )
    from OpenElectrophy import *
    print
    print "Successful import of OpenElectrophy 0.2"
    print
    oe_valid=True
except:
    print
    print "Failure to import OpenElectrophy 0.2"
    print
    oe_valid=False

def spectrum_analysis(sig,sr,burn_time=0.,plot_fig=False,verbose=False,return_full=False):
    """
    Compute the signal spectrum
    Assume sr is in Hz and burn time in second
    """
    
    # LFP analysis
    sig_mask=(1.*arange(sig.size)/sr)>burn_time
    NFFTpoints=int(.5*sr)
    pyplot.figure()
    Pxx,freqs=pyplot.psd(sig, NFFT=NFFTpoints, Fs=sr, detrend=pyplot.mlab.detrend_mean,
          window=pyplot.mlab.window_hanning, noverlap=NFFTpoints/2, pad_to=None)
    if plot_fig:
        pyplot.xlim(0,150)
    else:
        pyplot.close()
        
    if verbose: 
        print "LFP max (freq and value):"
        print freqs[Pxx.argmax()],Pxx.max()
    
    if return_full:
        return Pxx,freqs
    else:
        return freqs[Pxx.argmax()],Pxx.max()


def beta_gamma_detection(signal,sr,freq_cut=40,verbose=False,plot_fig=False,osc_th=None,burn_time=0.):
    """
    Detect beta and gamma oscillation (separation between both is given by freq_cut)
    
    Osc_th can be a float for a global threshold or a list of 2 float. In the latter case it assumed to be beta then gamma threshold
    Assume sr is in Hz and burn time in second
    """     

    if oe_valid:

        if type(osc_th)!=type([]):
            osc_th=[osc_th]

        anaSig=AnalogSignal(signal=signal-signal.mean(),sampling_rate=sr)

        # LineDetector has the same parameters as the UI:
        lineDetector = LineDetector(anaSig,
                            #scalogram
                            f_start=10.,
                            f_stop=100.,
                            # detection_zone
                            detection_zone = [ burn_time, inf , 10., 100.],
                            # threshold
                            manual_threshold = (osc_th[0]!=None),
                            abs_threshold= osc_th[0], 
                            std_relative_threshold = 3.,
                            reference_zone = [ -inf, 1, 10., 100.],
                            # clean
                            minimum_cycle_number= 2.0,
                            eliminate_simultaneous = True,
                            regroup_full_overlap = True , 
                            eliminate_partial_overlap = True,      
                            )

        if len(osc_th)==1:
            lineDetector.computeAllStep()
            
            # you want to inspect all detected oscillations:
            if verbose:
                for osci in lineDetector.list_oscillation:
                    print osci.time_start, osci.time_stop
                    print osci.freq_start, osci.freq_stop
                    print osci.time_line, osci.freq_line, osci.value_line
                
                    
            list_gamma=[]
            list_beta=[]
            for osci in lineDetector.list_oscillation:
                if osci.freq_max<freq_cut:
                    list_beta.append(osci)
                else:
                    list_gamma.append(osci)
        
        else:
            # First detect beta
            lineDetector.detection_zone=[ burn_time, inf , 10., freq_cut]
            lineDetector.computeAllStep()
            list_beta=copy(lineDetector.list_oscillation)
            # Second detect gamma
            lineDetector.detection_zone=[ burn_time, inf , freq_cut, 100.]
            lineDetector.abs_threshold=osc_th[1]
            lineDetector.computeAllStep()
            list_gamma=copy(lineDetector.list_oscillation)
            # Put back everything for plotting
            lineDetector.list_oscillation=r_[list_beta,list_gamma]
            lineDetector.detection_zone=[ burn_time, inf , 10., 100.]
            
        if 1: # This part can be useful only if 
            # (1) "eliminate_simultaneous" is set to False in lineDetector 
            # or (2) two threshold were used
            # It shortens simultaneous beta and gamma to avoid overlap
            # (this is slightly artificial)
            def recomp_osc_properties(osc):
                if osc.time_line.size>0:
                    osc.amplitude_max=float(abs(osc.value_line).max())
                    ind_max=abs(osc.value_line).argmax()
                    osc.time_max=float(osc.time_line[ind_max])
                    osc.freq_max= float(osc.freq_line[ind_max])
                    osc.time_start=float(osc.time_line[0])
                    osc.freq_start=float(osc.freq_line[0])
                    osc.time_stop=float(osc.time_line[-1])
                    osc.freq_stop=float(osc.freq_line[-1])
                return osc
            
            for gamma in list_gamma:
                for beta in list_beta:
                    if intersect1d(gamma.time_line,beta.time_line).size != 0 :
                        ind_gamma=in1d(gamma.time_line,beta.time_line)#in1d(gamma.time_line,beta.time_line)
                        ind_beta=in1d(beta.time_line,gamma.time_line)#in1d(beta.time_line,gamma.time_line)
                        compare_osc=(abs(gamma.value_line)[ind_gamma])>(abs(beta.value_line)[ind_beta])
                        #~ print "comp ",compare_osc
                        if compare_osc[0]:
                            osc1=gamma
                            ind1=ind_gamma
                            osc2=beta
                            ind2=ind_beta
                        else:
                            osc1=beta
                            ind1=ind_beta
                            osc2=gamma
                            ind2=ind_gamma
                        ind_osc1_keep=ind1.copy()
                        ind_osc1_keep=ind_osc1_keep[ind1]
                        if any(compare_osc)&any(~compare_osc):
                            cut=where(compare_osc!=compare_osc[0])[0][0]
                        elif not(compare_osc[0]):
                            cut=compare_osc.size
                        elif compare_osc[0]:
                            cut=0
                        ind_osc1_keep[cut:]=False
                        ind1_keep=r_[~(ind1[~ind1]),ind_osc1_keep]
                        ind2_keep=r_[~ind_osc1_keep,~(ind2[~ind2])]
                        osc1.time_line=osc1.time_line[ind1_keep]
                        osc1.freq_line=osc1.freq_line[ind1_keep]
                        osc1.value_line=osc1.value_line[ind1_keep]
                        osc1=recomp_osc_properties(osc1)
                        osc2.time_line=osc2.time_line[ind2_keep]
                        osc2.freq_line=osc2.freq_line[ind2_keep]
                        osc2.value_line=osc2.value_line[ind2_keep]
                        osc2=recomp_osc_properties(osc2)
            
            def f_test(osc): return osc.time_line.size>0
            if __name__=="__main__":
                list_filter=__builtins__.filter
            else:
                list_filter=__builtins__['filter']
            list_gamma=list_filter(f_test, list_gamma)
            list_beta=list_filter(f_test, list_beta)
            lineDetector.list_oscillation=list_filter(f_test, lineDetector.list_oscillation)
                    
        # for plotting PLotLineDetecor is based on matplotlib
        if plot_fig:
            fig = pyplot.figure()
            plotLineDetector = PlotLineDetector(figure = fig , 
                                        lineDetector = lineDetector,)
            plotLineDetector.reDrawAll()
            #~ pyplot.show()

                    

        return list_beta,list_gamma
        
    else: # if OpenElectrophy has not been properly imported

        print '# ### Oscillation analysis not possible ####'
        print '# ### OpenElectrophy 0.2 is required #######'
        print
        return [],[]
        
if __name__=="__main__":
    
    from scipy import *
    a=0.001*rand(10000)
    sr=1000.
    
    freq=25
    a[3000:4000]+=cos(2*pi*freq*arange(1000)/sr)
    freq=60
    a[3500:4500]+=cos(2*pi*freq*arange(1000)/sr)
    
    beta_gamma_detection(a,sr,freq_cut=40,plot_fig=True,osc_th=0.6)
    
    pyplot.show()
    
    