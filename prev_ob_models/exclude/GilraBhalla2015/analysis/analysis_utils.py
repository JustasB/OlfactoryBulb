import pickle
import sys
sys.path.extend(["..","../networks","../generators","../simulations"])

from sim_utils import * # has rebin() to alter binsize
import generate_firerates_odors as gen_frate_odors # has respiration_function()

## minimum error that must be present.
errcut = 1e-20

def set_min_err(errListOrig, errmin=None):
    ## put in a minimum error, else divide by zero problems.
    ## find the minimum error >= errcut
    errList = array(errListOrig) # make numpy array, else where() gives an empty list
    largeerrors = errList[where(errList>errcut)]
    ## if errmin is not provided, set it as min of errors > errcut
    ## if no errors > errcut, set errmin to errcut.
    if not errmin:
        if len(largeerrors)>0: errmin = largeerrors.min()
        else: errmin = errcut
    #print "Setting minumum error to ",errmin
    ## numpy where(), replace by errmin,
    ## all those elements in firingbinsList which are less than errmin 
    errList = where(errList>errcut, errList, errmin)
    return errList

def load_morph_responses_from_pulsefilename(pulsefilename,fitted_mitral,bindt,numresps):
    ## replace 'pulses' with 'morphs'/'morph' in pulses result filename to get morphs result filename
    fn_split = pulsefilename.split('pulses')
    morph_filename = fn_split[0]+'morphs'+fn_split[1]+'morph'+fn_split[2]
    f = open(morph_filename,'r')
    (morph_mitral_responses_list,morph_mitral_responses_binned_list) = pickle.load(f)
    f.close()
    morph_numavgs = len(morph_mitral_responses_list)
    ## rebin only takes the last resp period by default
    morph_mitral_responses_binned_list = \
        rebin(morph_mitral_responses_list, numbins=int(numresps*RESPIRATION/bindt),\
            bin_width_time=bindt, numresps=numresps)
    morph_mitral_responses_avg = mean(morph_mitral_responses_binned_list, axis=0)
    morph_mitral_responses_std = std(morph_mitral_responses_binned_list, axis=0)
    ## only fitted mitral
    morph_mitral_responses_avg = morph_mitral_responses_avg[:,fitted_mitral]
    morph_mitral_responses_std = morph_mitral_responses_std[:,fitted_mitral]
    ## min err for each mitral is set separately, over all morphs
    ## NOT SETTING MIN ERR FOR PULSE PREDICTIONS - will lead to worse fits
    #morph_mitral_responses_std = set_min_err(morph_mitral_responses_std)
    return morph_numavgs,morph_mitral_responses_avg,morph_mitral_responses_std

def read_morphfile(filename,fitted_mitral,numbins):
    f = open(filename,'r')
    #### mitral_responses_list[avgnum][odornum][mitralnum][binnum]
    #### mitral_responses_binned_list[avgnum][odornum][mitralnum][binnum]
    mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
    f.close()
    morph_numavgs = len(mitral_responses_list)

    ###################### Input conditioning
    ## smooths / overlapping bins.
    ## (for non-overlapping use bin_width_time= RESPIRATION/NUM_REBINS)
    ## Adil uses non-overlapping bins
    #bin_width_time = 2.0*RESPIRATION/numbins
    bin_width_time = RESPIRATION/numbins
    ## rebin only takes the last resp period by default
    mitral_responses_binned_list = \
        rebin(mitral_responses_list, numbins, bin_width_time)
    #### very important to convert to numpy array,
    #### else where() below returns empty list.
    mitral_responses_binned_list = array(mitral_responses_binned_list)
    mitral_responses_mean = mean(mitral_responses_binned_list, axis=0)
    mitral_responses_std = std(mitral_responses_binned_list, axis=0)
    ## since I fit the mean response,
    ## I must use standard error/deviation of the _mean_
    ## = standard deviation of a repeat / sqrt(num of repeats).
    ## NO! The model predicts the individual response not the mean.
    #NUMAVGs = len(mitral_responses_binned_list)
    #mitral_responses_se= mitral_responses_std/sqrt(NUMAVGs)
    ## take the odor responses of the mitral to be fitted
    firingbinsmeanList = mitral_responses_mean[:,fitted_mitral]
    firingbinserrList = mitral_responses_std[:,fitted_mitral]

    ## min err for each mitral is set separately, over all morphs
    firingbinserrList = set_min_err(firingbinserrList)
    
    return morph_numavgs,firingbinsmeanList,firingbinserrList

def read_pulsefile(filename):
    ######### load in the mitral pulse responses
    f = open(filename,'r')
    #### mitral_responses_list[avgnum][pulsenum][mitralnum][spikenum]
    #### mitral_responses_binned_list[avgnum][pulsenum][mitralnum][binnum]
    mitral_responses_list, mitral_responses_binned_list = pickle.load(f)
    f.close()
    return mitral_responses_list, mitral_responses_binned_list

def read_pulsestimuli(stimseed,filebase='../generators/firerates'):
    ########## load in the stimuli
    f = open(filebase+'/firerates_2sgm_'+str(stimseed)+'.pickle','r')
    frateOdorList,fratePulseList,randomPulseList,\
        randomPulseStepsList,randomResponseList,genORNkernels\
        = pickle.load(f)
    #del frateOdorList, fratePulseList
    pulseList = array(randomPulseList)
    pulseStepsList = array(randomPulseStepsList)
    ## randomResponseList[glomnum][frate,...]
    ORNfrateList = array(randomResponseList)[central_glom]
    ## frateOdorList[glomnum][inputnum]
    ORNrespList = array(frateOdorList)[central_glom]
    f.close()
    return pulseList,pulseStepsList,ORNfrateList,ORNrespList,genORNkernels

def read_pulseinresp_stimuli(stimseed,filebase='../generators/firerates'):
    ########## load in the stimuli
    f = open(filebase+'/firerates_pulseinresp_'+str(stimseed)+'.pickle','r')
    pulseInRespList,kernels = pickle.load(f)
    pulseInRespList = array(pulseInRespList)
    f.close()
    return pulseInRespList,kernels

def read_scaledpulses_stimuli(stimseed,width,filebase='../generators/firerates'):
    ########## load in the stimuli
    f = open(filebase+'/firerates_scaledpulses_width'+str(width)+'_'+str(stimseed)+'.pickle','r')
    scaledPulsesList,kernels = pickle.load(f)
    scaledPulsesList = array(scaledPulsesList)
    f.close()
    return scaledPulsesList,kernels

def rebin_mean(mitral_responses_list,pulserebindt,rebin_runtime=PULSE_RUNTIME):
    ############# rebin simulated responses into pulserebindt=50ms, then take mean over trials.
    ## rebin_pulses() from sim_utils.py
    ### for overlapping pulses (also need to uncomment in sim_utils.py)
    #mitral_responses_binned_list = \
    #    rebin_pulses(mitral_responses_list, int(rebin_runtime/pulserebindt), rebin_runtime, 0.0, bin_width_time)
    mitral_responses_binned_list = \
        rebin_pulses(mitral_responses_list, int(rebin_runtime/pulserebindt), rebin_runtime, 0.0)
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

    return numtrials,mitral_responses_mean,mitral_responses_std

def decimate(pulseList,pulsedt,genORNkernels,kerneldt,kernel_size):
    ############## decimate pulseList taking every pulsedt/FIRINGFILLDT bin value
    decimatefactor = int(round(pulsedt/FIRINGFILLDT)) # ensure that these are multiples
    ## decimating the pulses, i.e. taking one value after every decimatefactor
    ## but start from middle of bin i.e. from index decimatefactor/2,
    ## as the simulated response is avg of the bin
    pulserebinList = array( [ pulseList[pulsenum][decimatefactor/2::decimatefactor] \
        for pulsenum in range(len(pulseList)) ] )

    ### can also not decimate, rather take average of each bin, fits are not too good.
    ### binning the pulses: average of every decimatefactor number of values
    #pulserebinList = array( [ array( [pulseList[pulsenum][i:i+decimatefactor].mean() \
    #        for i in range(0,len(pulseList[pulsenum]),decimatefactor)] ) \
    #    for pulsenum in range(len(pulseList)) ] )
    
    ############## decimate ORNkernels every kerneldt/FIRINGFILLDT=50 values
    ## but start from middle of bin i.e. from index decimatefactorK/2,
    ## as the simulated response is avg of the bin
    decimatefactorK = int(kerneldt/FIRINGFILLDT)
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
    linORNkernelA = array(genORNkernels[1][decimatefactorK/2::decimatefactorK][:kernel_size])\
        *MitralORNfiringratio
    if len(linORNkernelA)<kernel_size:
        linORNkernelA = append(linORNkernelA,zeros(kernel_size-len(linORNkernelA)))
    linORNkernelB = array(genORNkernels[2][decimatefactorK/2::decimatefactorK][:kernel_size])\
        *MitralORNfiringratio
    if len(linORNkernelB)<kernel_size:
        linORNkernelB = append(linORNkernelB,zeros(kernel_size-len(linORNkernelB)))

    return pulserebinList,linORNkernelA,linORNkernelB

def myconvolve(seq1,seq2):
    """ Assume we want only len(seq1) output. seq1 is pulses, seq2 is kernel.
    Assume len(seq1) > len(seq2).
    Not periodic boundaries i.e. not circular convolution.
    Assume zeros when sequence is short.
    With above assumptions, this works same as using (of course maybe 10x slower!)
    numpy.convolve(seq1,seq2,mode='full'), then taking only [0:len1] values. """
    len1 = len(seq1)
    len2 = len(seq2)
    ## using generator inside sum (no []-s) for lazy evaluation
    ## tau only goes from max(0,t-len2+1) to t, so that seq2 is not indexed beyond end
    out = [ sum( seq1[tau]*seq2[t-tau] for tau in range(max(0,t-len2+1),t) )
            for t in range(len1) ]
    return array(out)
