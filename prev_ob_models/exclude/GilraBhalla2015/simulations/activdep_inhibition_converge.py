# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: nohup python2.6 activdep_inhibition_converge.py &> nohup_adiconv.out < /dev/null &
## if running multiple of these, change adi to adi1, adi2, etc above,
## change loop lists to result in largely non-overlapping runs or opposite order,
## change the ADIproc command to have adi1, adi2, etc unique str.
## also have different netseed below, to ensure that the output ADI result files do not have the same name.

import os,sys
import os.path
import pickle
import subprocess
cwd = os.getcwd() # current working directory
from pylab import *
from scipy import optimize

from lock_utils import *

IN_VIVO = False
directed = False
frac_directed = 0.05 # for activity dependent inhibition, only geometric connectivity
NONLINEAR_ORNS = False
reverse = False

mit_distance = 50.0
netseed = 200.0

mitral_granule_AMPA_Gbar_init = 0.35e-9 # Siemens ## ~=8mV EPSP, near 12mV EPSP of Trombley & Shepherd 1992 JNeurosci
granule_mitral_GABA_Gbar_init = 15e-9#15e-9#2e-9 # Siemens
self_mitral_GABA_Gbar_init = 50e-12 # Siemens
mitB_current_init = 1500e-12 # A

global iternum

def chisq_ADI(params):
    mitral_granule_AMPA_Gbar = params[0]
    granule_mitral_GABA_Gbar = params[1]
    self_mitral_GABA_Gbar = params[2]
    mitB_current = params[3]
    print "Trying params mitral_granule_AMPA_Gbar, granule_mitral_GABA_Gbar,"\
        "self_mitral_GABA_Gbar, mitB_current",params
    
    files_locked = True     # files have started to be opened for this iteration
                            # made False, once simulation is loaded and files closed
    ## a special lock file to keep track of locking,
    ## since portalocker didn't work properly with multiple files
    print "Acquiring Lock for ADI."
    sys.stdout.flush()
    #mylock('locksimfile.txt','ADI\n')
    lock_file = portalock_open('locksimfile.txt')
    print "Locked files for ADI."
    sys.stdout.flush()

    gen_file = open('../generators/stimuliConstantsMinimal.py','w') # blank file created
    gen_file.write('## This file is programmatically generated.\n')
    gen_file.write('\n')
    gen_file.write('## used by generate_firerates.py\n')
    gen_file.write('stim_rate_seednum = 1000.0#441.0#212.0#191.0\n')
    gen_file.write('## used by generate_neuroml.py\n')
    gen_file.write('stim_net_seed = '+str(netseed)+'\n')
    gen_file.write('## distance between 2 mitrals for activity dependent inhibition\n')
    gen_file.write('mit_distance = '+str(mit_distance)+' # microns\n')
    gen_file.write('## use thresholded erf() on ORN firing rate?\n')
    gen_file.write('NONLINEAR_ORNS = '+str(NONLINEAR_ORNS)+'\n')
    gen_file.close()

    net_file = open('../networks/networkConstantsMinimal.py','w') # blank file created    
    net_file.write('## actual number of modelled gloms could be 10 (for odor testing)\n')
    net_file.write('## or 2 (for inhibition testing) decided during neuroml generation.\n')
    net_file.write('## can set number of modelled glom to whatever you like.\n')
    net_file.write('## Randomly half of them will lie on central glom\'s mit0 or mit1.\n')
    net_file.write('## First half will receive odor A. Rest will receive odor B.\n')
    net_file.write('NUM_GLOMS = 2\n')
    net_file.write('\n')
    net_file.write('## Whether FRAC_DIRECTED of mits_per_syns will be\n')
    net_file.write('## connected between pairs listed in DIRECTED_CONNS.\n')
    net_file.write('## Keep directed True for simulating odors,\n')
    net_file.write('## Even for ADI, choose two connected mitrals.\n')
    net_file.write('directed = '+str(directed)+'\n')
    net_file.write('\n')
    net_file.write('## ensure that FRAC_DIRECTED * num of mitrals directed < 1.\n')
    net_file.write('## For NUM_GLOMS=10, 20mits all connected to mit0, FRAC_DIRECTED < 0.05.\n')
    net_file.write('## Can set FRAC_DIRECTED to 0.0 keeping DIRECTED=True. This will ensure that\n')
    net_file.write('## other mits lat dends are over directed centralmit\'s soma, if PROXIMAL_CONNECTION = True\n')
    net_file.write('frac_directed = '+str(frac_directed))
    net_file.write(' # I think you need to set this to 0.05 to get reasonable phase separation?\n')
    net_file.close()

    net_file = open('../synapses/synapseConstantsMinimal.py','w') # blank file created    
    net_file.write('## This file is programmatically generated for converging to best fit Activity Dependent Inhibition curve.\n\n')
    net_file.write('mitral_granule_AMPA_Gbar = '+str(mitral_granule_AMPA_Gbar)+' # Siemens\n')
    net_file.write('granule_mitral_GABA_Gbar =  '+str(granule_mitral_GABA_Gbar)+'# Siemens\n')
    net_file.write('self_mitral_GABA_Gbar = '+str(self_mitral_GABA_Gbar)+' # Siemens\n')
    net_file.close()

    OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed20_jointsclubbed1'\
        '_numgloms2_seed'+str(netseed)+"_mitdist"+str(mit_distance)
    if directed: OBNet_file += '_directed'+str(frac_directed)+'_proximal'
    OBNet_file += '_2GLOMS'
    if not IN_VIVO: OBNet_file += '_INVITRO.xml'
    else: OBNet_file += '.xml'
    if not os.path.exists(OBNet_file):
        print "Generating netfile",OBNet_file
        gen_command = 'python2.6 '+cwd+'/../generators/generate_neuroML.py 2GLOMS'
        if not IN_VIVO:
            gen_command += ' INVITRO'
        subprocess.check_call(gen_command,shell=True)
    else:
        print "Netfile",OBNet_file,"already exists."

    simset_file = open('simset_activinhibition_minimal.py','w') # blank file created
    simset_file.write('## This file is programmatically generated.\n')
    simset_file.write('\n')
    simset_file.write('netseedstr = "'+str(netseed)+'"\n')
    simset_file.write('mitdistance = '+str(mit_distance)+' # microns\n')
    simset_file.write('mitdistancestr = "_mitdist'+str(mit_distance)+'" # microns\n')
    simset_file.write('\n')
    simset_file.write('## When testing ADI (ASYM_TEST = False),'\
        ' fixed current in mitB to generate 80Hz. 1mM Mg++.\n')
    simset_file.write('## When testing asymmetry in inhibition (ASYM_TEST=True),'\
        ' same currents in mitA and mitB, and 0.2mM Mg++.\n')
    simset_file.write('ASYM_TEST = False\n')
    simset_file.write('## reverse roles of mitA and mitB in activity dependent inhibition\n')
    simset_file.write('REVERSED_ADI = '+str(reverse)+'\n')
    simset_file.write('IN_VIVO = '+str(IN_VIVO)+'\n')
    simset_file.write('oninject_ext = '+str(mitB_current)+' # A \n')
    simset_file.close()

    ## activdep_inhibition.py checks if the output files exists
    ## if there is already an output file, it quits.
    ## NOSHOW is for not showing plots, adi/adi2 is uniquestr
    ## for running multiple parallel activdep_inhibition_repeats.py.
    ADIproc = subprocess.Popen('mpiexec -machinefile ~/hostfile -n 61'\
        ' ~/Python-2.6.4/bin/python2.6 activdep_inhibition.py NOSHOW adiconv',\
        shell=True,stdout=subprocess.PIPE)
    while True:
        next_line = ADIproc.stdout.readline()
        if not next_line:
            break
        ## don't write each line!
        #sys.stdout.write(next_line)
        if files_locked and ('Loading' in next_line):
            ## now that the simulation has loaded,
            ## unlock files for the other process.
            ## only if files are locked still,
            ## else redundant since 'Loading' appears multiple times
            #myunlock('locksimfile.txt')
            portalocker.unlock(lock_file)
            lock_file.close()
            files_locked = False # files are closed now
            sys.stdout.write(next_line)
            print "UnLocked files for ADI."
        if 'Wrote' in next_line:
            sys.stdout.write(next_line)
            ADIfilename = next_line.split()[1]
            break
    ## read whatever remains of the process, but no need to print it
    ADIproc.communicate()
    ## unlock in case files are locked even after odor_morphs quits.
    if files_locked:
        #myunlock('locksimfile.txt')
        portalocker.unlock(lock_file)
        lock_file.close()
        print "UnLocked files for ADI after quit."

    f = open(ADIfilename,'r')
    Ainjectarray, both_firingratearrays = pickle.load(f)
    f.close()
    mit_alone = array(both_firingratearrays[0])
    mit_inhibited = array(both_firingratearrays[1])
    ## low end points
    low_pts = where(Ainjectarray<=0.5e-9)
    diff_frates = mit_alone[low_pts] - mit_inhibited[low_pts]
    avg_low_redux = mean(diff_frates)
    ## mid points1
    mid_pts = where((Ainjectarray>0.5e-9) & (Ainjectarray<=1.0e-9))
    diff_frates = mit_alone[mid_pts] - mit_inhibited[mid_pts]
    avg_mid1_redux = mean(diff_frates)
    ## mid points2
    mid_pts = where((Ainjectarray>1.0e-9) & (Ainjectarray<=1.5e-9))
    diff_frates = mit_alone[mid_pts] - mit_inhibited[mid_pts]
    avg_mid2_redux = mean(diff_frates)
    ## mid points3
    mid_pts = where((Ainjectarray>1.5e-9) & (Ainjectarray<=2.0e-9))
    diff_frates = mit_alone[mid_pts] - mit_inhibited[mid_pts]
    avg_mid3_redux = mean(diff_frates)
    ## high points
    high_pts = where(Ainjectarray>2.0e-9)
    diff_frates = mit_alone[high_pts] - mit_inhibited[high_pts]
    avg_high_redux = mean(diff_frates)
    ## remove the file since we'll iterate again, and ADI process will quit if resultfile exists.
    subprocess.check_call('rm '+ADIfilename,shell=True)

    global iternum
    iternum += 1
    print 'Iteration : ',iternum
    chisqarray = [avg_low_redux,avg_mid1_redux-3.0,avg_mid2_redux-10.0,avg_mid3_redux-3.0,avg_high_redux]
    print "Difference between fitted and desired inh reduxes =",chisqarray
    sys.stdout.flush()
    return chisqarray

#########
## leastsq() uses a modified version of Levenberg-Marquardt algorithm
## it optimizes M equations in N unknowns, where M>=N.
## Hence chisqfunc must return M numbers in an array
## which must be >= the number of params N in params0, 
## else: 'TypeError: Improper input parameters.'
iternum = 0
params_init = [ mitral_granule_AMPA_Gbar_init, granule_mitral_GABA_Gbar_init,\
    self_mitral_GABA_Gbar_init, mitB_current_init ]
params = optimize.leastsq( chisq_ADI, params_init,\
    full_output=1, maxfev=100000, epsfcn=1e-3 )
print params # print the status message
params = params[0] # take only fitted params; leastsq returns a tuple with errmsg, etc.
print "mitral_granule_AMPA_Gbar, granule_mitral_GABA_Gbar,"\
    "self_mitral_GABA_Gbar, mitB_current",params
print "Difference between fitted and desired inh reduxes =",chisq_ADI(params)
