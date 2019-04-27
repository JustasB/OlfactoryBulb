# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: nohup python2.6 activdep_inhibition_repeats.py &> nohup_adi.out < /dev/null &
## if running multiple of these, change adi to adi1, adi2, etc above,
## change loop lists to result in largely non-overlapping runs or opposite order,
## change the ADIproc command to have adi1, adi2, etc unique str.

import os,sys
import os.path
import pickle
import subprocess
cwd = os.getcwd() # current working directory

from pylab import *
from lock_utils import *

IN_VIVO = False#True
directed = True

## For STRONG_SYNAPSES i.e differential connectivity set frac_directed to 0.01,
## else, for random / uniform connectivity, set to 0.0.
#frac_directed = 0.0 # for activity dependent inhibition, only geometric connectivity
frac_directed = 0.01#0.003#0.03#0.05 # for activity dependent inhibition, only geometric connectivity
NONLINEAR_ORNS = False
mitB_current = 1200e-12#2000e-12 # A # current in mitB to cause lat inh

if IN_VIVO:
    reverse_list = [True,False]
    mit_distance_list = [50.0, 200.0, 400.0, 600.0, 800.0]
else:
    reverse_list = [False]
    mit_distance_list = [50.0]
netseeds = [100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0]

## inh_options = [ (no_singles,no_joints,no_PGs), ... ]
inh_options = [ (False,False,False) ]

ADI = []
for reverse in reverse_list:
    ADIdists = []
    for mit_distance in mit_distance_list:
        ADIfilenames = []
        for netseed in netseeds:
            for inh in inh_options:
                files_locked = True   # files have started to be opened for this iteration
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
                gen_file.write('scaledWidth = 0.2 # s # width of scaled pulses\n')
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

                OBNet_file = '../netfiles/syn_conn_array_10000_singlesclubbed100_jointsclubbed1'\
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
                ## inh = (no_singles,no_joints,no_PGs)
                simset_file.write('NO_SINGLES = '+str(inh[0])+'\n')
                simset_file.write('## spine inhibition and singles are self-inh\n')
                simset_file.write('## toggle them on/off together\n')
                simset_file.write('NO_SPINE_INH = NO_SINGLES\n')
                simset_file.write('NO_JOINTS = '+str(inh[1])+'\n')
                simset_file.write('NO_MULTIS = NO_JOINTS\n')
                simset_file.write('NO_PGS = '+str(inh[2])+'\n')
                simset_file.write('\n')
                simset_file.write('## When testing ADI (ASYM_TEST = False),'\
                    ' fixed current in mitB to generate 80Hz. 1mM Mg++.\n')
                simset_file.write('## When testing asymmetry in inhibition (ASYM_TEST=True),'\
                    ' same currents in mitA and mitB, and 0.2mM Mg++.\n')
                simset_file.write('ASYM_TEST = False\n')
                simset_file.write('## reverse roles of mitA and mitB in activity dependent inhibition\n')
                simset_file.write('REVERSED_ADI = '+str(reverse)+'\n')
                simset_file.write('IN_VIVO = '+str(IN_VIVO)+'\n')
                simset_file.write('## tuftinput: if ODORINH, use higher inputs to ORNs, higher gran bgnd not used.\n')
                simset_file.write('ODORINH = True\n')
                simset_file.write('oninject_ext = '+str(mitB_current)+' # A \n')
                simset_file.close()

                ## activdep_inhibition.py checks if the output files exists
                ## if there is already an output file, it quits.
                ## NOSHOW is for not showing plots, adi/adi2 is uniquestr
                ## for running multiple parallel activdep_inhibition_repeats.py.
                ADIproc = subprocess.Popen('mpiexec -machinefile ~/hostfile -n 61'\
                    ' ~/Python-2.6.4/bin/python2.6 activdep_inhibition.py NOSHOW adi1',\
                    shell=True,stdout=subprocess.PIPE)
                while True:
                    next_line = ADIproc.stdout.readline()
                    if not next_line:
                        break
                    sys.stdout.write(next_line)
                    if files_locked and ('Loading' in next_line):
                        ## now that the simulation has loaded,
                        ## unlock files for the other process.
                        ## only if files are locked still,
                        ## else redundant since 'Loading' appears multiple times
                        #myunlock('locksimfile.txt')
                        portalocker.unlock(lock_file)
                        lock_file.close()
                        files_locked = False # files are closed now
                        print "UnLocked files for ADI."
                    if 'Wrote' in next_line:
                        ADIfilename = next_line.split()[1]
                        ADIfilenames.append(ADIfilename)
                        break
                print ADIproc.communicate()[0]
                ## unlock in case files are locked even after odor_morphs quits.
                if files_locked:
                    #myunlock('locksimfile.txt')
                    portalocker.unlock(lock_file)
                    lock_file.close()
                    print "UnLocked files for ADI after quit."

        ADIdists.append(ADIfilenames)
    ADI.append(ADIdists)

print ADI
if IN_VIVO: invivo_str = '_invivo'
else: invivo_str = ''
if directed: directedstr = '_directed'+str(frac_directed)
else: directedstr = ''
fullfilename = '../results/ADI/ADI'+invivo_str+directedstr+'.pickle'
fullfile = open(fullfilename,'w')
pickle.dump(ADI, fullfile)
fullfile.close()
print "Wrote",fullfilename
