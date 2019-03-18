# -*- coding: utf-8 -*-
# ALL SI UNITS
# milliMolar is same as mol/m^3

## USAGE: nohup python2.6 odor_pulses_repeats.py &> nohup_pulses.out < /dev/null &
## if running multiple of these, change pulses to pulses1, pulses2, etc above,
## change loop lists to result in largely non-overlapping runs or opposite order,
## change the morphproc command far below to have pulses1, pulses2, etc unique str.

import os,sys
import os.path
import pickle
import subprocess
cwd = os.getcwd() # current working directory
from pylab import *

from lock_utils import *

IN_VIVO = True
directed = True
frac_directed = 0.01
NONLINEAR_ORNS = False
NONLINEAR_TYPE = 'P' # P for primary glom non-linear, L for lateral gloms non-linear
## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
## in order,below options are: all cells; no lat; no joints;
## no PGs; no singles + no joints; only mitrals
inh_options = [ \
    (False,False,False,False,False), (False,False,True,False,False), (False,True,False,False,False),\
    (False,False,False,True,False), (True,True,False,False,False), (True,True,False,True,False) ]
inh_options = [ (True,True,False,False,False) ]

FORCE_SALIENT = False#True
if FORCE_SALIENT:
    stim_seed_list = [-10]#,-10,-19,-28]#[-1,-2,-3,-4,-5,-6,-7,-8,-9]
else:
    #stim_seed_list = [802.0,853.0,753.0,751.0]#arange(800.0,805.0,1.0)#[157.0,160.0,190.0,191.0,212.0,441.0]
    stim_seed_list = arange(765.0,780.0,1.0)#[157.0,160.0,190.0,191.0,212.0,441.0]
net_seed_list = [200.0]#[100.0,200.0,300.0]

print "Starting automated pulses simulations ..."
morphs = []
for num_gloms in [3]:
    morphs_stimseeds = []
    for stimseed in stim_seed_list:
        morphs_netseeds = []
        ## for the non-salient stimuli, change the network also each time
        if stimseed>0: net_seed_list = [stimseed]
        for netseed in net_seed_list:
            morphs_inh = []
            for inh in inh_options:

                files_locked = True     # files have started to be opened for this iteration
                                        # made False, once simulation is loaded and files closed
                ## a special lock file to keep track of locking,
                ## since portalocker didn't work properly with multiple files
                print "Acquiring Lock for odor pulses."
                sys.stdout.flush()
                #mylock('locksimfile.txt','pulses\n')
                lock_file = portalock_open('locksimfile.txt')
                print "Locked files for odor pulses."
                sys.stdout.flush()

                ## generate the stimuli params file
                gen_file = open('../generators/stimuliConstantsMinimal.py','w')
                gen_file.write('## This file is programmatically generated.\n')
                gen_file.write('\n')
                gen_file.write('## used by generate_firerates.py\n')
                gen_file.write('stim_rate_seednum = '+str(stimseed)+'\n')
                gen_file.write('## used by generate_neuroml.py\n')
                gen_file.write('stim_net_seed = '+str(netseed)+'\n')
                gen_file.write('## distance between 2 mitrals for activity dependent inhibition\n')
                gen_file.write('mit_distance = 0 # microns ## irrelevant for odor responses\n')
                gen_file.write('## use thresholded erf() on ORN firing rate?\n')
                gen_file.write('NONLINEAR_ORNS = '+str(NONLINEAR_ORNS)+'\n')
                gen_file.write('NONLINEAR_TYPE = "'+ NONLINEAR_TYPE \
                    +'" # P for primary glom non-linear, L for lateral gloms non-linear\n')
                gen_file.write('scaledWidth = 0.2 # s # width of scaled pulses\n')
                gen_file.close()

                ## generate the network params file
                net_file = open('../networks/networkConstantsMinimal.py','w')
                net_file.write('## actual number of modelled gloms could be 10 (for odor testing)\n')
                net_file.write('## or 2 (for inhibition testing) decided during neuroml generation.\n')
                net_file.write('## can set number of modelled glom to whatever you like.\n')
                net_file.write('## Randomly half of them will lie on central glom\'s mit0 or mit1.\n')
                net_file.write('## First half will receive odor A. Rest will receive odor B.\n')
                net_file.write('NUM_GLOMS = '+str(num_gloms)+'\n')
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

                ## generate the neuroML netfile, don't if exists
                OBNet_file = '../netfiles/syn_conn_array_10000'
                OBNet_file += '_singlesclubbed100_jointsclubbed1'+\
                    '_numgloms'+str(num_gloms)+'_seed'+str(netseed)
                if directed: OBNet_file += '_directed'+str(frac_directed)+'_proximal'
                OBNet_file += '.xml'
                if not os.path.exists(OBNet_file):
                    print "Generating netfile",OBNet_file
                    gen_command = 'python2.6 '+cwd+'/../generators/generate_neuroML.py'
                    subprocess.check_call(gen_command,shell=True)
                else:
                    print "Netfile",OBNet_file,"already exists."

                ## Generate the odor frates, don't if exists
                stim_frate_filename = '../generators/firerates/firerates_2sgm_'+str(stimseed)
                stim_frate_filename += '.pickle'
                if not os.path.exists(stim_frate_filename):
                    print "Generating firerate file",stim_frate_filename
                    gen_command = 'python2.6 '+cwd+'/../generators/generate_firerates_odors.py NOSHOW'
                    subprocess.check_call(gen_command,shell=True)                
                else:
                    print "Firerate file",stim_frate_filename,"already exists."

                ## generate the odor morph firefiles (spiketimes),
                ## don't if an exemplary firefile exists
                stim_firefiles_dirname = '../firefiles/firefiles'+str(stimseed)
                if NONLINEAR_ORNS: stim_firefiles_dirname += '_NL'+NONLINEAR_TYPE
                if not os.path.exists(stim_firefiles_dirname):
                    print "Creating firefiles directory ",stim_firefiles_dirname
                    subprocess.check_call('mkdir '+stim_firefiles_dirname,shell=True)
                if not os.path.exists(stim_firefiles_dirname+\
                        '/firetimes_rndpulse_glom_'+str(num_gloms-1)+\
                        '_pulse_0_avgnum0.txt'):
                    print "Generating firefiles in",stim_firefiles_dirname
                    ## number of processes is numpulses*num_gloms+1
                    ## numpulses is RANDOM_PULSE_NUMS*2 = 6
                    gen_command = 'mpiexec -machinefile ~/hostfile -n '+str(6*num_gloms+1)+\
                        ' ~/Python-2.6.4/bin/python2.6 '+\
                        cwd+'/../generators/generate_firefiles_odors.py NOSHOW'
                    subprocess.check_call(gen_command,shell=True)                
                else:
                    print "Exemplary firefile in",stim_firefiles_dirname,"already exists."
                
                ## no need to presently generate baseline anew: I instead need to shuffle for trials!
                ## generate baseline firefiles also anew with this stimseed*netseed combo.
                #print "Generating baseline firefiles."
                #gen_command = 'python2.6 '+cwd+'/../generators/generate_firefiles_gran_baseline.py'
                #subprocess.check_call(gen_command,shell=True)
                
                ## generate the odor simulation params file
                simset_file = open('simset_odor_minimal.py','w')
                simset_file.write('## This file is programmatically generated.\n')
                simset_file.write('\n')
                simset_file.write('netseedstr = "'+str(netseed)+'"\n')
                simset_file.write('rateseedstr = "'+str(stimseed)+'"\n')
                simset_file.write('\n')
                simset_file.write('OBNet_file = "'+OBNet_file+'"\n')
                simset_file.write('ORNpathseedstr = "'+stim_firefiles_dirname+'/"\n')
                ## inh = (no_singles,no_joints,no_lat,no_PGs,varyRMP)
                simset_file.write('NO_SINGLES = '+str(inh[0])+'\n')
                simset_file.write('## spine inhibition and singles are self-inh\n')
                simset_file.write('## toggle them on/off together\n')
                simset_file.write('NO_SPINE_INH = NO_SINGLES\n')
                simset_file.write('NO_JOINTS = '+str(inh[1])+'\n')
                simset_file.write('NO_MULTIS = NO_JOINTS\n')
                simset_file.write('NO_PGS = '+str(inh[3])+'\n')
                simset_file.write('NO_LATERAL = '+str(inh[2])+'\n')
                simset_file.write('\n')
                simset_file.write('VARY_MITS_RMP = '+str(inh[4])+'\n')
                simset_file.close()

                ## odor_pulses.py checks if the output files exists
                ## if there is already an output file, it quits.
                morphproc = subprocess.Popen('mpiexec -machinefile ~/hostfile -n 55'\
                    ' ~/Python-2.6.4/bin/python2.6 odor_pulses.py NOSHOW pulses4',\
                    shell=True,stdout=subprocess.PIPE)
                while True:
                    next_line = morphproc.stdout.readline()
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
                        print "UnLocked files for odor pulses."
                    if 'Wrote' in next_line:
                        morphfilename = next_line.split()[1]
                        morphs_inh.append(morphfilename)
                        break
                print morphproc.communicate()[0]
                ## unlock in case files are locked even after odor_morphs quits.
                if files_locked:
                    #myunlock('locksimfile.txt')
                    portalocker.unlock(lock_file)
                    lock_file.close()
                    print "UnLocked files for odor pulses after quit."
            morphs_netseeds.append(morphs_inh)
        morphs_stimseeds.append(morphs_netseeds)
    morphs.append(morphs_stimseeds)

print morphs
fullfilename = '../results/odor_pulses/pulses_random'
if NONLINEAR_ORNS: fullfilename += 'NL'+NONLINEAR_TYPE
fullfilename += '.pickle'
fullfile = open(fullfilename,'w')
pickle.dump(morphs, fullfile)
fullfile.close()
print "Wrote",fullfilename
