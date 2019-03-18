# -*- coding: utf-8 -*-
"""
Default parameters for nrnproject.

AUTHORS:

- THOMAS MCTAVISH (2010-08-01): initial version
"""
from os.path import join

tstopval = 20000 # Used multiple times below

# SIMULATION PARAMETERS
sim_var = { \
'subsim' : False, #
'run_parallel' : False,
'tstop' : tstopval,  # Duration of simulation

# NETWORK PARAMETERS
'net_spatial_len' : 1000, # Lateral spread (width in microns) of the network
'secdenlen' : 1000,  # mitral lateral dendrite length of one branch
'num_granule' : 100,  # number of granule cells
'num_mitral' : 5,     # number of mitral cells
'ring' : False, # Whether or not the network wraps as a ring network.
'plasticity' : False, # Whether or not we are using plasticity

# CONNECTION PARAMETERS
'g2m_mean': 1, 
'g2m_var' : 0,
'g2m_ranseed' : 1,
'g2m_priden_seed_offset' : 0,  # 0 means all at 0.5

# WEIGHTS
'read_weights_normalized' : True, # Interpret the weights file as going from 0 to 1.
'write_weights_normalized' : True, # Write weights 0 to 1 from maximal conductances.
# Weight snapshot file to use as input.
# This should be something like "weights.dat". If None, then no 
# weight snapshot file is used. This weights file should be inside the weight_dir.
# If weight_dir is None, then the sim output directory will be used.
'weight_dir' : join('in','weights'),
'wt_input_file' : 'weight_input.dat',  # weight_zeros.dat
# The wt_input_file will be written dynamically if the wt_clusters variable is not
# None. Each cluster is expressed in 4 parameters. The first parameter is
# the location, centered by granule cell index. The second parameter is the
# width of the column, which should be an odd number. The third parameter is 
# the maximum normalized magnitude at that location. The fourth parameter
# is a list of mitral cell ids that should be connected with this magnitude.
'wt_clusters' : [(10,5,1,[0,4]), (90,5,1,[0,4]), \
                 (30,5,1,[1,3]),(70,5,1,[1,3])], # None #
'wt_cluster_seed' : None, # If specified, this will shuffle/permute the weight file.

# Weight snapshot file to write as output.
# If empty, this will write a file
# "weights_<stim_odors>_<stim_odor_mags>.dat"
'wt_output_file' : 'weight_snapshot.dat',
'ws_start' : tstopval, # Weight snapshot start time
'ws_interval' : tstopval, # Weight snapshot interval between snapshots.

'fi_gmax' : 0.005, # Fast inhibitory (granule-to-mitral) maximal conductance
'fi_tau1' : 0.1, # Fast inhibitory rise time
'fi_tau2' : 4., # Fast inhibitory decay time
'ampanmda_gmax' : 1.5, #Excitatory drive (mitral-to-granule) magnification.

# INPUT
'spike_input_file' : 'spikeout.spk', # Spike input file for subsimulations
'odor_dir' : 'in', # Directory relative to the project directory
'odor_file' : 'M5_Odors.txt',
'stim_odor_ids' : [0, 1, 2, 3, 4],   # Odor ID vector, zero-based.
'stim_odor_mags' : [.55, .55, .55, .55, .55],  # Odor magnitude vector.
#'stim_odor_mags' : [.0, .0, .0, .0, .0],  # Odor magnitude vector.
    # Individual odor vectors are multiplied by this value.
# End times of each odor. If using a mixture, duplicate the
# end times. For example, for a mixture of two odors followed by
# another mixture of two odors, each with a duration of 20000 ms,
# the vector would be "[20000, 20000, 40000, 40000]"
'stim_odors_start' : [2, 2, 2, 2, 2],
'stim_odors_end' : [tstopval, tstopval, tstopval, tstopval, tstopval],
'stim_odors_seed' : [100, 101, 102, 103, 104],
'stim_odor_max_delay' : 15, # Uniform random between 0 and this number that is
    # the onset of activation after the start of an inhilation.

'numodors' : 5, # Number of odors in the odors file. TODO: Deprecate
# Interval for time between "breaths", uniformly distributed.
'odorfreq' : 40, # Maximum stim rate at the beginning of odor onset.
'breath_interval' : [150, 250],
# Each breath can be thought of a occuring from 0 to 2*pi. During this interval,
# the cell will be driven with random synaptic inputs or weight breath_noise_mag.
# These events will occur as a Poisson process at rate breath_noise_freq.
# The selection of events can be thought of as drawing a random number between
# [0,1] at double the frequency. Those values > 0.5 will pass, and those below will
# not. We can also use this scheme to pass values more freqently at different times
# of the phase cycle.
'breath_seed' : 1, # Seed for breaths
#'breath_noise_mags': [.06, .06, .06, .08, .08], # On top of ampanmda_gmax
'breath_noise_mags': [.0, .0, .0, .0, .0], # On top of ampanmda_gmax
'breath_noise_freqs': [40, 40, 40, 30, 30],
'breath_noise_amps' : [0.0, 0.0, 0.0, 0.0, 0.0], # 0.5 +/- these values
                        # 0.2 at the bottom (better chance of firing), an
                        # 0.8 at the top (less of a chance).
'breath_noise_func' : 'cos', # Can be 'cos' or '-cos'. Simply set breath_noise_amp
                             # to zero to make evenly distributed.
'breath_events_file' : 'breathevents.txt',

# OUTPUT PARAMETERS
'sim_path' : '', # 'sims' # Main simulation output directory path relative to
        # project path. Individual simulations go in subdirectories of here.
'spike_file_name' : 'spikeout', # Name of spike output file (without any extension)
}
