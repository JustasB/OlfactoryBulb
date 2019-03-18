
# granules parameters
# This value changes the number of granules that will be generated
# the GCL volume is approximately 1910088333.383 (um^3)
# so the number of granules filling the GCL is volume/(grid_dim^3)

grid_dim = 25 # 122166 granules
# grid_dim = 49 # 16013 granules


# synapses interval on the mitral's dendrites
mean_synapse_interval = 10.

# synapses conductances
exc_gmax = 0.1 # nS
inh_gmax = 0.005 # uS

# sniffs parameters
# random ranges for syn weights
ods_wl = .7
ods_wh = 1.3

# random sniff frequency range
ods_freql = 2
ods_freqh = 10

# stream sniff
# you must change this to change the stimulation sequence
# of tuft weights or sniffs activation times
stream_ods_shift = 1

# odorsequence
# for each odor you must add a tuple in this manner
# the possible odors name are
# 'Apple', 'Banana', 'Basil', 'Black_Pepper',
# 'Cheese', 'Chocolate', 'Cinnamon',
# 'Cloves', 'Coffee', 'Garlic', 'Ginger',
# 'Lemongrass', 'Kiwi', 'Mint', 'Onion',
# 'Oregano', 'Pear', 'Pineapple'

# ('Mint', t init, t duration, rel. conc.), (...), (...)
odor_sequence = [ ('Mint', 50, 20000, 4e-3) ]
# sim. duration
tstop = 1000.

initial_weights = '' 

# sniff interval
# None is random
# the number was constant
sniff_invl = None

