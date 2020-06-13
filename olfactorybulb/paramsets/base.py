import os

class SilentNetwork:

    description = "A base, silent network with mild odor input and all synapses blocked"

    @property
    def name(self):
        return self.__class__.__name__

    rnd_seed = 0

    slice_dir = os.path.join('olfactorybulb', 'slices')
    slice_name = "DorsalColumnSlice"

    sim_dt = 1 / 10.0

    recording_period = 1 / 10.0  # ms

    tstop = 800.1  # ms

    # GJs disabled
    gap_juction_gmax = {
        "MC": 0,
        "TC": 0,
    }

    # Disable excitatory and inhibitory syns
    synapse_properties = {
        "AmpaNmdaSyn": {
            'gmax': 0,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        "GabaSyn": {
            'gmax': 0,
            'tau2': 100,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

    # Mild odor input spikes
    input_odors = {
        0:   {"name": "Apple", "rel_conc": 0.1},
        200: {"name": "Apple", "rel_conc": 0.2},
        400: {"name": "Apple", "rel_conc": 0.2},
        600: {"name": "Apple", "rel_conc": 0.2},
        800: {"name": "Apple", "rel_conc": 0.2},
        1000: {"name": "Apple", "rel_conc": 0.2},
        1200: {"name": "Apple", "rel_conc": 0.2},
        1400: {"name": "Apple", "rel_conc": 0.2},
        1600: {"name": "Apple", "rel_conc": 0.2},
        1800: {"name": "Apple", "rel_conc": 0.2},
    }

    # From Manabe & Mori (2013) fast: 100-150ms, slow: 150 ms
    inhale_duration = 125  # ms

    # ORN firing rate
    max_firing_rate = 150  # Hz from Duchamp-Viret et. al. (2000)

    # Gilra Bhalla (2016)
    input_syn_tau1 = 6
    input_syn_tau2 = 12

    # MC input disabled
    mc_input_delay = 50
    mc_input_weight = 0

    # TC input disabled
    tc_input_delay = 0
    tc_input_weight = 0

    # LFP electrode
    # Inside dorsal Granule Layer
    # Approximaly to where it was located in Manabe & Mori (2013)
    # In adult male Long-Evans rat:
    # 8.0 mm anterior to the bregma
    # 1.3 mm lateral to the midline
    # 2.5 mm from the skull surface
    lfp_electrode_location = [116, 1078, -61]

    record_from_somas = ['MC', 'TC', 'GC']

class ParameterSetBase(SilentNetwork):
    pass
