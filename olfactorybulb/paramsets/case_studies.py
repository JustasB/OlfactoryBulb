from olfactorybulb.paramsets.base import SilentNetwork

class PureMCs(SilentNetwork):

    description = "Pure MC input"
    mc_input_weight = 1.0

class PureTCs(SilentNetwork):

    description = "Pure TC input"
    tc_input_weight = 1.0

class PureMCsWithGJs(SilentNetwork):

    description = "Pure MC input and enabled gap junctions"

    gap_juction_gmax = {
        "MC": 32,
        # "TC": 0,
    }

    mc_input_weight = 1.0

class PureTCsWithGJs(SilentNetwork):

    description = "Pure TC input and enabled gap junctions"

    gap_juction_gmax = {
        # "MC": 0,
        "TC": 32,
    }

    tc_input_weight = 1.0

class MCsWithGJsGCs(SilentNetwork):

    gap_juction_gmax = {
        "MC": 32,
    }

    mc_input_weight = 1.0

    tstop = 400

    synapse_properties = {
        "AmpaNmdaSyn": {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        "GabaSyn": {
            'gmax': 2,
            'tau2': 16,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class TCsWithGJsGCs(SilentNetwork):

    gap_juction_gmax = {
        "TC": 32,
    }

    tc_input_weight = 1.0

    synapse_properties = {
        "AmpaNmdaSyn": {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        "GabaSyn": {
            'gmax': 2,
            'tau2': 16,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }



class MC_TC_Combined_Base(TCsWithGJsGCs):

    gap_juction_gmax = {
        "MC": 32,
        "TC": 32,
    }

    mc_input_weight = 1.0
    tc_input_weight = 1.0

class GammaSignature(MC_TC_Combined_Base):

    sniffs = 8
    tstop = (1+sniffs) * 200

    tc_input_weight = 0.8
    mc_input_weight = 0.2
    mc_input_delay = 0

    synapse_properties = {
        "AmpaNmdaSyn": {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        "GabaSyn": {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }


# Experiments

class GammaSignature_NoInhibition(GammaSignature):

    description = "Disabling inhibition should advance MC first spike times"

    def __init__(self):
        self.synapse_properties["GabaSyn"]["gmax"] = 0

class GammaSignature_NoTCGJs(GammaSignature):

    description = "Disabling TC GJs should abolish synchronized early TC firing"

    def __init__(self):
        self.gap_juction_gmax["TC"] = 0


class GammaSignature_NoMCGJs(GammaSignature):

    description = "Disabling MC GJs should abolish synchronized late MC firing"

    def __init__(self):
        self.gap_juction_gmax["MC"] = 0


class GammaSignature_EqualTCMCInputs(GammaSignature):

    description = "Increasing MC input weight to match TCs should advance MC first spike times"

    tc_input_weight = 0.8
    mc_input_weight = 0.8



class GammaSignature_DifferentOdor(GammaSignature):

    input_odors = {
        0:   {"name": "Apple", "rel_conc": 0.1},
        200: {"name": "Coffee", "rel_conc": 0.2},
        400: {"name": "Mint", "rel_conc": 0.2},
        600: {"name": "Apple", "rel_conc": 0.2},
        800: {"name": "Mint", "rel_conc": 0.2},
        1000: {"name": "Coffee", "rel_conc": 0.2},
        1200: {"name": "Apple", "rel_conc": 0.2},
        1400: {"name": "Mint", "rel_conc": 0.2},
        1600: {"name": "Mint", "rel_conc": 0.2},
        1800: {"name": "Apple", "rel_conc": 0.2}
    }


class GammaSignature_DifferentOdorConc(GammaSignature):

    input_odors = {
        0:   {"name": "Apple", "rel_conc": 0.1},
        200: {"name": "Apple", "rel_conc": 0.05},
        400: {"name": "Apple", "rel_conc": 0.1},
        600: {"name": "Apple", "rel_conc": 0.15},
        800: {"name": "Apple", "rel_conc": 0.20},
        1000: {"name": "Apple", "rel_conc": 0.25},
        1200: {"name": "Apple", "rel_conc": 0.30},
        1400: {"name": "Apple", "rel_conc": 0.35},
        1600: {"name": "Apple", "rel_conc": 0.4},
        1800: {"name": "Apple", "rel_conc": 0.45},
    }


class OneMsTest(GammaSignature):
    description = "Test of the simulation, for build testing only"
    tstop = 1