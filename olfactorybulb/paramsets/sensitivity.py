from olfactorybulb.paramsets.case_studies import GammaSignature

# class GammaSignature(GammaSignature):
#
#     sniffs = 8
#     tstop = (1+sniffs) * 200
#
#     tc_input_weight = 0.8
#     mc_input_weight = 0.2
#     mc_input_delay = 0
#
#     synapse_properties = {
#         "AmpaNmdaSyn": {
#             'gmax': 64,
#
#             'ltpinvl': 0,  # Disable plasticity
#             'ltdinvl': 0
#         },
#
#         "GabaSyn": {
#             'gmax': 2,
#             'tau2': 36,
#
#             'ltpinvl': 0,  # Disable plasticity
#             'ltdinvl': 0
#         }
#     }


class GammaSignature_GJ_0(GammaSignature):

    gap_juction_gmax = {
        'MC': 0,
        'TC': 0,
    }
class GammaSignature_GJ_1(GammaSignature):

    gap_juction_gmax = {
        'MC': 1,
        'TC': 1,
    }
class GammaSignature_GJ_2(GammaSignature):

    gap_juction_gmax = {
        'MC': 2,
        'TC': 2,
    }
class GammaSignature_GJ_4(GammaSignature):

    gap_juction_gmax = {
        'MC': 4,
        'TC': 4,
    }
class GammaSignature_GJ_8(GammaSignature):

    gap_juction_gmax = {
        'MC': 8,
        'TC': 8,
    }
class GammaSignature_GJ_16(GammaSignature):

    gap_juction_gmax = {
        'MC': 16,
        'TC': 16,
    }
class GammaSignature_GJ_32(GammaSignature):

    gap_juction_gmax = {
        'MC': 32,
        'TC': 32,
    }
class GammaSignature_GJ_64(GammaSignature):

    gap_juction_gmax = {
        'MC': 64,
        'TC': 64,
    }
class GammaSignature_GJ_128(GammaSignature):

    gap_juction_gmax = {
        'MC': 128,
        'TC': 128,
    }
    
    
# ------------------------------------------ #

class GammaSignature_AMPANMDA_0(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 0,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_1(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 1,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_2(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 2,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_4(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 4,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_8(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 8,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_16(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 16,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_32(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 32,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_64(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_128(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 128,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }

class GammaSignature_AMPANMDA_256(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 256,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }


# --------------------------------------------- #

class GammaSignature_GABA_0(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 0,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }


class GammaSignature_GABA_1(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 1,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }


class GammaSignature_GABA_2(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 2,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }


class GammaSignature_GABA_4(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 4,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }


class GammaSignature_GABA_8(GammaSignature):

    synapse_properties = {
        'AmpaNmdaSyn': {
            'gmax': 64,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        },

        'GabaSyn': {
            'gmax': 8,
            'tau2': 36,

            'ltpinvl': 0,  # Disable plasticity
            'ltdinvl': 0
        }
    }





# --------------------------------------------- #


class GammaSignature_TCWGHT_00(GammaSignature):
    tc_input_weight = 0



class GammaSignature_TCWGHT_02(GammaSignature):
    tc_input_weight = 0.2



class GammaSignature_TCWGHT_04(GammaSignature):
    tc_input_weight = 0.4



class GammaSignature_TCWGHT_06(GammaSignature):
    tc_input_weight = 0.6



class GammaSignature_TCWGHT_08(GammaSignature):
    tc_input_weight = 0.8



class GammaSignature_TCWGHT_10(GammaSignature):
    tc_input_weight = 1


# --------------------------------------------- #


class GammaSignature_MCWGHT_00(GammaSignature):
    mc_input_weight = 0



class GammaSignature_MCWGHT_01(GammaSignature):
    mc_input_weight = 0.1



class GammaSignature_MCWGHT_015(GammaSignature):
    mc_input_weight = 0.15



class GammaSignature_MCWGHT_02(GammaSignature):
    mc_input_weight = 0.2



class GammaSignature_MCWGHT_025(GammaSignature):
    mc_input_weight = 0.25



class GammaSignature_MCWGHT_03(GammaSignature):
    mc_input_weight = 0.3



class GammaSignature_MCWGHT_04(GammaSignature):
    mc_input_weight = 0.4



class GammaSignature_MCWGHT_06(GammaSignature):
    mc_input_weight = 0.6



class GammaSignature_MCWGHT_08(GammaSignature):
    mc_input_weight = 0.8



class GammaSignature_MCWGHT_10(GammaSignature):
    mc_input_weight = 1


