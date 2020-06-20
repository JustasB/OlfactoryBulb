# MOCKS for autodoc
import quantities as pq
if pq.__module__ == 'sphinx.ext.autodoc.mock':
    pq.pA = pq.nA = pq.mV = pq.ms = pq.Hz = 1
# END MOCKS



class BasePublication:
    '''
    A base class with common parameters. Subclasses of this class are used as mix-in classes to set the
    parameters of generic tests. E.g. A generic InputResistanceTest is combined with the publication Yu2015 by
    creating a specific "class InputResistanceTestYu2015(Yu2015, InputResistanceTest)"
    '''

    ss_delay = 1000 * pq.ms             # Delay before stimulus
    current_duration = 1000 * pq.ms     # Current injection time
    sag_window = 100 * pq.ms            # The value of membrane potential at this time after stimulus to use, if the voltage never drops below sag_testing_voltage
    temperature = 35                    # The temperature at which measurements were performed in the paper
    threshold_method = 'dv/dt=20'       # The definition of AP threshold 'dv/dt=20' or 'd3v/dt3'
    sag_testing_voltage = -90 * pq.mV   # The level from which sag voltage is measured
    rebound_rest_time = 1000 * pq.ms    # The amount of time to wait after stopping inhibitory current
    rebound_ap_method = "sag"           # Method to test for rebound spikes "sag" or "-300pA"
    spike_train_method = "target_freq"  # The method used to generate spike trains "target_freq" or "constant"

    required_capabilities = ()          # To meet interface requirements

class Angelo2012(BasePublication):
    ''' Angelo et. al. (2012) A biophysical signature of network affiliation and sensory processing in mitral cells '''
    current_duration = 1500 * pq.ms
    temperature = 36
    sag_testing_voltage = -103.5 * pq.mV  # LJP is not reported. Assuming uncorrected for 13.5mV, so cell was at 103.5mV

class BurtonUrban2014(BasePublication):
    ''' Burton and Urban (2014) Greater excitability and firing irregularity of tufted cells underlies distinct afferent-evoked activity of olfactory bulb mitral and tufted cells '''

    current_duration = 2000 * pq.ms
    temperature = 37
    sag_testing_voltage = -103 * pq.mV  # Paper used -90mV but with unccorected 13mV LJP so cell was at -103 mV
    ahp_time_method = 'threshold2amplitude50%'
    ahp_amplitude_method = 'threshold2minWithin10ms'
    rebound_ap_method = "sag"
    spike_train_target_freq = 20*pq.Hz

class BurtonUrban2015(BurtonUrban2014):
    ''' Burton & Urban (2015) Rapid Feedforward Inhibition and Asynchronous Excitation Regulate Granule Cell Activity in the Mammalian Main Olfactory Bulb'''
    temperature = 32


class Yu2015(BasePublication):
    ''' Yu et. al. (2015) Postnatal development attunes olfactory bulb mitral cells to high-frequency signaling '''

    current_duration = 2000 * pq.ms  # Current injection time
    sag_testing_voltage = -103.5 * pq.mV  # LJP is not corrected. Assuming 13.5mV, so cell was at 103.5mV
    sag_window = 120 * pq.ms
    threshold_method = 'd3v/dt3'
    ahp_time_method = 'threshold2min'
    ahp_amplitude_method = 'threshold2min'
    spike_train_target_freq = 30 * pq.Hz


class Hu2016(BasePublication):
    ''' Hu et. al. (2016) Hyperpolarization-Activated Currents and Subthreshold Resonance in Granule Cells of the Olfactory Bulb '''
    temperature = 32
    sag_testing_voltage = -103.5 * pq.mV  # LJP is not corrected. Assuming 13.5mV, so cell was at 103.5mV

class JohnsonDelaney2010(BasePublication):
    ''' Johnson and Delaney (2010)	Synaptic Activation of T-Type Ca2+ Channels Via mGluR Activation in the Primary Dendrite of Mitral Cells '''
    temperature = 35
    current_duration = 4000 * pq.ms
    rebound_ap_method = "-300pA"

class Zibman2011(BasePublication):
    ''' Zibman et. al. (2011)	DISTINCT INTRINSIC MEMBRANE PROPERTIES DETERMINE DIFFERENTIAL INFORMATION PROCESSING BETWEEN MAIN AND ACCESSORY OLFACTORY BULB MITRAL CELLS '''
    temperature = 22
    current_duration = 500 * pq.ms
    spike_train_method = "constant_current"
    spike_train_current = 0.100*pq.nA


class Stroh2012(BasePublication):
    ''' Stroh et. al. (2012)	NMDA Receptor-Dependent Synaptic Activation of TRPC Channels in Olfactory Bulb Granule Cells '''
    temperature = 21
    adp_current_duration = 1*pq.ms
    adp_current_amplitude = 1000*pq.pA


class Abraham2010(BasePublication):
    temperature = 24.5


class Hovis2010(BasePublication):
    temperature = 35


class Shpak2012(BasePublication):
    temperature = 22


class Christie2005(BasePublication):
    temperature = 33.5


class Fukunaga2012(BasePublication):
    temperature = 36


