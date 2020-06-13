debug = False
plot = False

from matplotlib import pyplot as plt
from neuronunit import capabilities as ncap
from neuronunit.tests.base import scores
from sciunit import capabilities as scap

from olfactorybulb.neuronunit import capabilities as obncap
from olfactorybulb.neuronunit.tests import OlfactoryBulbCellTest, OlfactoryBulbCellSpikeTest
from olfactorybulb.neuronunit.tests.utilities import *


class PassivePropertyTestHelper(OlfactoryBulbCellTest):

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime)

    required_properties = [
        "temperature",
        "ss_delay",
        "current_duration",
    ]

    def generate_prediction_nocache(self, model):

        model.set_temperature(self.temperature)

        # Find the current needed to set the membrane to a negative voltage
        # Using the sag test voltage as it's also needed for later tests
        model.set_stop_time(self.ss_delay + self.current_duration)
        current = model.clamp_voltage([self.sag_testing_voltage, self.sag_testing_voltage, 0*pq.mV],
                                      [self.ss_delay, self.current_duration, 0])
        inhibitory_current = np.median(current.magnitude[-10:])*pq.nA

        # Get change in voltage in response to negative current
        model.set_stop_time(self.ss_delay + self.current_duration + self.current_duration)
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": inhibitory_current})

        if plot:
            plt.plot(voltage.times, voltage)
            plt.title(str(self))
            plt.show()

        return voltage, inhibitory_current


class RestingVoltageTest(OlfactoryBulbCellTest):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the resting voltage of a cell"
    name = "Resting potential test"

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,)

    def generate_prediction_nocache(self, model):
        voltage, _ = self.get_dependent_prediction(PassivePropertyTestHelper, model)

        # Get the voltage right before the current step
        voltage = voltage.magnitude[np.where(voltage.times < self.ss_delay)]
        resting_v = np.median(voltage[-10:]) * pq.mV

        if plot:
            plt.plot(voltage)
            plt.plot([0, len(voltage)], [resting_v, resting_v])
            plt.title(str(self) + " " + str(resting_v))
            plt.show()

        return resting_v


class InputResistanceTest(OlfactoryBulbCellTest):
    units = pq.MOhm
    score_type = scores.ZScore
    description = "A test of the input resistance of a cell"
    name = "Input resistance test"

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime)

    def generate_prediction_nocache(self, model):
        # Get the voltage just before the start of the current
        resting_v = self.get_dependent_prediction(RestingVoltageTest, model)

        # Get the voltage right before the end of current step
        voltage, test_current = self.get_dependent_prediction(PassivePropertyTestHelper, model)
        voltage = voltage.magnitude[np.where(voltage.times < self.ss_delay + self.current_duration)]
        response_v = np.median(voltage[-10:]) * pq.mV


        # V=IR -> R=V/I -> R=deltaV/deltaI
        delta_v = response_v - resting_v
        delta_i = test_current
        Rin = delta_v/delta_i
        Rin.units = pq.MOhm

        if plot:
            plt.plot(voltage)
            plt.plot([0, len(voltage)], [resting_v, resting_v])
            plt.plot([0, len(voltage)], [response_v, response_v])
            plt.title(str(self) + " " + str(Rin))
            plt.show()

        return Rin


class MembraneTimeConstantTest(OlfactoryBulbCellTest):
    units = pq.ms
    score_type = scores.ZScore
    description = "A test of membrane time constant of a cell"
    name = "Membrane time constant test"

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime)

    def generate_prediction_nocache(self, model):
        # import pydevd
        # pydevd.settrace('192.168.0.100', port=4200)

        voltage, test_current = self.get_dependent_prediction(PassivePropertyTestHelper, model)
        sampling_rate = voltage.sampling_rate

        # Get the voltage right after the end of current step
        voltage = voltage.magnitude[np.where(voltage.times > self.ss_delay + self.current_duration)]
        start_v = voltage[0]
        end_v = voltage[-1]
        delta_v = end_v - start_v

        t_above_63pct = np.where(voltage > start_v + 0.6321206*delta_v)
        tau = t_above_63pct[0][0] / sampling_rate


        if plot:
            plt.plot(voltage)
            plt.axhline(y=start_v + 0.6321206*delta_v)
            plt.axvline(x=t_above_63pct[0][0])
            plt.title(str(self) + " " + str(tau))
            plt.xlim((0, tau* sampling_rate *5))
            plt.show()

        return tau


class CellCapacitanceTest(OlfactoryBulbCellTest):
    units = pq.pF
    score_type = scores.ZScore
    description = "A test of membrane capacitance of a cell"
    name = "Membrane capacitance test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime)

    def generate_prediction_nocache(self, model):
        tau = self.get_dependent_prediction(MembraneTimeConstantTest, model)
        rin = self.get_dependent_prediction(InputResistanceTest, model)
        cap = tau / rin
        cap.units = pq.pF
        return cap

class RheobaseTest(OlfactoryBulbCellTest):
    units = pq.pA
    score_type = scores.ZScore
    description = "A test of the rheobase current of a cell"
    name = "Rheobase current test"
    max_iterations = 20
    min_precision = 0.01 # find the RB within this precision (e.g. 0.01 is 1%)

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,)

    def generate_prediction_nocache(self, model):

        resting_v = self.get_dependent_prediction(RestingVoltageTest, model)

        Rin = self.get_dependent_prediction(InputResistanceTest, model)

        # Estimate the current necessary to generate APs from passive response
        nA_to_0mV = -resting_v / Rin
        nA_to_0mV.units = pq.nA

        upper_bound = nA_to_0mV
        lower_bound = 0 * pq.nA
        trial_current = lower_bound
        spikes_found = False
        iteration = 0

        if plot:
            trace_spike = None
            trace_nospike = None

        # Binary search to find the boundaries that enclose the true rheobase
        # Terminate when: spikes under no stim, or narrowed to 1%
        while iteration < self.max_iterations and \
              upper_bound > 0*pq.nA and \
              (upper_bound-lower_bound)/upper_bound > self.min_precision:

            model.set_temperature(self.temperature)
            model.set_stop_time(self.ss_delay + self.current_duration)

            voltage = model.inject_square_current({"delay":     self.ss_delay,
                                                   "duration":  self.current_duration,
                                                   "amplitude": trial_current},
                                                  stop_on_spike=True)

            # Quickly check for APs without extracting other AP properties
            crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

            if len(crossings) >= 1:
                upper_bound = trial_current
                spikes_found = True

                if plot:
                    trace_spike = voltage
            else:
                lower_bound = trial_current

                if plot:
                    trace_nospike = voltage

            # Initially test if APs are produced within the maximum bounds
            if iteration > 0:
                # If fails to produce a spike at expected 0mV current, try it's double
                if trial_current == nA_to_0mV and len(crossings) < 1:
                    upper_bound = trial_current = nA_to_0mV * 2
                else:
                    trial_current = (upper_bound + lower_bound) / 2.0
            else:
                trial_current = upper_bound

            iteration += 1

        if plot:
            if trace_spike is not None:
                plt.plot(trace_spike.times, trace_spike)

            if trace_nospike is not None:
                plt.plot(trace_nospike.times, trace_nospike)

            plt.title(str(self) + " " + str(upper_bound))
            plt.show()

        if upper_bound == 0 * pq.nA:
            raise Exception("Negative rheobase: Model produces spikes without stimulation")

        if not spikes_found:
            raise Exception("Rheobase was not found between 0 and " + str(upper_bound))

        # If no rounding error will result from unit conversion
        if upper_bound == upper_bound.rescale(self.units):
            upper_bound.units = self.units

        # Otherwise correct the rounding error
        else:
            rescaled = upper_bound.rescale(self.units)
            unit_ratio = np.round((upper_bound.units / self.units).simplified.magnitude)
            rounding_error = upper_bound.magnitude * unit_ratio - rescaled.magnitude
            corrected = rescaled + rounding_error * self.units
            upper_bound = corrected

        return upper_bound


class RheobaseResponseTestHelper(OlfactoryBulbCellSpikeTest):
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,)

    def generate_prediction_nocache(self, model):
        if hasattr(self, "rheobase"):
            rheobase = self.rheobase
        else:
            rheobase = self.get_dependent_prediction(RheobaseTest, model)

        model.set_temperature(self.temperature)

        # Include a cool-off period to ensure any APs complete
        model.set_stop_time(self.ss_delay + self.current_duration + 100*pq.ms)

        # Inject rheobase - sampling at finer rate to capture AP shape better
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": rheobase,
                                               "sampling_period": 0.125})

        if plot:
            plt.plot(voltage)
            plt.title(str(self))
            plt.show()

        return voltage

class RheobaseSpikesTestHelper(OlfactoryBulbCellSpikeTest):
    units = pq.dimensionless
    score_type = scores.ZScore

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,)

    def generate_prediction_nocache(self, model):

        # Get voltage at rheobase
        voltage = self.get_dependent_prediction(RheobaseResponseTestHelper, model)

        # Quickly check for APs without extracting other AP properties
        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

        return len(crossings)

class FirstSpikeTestHelper(OlfactoryBulbCellSpikeTest):
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,)

    def get_first_ap(self, model):
        # Get voltage at rheobase
        voltage = self.get_dependent_prediction(RheobaseResponseTestHelper, model)

        aps = self.get_aps(voltage)

        if plot:
            plt.plot(voltage)
            plt.title(str(self))
            plt.show()

            for i, ap in enumerate(aps):
                plt.plot(ap["voltage"].times, ap["voltage"])
                plt.plot([ap["threshold_t"].rescale(pq.sec)], [ap["threshold_v"]], "o")
                plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["peak_v"], ap["peak_v"]])
                plt.title(str(self) + " AP " + str(i+1) + " of " + str(len(aps)))
                plt.xlim((ap["voltage"].times[0], ap["voltage"].times[0] + 30 * pq.ms))
                plt.show()

        if len(aps) < 1:
            raise Exception("The voltage trace does not contain any detectable action potentials using method: " + self.threshold_method)

        return aps[0]

    def generate_prediction_nocache(self, model):
        return self.get_first_ap(model)


class SpikeThresholdTest(FirstSpikeTestHelper):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the spike threshold"
    name = "Spike threshold test"

    def generate_prediction_nocache(self, model):
        ap = self.get_dependent_prediction(FirstSpikeTestHelper, model)

        if type(ap) == str and "Exception" in ap:
            return ap

        return ap["threshold_v"]


class SpikePeakTest(FirstSpikeTestHelper):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the spike peak voltage"
    name = "Spike peak test"

    def generate_prediction_nocache(self, model):
        ap = self.get_dependent_prediction(FirstSpikeTestHelper, model)

        if type(ap) == str and "Exception" in ap:
            return ap

        return ap["peak_v"]


class SpikeAmplitudeTest(FirstSpikeTestHelper):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the spike amplitude"
    name = "Spike amplitude test"

    def generate_prediction_nocache(self, model):
        ap = self.get_dependent_prediction(FirstSpikeTestHelper, model)

        if type(ap) == str and "Exception" in ap:
            return ap

        return ap["amplitude"]


class SpikeHalfWidthTest(FirstSpikeTestHelper):
    units = pq.ms
    score_type = scores.ZScore
    description = "A test of the spike width at half-amplitude"
    name = "Spike half-width test"

    def generate_prediction_nocache(self, model):
        ap = self.get_dependent_prediction(FirstSpikeTestHelper, model)

        if type(ap) == str and "Exception" in ap:
            return ap

        if plot:
            upper = np.where(ap["voltage"].magnitude > ap["threshold_v"] + 0.5 * ap["amplitude"])
            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.fill_between(ap["voltage"].times[upper[0]], ap["voltage"].magnitude[upper[0]].flatten(),
                             np.min(ap["voltage"]),
                             facecolor="orange",  # The fill color
                             color='blue',  # The outline color
                             alpha=0.2)
            plt.title(str(self) + " " + str(ap["half_width"]))
            plt.xlim((ap["voltage"].times[0], ap["voltage"].times[0] + 10 * pq.ms))
            plt.show()

        return ap["half_width"]


class SagVoltageTest(OlfactoryBulbCellTest):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the sag voltage of a cell"
    name = "Sag voltage test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,
                             obncap.SupportsVoltageClamp,)

    required_properties = [
        "sag_testing_voltage",
        "sag_window",
    ]

    def generate_prediction_nocache(self, model):

        # Get the current needed to reach the target sag testing voltage
        voltage, _ = self.get_dependent_prediction(PassivePropertyTestHelper, model)

        times = voltage.times

        # Find the minimum v within the ROI after SS and before end of injection
        t_roi = np.where((times > self.ss_delay) & (times <= self.ss_delay + self.current_duration))
        v_roi = voltage.base[t_roi]

        if plot:
            plt.plot(voltage.times, voltage)

        v_injection_ss = v_roi[-1]
        v_min_index = np.argmin(v_roi)
        v_min = v_roi[v_min_index]

        if plot:
            plt.plot([times[t_roi][-1]], [v_injection_ss], 'o', label="Current SS: " + str(v_injection_ss))
            plt.plot([times[t_roi][v_min_index]], [v_min], 'o', label="Min: " + str(v_min))
            plt.ylim((v_roi.min().magnitude - 5, v_roi.max().magnitude + 5))

        if v_min < v_injection_ss:
            sag = v_injection_ss - v_min

        # If there is no sag within the window, use the voltage at specific time point e.g. 100ms after current onset
        else:
            t_roi = np.where((times > self.ss_delay + self.sag_window - 1*pq.ms) & (times < self.ss_delay + self.sag_window + 1*pq.ms))
            v_roi = voltage.base[t_roi]
            v_at_sag_window = np.median(v_roi.magnitude) * pq.mV

            sag = v_injection_ss - v_at_sag_window

            if plot:
                plt.plot([np.median(times[t_roi])], [v_at_sag_window], 'o', label="Alt min: "+str(v_at_sag_window))


        if plot:
            plt.legend()
            plt.title(str(self) + " " + str(sag))
            plt.show()

        return sag


class ReboundSpikingTest(OlfactoryBulbCellSpikeTest):
    units = pq.dimensionless
    score_type = scores.BooleanScore
    description = "A test of the presence of rebound spikes"
    name = "Rebound spiking test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,
                             obncap.SupportsVoltageClamp,)

    required_properties = [
        "temperature",
        "ss_delay",
        "current_duration",
        "rebound_ap_method",
        "sag_testing_voltage",
        "rebound_rest_time",
    ]

    def generate_prediction_nocache(self, model):

        model.set_temperature(self.temperature)

        if self.rebound_ap_method == "sag":
            voltage, inhibitory_current = self.get_dependent_prediction(PassivePropertyTestHelper, model)

        elif self.rebound_ap_method == "-300pA":
            inhibitory_current = -0.3*pq.nA

            # Inject that current - but wait for additional time after end
            model.set_stop_time(self.ss_delay + self.current_duration + self.rebound_rest_time)

            voltage = model.inject_square_current({"delay":self.ss_delay,
                                                    "duration": self.current_duration,
                                                    "amplitude": inhibitory_current})

        # Look for APs after current is withdrawn
        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay + self.current_duration)

        if plot:
            plt.plot(voltage.times, voltage)
            plt.title(str(self) + " rebound APs: " + str(len(crossings)) + " at: " + str(inhibitory_current))
            plt.show()

        return len(crossings) > 0


class AfterHyperpolarizationTest(FirstSpikeTestHelper):

    required_properties = [
        "ahp_amplitude_method",
    ]

    def compute_amplitude(self, ap):
        if self.ahp_amplitude_method == 'threshold2min':
            roi_v = ap["voltage"].magnitude[np.where(ap["voltage"].times < ap["threshold_t"] + 100 * pq.ms)]
            i_min_v = np.argmin(roi_v)

        elif self.ahp_amplitude_method == 'threshold2minWithin10ms':
            within10ms = np.where(ap["voltage"].times < ap["threshold_t"] + 10*pq.ms)
            i_min_v = np.argmin(ap["voltage"].magnitude[within10ms])

        else:
            raise Exception("Unrecognized AHP Amplitude method: " + str(self.ahp_amplitude_method))

        min_v = ap["voltage"][i_min_v]
        ahp_amplitude = ap["threshold_v"] - min_v

        return ahp_amplitude


class AfterHyperpolarizationAmplitudeTest(AfterHyperpolarizationTest):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the AHP amplitude of a cell"
    name = "AHP amplitude"

    def generate_prediction_nocache(self, model):
        ap = self.get_dependent_prediction(FirstSpikeTestHelper, model)

        if type(ap) == str:
            return 0 * pq.mV

        amp = self.compute_amplitude(ap)

        if plot:
            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["threshold_v"], ap["threshold_v"]])
            plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["threshold_v"]-amp, ap["threshold_v"]-amp])
            plt.title(str(self) + " " + str(amp))
            plt.xlim((ap["voltage"].times[0], ap["voltage"].times[0] + 10 * pq.ms))
            plt.show()

            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["threshold_v"], ap["threshold_v"]])
            plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["threshold_v"]-amp, ap["threshold_v"]-amp])
            plt.title(str(self) + " " + str(amp))
            plt.xlim((ap["voltage"].times[0], ap["voltage"].times[0] + 110 * pq.ms))
            plt.show()

        return amp


class AfterHyperpolarizationTimeTest(AfterHyperpolarizationTest):
    units = pq.ms
    score_type = scores.ZScore
    description = "A test of the AHP duration of a cell"
    name = "AHP duration"

    required_properties = [
        "ahp_amplitude_method",
        "ahp_time_method",
    ]

    def generate_prediction_nocache(self, model):
        _, _, voltage = self.get_dependent_prediction(TargetFreqTestHelper, model)
        aps = self.get_aps(voltage)

        if len(aps) == 0:
            return 0 * self.units

        durations = np.vectorize(self.compute_duration)(aps)

        return np.mean(durations) * self.units

    def compute_duration(self, ap):
        if type(ap) == str:
            return 0 * pq.ms

        # Note - using crossing below threshold is prone to under-measuring the AHP
        # when the APs are very wide. The below code is replaced with measuring from the time of AP peak.

        # When the voltage droppes below the threshold voltage
        # crossings = get_zero_crossings_pos2neg(ap["voltage"] - ap["threshold_v"])
        # crossing = crossings[0]
        # crossing_time = ap["voltage"].times[crossing]

        # Starting at the time of AP peak
        crossing_time = ap["voltage"].times[np.argmax(ap["voltage"])]


        if self.ahp_time_method == 'threshold2min':
            # Look within 100ms of the AP onset
            roi_v = ap["voltage"].magnitude[np.where(ap["voltage"].times < ap["threshold_t"] + 100*pq.ms)]

            # Find the minimum within the ROI
            i_min_v = np.argmin(roi_v)

            # Get the time of the minimum
            min_v_time = ap["voltage"].times[i_min_v]

            # Compute the duration from AP onset to minimum
            ahp_time = (min_v_time - crossing_time).rescale(pq.ms)

        elif self.ahp_time_method == 'threshold2amplitude50%':
            # Compute the AHP amplitude
            amplitude = self.compute_amplitude(ap)

            # Compute the half amplitude
            amp50 = ap["threshold_v"] - amplitude / 2.0

            # Find where the voltage first crosses below it
            i_amp50 = np.where(ap["voltage"].magnitude[:,0] <= amp50)[0][0]

            # Get the time of the crossing
            t_amp50 = ap["voltage"].times[i_amp50]

            # Compute tye duration from AP onset
            ahp_time = (t_amp50 - crossing_time).rescale(pq.ms)

            # Make sure it's not negative
            ahp_time = max(0*pq.ms, ahp_time)

        else:
            raise Exception("Unrecognized AHP Time method: " + str(self.ahp_time_method))



        if plot:
            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["threshold_v"], ap["threshold_v"]])
            plt.plot([crossing_time], [ap["threshold_v"]], 'o')
            plt.axvline(x=crossing_time + ahp_time)
            plt.title(str(self) + " " + str(ahp_time) + " " + self.ahp_time_method)
            plt.xlim((ap["voltage"].times[0], ap["voltage"].times[0] + 10 * pq.ms))
            plt.show()

            plt.plot(ap["voltage"].times, ap["voltage"])
            plt.plot([ap["voltage"].times[0], ap["voltage"].times[-1]], [ap["threshold_v"], ap["threshold_v"]])
            plt.plot([crossing_time], [ap["threshold_v"]], 'o')
            plt.axvline(x=crossing_time + ahp_time)
            plt.title(str(self) + " " + str(ahp_time) + " " + self.ahp_time_method)
            plt.xlim((ap["voltage"].times[0], ap["voltage"].times[0] + 110 * pq.ms))
            plt.show()


        return ahp_time


class AfterDepolarizationResponseHelper(OlfactoryBulbCellSpikeTest):

    required_properties = [
        "adp_current_duration",
        "adp_current_amplitude",
    ]

    def generate_prediction_nocache(self, model):
        model.set_temperature(self.temperature)

        # Include a cool-off period to ensure any APs complete
        model.set_stop_time(self.ss_delay + self.adp_current_duration + 100*pq.ms)

        # Inject ADP current - sampling at finer rate to capture AP shape better
        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.adp_current_duration,
                                               "amplitude": self.adp_current_amplitude,
                                               "sampling_period": 0.125})

        return voltage

class AfterDepolarizationTimeTest(OlfactoryBulbCellSpikeTest):
    units = pq.ms
    score_type = scores.ZScore
    description = "A test of the ADP duration of a cell"
    name = "ADP duration"

    def generate_prediction_nocache(self, model):
        resting_v = self.get_dependent_prediction(RestingVoltageTest, model)

        voltage = self.get_dependent_prediction(AfterDepolarizationResponseHelper, model)

        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

        if len(crossings) != 1:
            raise Exception("Could not compute ADP duration because no APs were detected in the waveform")

        roi = np.where(voltage.times >= self.ss_delay.rescale(pq.ms))
        times = voltage.times[roi]
        voltage = voltage.magnitude[roi].flatten()

        max_v_i = np.argmax(voltage)
        max_v_t = times[max_v_i]

        crossing_t = times[crossings[0]]
        max_v = voltage[max_v_i] * pq.mV
        v_half_amp = (max_v - resting_v) / 2.0 + resting_v
        half_of_v_half = resting_v + (v_half_amp - resting_v) / 2.0
        half_of_v_half_t = times[np.where((times > crossing_t) & (voltage < half_of_v_half))][0]

        adp_duration = (half_of_v_half_t - crossing_t).rescale(pq.ms)

        if plot:
            plt.plot(times, voltage)
            plt.axhline(y=resting_v)
            plt.axhline(y=v_half_amp)
            plt.axhline(y=half_of_v_half)
            plt.axvline(x=half_of_v_half_t)
            plt.title(str(self) + " " + str(adp_duration))
            plt.xlim((crossing_t, half_of_v_half_t + 10 * pq.ms))
            plt.show()

        return adp_duration

class AfterDepolarizationDepthTest(OlfactoryBulbCellSpikeTest):
    units = pq.mV
    score_type = scores.ZScore
    description = "A test of the ADP depth of a cell"
    name = "ADP depth"

    def generate_prediction_nocache(self, model):
        '''
        For a post-AP waveform to be considered an "ADP", the AP should not have an AHP.
        In other words, the spike forms and then slowly descends back to resting Vm without
        crossing below the resting Vm (which is what happens with AHP)
        '''

        resting_v = self.get_dependent_prediction(RestingVoltageTest, model)

        voltage = self.get_dependent_prediction(AfterDepolarizationResponseHelper, model)

        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

        if len(crossings) != 1:
            raise Exception("Could not compute ADP duration because no APs were detected in the waveform")

        roi = np.where(voltage.times >= self.ss_delay.rescale(pq.ms))
        times = voltage.times[roi]
        voltage = voltage.magnitude[roi].flatten()

        min_v = np.min(voltage) * pq.mV

        adp_depth = resting_v - min_v

        if plot:
            crossing_t = times[crossings[0]]

            plt.plot(times, voltage)
            plt.axhline(y=resting_v)
            plt.axhline(y=min_v)
            plt.title(str(self) + " " + str(adp_depth))
            plt.xlim((crossing_t, crossing_t + 20 * pq.ms))
            plt.show()

        return adp_depth

class FISlopeTest(OlfactoryBulbCellSpikeTest):
    units = pq.Hz/pq.nA
    score_type = scores.ZScore
    description = "A test of the FI curve slope/gain"
    name = "FI slope test"
    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime,)

    required_properties = [
        "temperature",
        "ss_delay",
        "current_duration",
    ]

    def generate_prediction_nocache(self, model):

        model.set_temperature(self.temperature)
        model.set_stop_time(self.ss_delay + self.current_duration)

        # Current steps 50-300 pA in 50pA steps
        currents = np.concatenate((np.arange(50, 300, 50), [300])) * pq.pA
        frequencies = []

        for i, trial_current in enumerate(currents):
            # Inject current
            voltage = model.inject_square_current({"delay":self.ss_delay,
                                                   "duration": self.current_duration,
                                                   "amplitude": trial_current})

            # Count APs
            crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

            freq = len(crossings) / self.current_duration.rescale(pq.sec)

            frequencies.append(freq)

        fi_gain = np.max(np.diff(frequencies))
        slope = fi_gain / 50 * 1000 * pq.Hz / pq.nA # 50 is the step size

        if plot:
            plt.plot(currents, frequencies, 'o')
            plt.title(str(self) + " " + str(slope))
            plt.show()

        return slope


class SpikeTrainTest(OlfactoryBulbCellSpikeTest):

    def current_for_target_freq(self, model, rheobase, target_freq):
        model.set_stop_time(self.ss_delay + self.current_duration)

        freq_rb, _, voltage = self.freq_at(model, rheobase)

        if plot:
            plt.plot(voltage.times, voltage, label="1X RB")

        freq_2rb, _, voltage = self.freq_at(model, rheobase * 2)

        if plot:
            plt.plot(voltage.times, voltage, label="2X RB")
            plt.legend()
            plt.show()

        freq_slope = (freq_2rb - freq_rb) / rheobase

        if freq_rb != freq_2rb:
            current_targetHz = rheobase + (target_freq - freq_rb) / freq_slope
            return current_targetHz
        else:
            return None

    def freq_at(self, model, current):
        rheobase = self.get_dependent_prediction(RheobaseTest, model)

        # Reuse the rheobase response
        if current == rheobase:
            voltage = self.get_dependent_prediction(RheobaseResponseTestHelper, model)

        else:
            # Inject current
            voltage = model.inject_square_current({"delay": self.ss_delay,
                                                   "duration": self.current_duration,
                                                   "amplitude": current})

        # Count APs
        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)
        freq = len(crossings) / self.current_duration.rescale(pq.sec)

        return freq, crossings, voltage



class TargetFreqTestHelper(SpikeTrainTest):

    required_properties = [
        "spike_train_target_freq",
    ]

    def generate_prediction_nocache(self, model):

        rheobase = self.get_dependent_prediction(RheobaseTest, model)

        model.set_temperature(self.temperature)

        current_targetHz = self.current_for_target_freq(model,
                                                        rheobase,
                                                        self.spike_train_target_freq)

        if current_targetHz is not None:
            freq_targetHz, crossings, voltage = self.freq_at(model, current_targetHz)

            if plot:
                plt.plot(voltage.times, voltage)
                plt.title(str(self) + " Hz at target Hz current: " + str(freq_targetHz))
                plt.show()

            return freq_targetHz, crossings, voltage

        else:
            return None, None, None



class ISICVTest(SpikeTrainTest):
    units = pq.dimensionless
    score_type = scores.ZScore
    description = "A test of the interspike interval (ISI) coefficient of variation (CV)"
    name = "ISI CV test"

    required_properties = [
        "spike_train_target_freq",
    ]

    def generate_prediction_nocache(self, model):
        if self.spike_train_method == "target_freq":
            freq_targetHz, crossings, voltage = self.get_dependent_prediction(TargetFreqTestHelper, model)

            if freq_targetHz is None:
                return 0

        elif self.spike_train_method == "constant_current":
            crossings, voltage = self.get_dependent_prediction(SpikeAccommodationHelperConstantCurrentTest, model)

        if len(crossings) < 2:
            return 0

        isis = np.diff(crossings / voltage.sampling_rate)
        cv = np.std(isis) / np.mean(isis)

        if plot:
            plt.plot(voltage.times, voltage)
            plt.title(str(self))
            plt.show()

            plt.plot(isis, 'o')
            plt.title(str(self) + " " + str(cv))
            plt.show()

        return cv

class SpikeAccommodationHelperConstantCurrentTest(SpikeTrainTest):

    required_properties = [
        "temperature",
        "ss_delay",
        "current_duration",
        "spike_train_current",
    ]

    def generate_prediction_nocache(self, model):
        model.set_temperature(self.temperature)

        model.set_stop_time(self.ss_delay + self.current_duration)

        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": self.spike_train_current})

        # Count APs
        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

        if plot:
            plt.plot(voltage.times, voltage)
            plt.title(str(self))
            plt.show()

        return crossings, voltage

class SpikesAtConstantTrainCurrentTest(SpikeTrainTest):
    units = pq.dimensionless
    score_type = scores.ZScore

    def generate_prediction_nocache(self, model):
        crossings, _ = self.get_dependent_prediction(SpikeAccommodationHelperConstantCurrentTest, model)
        return len(crossings)

class SpikeAccommodationTest(SpikeTrainTest):
    units = pq.Hz
    score_type = scores.ZScore
    description = "A test of the firing rate accommodation"
    name = "Spike Accommodation Test"

    required_properties = [
        "spike_train_method",
    ]

    def generate_prediction_nocache(self, model):
        if self.spike_train_method == "target_freq":
            freq_targetHz, crossings, voltage = self.get_dependent_prediction(TargetFreqTestHelper, model)

            if freq_targetHz is None:
                return 0 * pq.Hz

        elif self.spike_train_method == "constant_current":
            crossings, voltage = self.get_dependent_prediction(SpikeAccommodationHelperConstantCurrentTest, model)

        else:
            raise Exception("Unrecognized spike_train_method: " + self.spike_train_method)

        isis = np.diff(crossings / voltage.sampling_rate)
        isis.units = pq.sec

        if len(isis) < 2:
            if plot:
                plt.plot(voltage.times, voltage)
                plt.title(str(self) + " Can't perform test: fewer than 2 APs")
                plt.show()

            return 0 * pq.Hz

        ifr_first = 1.0/isis[0]
        ifr_last =  1.0/isis[-1]

        ifr_change = ifr_last - ifr_first
        ifr_change.units = pq.Hz

        if plot:
            plt.plot(isis.rescale(pq.ms), 'o')
            plt.title(str(self) + " ISI plot")
            plt.show()

        return ifr_change


class SpikeAccommodationTimeConstantTest(SpikeTrainTest):
    units = pq.ms
    score_type = scores.ZScore
    description = "A test of the firing rate accommodation time constant"
    name = "Spike Accommodation Time Constant Test"

    def generate_prediction_nocache(self, model):

        if self.spike_train_method == "target_freq":
            freq_targetHz, crossings, voltage = self.get_dependent_prediction(TargetFreqTestHelper, model)

        elif self.spike_train_method == "constant_current":
            crossings, voltage = self.get_dependent_prediction(SpikeAccommodationHelperConstantCurrentTest, model)

        if crossings is None or len(crossings) < 4:
            return self.current_duration

        crossing_times = crossings / voltage.sampling_rate
        crossing_times -= crossing_times[0]

        isis = np.diff(crossing_times)
        isis.units = pq.sec

        ifrs = 1.0 / isis
        ifrs.units = pq.Hz
        ifrs = ifrs.magnitude
        ifr_times = (crossing_times - crossing_times[1])[1:]
        ifr_times = ifr_times.magnitude

        def ifr_func(t, start, finish, tau):
            return (start - finish) * np.exp(-t / tau) + finish

        from lmfit import Model

        model = Model(ifr_func)
        params = model.make_params(start=ifrs[0], finish=ifrs[-1], tau=10.0)
        params['tau'].min = 0
        result = model.fit(ifrs, t=ifr_times, params=params)

        start = result.best_values["start"]
        finish = result.best_values["finish"]
        tau = result.best_values["tau"] * pq.ms

        if plot:
            print(result.fit_report())

            plt.plot(crossing_times[1:], ifrs, 'bo')
            plt.plot(crossing_times[1:], result.best_fit, 'r-')
            plt.title(str(self) + " tau: " + str(tau))
            plt.show()

        return tau

class SpikesAtCurrentTest(OlfactoryBulbCellTest):
    units = pq.dimensionless
    score_type = scores.ZScore

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime)

    required_properties = [
        "temperature",
        "ss_delay",
        "current_duration",
        "current"
    ]

    def generate_prediction_nocache(self, model):

        model.set_temperature(self.temperature)

        model.set_stop_time(self.ss_delay + self.current_duration)

        voltage = model.inject_square_current({"delay": self.ss_delay,
                                               "duration": self.current_duration,
                                               "amplitude": self.current})

        crossings = get_zero_crossings_neg2pos(voltage, self.ss_delay)

        if plot:
            plt.plot(voltage.times, voltage)
            plt.title(str(self))
            plt.show()

        return len(crossings)