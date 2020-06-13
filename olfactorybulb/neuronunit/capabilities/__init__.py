# MOCKS for autodoc
import quantities as pq
if pq.mV.__class__.__module__ == 'sphinx.ext.autodoc.importer':
    pq.mV = pq.ms = pq.Hz = 1
# END MOCKS

import quantities as pq
import sciunit

class SupportsVoltageClamp(sciunit.Capability):
    """Indicates that the model can be held at three levels of voltages using a voltage clamp"""

    def clamp_voltage(self, voltages=[0*pq.mV]*3, durations=[0*pq.ms]*3):
        '''
        Maintains the model membrane potential for the specified durations at the specified voltages

        :param voltages: a 3-element array of voltages to clamp to
        :param durations: a 3-element array of durations to maintain the corresponding voltage levels
        :return: neo.core.AnalogSignal of the current required to keep the model at the specified voltages
        '''
        raise NotImplementedError()


class SupportsSettingTemperature(sciunit.Capability):
    """Indicates that the model can be executed using a specific temperature in Celsius"""

    def set_temperature(self, temperature=6.3):
        '''
        Specifies the simulator temperature

        :param temperature: the simulator temperature in degrees Celsius
        :return: Nothing
        '''
        raise NotImplementedError()


class SupportsSettingStopTime(sciunit.Capability):
    """Indicates that the model's simulation stop time can be specified"""

    def set_stop_time(self, tstop):
        '''
        Specifies the simulator stop time

        :param temperature: the simulator stop time in ms
        :return: Nothing
        '''
        raise NotImplementedError()
