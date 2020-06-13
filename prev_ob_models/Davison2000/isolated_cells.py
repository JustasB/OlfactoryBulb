from prev_ob_models.utils import RunInClassDirectory

class MC():
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui

            h.xopen("mit4.hoc")

            h.IClamp[0].amp = 0
            h.IClamp[1].amp = 0

            h.celsius = 6.3

            h.cvode_active(1)

            self.h = h
            self.soma = h.soma
