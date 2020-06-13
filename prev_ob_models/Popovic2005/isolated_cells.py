from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h,gui

            h.load_file("nrngui.hoc")
            h.load_file("mulfit.hoc")
            h.load_file("parameters.hoc")
            h.load_file("morphology.hoc")
            h.load_file("membrane.hoc")
            h.celsius = 36

            self.h = h
            self.soma = self.h.soma_sc

            h.cvode_active(1)
