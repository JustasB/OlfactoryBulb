from prev_ob_models.utils import RunInClassDirectory

class MC():
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui

            h.load_file('init.hoc')

            h.celsius = 23

            h.cvode_active(1)

            self.h = h
            self.soma = h.soma
