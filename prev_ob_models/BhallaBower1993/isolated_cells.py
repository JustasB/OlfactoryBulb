from prev_ob_models.utils import RunInClassDirectory

class MC():
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui

            h.xopen("mit_param.hoc")
            h.xopen("mit_morph.hoc")
            h.xopen("mit_memb.hoc")
            h.electrode_leak()

            h.celsius = 6.3

            h.cvode_active(1)

            self.h = h
            self.soma = h.soma
