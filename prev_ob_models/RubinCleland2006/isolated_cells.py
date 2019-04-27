from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h,gui

            h.load_file("mitral.hoc")
            h.celsius = 6.3

            self.h = h
            self.soma = self.h.soma

            h.cvode_active(0)
