from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h,gui

            h.xopen("mitral.hoc")
            h.xopen("memb.hoc")

            h.celsius = 23

            self.h = h
            self.soma = self.h.soma

            h.cvode_active(1)




