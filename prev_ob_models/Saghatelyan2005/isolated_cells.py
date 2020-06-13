from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h,gui

            h.load_file("mitral-occ.hoc")

            h.celsius=25

            self.cell = h.Mitral()

            self.h = h
            self.soma = self.cell.soma

            h.cvode_active(1)

class GC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h,gui

            h.load_file("gc-occ.hoc")

            h.celsius=25

            self.cell = h.GC()

            self.h = h
            self.soma = self.cell.somagc

            h.cvode_active(1)




