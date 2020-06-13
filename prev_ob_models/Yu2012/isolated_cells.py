from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui
            h.load_file('mitral.hoc')

            self.h = h
            self.cell = h.Mitral()
            self.soma = self.cell.soma

            h.celsius = 35
            h.cvode_active(1)

class GC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h, gui
            h.load_file('granule.hoc')

            self.h = h
            self.cell = h.Granule()
            self.soma = self.cell.soma

            h.celsius = 35
            h.cvode_active(1)
