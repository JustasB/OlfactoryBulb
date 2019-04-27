from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h,gui

            from mkmitral import mkmitral
            self.cell = mkmitral(0)

            h.celsius=35

            self.h = h
            self.soma = self.cell.soma

            h.cvode_active(1)

class GC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h,gui

            from net_mitral_centric import mkgranule
            import granules

            self.cell = mkgranule(0) 

            h.celsius=35

            self.h = h
            self.soma = self.cell.soma

            h.cvode_active(1)
