from prev_ob_models.utils import IsolatedCell, RunInClassDirectory

class MC():
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui

            h.xopen("tabchannels.hoc")

            h("thresh = 10")

            h.load_file('mitral.tem')
            self.cell = h.Mit(0)

            h.celsius = 36

            h.cvode_active(1)

            self.h = h
            self.soma = self.cell.soma


class GC():
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h, gui

            h.xopen("tabchannels.hoc")

            h("thresh = 10")

            h.load_file('granule.tem')
            self.cell = h.Gran(0)

            h.celsius = 36

            h.cvode_active(1)

            self.h = h
            self.soma = self.cell.soma