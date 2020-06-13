from prev_ob_models.utils import RunInClassDirectory

class GC():
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h, gui

            h.xopen("mathslib.hoc")
            h.xopen("tabchannels.hoc")
            h.load_file('granule.tem')

            self.cell = h.Gran()

            h.celsius = 6.3

            h.cvode_active(1)

            self.h = h
            self.soma = self.cell.soma
