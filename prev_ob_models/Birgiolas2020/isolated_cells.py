from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h #,gui
            h.load_file("stdrun.hoc")

            h.celsius = 35
            h.cvode_active(1)

            h.load_file("MC1.hoc")
            self.cell = h.MC1()

            self.h = h
            self.soma = self.cell.soma




