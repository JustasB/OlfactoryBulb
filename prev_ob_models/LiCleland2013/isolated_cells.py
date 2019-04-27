from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class PGC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(PGC):
            from neuron import h,gui
            h.load_file('PG_Stim.hoc')

            self.close_window("Graph")

            h.celsius = 35

            self.h = h
            self.cell = h.pg
            self.soma = self.cell.soma

            h.cvode_active(1)



class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h,gui
            h.load_file('MC_Stim.hoc')

            # Close various windows
            self.close_window("Graph")
            self.close_window("RunControl")

            h.celsius = 35

            self.h = h
            self.cell = h.mit
            self.soma = self.cell.soma

            h.cvode_active(1)

class GC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h,gui
            h.load_file('GC_Stim.hoc')

            self.close_window("Graph")

            h.celsius = 35

            self.h = h
            self.cell = h.Gran
            self.soma = self.cell.soma

            h.cvode_active(1)
