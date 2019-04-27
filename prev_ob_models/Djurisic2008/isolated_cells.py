from prev_ob_models.utils import IsolatedCell, RunInClassDirectory

class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h

            h.load_file("init_active.hoc")

            # Close cell builder window
            self.close_window("CellBuild")

            h.celsius = 6.3

            h.cvode_active(1)

            self.h = h
            self.soma = h.soma
