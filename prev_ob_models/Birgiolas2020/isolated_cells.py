from prev_ob_models.utils import RunInClassDirectory, IsolatedCell
import os, sys

class MC(IsolatedCell):
    def __init__(self, mc_id):
        mc_id = str(mc_id)

        with RunInClassDirectory(MC):
            # Load the channels
            os.chdir("Mechanisms")
            from neuron import h#, gui
            os.chdir("..")

            h.load_file("stdrun.hoc")
            h.celsius = 35
            h.cvode_active(1)

            # Load the cell HOC file
            os.chdir("Cells")
            h.load_file("MC"+mc_id+".hoc")
            os.chdir("..")

            # Build the cell
            self.cell = getattr(h,"MC"+mc_id)()
            self.h = h
            self.soma = self.cell.soma

class MC1(MC):
    def __init__(self):
        super(MC1, self).__init__(mc_id=1)

class MC2(MC):
    def __init__(self):
        super(MC2, self).__init__(mc_id=2)

class MC3(MC):
    def __init__(self):
        super(MC3, self).__init__(mc_id=3)

class MC4(MC):
    def __init__(self):
        super(MC4, self).__init__(mc_id=4)

class MC5(MC):
    def __init__(self):
        super(MC5, self).__init__(mc_id=5)