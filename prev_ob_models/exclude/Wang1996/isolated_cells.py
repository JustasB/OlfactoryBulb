from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class PGC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(PGC):
            from neuron import h, gui
            h.load_file('PG_def.hoc')

            self.h = h
            self.cell = h.PGcell(0) # 0 is nicotinic current
            self.soma = self.cell.soma

            h.celsius = 6.3
            h.cvode_active(1)


class MC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui
            h.load_file('mitral.hoc')

            self.h = h
            self.cell = h.Mitral()
            self.soma = self.cell.soma

            h.celsius = 6.3
            h.cvode_active(1)

class GC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC):
            from neuron import h, gui
            h.load_file('granule.hoc')

            self.h = h
            self.cell = h.Granule()
            self.soma = self.cell.soma

            h.celsius = 6.3
            h.cvode_active(1)

class ETC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(ETC):
            from neuron import h, gui
            h.load_file('et.hoc')

            self.h = h
            self.cell = h.ET()
            self.soma = self.cell.soma

            h.celsius = 6.3
            h.cvode_active(1)

