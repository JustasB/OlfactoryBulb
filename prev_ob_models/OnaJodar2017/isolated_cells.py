from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class GC1(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC1):
            from neuron import h,gui

            import GC_gna

            roots = h.SectionList()
            roots.allroots()
            roots = [r for r in roots]
            self.cell = roots[0]

            h.celsius=22

            self.h = h
            self.soma = self.cell

            h.cvode_active(1)

class GC2(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC2):
            from neuron import h,gui

            import GC_gna

            roots = h.SectionList()
            roots.allroots()
            roots = [r for r in roots]
            self.cell = roots[1]

            h.celsius=22

            self.h = h
            self.soma = self.cell

            h.cvode_active(1)

class GC3(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC3):
            from neuron import h,gui

            import GC_gna

            roots = h.SectionList()
            roots.allroots()
            roots = [r for r in roots]
            self.cell = roots[2]

            h.celsius=22

            self.h = h
            self.soma = self.cell

            h.cvode_active(1)

class GC4(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(GC4):
            from neuron import h,gui

            import GC_gna

            roots = h.SectionList()
            roots.allroots()
            roots = [r for r in roots]
            self.cell = roots[3]

            h.celsius=22

            self.h = h
            self.soma = self.cell

            h.cvode_active(1)
