from neuron import h

class PGC():
    def __init__(self):
        h.load_file('PG_def.hoc')

        self.h = h
        self.cell = h.PGcell(0)
        self.soma = self.cell.soma


class MC():
    def __init__(self):
        h.load_file('mitral.hoc')

        self.h = h
        self.cell = h.Mitral()
        self.soma = self.cell.soma

class GC():
    def __init__(self):
        h.load_file('granule.hoc')

        self.h = h
        self.cell = h.Granule()
        self.soma = self.cell.soma

class ETC():
    def __init__(self):
        h.load_file('et.hoc')

        self.h = h
        self.cell = h.ET()
        self.soma = self.cell.soma
