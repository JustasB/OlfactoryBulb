from neuron import h

class PGC():
    def __init__(self):
        h.load_file('PG_Stim.hoc')

        self.h = h
        self.cell = h.pg
        self.soma = self.cell.soma

class MC():
    def __init__(self):
        h.load_file('MC_Stim.hoc')

        self.h = h
        self.cell = h.mit
        self.soma = self.cell.soma

class GC():
    def __init__(self):
        h.load_file('GC_Stim.hoc')

        self.h = h
        self.cell = h.Gran
        self.soma = self.cell.soma
