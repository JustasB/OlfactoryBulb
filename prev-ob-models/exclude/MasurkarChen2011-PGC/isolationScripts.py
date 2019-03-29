from neuron import h, gui

class PGC():
    def __init__(self):
        h.load_file('mosinit.hoc')
        h.IClamp[0].amp = 0

        self.h = h
        self.cell = h.soma
        self.soma = h.soma
