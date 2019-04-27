from neuron import h, gui

h.load_file('mitral.hoc')
h.load_file('granule.hoc')

gc = h.Granule()
mc = h.Mitral()

from hoc2swc import neuron2swc
neuron2swc("cell.swc")