from neuron import h,gui

h.load_file('MC_def.hoc')
h.load_file('GC_def.hoc')

gc = h.Granule(0)
mc = h.Mitral(0)

from hoc2swc import neuron2swc
neuron2swc("cell.swc")