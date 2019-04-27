from neuron import h,gui

h.load_file('granule.tem')

gc = h.Gran()

from hoc2swc import neuron2swc
neuron2swc('granule.swc')
