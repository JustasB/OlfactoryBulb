from neuron import h,gui

h("thresh = 10")

h.load_file('granule.tem')
h.load_file('mitral.tem')

gc = h.Gran(0)
mc = h.Mit(0)

from hoc2swc import neuron2swc
neuron2swc("cell.swc")
