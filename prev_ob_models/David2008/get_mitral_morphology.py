from neuron import h, gui

# Load a multi-compartment cell template
h.load_file('mitral.tem')

# Create an instance of the cell
cell1 = h.Mit(100)

h.define_shape()

from hoc2swc import neuron2swc
neuron2swc('mitral.swc')
quit()
