from neuron import h
h.load_file("loadbal.hoc")
try:
    lb = h.LoadBalance()
except AttributeError:
    pass

