from neuron import h
startsw = h.startsw()
h.load_file("stdgui.hoc")
import params
pc = h.ParallelContext()
rank = int(pc.id())
from modeldata import getmodel

nhost = int(pc.nhost())
nmitral = params.Nmitral

import granules

ncell = nmitral + granules.ngranule


